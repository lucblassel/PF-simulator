//! Tree priors delegated to phylotree, with the Python dataset's rate/diameter pipeline.
use crate::tree::Tree;
use clap::{Args, ValueEnum};
use phylotree::{
    simulation::{self, BirthDeath, BranchLengths, Topology},
    tree::{NodeId, Tree as PhyloTree},
};
use rand::{Rng, SeedableRng, rngs::StdRng};
use rand_distr::{Distribution, Exp, LogNormal, Normal};
use std::{collections::HashMap, path::PathBuf};

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum TreePrior {
    Uniform,
    UniformRooted,
    Ete3,
    Yule,
    BirthDeath,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum BranchPrior {
    Exponential,
    Gamma,
    ProcessTime,
}

/// Named tree-generation pipelines. The CLI also applies alignment defaults.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum)]
pub enum Preset {
    /// Uniform unrooted topology and exponential branches with mean 0.1.
    #[default]
    Default,
    /// Python birth-death pipeline with inherited rates and empirical diameters.
    Python,
}

#[derive(Clone, Debug, Args)]
pub struct TreeArgs {
    /// Tree simulation strategy: default or python.
    #[arg(
        long,
        value_enum,
        default_value_t = Preset::Default,
        hide_possible_values = true,
        long_help = "Tree simulation strategy. In the CLI, python also defaults to LG + empirical G8 unless --model or --gamma overrides it.\n\n  default: uniform unrooted topology with exponential branch lengths (mean 0.1).\n  python: birth-death process times (birth 1, death 0.5), inherited rate changes, bundled empirical diameter scaling, then terminal repair below 0.001.\n\nNumeric settings and the empirical data file remain configurable. The python preset requires birth-death topology and process-time branches."
    )]
    pub preset: Preset,
    /// Topology prior (default: uniform unrooted; the python preset selects birth-death).
    #[arg(long, value_enum, hide_short_help = true, help_heading = "Tree priors")]
    pub tree_prior: Option<TreePrior>,
    /// Edge prior (default: process-time for birth-death, exponential otherwise).
    #[arg(long, value_enum, hide_short_help = true, help_heading = "Tree priors")]
    pub branch_prior: Option<BranchPrior>,
    /// Exponential mean (not rate), in substitutions per site.
    #[arg(
        long,
        default_value_t = 0.1,
        hide_short_help = true,
        help_heading = "Tree priors"
    )]
    pub branch_mean: f64,
    /// Gamma shape; its mean is shape * scale.
    #[arg(
        long,
        default_value_t = 2.0,
        hide_short_help = true,
        help_heading = "Tree priors"
    )]
    pub gamma_shape: f64,
    /// Gamma scale (not rate).
    #[arg(
        long,
        default_value_t = 0.05,
        hide_short_help = true,
        help_heading = "Tree priors"
    )]
    pub gamma_scale: f64,
    /// Birth rate per extant lineage.
    #[arg(
        long,
        default_value_t = 1.0,
        hide_short_help = true,
        help_heading = "Tree priors"
    )]
    pub birth_rate: f64,
    /// Death rate per extant lineage.
    #[arg(
        long,
        default_value_t = 0.5,
        hide_short_help = true,
        help_heading = "Tree priors"
    )]
    pub death_rate: f64,
    /// Compatibility shorthand for --preset python.
    #[arg(long, hide_short_help = true, help_heading = "Python rate process")]
    pub python_bd: bool,
    /// Apply Python's inherited lognormal rate-change process before diameter scaling.
    #[arg(long, hide_short_help = true, help_heading = "Python rate process")]
    pub python_rescaling: bool,
    /// Use each edge's actual time for rate changes (Python use_bl=True; default False).
    #[arg(long, hide_short_help = true, help_heading = "Python rate process")]
    pub use_branch_times: bool,
    /// Mean waiting time for small rate events (Python rate_small is a scale).
    #[arg(
        long,
        default_value_t = 0.03,
        hide_short_help = true,
        help_heading = "Python rate process"
    )]
    pub small_event_mean: f64,
    /// Lognormal sigma for small events (standard deviation of log multiplier).
    #[arg(
        long,
        default_value_t = 0.1,
        hide_short_help = true,
        help_heading = "Python rate process"
    )]
    pub small_log_sd: f64,
    /// Mean waiting time for large rate events.
    #[arg(
        long,
        default_value_t = 1.0,
        hide_short_help = true,
        help_heading = "Python rate process"
    )]
    pub large_event_mean: f64,
    /// Lognormal sigma for large events.
    #[arg(
        long,
        default_value_t = 0.2,
        hide_short_help = true,
        help_heading = "Python rate process"
    )]
    pub large_log_sd: f64,
    /// Draw jittered target diameters from the bundled empirical mixture.
    #[arg(long, hide_short_help = true, help_heading = "Diameter scaling")]
    pub empirical_diameter: bool,
    /// Override empirical data with a text file: diameter [integer weight] per line.
    #[arg(
        long,
        value_name = "FILE",
        hide_short_help = true,
        help_heading = "Diameter scaling"
    )]
    pub empirical_diameters: Option<PathBuf>,
    /// After diameter scaling, replace shorter terminal edges by N(min,0.005) truncated below min.
    /// Default is 0.001 for the python preset/python-rescaling, otherwise disabled.
    #[arg(long, hide_short_help = true, help_heading = "Diameter scaling")]
    pub minimum_terminal: Option<f64>,
}
impl Default for TreeArgs {
    fn default() -> Self {
        Self {
            preset: Preset::Default,
            tree_prior: None,
            branch_prior: None,
            branch_mean: 0.1,
            gamma_shape: 2.0,
            gamma_scale: 0.05,
            birth_rate: 1.0,
            death_rate: 0.5,
            python_bd: false,
            python_rescaling: false,
            use_branch_times: false,
            small_event_mean: 0.03,
            small_log_sd: 0.1,
            large_event_mean: 1.0,
            large_log_sd: 0.2,
            empirical_diameter: false,
            empirical_diameters: None,
            minimum_terminal: None,
        }
    }
}

#[derive(Debug)]
pub struct EmpiricalDiameters {
    values: Vec<f64>,
    cumulative: Vec<u64>,
    total: u64,
}
impl EmpiricalDiameters {
    pub fn parse(text: &str) -> Result<Self, String> {
        let mut values = Vec::new();
        let mut cumulative = Vec::new();
        let mut total = 0_u64;
        for (i, line) in text.lines().enumerate() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let parts: Vec<_> = line.split_whitespace().collect();
            let fail = || format!("Invalid empirical diameter at line {}", i + 1);
            if parts.is_empty() || parts.len() > 2 {
                return Err(fail());
            }
            let value = parts[0].parse::<f64>().map_err(|_| fail())?;
            let weight = if parts.len() == 2 {
                parts[1].parse::<u64>().map_err(|_| fail())?
            } else {
                1
            };
            if !value.is_finite() || value <= 0.02 || value >= 15.0 || weight == 0 {
                return Err(fail());
            }
            total = total.checked_add(weight).ok_or_else(fail)?;
            values.push(value);
            cumulative.push(total);
        }
        if values.is_empty() {
            return Err("Empty empirical diameter prior".into());
        }
        Ok(Self {
            values,
            cumulative,
            total,
        })
    }
    pub fn bundled() -> Self {
        Self::parse(include_str!("../data/empirical_diameters.tsv"))
            .expect("validated empirical data")
    }
    pub fn total_weight(&self) -> u64 {
        self.total
    }
    pub fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> f64 {
        // Draw once over the replicated empirical list, represented compactly by weights.
        let draw = rng.gen_range(0..self.total);
        let index = self.cumulative.partition_point(|&end| end <= draw);
        let mean = self.values[index];
        Normal::new(mean, mean / 10.0)
            .unwrap()
            .sample(rng)
            .max(0.02)
    }
}

#[derive(Debug)]
pub struct TreeSimulation {
    topology: Topology,
    branches: BranchLengths,
    rates: Option<RateVariation>,
    empirical: Option<EmpiricalDiameters>,
    minimum: f64,
}
impl TreeSimulation {
    pub fn new(args: &TreeArgs) -> Result<Self, String> {
        let python = args.python_bd || args.preset == Preset::Python;
        let python_option = if args.python_bd {
            "--python-bd"
        } else {
            "--preset python"
        };
        if python && args.tree_prior.is_some_and(|p| p != TreePrior::BirthDeath) {
            return Err(format!("{python_option} requires the birth-death topology"));
        }
        if python
            && args
                .branch_prior
                .is_some_and(|p| p != BranchPrior::ProcessTime)
        {
            return Err(format!("{python_option} requires process-time branches"));
        }
        let prior = args.tree_prior.unwrap_or(if python {
            TreePrior::BirthDeath
        } else {
            TreePrior::Uniform
        });
        let topology = match prior {
            TreePrior::Uniform => Topology::UniformUnrooted,
            TreePrior::UniformRooted => Topology::UniformRooted,
            TreePrior::Ete3 => Topology::Ete3,
            TreePrior::Yule => Topology::Yule,
            TreePrior::BirthDeath => Topology::BirthDeath(BirthDeath {
                birth_rate: args.birth_rate,
                death_rate: args.death_rate,
                ..Default::default()
            }),
        };
        let bp = args
            .branch_prior
            .unwrap_or(if prior == TreePrior::BirthDeath {
                BranchPrior::ProcessTime
            } else {
                BranchPrior::Exponential
            });
        let positive = |x: f64| x.is_finite() && x > 0.0;
        if !positive(args.branch_mean) || !(1e-8..=100.0).contains(&args.branch_mean) {
            return Err("Branch mean must be finite and between 1e-8 and 100".into());
        }
        if !positive(args.gamma_shape) || !positive(args.gamma_scale) {
            return Err("Gamma shape and scale must be finite and positive".into());
        }
        if !positive(args.birth_rate) || !args.death_rate.is_finite() || args.death_rate < 0.0 {
            return Err("Birth rate must be positive and death rate nonnegative".into());
        }
        if bp == BranchPrior::ProcessTime && prior != TreePrior::BirthDeath {
            return Err("Process-time branches require birth-death topology".into());
        }
        let branches = match bp {
            BranchPrior::Exponential => BranchLengths::Exponential {
                mean: args.branch_mean,
            },
            BranchPrior::Gamma => BranchLengths::Gamma {
                shape: args.gamma_shape,
                scale: args.gamma_scale,
            },
            BranchPrior::ProcessTime => BranchLengths::ProcessTime,
        };
        let alter = python || args.python_rescaling;
        let rates = alter.then_some(RateVariation {
            use_branch_times: args.use_branch_times,
            small_event_mean: args.small_event_mean,
            small_log_sd: args.small_log_sd,
            large_event_mean: args.large_event_mean,
            large_log_sd: args.large_log_sd,
        });
        if let Some(r) = &rates {
            r.validate()?;
        }
        if args.use_branch_times && !alter {
            return Err(
                "--use-branch-times requires --python-rescaling, --preset python or --python-bd"
                    .into(),
            );
        }
        let minimum = args
            .minimum_terminal
            .unwrap_or(if alter { 0.001 } else { 0.0 });
        if !minimum.is_finite() || minimum < 0.0 {
            return Err("Minimum terminal length must be finite and nonnegative".into());
        }
        let empirical = if let Some(file) = &args.empirical_diameters {
            Some(EmpiricalDiameters::parse(
                &std::fs::read_to_string(file).map_err(|e| format!("{}: {e}", file.display()))?,
            )?)
        } else if args.empirical_diameter || python {
            Some(EmpiricalDiameters::bundled())
        } else {
            None
        };
        Ok(Self {
            topology,
            branches,
            rates,
            empirical,
            minimum,
        })
    }
    pub fn generate(&self, tips: usize, seed: u64) -> Result<Tree, String> {
        self.generate_with_target(tips, seed).map(|(tree, _)| tree)
    }
    pub fn generate_with_target(
        &self,
        tips: usize,
        seed: u64,
    ) -> Result<(Tree, Option<f64>), String> {
        let mut rng = StdRng::seed_from_u64(seed);
        let target = self.empirical.as_ref().map(|e| e.sample(&mut rng));
        let mut tree = simulation::simulate(tips, self.topology, self.branches, &mut rng)
            .map_err(|e| e.to_string())?;
        if let Some(rates) = &self.rates {
            rates.apply(&mut tree, &mut rng)?;
        }
        if let Some(target) = target {
            tree.rescale_diameter(target).map_err(|e| e.to_string())?;
        }
        if self.minimum > 0.0 {
            let dist = Normal::new(self.minimum, 0.005).map_err(|e| e.to_string())?;
            for leaf in tree.leaves() {
                let mut length = tree
                    .get(leaf)
                    .map_err(|e| e.to_string())?
                    .parent_edge
                    .ok_or("Missing terminal branch")?;
                while length < self.minimum {
                    length = dist.sample(&mut rng);
                }
                set_edge(&mut tree, leaf, length)?;
            }
        }
        Ok((Tree::from_phylogeny(tree)?, target))
    }
}

#[derive(Clone, Copy, Debug)]
pub struct RateVariation {
    pub use_branch_times: bool,
    pub small_event_mean: f64,
    pub small_log_sd: f64,
    pub large_event_mean: f64,
    pub large_log_sd: f64,
}
impl RateVariation {
    fn validate(&self) -> Result<(), String> {
        if [self.small_event_mean, self.large_event_mean]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0)
            || [self.small_log_sd, self.large_log_sd]
                .iter()
                .any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err(
                "Event means must be positive; lognormal sigmas nonnegative; all finite".into(),
            );
        }
        Ok(())
    }
    pub fn apply<R: Rng + ?Sized>(&self, tree: &mut PhyloTree, rng: &mut R) -> Result<(), String> {
        self.validate()?;
        let root = tree.root().map_err(|e| e.to_string())?;
        let order = crate::tree::preorder(tree)?;
        let lengths = order
            .iter()
            .skip(1)
            .map(|&id| {
                let length = tree
                    .get(id)
                    .map_err(|e| e.to_string())?
                    .parent_edge
                    .ok_or("Missing branch time")?;
                if !length.is_finite() || length < 0.0 {
                    return Err("Invalid branch time".to_string());
                }
                Ok(length)
            })
            .collect::<Result<Vec<_>, String>>()?;
        let average = lengths.iter().sum::<f64>() / lengths.len() as f64;
        if !average.is_finite() || average <= 0.0 {
            return Err("Rate changes require a positive mean branch time".into());
        }
        let small = Exp::new(1.0 / self.small_event_mean).map_err(|e| e.to_string())?;
        let big = Exp::new(1.0 / self.large_event_mean).map_err(|e| e.to_string())?;
        let small_mult = LogNormal::new(0.0, self.small_log_sd).map_err(|e| e.to_string())?;
        let big_mult = LogNormal::new(0.0, self.large_log_sd).map_err(|e| e.to_string())?;
        let mut inherited = HashMap::from([(root, 1.0)]);
        let mut events = 0_u64;
        for id in order.into_iter().skip(1) {
            let n = tree.get(id).map_err(|e| e.to_string())?;
            let original = n.parent_edge.ok_or("Missing branch time")?;
            let d = if self.use_branch_times {
                original
            } else {
                average
            };
            let normalizer = if self.use_branch_times {
                1.0
            } else {
                original / average
            };
            let mut rate = inherited[&n.parent.ok_or("Missing parent")?];
            let mut last = 0.0;
            let mut integral = 0.0;
            while last < d {
                events += 1;
                if events > 10_000_000 {
                    return Err("Rate-change process exceeded its event limit".into());
                }
                let ts = small.sample(rng);
                let tb = big.sample(rng);
                let (wait, multiplier) = if ts < tb {
                    (ts, small_mult.sample(rng))
                } else {
                    (tb, big_mult.sample(rng))
                };
                let time = last + wait;
                integral += rate * (time.min(d) - last);
                if time >= d {
                    break;
                }
                rate *= multiplier;
                last = time;
                if !rate.is_finite() {
                    return Err("Rate process overflowed".into());
                }
            }
            let length = integral * normalizer;
            if !length.is_finite() || length < 0.0 {
                return Err("Invalid altered branch length".into());
            }
            inherited.insert(id, rate);
            set_edge(tree, id, length)?;
        }
        Ok(())
    }
}
fn set_edge(tree: &mut PhyloTree, id: NodeId, length: f64) -> Result<(), String> {
    let parent = tree
        .get(id)
        .map_err(|e| e.to_string())?
        .parent
        .ok_or("Cannot assign root stem")?;
    tree.get_mut(id).map_err(|e| e.to_string())?.parent_edge = Some(length);
    tree.get_mut(parent)
        .map_err(|e| e.to_string())?
        .set_child_edge(id, Some(length));
    Ok(())
}
