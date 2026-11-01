//! Export independent tree-level statistics from the same generator as pf-sim.
use clap::Parser;
use pf_simulator::{
    tree::sample_seed,
    tree_simulation::{TreeArgs, TreeSimulation},
};
use std::{
    collections::HashMap,
    error::Error,
    fs::File,
    io::{BufWriter, Write},
    path::PathBuf,
};
#[derive(Parser)]
struct Args {
    #[arg(long, default_value_t = 3000)]
    samples: u64,
    #[arg(long, default_value_t = 50)]
    tips: usize,
    #[arg(long, default_value_t = 271828)]
    seed: u64,
    #[arg(long)]
    output: PathBuf,
    #[command(flatten)]
    tree: TreeArgs,
}
fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();
    let simulator = TreeSimulation::new(&args.tree)?;
    let mut out = BufWriter::new(File::create(args.output)?);
    writeln!(
        out,
        "sample,diameter,total_length,crown_height,height_span,terminal_mean,terminal_cv,internal_mean,min_terminal,zero_terminals,cherries,sackin,colless,root_balance,parent_child_correlation,diameter_target,diameter_increase"
    )?;
    for i in 0..args.samples {
        let (wrapper, target) =
            simulator.generate_with_target(args.tips, sample_seed(args.seed, i))?;
        let t = wrapper.phylogeny();
        let root = t.root()?;
        let order = t.preorder(root)?;
        let mut heights = HashMap::from([(root, 0.0)]);
        let mut terminal = Vec::new();
        let mut internal = Vec::new();
        let mut pairs = Vec::new();
        let mut leaf_heights = Vec::new();
        let mut sackin = 0;
        let mut cherries = 0;
        for &id in &order {
            let n = t.get(id)?;
            if let Some(p) = n.parent {
                let edge = n.parent_edge.unwrap();
                heights.insert(id, heights[&p] + edge);
                if p != root {
                    pairs.push((edge, t.get(p)?.parent_edge.unwrap()));
                }
            }
            if n.children.is_empty() {
                terminal.push(n.parent_edge.unwrap());
                leaf_heights.push(heights[&id]);
                sackin += n.get_depth();
            } else {
                if id != root {
                    internal.push(n.parent_edge.unwrap());
                }
                if n.children.len() == 2
                    && n.children
                        .iter()
                        .all(|&id| t.get(id).unwrap().children.is_empty())
                {
                    cherries += 1;
                }
            }
        }
        let mut counts = HashMap::new();
        let mut colless = 0;
        for &id in order.iter().rev() {
            let n = t.get(id)?;
            let sizes: Vec<usize> = n.children.iter().map(|c| counts[c]).collect();
            if sizes.len() == 2 {
                colless += sizes[0].abs_diff(sizes[1]);
            }
            counts.insert(
                id,
                if sizes.is_empty() {
                    1
                } else {
                    sizes.iter().sum()
                },
            );
        }
        let root_balance = t
            .get(root)?
            .children
            .iter()
            .map(|c| counts[c])
            .min()
            .unwrap() as f64
            / args.tips as f64;
        let mean = |v: &[f64]| v.iter().sum::<f64>() / v.len() as f64;
        let tm = mean(&terminal);
        let cv = if tm > 0.0 {
            (terminal.iter().map(|v| (v - tm).powi(2)).sum::<f64>() / terminal.len() as f64).sqrt()
                / tm
        } else {
            0.0
        };
        let high = leaf_heights.iter().copied().fold(0.0, f64::max);
        let low = leaf_heights.iter().copied().fold(f64::INFINITY, f64::min);
        let span = if (high - low) < 1e-10 * high.max(1.0) {
            0.0
        } else {
            high - low
        };
        let x = pairs.iter().map(|p| p.0).sum::<f64>() / pairs.len() as f64;
        let y = pairs.iter().map(|p| p.1).sum::<f64>() / pairs.len() as f64;
        let covariance = pairs.iter().map(|p| (p.0 - x) * (p.1 - y)).sum::<f64>();
        let vx = pairs.iter().map(|p| (p.0 - x).powi(2)).sum::<f64>();
        let vy = pairs.iter().map(|p| (p.1 - y).powi(2)).sum::<f64>();
        let corr = if vx * vy > 0.0 {
            covariance / (vx * vy).sqrt()
        } else {
            0.0
        };
        let diam = t.diameter()?;
        let increase = target.map_or(f64::NAN, |target| {
            let delta = diam - target;
            if delta.abs() < 1e-10 * diam.max(1.0) {
                0.0
            } else {
                delta
            }
        });
        let minimum = terminal.iter().copied().fold(f64::INFINITY, f64::min);
        let zeros = terminal.iter().filter(|&&v| v == 0.0).count();
        writeln!(
            out,
            "{i},{diam},{},{high},{span},{tm},{cv},{},{minimum},{zeros},{cherries},{sackin},{colless},{root_balance},{corr},{},{}",
            terminal.iter().chain(&internal).sum::<f64>(),
            mean(&internal),
            target.unwrap_or(f64::NAN),
            increase
        )?;
    }
    out.flush()?;
    Ok(())
}
