use pf_simulator::tree_simulation::{
    BranchPrior, EmpiricalDiameters, Preset, RateVariation, TreeArgs, TreePrior, TreeSimulation,
};
use phylotree::simulation::{self, BirthDeath, BranchLengths, Topology};
use rand::{SeedableRng, rngs::StdRng};

#[test]
fn presets_preserve_existing_pipelines_and_numeric_customization() {
    for customize in [false, true] {
        let mut args = TreeArgs {
            preset: Preset::Python,
            ..Default::default()
        };
        if customize {
            args.birth_rate = 2.0;
            args.death_rate = 0.25;
            args.use_branch_times = true;
            args.small_log_sd = 0.05;
            args.minimum_terminal = Some(0.0);
        }
        let preset = TreeSimulation::new(&args).unwrap();
        let legacy = TreeSimulation::new(&TreeArgs {
            preset: Preset::Default,
            python_bd: true,
            ..args
        })
        .unwrap();
        for seed in 0..10 {
            let (a, a_target) = preset.generate_with_target(30, seed).unwrap();
            let (b, b_target) = legacy.generate_with_target(30, seed).unwrap();
            assert_eq!(a.newick(), b.newick());
            assert_eq!(a_target, b_target);
        }
    }
    let baseline = TreeSimulation::new(&TreeArgs::default()).unwrap();
    let changed_mean = TreeSimulation::new(&TreeArgs {
        preset: Preset::Default,
        branch_mean: 0.2,
        ..Default::default()
    })
    .unwrap();
    let a = baseline.generate(20, 18).unwrap();
    let b = changed_mean.generate(20, 18).unwrap();
    let a = a.phylogeny();
    let b = b.phylogeny();
    for id in a.preorder(a.root().unwrap()).unwrap() {
        assert_eq!(
            a.get(id).unwrap().parent_edge.map(|length| length * 2.0),
            b.get(id).unwrap().parent_edge
        );
    }
    for args in [
        TreeArgs {
            tree_prior: Some(TreePrior::Yule),
            ..Default::default()
        },
        TreeArgs {
            branch_prior: Some(BranchPrior::Gamma),
            ..Default::default()
        },
    ] {
        let error = TreeSimulation::new(&TreeArgs {
            preset: Preset::Python,
            ..args
        })
        .unwrap_err();
        assert!(error.contains("--preset python requires"));
    }
}

#[test]
fn empirical_mixture_and_normal_jitter_match_the_python_definition() {
    assert_eq!(EmpiricalDiameters::bundled().total_weight(), 103120);
    let prior = EmpiricalDiameters::parse("1 1\n10 3\n").unwrap();
    let mut rng = StdRng::seed_from_u64(342);
    let n = 50_000;
    let big = (0..n).filter(|_| prior.sample(&mut rng) > 5.0).count();
    assert!((big as f64 / n as f64 - 0.75).abs() < 0.01);
    let one = EmpiricalDiameters::parse("1\n").unwrap();
    let x: Vec<_> = (0..n).map(|_| one.sample(&mut rng)).collect();
    let mean = x.iter().sum::<f64>() / n as f64;
    let var = x.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n as f64;
    assert!((mean - 1.0).abs() < 0.002);
    assert!((var - 0.01).abs() < 0.0005);
    for bad in [
        "",
        "NaN",
        "0.02",
        "15",
        "1 0",
        "1 -1",
        "1 18446744073709551615\n2 1",
    ] {
        assert!(EmpiricalDiameters::parse(bad).is_err());
    }
}
#[test]
fn zero_log_variance_preserves_times_and_missing_times_fail_cleanly() {
    for use_branch_times in [true, false] {
        let mut rng = StdRng::seed_from_u64(41);
        let mut t = simulation::simulate(
            20,
            Topology::BirthDeath(BirthDeath::default()),
            BranchLengths::ProcessTime,
            &mut rng,
        )
        .unwrap();
        let original = t.clone();
        let rates = RateVariation {
            use_branch_times,
            small_event_mean: 0.03,
            small_log_sd: 0.0,
            large_event_mean: 1.0,
            large_log_sd: 0.0,
        };
        rates.apply(&mut t, &mut rng).unwrap();
        for id in original.preorder(original.root().unwrap()).unwrap() {
            assert!(
                (original.get(id).unwrap().parent_edge.unwrap_or(0.0)
                    - t.get(id).unwrap().parent_edge.unwrap_or(0.0))
                .abs()
                    < 1e-12
            );
        }
        let mut missing = phylotree::tree::Tree::from_newick("(A,B,C);").unwrap();
        assert!(rates.apply(&mut missing, &mut rng).is_err());
    }
}
#[test]
fn full_pipeline_is_reproducible_preserves_floor_order_and_diameter_without_floor() {
    for prior in [
        TreePrior::Uniform,
        TreePrior::UniformRooted,
        TreePrior::Ete3,
        TreePrior::Yule,
        TreePrior::BirthDeath,
    ] {
        let generator = TreeSimulation::new(&TreeArgs {
            tree_prior: Some(prior),
            empirical_diameter: true,
            minimum_terminal: Some(0.0),
            ..Default::default()
        })
        .unwrap();
        let (a, target) = generator.generate_with_target(30, 72).unwrap();
        let b = generator.generate(30, 72).unwrap();
        assert_eq!(a.newick(), b.newick());
        assert!((a.phylogeny().diameter().unwrap() - target.unwrap()).abs() < 1e-12);
    }
    let generator = TreeSimulation::new(&TreeArgs {
        python_bd: true,
        ..Default::default()
    })
    .unwrap();
    for seed in 0..20 {
        let (t, target) = generator.generate_with_target(20, seed).unwrap();
        let p = t.phylogeny();
        for id in p.leaves() {
            assert!(p.get(id).unwrap().parent_edge.unwrap() >= 0.001);
        }
        assert!(p.diameter().unwrap() + 1e-12 >= target.unwrap());
    }
    assert!(
        TreeSimulation::new(&TreeArgs {
            python_bd: true,
            tree_prior: Some(TreePrior::Ete3),
            ..Default::default()
        })
        .is_err()
    );
    assert!(
        TreeSimulation::new(&TreeArgs {
            python_bd: true,
            branch_prior: Some(BranchPrior::Gamma),
            ..Default::default()
        })
        .is_err()
    );
}
