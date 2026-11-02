//! Controlled single-thread ablation; does not change production simulation behavior.
use clap::Parser;
use pf_simulator::{
    alignment::{AlignmentConfig, simulate_alignment_from_tree},
    selection::{TipChecks, WindowPolicy, select},
    site_rates::GammaMode,
    tree_simulation::{Preset, TreeArgs, TreeSimulation},
};
use std::{
    error::Error,
    fs,
    hint::black_box,
    io::Write,
    path::PathBuf,
    process::{Command, Stdio},
    time::Instant,
};

#[derive(Parser)]
struct Args {
    #[arg(long)]
    output: PathBuf,
    #[arg(long)]
    reference: PathBuf,
    #[arg(long, default_value_t = 3)]
    rounds: u32,
    /// Site-rate modes to compare, using paired empirical shape draws.
    #[arg(long, value_delimiter = ',', default_value = "none,GC")]
    gammas: Vec<GammaMode>,
}
fn sample_seed(seed: u64, index: u64) -> u64 {
    let mut z = seed.wrapping_add(index.wrapping_add(1).wrapping_mul(0x9e3779b97f4a7c15));
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
    z ^ (z >> 31)
}
fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();
    fs::create_dir_all(&args.output)?;
    let work = args.output.join("work");
    fs::create_dir(&work)?;
    let mut csv = fs::File::create(args.output.join("runs.csv"))?;
    writeln!(
        csv,
        "tips,sites,round,candidate,gamma,ancestral,model,tree_s,engine_s,gapped_s,window_s,full_format_s,extra_format_s,write_all_s,standalone_s,accepted_gapped,accepted_window,columns,rows,bytes"
    )?;
    let generator = TreeSimulation::new(&TreeArgs {
        preset: Preset::Python,
        ..Default::default()
    })?;
    for (tips, sites, samples) in [(50, 500, 16), (200, 5000, 4)] {
        for round in 0..args.rounds {
            let master = 20261102 + u64::from(round);
            let mut variants: Vec<_> = args
                .gammas
                .iter()
                .flat_map(|&gamma| [(gamma, false), (gamma, true)])
                .collect();
            if round % 2 == 1 {
                variants.reverse();
            }
            for (gamma, ancestral) in variants {
                eprintln!(
                    "{tips} tips × {sites} sites, round {round}, {gamma}, ancestors={ancestral}"
                );
                for candidate in 0..samples {
                    let start = Instant::now();
                    let tree_seed = sample_seed(master, candidate);
                    let tree = generator.generate(tips, tree_seed)?;
                    let tree_s = start.elapsed().as_secs_f64();
                    let seed = 1 + (sample_seed(tree_seed, 0) % 2147483646) as i32;
                    let model = gamma.model("LG", master, candidate, 0);
                    let config = AlignmentConfig {
                        length: sites,
                        model: model.clone(),
                        ancestral,
                        indel: Some("0.002,0.002".into()),
                        indel_size: Some("POW{1.7/100},POW{1.7/100}".into()),
                        ..Default::default()
                    };
                    fs::write(work.join("input.nwk"), tree.newick())?;
                    let prefix = work.join("native");
                    let mut native = Command::new(&args.reference);
                    native
                        .args(["--alisim"])
                        .arg(&prefix)
                        .args(["-t"])
                        .arg(work.join("input.nwk"))
                        .args([
                            "-m",
                            &model,
                            "--seqtype",
                            "AA",
                            "--length",
                            &sites.to_string(),
                            "--seed",
                            &seed.to_string(),
                            "--indel",
                            "0.002,0.002",
                            "--indel-size",
                            "POW{1.7/100},POW{1.7/100}",
                            "--out-format",
                            "fasta",
                            "-nt",
                            "1",
                            "-redo",
                        ])
                        .stdout(Stdio::null())
                        .stderr(Stdio::null());
                    if ancestral {
                        native.arg("--write-all");
                    }
                    // Alternate engine order; standalone includes process startup and native FASTA I/O.
                    let mut standalone_s = 0.0;
                    if round % 2 == 0 {
                        let start = Instant::now();
                        if !native.status()?.success() {
                            return Err("standalone failed".into());
                        }
                        standalone_s = start.elapsed().as_secs_f64();
                    }
                    let start = Instant::now();
                    let raw = simulate_alignment_from_tree(
                        &tree,
                        &work.join("numeric"),
                        &config,
                        seed,
                        None,
                    )?;
                    let engine_s = start.elapsed().as_secs_f64();
                    if round % 2 == 1 {
                        let start = Instant::now();
                        if !native.status()?.success() {
                            return Err("standalone failed".into());
                        }
                        standalone_s = start.elapsed().as_secs_f64();
                    }
                    let checks = TipChecks {
                        gapped: true,
                        ungapped: false,
                        window: false,
                        window_policy: WindowPolicy::Python,
                    };
                    // Recheck the identical numeric data ten times to stabilize sub-ms timings.
                    let start = Instant::now();
                    for _ in 0..10 {
                        let _ = black_box(select(&raw, tips, sites as usize, checks)?);
                    }
                    let gapped_s = start.elapsed().as_secs_f64() / 10.0;
                    let start = Instant::now();
                    for _ in 0..10 {
                        let _ = black_box(select(
                            &raw,
                            tips,
                            sites as usize,
                            TipChecks {
                                window: true,
                                ..checks
                            },
                        )?);
                    }
                    let window_s = start.elapsed().as_secs_f64() / 10.0;
                    let accepted_gapped = select(&raw, tips, sites as usize, checks)?.is_ok();
                    let window = select(
                        &raw,
                        tips,
                        sites as usize,
                        TipChecks {
                            window: true,
                            ..checks
                        },
                    )?;
                    // Format/write all candidates even if rejected: keep I/O amount independent of checks.
                    let start = Instant::now();
                    let full = raw.to_fasta(false);
                    let full_format_s = start.elapsed().as_secs_f64();
                    let start = Instant::now();
                    let ungapped = raw.to_fasta(true);
                    let start_col = window.as_ref().map_or(0, |w| w.start);
                    let trimmed = raw.window_fasta(
                        start_col,
                        (start_col + sites as usize).min(raw.rows[0].len()),
                    );
                    let extra_format_s = start.elapsed().as_secs_f64();
                    let start = Instant::now();
                    fs::write(work.join("full.fa"), &full)?;
                    fs::write(work.join("trimmed.fa"), &trimmed)?;
                    fs::write(work.join("ungapped.fa"), &ungapped)?;
                    let write_all_s = start.elapsed().as_secs_f64();
                    writeln!(
                        csv,
                        "{tips},{sites},{round},{candidate},{gamma},{ancestral},{model},{tree_s},{engine_s},{gapped_s},{window_s},{full_format_s},{extra_format_s},{write_all_s},{standalone_s},{accepted_gapped},{},{},{},{}",
                        window.is_ok(),
                        raw.rows[0].len(),
                        raw.rows.len(),
                        full.len() + trimmed.len() + ungapped.len()
                    )?;
                    csv.flush()?;
                }
            }
        }
    }
    fs::remove_dir_all(work)?;
    Ok(())
}
