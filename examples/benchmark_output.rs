//! Compare the reference FASTA transport with the default numeric transport.
use std::collections::BTreeMap;
use std::error::Error;
use std::fs;
use std::hint::black_box;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

use clap::Parser;
use pf_simulator::alignment::{
    AlignmentConfig, SequenceType, has_duplicate_sequences, simulate_alignment_in_memory,
    simulate_alignment_raw_benchmark,
};
use pf_simulator::tree::Tree;

#[derive(Parser)]
struct Args {
    #[arg(long, default_value_t = 7, value_parser = clap::value_parser!(u32).range(3..))]
    rounds: u32,
    #[arg(long, default_value_t = 120)]
    target_ms: u64,
    #[arg(long)]
    output: Option<PathBuf>,
    /// Only run the 10-tip smoke-test workloads.
    #[arg(long)]
    quick: bool,
    /// Run selected cases: LG, LG+GC{0.5}, JC, indel, or indel_LG.
    #[arg(long, value_delimiter = ',')]
    cases: Vec<String>,
    /// Taxa x root sites, separated by commas.
    #[arg(long, value_delimiter = ',', default_value = "10x500,50x500,200x5000")]
    workloads: Vec<String>,
    #[arg(long, default_value_t = 0.1)]
    branch_mean: f64,
}

#[derive(Default)]
struct Timing {
    simulation: f64,
    check: f64,
    serialization: f64,
    pipeline: f64,
    duplicates: usize,
    flags: Vec<bool>,
}

fn records(bytes: &[u8]) -> BTreeMap<String, Vec<u8>> {
    let mut result = BTreeMap::new();
    let mut name = String::new();
    for line in bytes.split(|&byte| byte == b'\n').map(<[u8]>::trim_ascii) {
        if let Some(header) = line.strip_prefix(b">") {
            name = String::from_utf8(header.trim_ascii().to_vec()).unwrap();
            result.insert(name.clone(), Vec::new());
        } else if !line.is_empty() {
            result.get_mut(&name).unwrap().extend_from_slice(line);
        }
    }
    result
}

fn verify(
    tree: &Path,
    prefix: &Path,
    config: &AlignmentConfig,
    tips: usize,
) -> Result<(), Box<dyn Error>> {
    let fasta = simulate_alignment_in_memory(tree, prefix, config, 123, None)?;
    let raw = simulate_alignment_raw_benchmark(tree, prefix, config, 123)?;
    assert_eq!(raw.rows.len(), tips);
    assert_eq!(raw.names.len(), tips);
    assert!(raw.rows.iter().all(|row| row.len() == raw.rows[0].len()));
    assert_eq!(records(&fasta.aligned_fasta), records(&raw.to_fasta(false)));
    assert_eq!(
        has_duplicate_sequences(&fasta.aligned_fasta)?,
        raw.has_duplicates()
    );
    if let Some(unaligned) = fasta.unaligned_fasta {
        assert_eq!(records(&unaligned), records(&raw.to_fasta(true)));
    }
    Ok(())
}

fn measure(
    raw: bool,
    tree: &Path,
    prefix: &Path,
    config: &AlignmentConfig,
    count: usize,
    first_seed: i32,
) -> Result<Timing, Box<dyn Error>> {
    let mut timing = Timing::default();
    let all = Instant::now();
    for index in 0..count {
        let start = Instant::now();
        let duplicate = if raw {
            let data =
                simulate_alignment_raw_benchmark(tree, prefix, config, first_seed + index as i32)?;
            timing.simulation += start.elapsed().as_secs_f64();
            let start = Instant::now();
            let duplicate = black_box(data.has_duplicates());
            timing.check += start.elapsed().as_secs_f64();
            if !duplicate {
                let start = Instant::now();
                black_box(data.to_fasta(false));
                if config.indel.is_some() && !config.no_unaligned {
                    black_box(data.to_fasta(true));
                }
                timing.serialization += start.elapsed().as_secs_f64();
            }
            duplicate
        } else {
            let data = simulate_alignment_in_memory(
                tree,
                prefix,
                config,
                first_seed + index as i32,
                None,
            )?;
            timing.simulation += start.elapsed().as_secs_f64();
            let start = Instant::now();
            let duplicate = black_box(has_duplicate_sequences(&data.aligned_fasta)?);
            timing.check += start.elapsed().as_secs_f64();
            black_box(&data);
            duplicate
        };
        timing.duplicates += usize::from(duplicate);
        timing.flags.push(duplicate);
    }
    timing.pipeline = all.elapsed().as_secs_f64();
    Ok(timing)
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();
    let models = ["LG", "LG+GC{0.5}", "JC", "indel", "indel_LG"];
    if args
        .cases
        .iter()
        .any(|case| !models.contains(&case.as_str()))
    {
        return Err("Unknown benchmark case".into());
    }
    let timestamp = SystemTime::now().duration_since(UNIX_EPOCH)?.as_secs();
    let output = args
        .output
        .unwrap_or_else(|| format!("benchmarks/output-comparison-{timestamp}.csv").into());
    if let Some(parent) = output.parent().filter(|path| !path.as_os_str().is_empty()) {
        fs::create_dir_all(parent)?;
    }
    let mut csv = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&output)?;
    writeln!(
        csv,
        "case,tips,root_sites,branch_mean,mode,round,iterations,simulation_ms,check_ms,accepted_serialization_ms,pipeline_ms,rejections"
    )?;
    let directory = std::env::temp_dir().join(format!(
        "pf-output-bench-{}-{timestamp}",
        std::process::id()
    ));
    fs::create_dir(&directory)?;
    let sizes = if args.quick {
        vec![(10, 500)]
    } else {
        args.workloads
            .iter()
            .map(|workload| -> Result<(usize, u32), Box<dyn Error>> {
                let (tips, sites) = workload
                    .split_once('x')
                    .ok_or("Expected workload like 50x500")?;
                Ok((tips.parse()?, sites.parse()?))
            })
            .collect::<Result<Vec<_>, _>>()?
    };
    for (tips, length) in sizes {
        let tree = directory.join("tree.nwk");
        fs::write(&tree, Tree::uniform(tips, args.branch_mean, 42)?.newick())?;
        let prefix = directory.join("virtual");
        for model in models {
            if !args.cases.is_empty() && !args.cases.iter().any(|case| case == model) {
                continue;
            }
            let config = AlignmentConfig {
                length,
                model: match model {
                    "indel" => "LG+GC{0.5}",
                    "indel_LG" => "LG",
                    _ => model,
                }
                .into(),
                sequence_type: if model == "JC" {
                    SequenceType::Dna
                } else {
                    SequenceType::AminoAcid
                },
                indel: model.starts_with("indel").then(|| "0.002,0.002".into()),
                indel_size: model
                    .starts_with("indel")
                    .then(|| "POW{1.7/100},POW{1.7/100}".into()),
                no_unaligned: false,
            };
            verify(&tree, &prefix, &config, tips)?;
            let estimate = measure(false, &tree, &prefix, &config, 1, 123)?.pipeline;
            let count = ((args.target_ms as f64 / 1000.0 / estimate).ceil() as usize).clamp(1, 128);
            eprintln!(
                "{tips} tips, {length} sites, {model}: {} paired rounds of {count} attempts",
                args.rounds
            );
            for round in 0..args.rounds {
                // Alternate order to reduce bias from warm-up and thermal drift.
                let first_raw = round % 2 == 1;
                let seed = 1000 + round as i32 * 1000;
                let first = measure(first_raw, &tree, &prefix, &config, count, seed)?;
                let second = measure(!first_raw, &tree, &prefix, &config, count, seed)?;
                assert_eq!(
                    first.flags, second.flags,
                    "Duplicate decisions differ between transports"
                );
                for (raw, timing) in [(first_raw, first), (!first_raw, second)] {
                    let scale = 1000.0 / count as f64;
                    writeln!(
                        csv,
                        "{model},{tips},{length},{},{},{round},{count},{:.9},{:.9},{:.9},{:.9},{}",
                        args.branch_mean,
                        if raw { "numeric" } else { "fasta" },
                        timing.simulation * scale,
                        timing.check * scale,
                        timing.serialization * scale,
                        timing.pipeline * scale,
                        timing.duplicates
                    )?;
                }
            }
            csv.flush()?;
        }
    }
    fs::remove_dir_all(directory)?;
    println!("Results: {}", output.display());
    Ok(())
}
