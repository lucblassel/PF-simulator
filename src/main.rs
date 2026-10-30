use std::error::Error;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::Duration;

use clap::{Args, Parser, ValueEnum};
use indicatif::{ProgressBar, ProgressStyle};
use pf_simulator::alignment::{AlignmentConfig, SequenceType, simulate_alignment_numeric};
use pf_simulator::tree::{Tree, sample_seed};

/// Generate uniform unrooted trees and tip MSAs using linked AliSim.
#[derive(Debug, Parser)]
#[command(name = "pf-sim", version)]
struct Config {
    /// Output directory (must not already exist).
    #[arg(long = "out", default_value = "simulated_data", value_name = "DIR")]
    output: PathBuf,

    /// Number of tree/alignment pairs.
    #[arg(long, default_value_t = 1, value_parser = clap::value_parser!(u64).range(1..))]
    samples: u64,

    /// Number of tips per tree.
    #[arg(long, default_value_t = 10, value_parser = clap::value_parser!(u32).range(3..=10_000))]
    tips: u32,

    /// Master random seed.
    #[arg(long, default_value_t = 1)]
    seed: u64,

    /// Mean exponential branch length, in substitutions per site.
    #[arg(long, default_value_t = 0.1, value_parser = parse_branch_mean)]
    branch_mean: f64,

    /// Save AliSim's informational output to a new file.
    #[arg(long, value_name = "FILE")]
    alisim_log: Option<PathBuf>,

    /// Maximum MSA attempts per fixed tree before failing (includes the first attempt).
    #[arg(long, default_value_t = 20, value_parser = clap::value_parser!(u32).range(1..=2_147_483_646))]
    max_attempts: u32,

    /// Accept identical aligned tip sequences without retrying.
    #[arg(long)]
    allow_duplicate_sequences: bool,

    #[command(flatten)]
    alignment: AlignmentArgs,
}

#[derive(Debug, Args)]
struct AlignmentArgs {
    /// Root sequence length in sites (indels can change the output alignment length).
    #[arg(long, default_value_t = 500, value_parser = clap::value_parser!(u32).range(1..=i32::MAX as i64))]
    length: u32,

    /// AliSim substitution model, optionally including rate variation.
    #[arg(long, default_value = "LG", value_parser = clap::builder::NonEmptyStringValueParser::new())]
    model: String,

    /// Sequence alphabet (must agree with the model).
    #[arg(long = "seqtype", value_enum, default_value = "AA")]
    sequence_type: SequenceTypeArg,

    /// AliSim insertion,deletion rates relative to substitutions (e.g. 0.002,0.002).
    #[arg(long, value_name = "INS,DEL", value_parser = clap::builder::NonEmptyStringValueParser::new())]
    indel: Option<String>,

    /// AliSim indel size distributions (e.g. 'POW{1.7/100},POW{1.7/100}').
    #[arg(long, value_name = "INS_DIST,DEL_DIST", requires = "indel", value_parser = clap::builder::NonEmptyStringValueParser::new())]
    indel_size: Option<String>,

    /// Disable AliSim's additional gap-free FASTA output when using indels.
    #[arg(long)]
    no_unaligned: bool,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum SequenceTypeArg {
    #[value(name = "AA")]
    AminoAcid,
    #[value(name = "DNA")]
    Dna,
}

impl From<AlignmentArgs> for AlignmentConfig {
    fn from(arguments: AlignmentArgs) -> Self {
        Self {
            length: arguments.length,
            model: arguments.model,
            sequence_type: match arguments.sequence_type {
                SequenceTypeArg::AminoAcid => SequenceType::AminoAcid,
                SequenceTypeArg::Dna => SequenceType::Dna,
            },
            indel: arguments.indel,
            indel_size: arguments.indel_size,
            no_unaligned: arguments.no_unaligned,
        }
    }
}

fn parse_branch_mean(value: &str) -> Result<f64, String> {
    let mean: f64 = value
        .parse()
        .map_err(|error| format!("Invalid branch mean: {error}"))?;
    if !mean.is_finite() || !(1e-8..=100.0).contains(&mean) {
        return Err("Branch mean must be finite and between 1e-8 and 100".into());
    }
    Ok(mean)
}

fn log_attempt(path: Option<&Path>, message: &str) -> std::io::Result<()> {
    if let Some(path) = path {
        writeln!(fs::OpenOptions::new().append(true).open(path)?, "{message}")?;
    }
    Ok(())
}

fn run() -> Result<(), Box<dyn Error>> {
    let config = Config::parse();
    if fs::symlink_metadata(&config.output).is_ok() {
        return Err(format!("Output already exists: {}", config.output.display()).into());
    }
    if let Some(parent) = config
        .output
        .parent()
        .filter(|path| !path.as_os_str().is_empty())
    {
        fs::create_dir_all(parent)?;
    }
    fs::create_dir(&config.output)?;
    let trees = config.output.join("trees");
    let sequences = config.output.join("sequences");
    let unaligned_sequences = config.output.join("unaligned_sequences");
    fs::create_dir(&trees)?;
    fs::create_dir(&sequences)?;

    if let Some(log) = &config.alisim_log {
        let parent = log
            .parent()
            .filter(|path| !path.as_os_str().is_empty())
            .unwrap_or_else(|| std::path::Path::new("."));
        let parent = parent.canonicalize()?;
        if parent == trees.canonicalize()? || parent == sequences.canonicalize()? {
            return Err("AliSim log must be outside the trees and sequences directories".into());
        }
        fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(log)?;
    }

    let alignment = AlignmentConfig::from(config.alignment);
    let progress = ProgressBar::new(config.samples).with_style(
        ProgressStyle::with_template(
            "{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} pairs (ETA {eta}) {msg}",
        )?,
    );
    progress.enable_steady_tick(Duration::from_millis(100));

    for sample in 0..config.samples {
        let seed = sample_seed(config.seed, sample);
        let tree = Tree::uniform(config.tips as usize, config.branch_mean, seed)?;
        let name = format!("sample_{sample:06}");
        let tree_file = trees.join(format!("{name}.nwk"));
        fs::write(&tree_file, tree.newick())?;
        let data = (|| -> Result<_, Box<dyn Error>> {
            let attempts = if config.allow_duplicate_sequences {
                1
            } else {
                config.max_attempts
            };
            let initial_seed = sample_seed(seed, 0) % 2_147_483_646;
            for attempt in 0..attempts {
                // Preserve the original first-attempt seed; later attempts use
                // distinct seeds without changing this tree or other samples.
                let alignment_seed =
                    (1 + (initial_seed + u64::from(attempt)) % 2_147_483_646) as i32;
                progress.set_message(format!("{name}, attempt {}/{}", attempt + 1, attempts));
                log_attempt(
                    config.alisim_log.as_deref(),
                    &format!(
                        "Sample {name}, attempt {}/{attempts}, seed {alignment_seed}",
                        attempt + 1
                    ),
                )?;
                let data = simulate_alignment_numeric(
                    &tree_file,
                    &sequences.join(&name),
                    &alignment,
                    alignment_seed,
                    config.alisim_log.as_deref(),
                )?;
                if config.allow_duplicate_sequences || !data.has_duplicates() {
                    log_attempt(config.alisim_log.as_deref(), "Accepted alignment")?;
                    return Ok(data);
                }
                log_attempt(
                    config.alisim_log.as_deref(),
                    "Rejected alignment: duplicate tip sequences",
                )?;
            }
            Err(format!(
                "No duplicate-free MSA for {name} after {attempts} attempts; \
                increase --max-attempts or use --allow-duplicate-sequences"
            )
            .into())
        })();
        let data = match data {
            Ok(data) => data,
            Err(error) => {
                progress.abandon_with_message("Simulation failed");
                fs::remove_file(&tree_file)?;
                return Err(error);
            }
        };

        // Format only accepted numeric rows.
        let data = data.into_fasta();
        // Only accepted alignment bytes reach the filesystem.
        fs::write(sequences.join(format!("{name}.fa")), data.aligned_fasta)?;
        if let Some(unaligned) = data.unaligned_fasta {
            fs::create_dir_all(&unaligned_sequences)?;
            fs::write(unaligned_sequences.join(format!("{name}.fa")), unaligned)?;
        }
        progress.inc(1);
    }
    progress.finish_and_clear();
    println!(
        "Wrote {} tree/FASTA pairs to {}",
        config.samples,
        config.output.display()
    );
    Ok(())
}

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("pf-sim: {error}");
            ExitCode::FAILURE
        }
    }
}
