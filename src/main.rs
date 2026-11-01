use std::error::Error;
use std::fs;
use std::io::{BufRead, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::{Duration, Instant, SystemTime};

use clap::{Args, Parser, ValueEnum};
use indicatif::{ProgressBar, ProgressStyle};
use pf_simulator::alignment::{AlignmentConfig, SequenceType, simulate_alignment_from_tree};
use pf_simulator::error::SimulationError;
use pf_simulator::tree::sample_seed;
use pf_simulator::tree_simulation::{TreeArgs, TreeSimulation};

mod run_summary;
mod workers;
use workers::{Counts, Outcome};

/// Generate phylogenetic trees and MSAs.
#[derive(Debug, Parser)]
#[command(
    name = "pf-sim",
    version,
    long_about = "Generate phylogenetic trees with phylotree and tip MSAs with linked AliSim.\n\nThe default strategy uses uniform unrooted topologies and exponential branches with mean 0.1. Choose --preset python for the Python birth-death pipeline. Alignments default to 500 amino acids under LG, with duplicate tip sequences rejected.",
    after_help = "Use --help for all options and preset details.",
    after_long_help = "Examples:\n  pf-sim --out data --samples 100 --tips 50\n  pf-sim --out python_data --preset python --samples 100 --tips 50\n  pf-sim --out dna_data --seqtype DNA --model JC --length 1000\n  pf-sim --out gamma_data --branch-prior gamma --gamma-shape 2 --gamma-scale 0.05"
)]
struct Config {
    /// Output directory (must not already exist).
    #[arg(long = "out", default_value = "simulated_data", value_name = "DIR")]
    output: PathBuf,

    /// Number of simulations to attempt.
    ///
    /// Samples that exhaust their duplicate retries are skipped and counted;
    /// the run may write fewer tree/alignment pairs than requested.
    #[arg(long, default_value_t = 1, value_parser = clap::value_parser!(u64).range(1..))]
    samples: u64,

    /// Number of tips per tree.
    #[arg(long, default_value_t = 10, value_parser = clap::value_parser!(u32).range(3..=10_000))]
    tips: u32,

    /// Master seed (random when omitted).
    ///
    /// The resolved seed and a replay command are saved to simulation_summary.json.
    #[arg(long)]
    seed: Option<u64>,

    /// Number of worker processes (1 runs directly in this process).
    ///
    /// Workers run independent AliSim instances and write separate sample files.
    /// Sample IDs and seeds are identical for every job count. Fatal errors stop
    /// the batch; already completed samples remain. Requested jobs are capped
    /// at the number of samples. Log fragments are merged after workers stop.
    #[arg(short = 'j', long, default_value_t = 1, value_parser = clap::value_parser!(u16).range(1..))]
    jobs: u16,

    /// Internal process-worker protocol over stdin/stdout.
    #[arg(long, hide = true)]
    worker: bool,

    #[command(flatten)]
    tree: TreeArgs,

    /// Generate tree files without simulating alignments.
    #[arg(long, hide_short_help = true, help_heading = "Output and retries")]
    trees_only: bool,

    /// Save AliSim's informational output to a new file.
    #[arg(
        long,
        value_name = "FILE",
        hide_short_help = true,
        help_heading = "Output and retries"
    )]
    alisim_log: Option<PathBuf>,

    /// Maximum MSA attempts per fixed tree before skipping it (includes the first attempt).
    #[arg(long, default_value_t = 20, value_parser = clap::value_parser!(u32).range(1..=2_147_483_646), hide_short_help = true, help_heading = "Output and retries")]
    max_attempts: u32,

    /// Accept identical aligned tip sequences without retrying.
    #[arg(long, hide_short_help = true, help_heading = "Output and retries")]
    allow_duplicate_sequences: bool,

    #[command(flatten)]
    alignment: AlignmentArgs,
}

#[derive(Debug, Args)]
struct AlignmentArgs {
    /// Root sequence length in sites.
    ///
    /// Indels can change the output alignment length.
    #[arg(long, default_value_t = 500, value_parser = clap::value_parser!(u32).range(1..=i32::MAX as i64))]
    length: u32,

    /// AliSim substitution model.
    ///
    /// Rate variation may be included, for example 'LG+G4{0.5}'.
    #[arg(long, default_value = "LG", value_parser = clap::builder::NonEmptyStringValueParser::new())]
    model: String,

    /// Sequence alphabet: AA or DNA.
    #[arg(
        long = "seqtype",
        value_enum,
        default_value = "AA",
        hide_possible_values = true,
        long_help = "Sequence alphabet: AA or DNA. Must agree with the substitution model."
    )]
    sequence_type: SequenceTypeArg,

    /// AliSim insertion,deletion rates.
    ///
    /// Rates are relative to substitutions, for example 0.002,0.002.
    #[arg(long, value_name = "INS,DEL", value_parser = clap::builder::NonEmptyStringValueParser::new())]
    indel: Option<String>,

    /// AliSim indel size distributions (e.g. 'POW{1.7/100},POW{1.7/100}').
    #[arg(long, value_name = "INS_DIST,DEL_DIST", requires = "indel", value_parser = clap::builder::NonEmptyStringValueParser::new(), hide_short_help = true, help_heading = "Indels")]
    indel_size: Option<String>,

    /// Disable AliSim's additional gap-free FASTA output when using indels.
    #[arg(long, hide_short_help = true, help_heading = "Indels")]
    no_unaligned: bool,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum SequenceTypeArg {
    #[value(name = "AA")]
    AminoAcid,
    #[value(name = "DNA")]
    Dna,
}

impl From<&AlignmentArgs> for AlignmentConfig {
    fn from(arguments: &AlignmentArgs) -> Self {
        Self {
            length: arguments.length,
            model: arguments.model.clone(),
            sequence_type: match arguments.sequence_type {
                SequenceTypeArg::AminoAcid => SequenceType::AminoAcid,
                SequenceTypeArg::Dna => SequenceType::Dna,
            },
            indel: arguments.indel.clone(),
            indel_size: arguments.indel_size.clone(),
            no_unaligned: arguments.no_unaligned,
        }
    }
}

fn log_attempt(path: Option<&Path>, message: &str) -> std::io::Result<()> {
    if let Some(path) = path {
        writeln!(fs::OpenOptions::new().append(true).open(path)?, "{message}")?;
    }
    Ok(())
}

/// One fixed-tree simulation. Only accepted outputs reach the filesystem.
fn simulate_sample(
    config: &Config,
    tree_simulation: &TreeSimulation,
    alignment: &AlignmentConfig,
    sample: u64,
    mut on_attempt: impl FnMut(u32, u32),
) -> Result<(), SimulationError> {
    let seed = sample_seed(
        config.seed.expect("master seed resolved before simulation"),
        sample,
    );
    let tree = tree_simulation.generate(config.tips as usize, seed)?;
    let name = format!("sample_{sample:06}");
    let tree_file = config.output.join("trees").join(format!("{name}.nwk"));
    if config.trees_only {
        fs::write(tree_file, tree.newick())?;
        return Ok(());
    }
    let sequences = config.output.join("sequences");
    let attempts = if config.allow_duplicate_sequences {
        1
    } else {
        config.max_attempts
    };
    let initial_seed = sample_seed(seed, 0) % 2_147_483_646;
    for attempt in 0..attempts {
        let alignment_seed = (1 + (initial_seed + u64::from(attempt)) % 2_147_483_646) as i32;
        on_attempt(attempt + 1, attempts);
        log_attempt(
            config.alisim_log.as_deref(),
            &format!(
                "Sample {name}, attempt {}/{attempts}, seed {alignment_seed}",
                attempt + 1
            ),
        )?;
        let data = simulate_alignment_from_tree(
            &tree,
            &sequences.join(&name),
            alignment,
            alignment_seed,
            config.alisim_log.as_deref(),
        )?;
        if !config.allow_duplicate_sequences && data.has_duplicates() {
            log_attempt(
                config.alisim_log.as_deref(),
                "Rejected alignment: duplicate tip sequences",
            )?;
            continue;
        }
        log_attempt(config.alisim_log.as_deref(), "Accepted alignment")?;
        let data = data.into_fasta();
        fs::write(tree_file, tree.newick())?;
        fs::write(sequences.join(format!("{name}.fa")), data.aligned_fasta)?;
        if let Some(unaligned) = data.unaligned_fasta {
            let directory = config.output.join("unaligned_sequences");
            fs::create_dir_all(&directory)?;
            fs::write(directory.join(format!("{name}.fa")), unaligned)?;
        }
        return Ok(());
    }
    Err(SimulationError::DuplicateSequences {
        sample: name,
        attempts,
    })
}

fn process_sample(
    config: &Config,
    trees: &TreeSimulation,
    alignment: &AlignmentConfig,
    sample: u64,
    on_attempt: impl FnMut(u32, u32),
) -> Result<Outcome, SimulationError> {
    match simulate_sample(config, trees, alignment, sample, on_attempt) {
        Ok(()) => Ok(Outcome::Written),
        Err(error @ SimulationError::DuplicateSequences { .. }) => {
            log_attempt(
                config.alisim_log.as_deref(),
                &format!("Skipped sample: {error}"),
            )?;
            Ok(Outcome::Skipped)
        }
        Err(error) => Err(error),
    }
}

fn run_worker(
    config: &Config,
    trees: &TreeSimulation,
    alignment: &AlignmentConfig,
) -> Result<(), Box<dyn Error>> {
    if !config.output.join("trees").is_dir() || !config.output.join("sequences").is_dir() {
        return Err("Worker output directories must be created by the coordinator".into());
    }
    let stdin = std::io::stdin();
    let stdout = std::io::stdout();
    let mut output = stdout.lock();
    for line in stdin.lock().lines() {
        let sample: u64 = line?.parse().map_err(|_| "Invalid worker sample ID")?;
        if sample >= config.samples {
            return Err("Worker sample ID is outside this batch".into());
        }
        let outcome = process_sample(config, trees, alignment, sample, |_, _| {})?;
        let status = match outcome {
            Outcome::Written => "written",
            Outcome::Skipped => "skipped",
        };
        writeln!(output, "{sample}\t{status}")?;
        output.flush()?;
    }
    Ok(())
}

fn run_serial(
    config: &Config,
    trees: &TreeSimulation,
    alignment: &AlignmentConfig,
    progress: &ProgressBar,
    counts: &mut Counts,
) -> Result<(), Box<dyn Error>> {
    for sample in 0..config.samples {
        let outcome = process_sample(config, trees, alignment, sample, |attempt, attempts| {
            progress.set_message(format!(
                "sample_{sample:06}, attempt {attempt}/{attempts} ({} failed)",
                counts.failed
            ));
        })?;
        counts.record(outcome);
        progress.set_message(format!(
            "{} written, {} failed",
            counts.written, counts.failed
        ));
        progress.inc(1);
    }
    Ok(())
}

fn run() -> Result<(), Box<dyn Error>> {
    let started = Instant::now();
    let started_at = SystemTime::now();
    let mut config = Config::parse();
    let provided_seed = config.seed.is_some();
    if config.worker && !provided_seed {
        return Err("Worker requires an explicit seed from its coordinator".into());
    }
    config.seed.get_or_insert_with(rand::random);
    let tree_simulation = TreeSimulation::new(&config.tree)?;
    let alignment = AlignmentConfig::from(&config.alignment);
    if config.worker {
        return run_worker(&config, &tree_simulation, &alignment);
    }
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
    let mut summary = run_summary::RunSummary::new(&config, provided_seed, started, started_at)?;
    println!(
        "Seed: {} (saved in {})",
        config.seed.unwrap(),
        config.output.join(run_summary::FILE_NAME).display()
    );
    let mut counts = Counts::default();
    let result = (|| -> Result<(), Box<dyn Error>> {
        let trees = config.output.join("trees");
        let sequences = config.output.join("sequences");
        fs::create_dir(&trees)?;
        fs::create_dir(&sequences)?;
        if let Some(log) = &config.alisim_log {
            let parent = log
                .parent()
                .filter(|path| !path.as_os_str().is_empty())
                .unwrap_or_else(|| Path::new("."));
            let parent = parent.canonicalize()?;
            if parent == config.output.canonicalize()?
                && log.file_name().is_some_and(|name| {
                    name == run_summary::FILE_NAME || name == run_summary::TEMP_FILE_NAME
                })
            {
                return Err("AliSim log must not use a simulation summary filename".into());
            }
            if parent == trees.canonicalize()? || parent == sequences.canonicalize()? {
                return Err(
                    "AliSim log must be outside the trees and sequences directories".into(),
                );
            }
            fs::OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(log)?;
        }
        let progress = ProgressBar::new(config.samples).with_style(ProgressStyle::with_template(
            "{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} samples (ETA {eta}) {msg}",
        )?);
        progress.enable_steady_tick(Duration::from_millis(100));
        let result = if config.jobs > 1 && config.samples > 1 {
            workers::run(&config, &progress, &mut counts)
        } else {
            run_serial(
                &config,
                &tree_simulation,
                &alignment,
                &progress,
                &mut counts,
            )
        };
        if result.is_ok() {
            progress.finish_and_clear();
        } else {
            progress.abandon_with_message("Simulation failed");
        }
        result
    })();
    let recorded = summary.finish(&counts, result.as_ref().err().map(|error| error.as_ref()));
    if let Err(error) = result {
        if let Err(recording) = recorded {
            return Err(
                format!("{error}; also failed to update simulation summary: {recording}").into(),
            );
        }
        return Err(error);
    }
    recorded?;
    println!(
        "Wrote {} {} to {} ({} failed due to duplicate sequences; {} requested)",
        counts.written,
        if config.trees_only {
            "trees"
        } else {
            "tree/FASTA pairs"
        },
        config.output.display(),
        counts.failed,
        config.samples,
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
