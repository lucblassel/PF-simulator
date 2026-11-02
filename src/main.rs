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

use pf_simulator::selection::{self, Rejection, TipChecks, WindowPolicy};
use pf_simulator::site_rates::GammaMode;

mod packed;
mod run_summary;
mod workers;
use workers::{Counts, Outcome};

/// Generate phylogenetic trees and MSAs.
#[derive(Debug, Parser)]
#[command(
    name = "pf-sim",
    version,
    long_about = "Generate phylogenetic trees with phylotree and tip MSAs with linked AliSim.\n\nThe default strategy uses uniform unrooted topologies and exponential branches with mean 0.1. Choose --preset python for Python birth-death trees and LG with empirical eight-category Gamma (G8) site rates. Alignments default to 500 amino acids under LG, with duplicate tip sequences rejected.",
    after_help = "Use --help for all options and preset details.",
    after_long_help = "Examples:\n  pf-sim --out data --samples 100 --tips 50\n  pf-sim --out python_data --preset python --samples 100 --tips 50\n  pf-sim --out dna_data --seqtype DNA --model JC --length 1000\n  pf-sim --out gamma_data --branch-prior gamma --gamma-shape 2 --gamma-scale 0.05"
)]
struct Config {
    /// Output directory (must not already exist).
    #[arg(long = "out", default_value = "simulated_data", value_name = "DIR")]
    output: PathBuf,

    /// Store examples as individual files or worker-owned tar shards plus SQLite.
    #[arg(long, value_enum, default_value_t = OutputFormat::Files, hide_short_help = true, help_heading = "Packed output")]
    output_format: OutputFormat,

    /// Target tar size in bytes (1 GiB); pairs stay together and may exceed it.
    #[arg(long, default_value_t = 1_073_741_824, value_parser = clap::value_parser!(u64).range(1024..), hide_short_help = true, help_heading = "Packed output")]
    shard_size_bytes: u64,

    /// Optional identity for merging independent batches (does not change the seed).
    #[arg(long, hide_short_help = true, help_heading = "Packed output")]
    batch_id: Option<String>,

    /// Prefix sample and shard filenames with PREFIX- (index/summary names stay fixed).
    ///
    /// Use a different prefix per cluster batch for collision-free tar extraction.
    /// Accepts 1–40 ASCII letters, digits, dots, underscores or hyphens.
    #[arg(
        long,
        value_name = "PREFIX",
        hide_short_help = true,
        help_heading = "Output and retries"
    )]
    file_prefix: Option<String>,

    #[arg(long, hide = true)]
    worker_id: Option<u16>,

    /// Number of simulations to attempt.
    ///
    /// Samples that exhaust their duplicate retries are skipped and counted;
    /// the run may write fewer tree/alignment pairs than requested unless
    /// --ensure-samples is enabled.
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
    /// Workers run independent AliSim instances and own their output files/shards.
    /// Sample IDs and seeds are identical for every job count; packed shard
    /// assignment may vary. Fatal errors stop the batch; already completed
    /// samples remain. Requested jobs are capped
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

    /// Generate fresh candidates until --samples outputs have been written.
    ///
    /// WARNING: this can significantly increase runtime and may never finish
    /// if acceptable outputs cannot be generated. Each fixed tree still
    /// uses at most --max-attempts MSA attempts before being replaced.
    #[arg(long, hide_short_help = true, help_heading = "Output and retries")]
    ensure_samples: bool,

    /// Bypass all three tip checks and accept the first attempt.
    ///
    /// Window metadata and optional trimmed FASTA use the first L columns.
    #[arg(long, hide_short_help = true, help_heading = "Output and retries")]
    allow_duplicate_sequences: bool,

    #[command(flatten)]
    alignment: AlignmentArgs,

    #[command(flatten)]
    selection: SelectionArgs,

    /// Additional FASTAs to write: trimmed, unaligned, or all (comma-separated).
    ///
    /// Full MSAs and per-sample metadata are always saved. Unaligned exports
    /// remove gaps from full rows; trimmed exports use alignment columns.
    #[arg(
        long,
        value_enum,
        value_delimiter = ',',
        hide_short_help = true,
        help_heading = "Output and retries"
    )]
    materialize: Vec<Materialize>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
enum OutputFormat {
    Files,
    Packed,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
enum Materialize {
    Trimmed,
    Unaligned,
    All,
}

#[derive(Debug, Args)]
struct SelectionArgs {
    /// Require distinct full gapped tip rows (set =false to disable).
    #[arg(long, num_args = 0..=1, require_equals = true, default_missing_value = "true", default_value_t = true, hide_short_help = true, help_heading = "Tip selection")]
    require_distinct_gapped: bool,
    /// Require distinct full ungapped tip rows (set =false to disable).
    #[arg(long, num_args = 0..=1, require_equals = true, default_missing_value = "true", default_value_t = false, hide_short_help = true, help_heading = "Tip selection")]
    require_distinct_ungapped: bool,
    /// Require a length-L window with distinct aligned tip rows.
    #[arg(long, num_args = 0..=1, require_equals = true, default_missing_value = "true", default_value_t = false, hide_short_help = true, help_heading = "Tip selection")]
    require_unique_window: bool,
    /// Sliding-window starts: inclusive (all starts), python (exclude last for M > L).
    ///
    /// Defaults to python for --preset python, otherwise inclusive. For M == L
    /// both test start zero; M < L fails the window check. Only tips are compared.
    #[arg(
        long,
        value_enum,
        hide_short_help = true,
        help_heading = "Tip selection"
    )]
    window_policy: Option<WindowPolicy>,
}

impl Config {
    fn file_prefix(&self) -> String {
        self.file_prefix
            .as_ref()
            .map(|prefix| format!("{prefix}-"))
            .unwrap_or_default()
    }

    fn materializes(&self, kind: Materialize) -> bool {
        self.materialize.contains(&Materialize::All) || self.materialize.contains(&kind)
    }
    fn checks(&self) -> TipChecks {
        let enabled = !self.allow_duplicate_sequences;
        TipChecks {
            gapped: enabled && self.selection.require_distinct_gapped,
            ungapped: enabled && self.selection.require_distinct_ungapped,
            window: enabled && self.selection.require_unique_window,
            window_policy: self
                .selection
                .window_policy
                .expect("resolved window policy"),
        }
    }
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
    /// Defaults to LG. Explicit --model overrides the preset's empirical Gamma
    /// default. Use --gamma G8 or GC to add empirical rates to an explicit model.
    /// Fixed rates can be included directly, for example 'LG+G4{0.5}'.
    #[arg(long, value_parser = clap::builder::NonEmptyStringValueParser::new())]
    model: Option<String>,

    /// Empirical site-rate Gamma: none, GC, or G[n] (e.g. G8).
    ///
    /// Each attempt draws Normal(a, a/10), floored at 0.05, with a uniformly
    /// sampled from bundled HOGENOM observations. Default: G8 for the Python
    /// preset without an explicit --model, otherwise none. This flag overrides
    /// that default; do not combine it with an existing +G/+R rate modifier.
    ///
    /// WARNING: GC (continuous Gamma) can be much slower than G8, especially
    /// for large trees and alignments. Use --gamma GC explicitly to reproduce
    /// the Python pipeline's continuous site-rate prior.
    #[arg(
        long,
        value_name = "MODE",
        hide_short_help = true,
        help_heading = "Site rates"
    )]
    gamma: Option<GammaMode>,

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

    /// Legacy opt-out for gap-free output (already off unless materialized).
    ///
    /// Conflicts with --materialize unaligned or all.
    #[arg(long, hide_short_help = true, help_heading = "Indels")]
    no_unaligned: bool,

    /// Include simulated ancestral rows, labelled to match internal tree nodes.
    #[arg(long, hide_short_help = true, help_heading = "Output and retries")]
    ancestral: bool,
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
            model: arguments.model.clone().unwrap_or_else(|| "LG".into()),
            sequence_type: match arguments.sequence_type {
                SequenceTypeArg::AminoAcid => SequenceType::AminoAcid,
                SequenceTypeArg::Dna => SequenceType::Dna,
            },
            indel: arguments.indel.clone(),
            indel_size: arguments.indel_size.clone(),
            no_unaligned: arguments.no_unaligned,
            ancestral: arguments.ancestral,
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
    packed: &mut Option<packed::PackedWriter>,
    mut on_attempt: impl FnMut(u32, u32),
) -> Result<(), SimulationError> {
    let seed = sample_seed(
        config.seed.expect("master seed resolved before simulation"),
        sample,
    );
    let tree = tree_simulation.generate(config.tips as usize, seed)?;
    let name = format!("{}sample_{sample:06}", config.file_prefix());
    let tree_file = config.output.join("trees").join(format!("{name}.nwk"));
    if config.trees_only {
        fs::write(
            tree_file,
            if alignment.ancestral {
                tree.newick_with_internal_names()?
            } else {
                tree.newick()
            },
        )?;
        return Ok(());
    }
    let sequences = config.output.join("msa_untrimmed");
    let attempts = if config.allow_duplicate_sequences {
        1
    } else {
        config.max_attempts
    };
    let initial_seed = sample_seed(seed, 0) % 2_147_483_646;
    let mut last_rejection = Rejection::Gapped;
    let mut rejections = [0_u32; 3];
    for attempt in 0..attempts {
        let alignment_seed = (1 + (initial_seed + u64::from(attempt)) % 2_147_483_646) as i32;
        let mut attempt_alignment = alignment.clone();
        attempt_alignment.model = config.alignment.gamma.expect("resolved rate mode").model(
            &alignment.model,
            config.seed.unwrap(),
            sample,
            attempt,
        );
        on_attempt(attempt + 1, attempts);
        log_attempt(
            config.alisim_log.as_deref(),
            &format!(
                "Sample {name}, attempt {}/{attempts}, seed {alignment_seed}, model {}",
                attempt + 1,
                attempt_alignment.model
            ),
        )?;
        let data = simulate_alignment_from_tree(
            &tree,
            &sequences.join(&name),
            &attempt_alignment,
            alignment_seed,
            config.alisim_log.as_deref(),
        )?;
        let window = match selection::select(
            &data,
            tree.tips(),
            alignment.length as usize,
            config.checks(),
        )? {
            Ok(window) => window,
            Err(reason) => {
                last_rejection = reason;
                rejections[match reason {
                    Rejection::Gapped => 0,
                    Rejection::Ungapped => 1,
                    Rejection::Window => 2,
                }] += 1;
                log_attempt(
                    config.alisim_log.as_deref(),
                    &format!(
                        "Rejected alignment: {}",
                        match reason {
                            Rejection::Gapped => "duplicate tip sequences (gapped)",
                            Rejection::Ungapped => "duplicate tip sequences (ungapped)",
                            Rejection::Window => "no unique length-L tip window",
                        }
                    ),
                )?;
                continue;
            }
        };
        log_attempt(config.alisim_log.as_deref(), "Accepted alignment")?;
        let checks = config.checks();
        let metadata = serde_json::json!({
            "schema_version": 1, "sample": name, "candidate_id": sample, "batch_id": config.batch_id,
            "master_seed": config.seed.unwrap(), "tree_seed": seed,
            "attempt_index": attempt, "attempt_number": attempt + 1,
            "alisim_seed": alignment_seed, "model": attempt_alignment.model,
            "alignment_columns": data.rows[0].len(), "requested_length": alignment.length,
            "tip_names": (1..=tree.tips()).map(|i| format!("T{i}")).collect::<Vec<_>>(),
            "window": { "start": window.start, "end": window.end,
                "coordinates": "zero-based half-open full-alignment columns",
                "uniqueness_checked": window.checked,
                "selection": if window.checked { "first_unique_tip_window" } else { "prefix" },
                "endpoint_policy": format!("{:?}", checks.window_policy).to_lowercase() },
            "checks": { "distinct_gapped": checks.gapped, "distinct_ungapped": checks.ungapped, "unique_window": checks.window },
            "allow_duplicate_sequences": config.allow_duplicate_sequences,
            "prior_attempt_rejections": { "gapped": rejections[0], "ungapped": rejections[1], "window": rejections[2] },
            "unaligned_source": "full alignment, not selected window",
            "ancestral": alignment.ancestral,
            "ancestral_names": data.names.iter().filter(|name| !name.starts_with('T')).collect::<Vec<_>>(),
        });
        if let Some(writer) = packed {
            let newick = if alignment.ancestral {
                tree.newick_with_internal_names()?
            } else {
                tree.newick()
            };
            writer.write(
                sample,
                &name,
                newick.as_bytes(),
                &data.to_fasta(false),
                &metadata,
            )?;
        } else {
            fs::write(
                tree_file,
                if alignment.ancestral {
                    tree.newick_with_internal_names()?
                } else {
                    tree.newick()
                },
            )?;
            fs::write(sequences.join(format!("{name}.fa")), data.to_fasta(false))?;
            if config.materializes(Materialize::Unaligned) {
                fs::write(
                    config.output.join("unaligned").join(format!("{name}.fa")),
                    data.to_fasta(true),
                )?;
            }
            if config.materializes(Materialize::Trimmed) {
                fs::write(
                    config.output.join("msa_trimmed").join(format!("{name}.fa")),
                    data.window_fasta(window.start, window.end),
                )?;
            }
            fs::write(
                config.output.join("metadata").join(format!("{name}.json")),
                serde_json::to_vec_pretty(&metadata).map_err(std::io::Error::other)?,
            )?;
        }
        return Ok(());
    }
    Err(if last_rejection == Rejection::Gapped {
        SimulationError::DuplicateSequences {
            sample: name,
            attempts,
        }
    } else {
        SimulationError::SelectionFailed {
            sample: name,
            attempts,
            reason: last_rejection,
        }
    })
}

fn process_sample(
    config: &Config,
    trees: &TreeSimulation,
    alignment: &AlignmentConfig,
    sample: u64,
    packed: &mut Option<packed::PackedWriter>,
    on_attempt: impl FnMut(u32, u32),
) -> Result<Outcome, SimulationError> {
    match simulate_sample(config, trees, alignment, sample, packed, on_attempt) {
        Ok(()) => Ok(Outcome::Written),
        Err(error @ SimulationError::DuplicateSequences { .. }) => {
            log_attempt(
                config.alisim_log.as_deref(),
                &format!("Skipped sample: {error}"),
            )?;
            Ok(Outcome::Skipped(Rejection::Gapped))
        }
        Err(error @ SimulationError::SelectionFailed { .. }) => {
            let reason = match &error {
                SimulationError::SelectionFailed { reason, .. } => *reason,
                _ => unreachable!(),
            };
            log_attempt(
                config.alisim_log.as_deref(),
                &format!("Skipped sample: {error}"),
            )?;
            Ok(Outcome::Skipped(reason))
        }
        Err(error) => Err(error),
    }
}

fn run_worker(
    config: &Config,
    trees: &TreeSimulation,
    alignment: &AlignmentConfig,
) -> Result<(), Box<dyn Error>> {
    if config.output_format == OutputFormat::Files
        && (!config.output.join("trees").is_dir() || !config.output.join("msa_untrimmed").is_dir())
    {
        return Err("Worker output directories must be created by the coordinator".into());
    }
    let mut packed = new_packed_writer(config, config.worker_id.ok_or("Missing worker identity")?)?;
    let stdin = std::io::stdin();
    let stdout = std::io::stdout();
    let mut output = stdout.lock();
    for line in stdin.lock().lines() {
        let sample: u64 = line?.parse().map_err(|_| "Invalid worker sample ID")?;
        if !config.ensure_samples && sample >= config.samples {
            return Err("Worker sample ID is outside this batch".into());
        }
        let outcome = process_sample(config, trees, alignment, sample, &mut packed, |_, _| {})?;
        let status = match outcome {
            Outcome::Written => "written",
            Outcome::Skipped(reason) => reason.name(),
        };
        writeln!(output, "{sample}\t{status}")?;
        output.flush()?;
    }
    if let Some(writer) = packed {
        writer.finish()?;
    }
    Ok(())
}

fn new_packed_writer(config: &Config, id: u16) -> std::io::Result<Option<packed::PackedWriter>> {
    if config.output_format == OutputFormat::Packed {
        packed::PackedWriter::new(
            config.output.join("shards"),
            config.shard_size_bytes,
            id,
            config.file_prefix(),
        )
        .map(Some)
    } else {
        Ok(None)
    }
}

fn run_serial(
    config: &Config,
    trees: &TreeSimulation,
    alignment: &AlignmentConfig,
    progress: &ProgressBar,
    counts: &mut Counts,
) -> Result<(), Box<dyn Error>> {
    let mut packed = new_packed_writer(config, 0)?;
    let mut next = 0_u64;
    while !counts.finished(config) {
        let sample = next;
        next = next
            .checked_add(1)
            .ok_or("Candidate sample IDs exhausted")?;
        let outcome = process_sample(
            config,
            trees,
            alignment,
            sample,
            &mut packed,
            |attempt, attempts| {
                progress.set_message(format!(
                    "sample_{sample:06}, attempt {attempt}/{attempts} ({} failed)",
                    counts.failed
                ));
            },
        )?;
        counts.record(outcome);
        counts.update_progress(config, progress);
    }
    if let Some(writer) = packed {
        writer.finish()?;
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
    // Bound generated basenames so each tar member fits one ustar header,
    // including the largest worker/candidate IDs and the fixed directory name.
    if let Some(prefix) = &config.file_prefix
        && (prefix.is_empty()
            || prefix.len() > 40
            || !prefix
                .bytes()
                .all(|c| c.is_ascii_alphanumeric() || b"._-".contains(&c)))
    {
        return Err(
            "--file-prefix must contain 1–40 ASCII letters, digits, dots, underscores or hyphens"
                .into(),
        );
    }
    config.seed.get_or_insert_with(rand::random);
    let default_batch_id = format!("seed-{}", config.seed.unwrap());
    let batch_id = config.batch_id.get_or_insert(default_batch_id);
    if batch_id.is_empty()
        || batch_id.len() > 128
        || !batch_id
            .bytes()
            .all(|c| c.is_ascii_alphanumeric() || b"_-".contains(&c))
    {
        return Err(
            "--batch-id must contain 1–128 ASCII letters, digits, underscores or hyphens".into(),
        );
    }
    if config.output_format == OutputFormat::Packed
        && (config.trees_only || !config.materialize.is_empty())
    {
        return Err("Packed output stores one full MSA/tree pair; --trees-only and --materialize require --output-format files".into());
    }
    let python = config.tree.python_bd
        || config.tree.preset == pf_simulator::tree_simulation::Preset::Python;
    let gamma =
        *config
            .alignment
            .gamma
            .get_or_insert(if python && config.alignment.model.is_none() {
                GammaMode::Discrete(8)
            } else {
                GammaMode::None
            });
    if gamma != GammaMode::None
        && config
            .alignment
            .model
            .as_deref()
            .unwrap_or("LG")
            .split('+')
            .skip(1)
            .any(|part| part.starts_with('G') || part.starts_with('R'))
    {
        return Err("--gamma cannot be combined with an existing +G/+R model modifier; use --gamma none for a fixed model".into());
    }
    config.selection.window_policy.get_or_insert(if python {
        WindowPolicy::Python
    } else {
        WindowPolicy::Inclusive
    });
    if config.alignment.no_unaligned && config.materializes(Materialize::Unaligned) {
        return Err("--no-unaligned conflicts with materializing unaligned FASTA".into());
    }
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
    if config.ensure_samples {
        eprintln!(
            "Warning: --ensure-samples can significantly increase runtime and may never finish if acceptable outputs cannot be generated."
        );
    }
    let mut counts = Counts::default();
    let result = (|| -> Result<(), Box<dyn Error>> {
        let mut output_directories = Vec::new();
        if config.output_format == OutputFormat::Packed {
            let directory = config.output.join("shards");
            fs::create_dir(&directory)?;
            output_directories.push(directory);
        } else {
            let trees = config.output.join("trees");
            let sequences = config.output.join("msa_untrimmed");
            fs::create_dir(&trees)?;
            fs::create_dir(&sequences)?;
            let metadata = config.output.join("metadata");
            fs::create_dir(&metadata)?;
            output_directories.extend([trees.clone(), sequences.clone(), metadata]);
            for (kind, name) in [
                (Materialize::Trimmed, "msa_trimmed"),
                (Materialize::Unaligned, "unaligned"),
            ] {
                if config.materializes(kind) {
                    let directory = config.output.join(name);
                    fs::create_dir(&directory)?;
                    output_directories.push(directory);
                }
            }
        }
        if let Some(log) = &config.alisim_log {
            let parent = log
                .parent()
                .filter(|path| !path.as_os_str().is_empty())
                .unwrap_or_else(|| Path::new("."));
            let parent = parent.canonicalize()?;
            if parent == config.output.canonicalize()?
                && log.file_name().is_some_and(|name| {
                    name == run_summary::FILE_NAME
                        || name == run_summary::TEMP_FILE_NAME
                        || name == packed::INDEX
                        || name == packed::PARTIAL_INDEX
                })
            {
                return Err("AliSim log must not use a simulation summary filename or packed index filename".into());
            }
            if output_directories
                .iter()
                .any(|dir| dir.canonicalize().is_ok_and(|path| path == parent))
            {
                return Err("AliSim log must be outside the sample output directories".into());
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
        let result = result.and_then(|()| {
            if config.output_format == OutputFormat::Packed {
                progress.set_message("Finalizing packed index");
                packed::merge_workers(
                    &config.output,
                    u64::from(config.jobs).min(config.samples),
                    counts.written,
                )?;
            }
            Ok(())
        });
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
        "Wrote {} {} to {} ({} failed {}; {} requested)",
        counts.written,
        if config.trees_only {
            "trees"
        } else {
            "tree/FASTA pairs"
        },
        config.output.display(),
        counts.failed,
        if config.checks().window {
            "selection checks"
        } else {
            "due to duplicate sequences"
        },
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
