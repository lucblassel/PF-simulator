//! Machine-readable provenance and timings for one coordinator run.
use std::ffi::{OsStr, OsString};
use std::io;
use std::path::PathBuf;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

use clap::ValueEnum;
use pf_simulator::tree_simulation::{BranchPrior, Preset, TreePrior};
use serde_json::{Value, json};

use crate::workers::Counts;
use crate::{Config, SequenceTypeArg};

pub(crate) const FILE_NAME: &str = "simulation_summary.json";
pub(crate) const TEMP_FILE_NAME: &str = "simulation_summary.json.tmp";

pub(crate) struct RunSummary {
    path: PathBuf,
    started: Instant,
    document: Value,
}

/// Bash quoting also supports non-UTF-8 Unix paths without lossy conversion.
fn quote(argument: &OsStr) -> String {
    match argument.to_str() {
        Some(text)
            if !text.is_empty()
                && text
                    .bytes()
                    .all(|byte| byte.is_ascii_alphanumeric() || b"_@%+=:,./-".contains(&byte)) =>
        {
            text.to_owned()
        }
        Some(text) => format!("'{}'", text.replace('\'', "'\"'\"'")),
        None => {
            let bytes: String = argument
                .as_encoded_bytes()
                .iter()
                .map(|byte| format!("\\x{byte:02x}"))
                .collect();
            format!("$'{bytes}'")
        }
    }
}

fn value<T: ValueEnum>(item: T) -> String {
    item.to_possible_value()
        .expect("CLI enum variant")
        .get_name()
        .to_owned()
}

/// Spell out defaults too, so the saved command does not rely on CLI defaults.
fn resolved_arguments(config: &Config) -> io::Result<Vec<OsString>> {
    let tree = &config.tree;
    let python = tree.python_bd || tree.preset == Preset::Python;
    let topology = tree.tree_prior.unwrap_or(if python {
        TreePrior::BirthDeath
    } else {
        TreePrior::Uniform
    });
    let branches = tree
        .branch_prior
        .unwrap_or(if topology == TreePrior::BirthDeath {
            BranchPrior::ProcessTime
        } else {
            BranchPrior::Exponential
        });
    let mut arguments = Vec::new();
    let mut add = |name: &str, text: String| {
        arguments.push(OsString::from(format!("{name}={text}")));
    };
    for (name, text) in [
        ("--samples", config.samples.to_string()),
        ("--tips", config.tips.to_string()),
        ("--seed", config.seed.expect("resolved seed").to_string()),
        ("--jobs", config.jobs.to_string()),
        ("--preset", value(tree.preset)),
        ("--tree-prior", value(topology)),
        ("--branch-prior", value(branches)),
        ("--branch-mean", tree.branch_mean.to_string()),
        ("--gamma-shape", tree.gamma_shape.to_string()),
        ("--gamma-scale", tree.gamma_scale.to_string()),
        ("--birth-rate", tree.birth_rate.to_string()),
        ("--death-rate", tree.death_rate.to_string()),
        ("--small-event-mean", tree.small_event_mean.to_string()),
        ("--small-log-sd", tree.small_log_sd.to_string()),
        ("--large-event-mean", tree.large_event_mean.to_string()),
        ("--large-log-sd", tree.large_log_sd.to_string()),
        (
            "--minimum-terminal",
            tree.minimum_terminal
                .unwrap_or(if python || tree.python_rescaling {
                    0.001
                } else {
                    0.0
                })
                .to_string(),
        ),
        (
            "--require-distinct-gapped",
            config.selection.require_distinct_gapped.to_string(),
        ),
        (
            "--require-distinct-ungapped",
            config.selection.require_distinct_ungapped.to_string(),
        ),
        (
            "--require-unique-window",
            config.selection.require_unique_window.to_string(),
        ),
        (
            "--window-policy",
            value(
                config
                    .selection
                    .window_policy
                    .expect("resolved window policy"),
            ),
        ),
        ("--max-attempts", config.max_attempts.to_string()),
        ("--length", config.alignment.length.to_string()),
        (
            "--model",
            config
                .alignment
                .model
                .clone()
                .unwrap_or_else(|| "LG".into()),
        ),
        (
            "--gamma",
            config
                .alignment
                .gamma
                .expect("resolved rate mode")
                .to_string(),
        ),
        (
            "--seqtype",
            match config.alignment.sequence_type {
                SequenceTypeArg::AminoAcid => "AA",
                SequenceTypeArg::Dna => "DNA",
            }
            .into(),
        ),
    ] {
        add(name, text);
    }
    if !config.materialize.is_empty() {
        add(
            "--materialize",
            config
                .materialize
                .iter()
                .map(|&kind| value(kind))
                .collect::<Vec<_>>()
                .join(","),
        );
    }
    if let Some(indel) = &config.alignment.indel {
        add("--indel", indel.clone());
    }
    if let Some(size) = &config.alignment.indel_size {
        add("--indel-size", size.clone());
    }
    for (name, enabled) in [
        ("--trees-only", config.trees_only),
        ("--ensure-samples", config.ensure_samples),
        (
            "--allow-duplicate-sequences",
            config.allow_duplicate_sequences,
        ),
        ("--no-unaligned", config.alignment.no_unaligned),
        ("--ancestral", config.alignment.ancestral),
        ("--python-bd", tree.python_bd),
        ("--python-rescaling", tree.python_rescaling),
        ("--use-branch-times", tree.use_branch_times),
        ("--empirical-diameter", tree.empirical_diameter || python),
    ] {
        if enabled {
            arguments.push(name.into());
        }
    }
    if let Some(path) = &tree.empirical_diameters {
        let mut argument = OsString::from("--empirical-diameters=");
        argument.push(path.canonicalize()?);
        arguments.push(argument);
    }
    Ok(arguments)
}

impl RunSummary {
    pub fn new(
        config: &Config,
        provided_seed: bool,
        started: Instant,
        started_at: SystemTime,
    ) -> io::Result<Self> {
        let directory = std::env::current_dir()?;
        let executable = std::env::current_exe()?;
        let arguments = resolved_arguments(config)?;
        // A fresh parent directory avoids overwriting the original data or logs,
        // and makes the recorded command safe to run more than once.
        let mut command = format!(
            "replay_dir=$(mktemp -d) &&\ncd {} &&\n{} --out \"$replay_dir/data\"",
            quote(directory.as_os_str()),
            quote(executable.as_os_str())
        );
        for argument in &arguments {
            command.push(' ');
            command.push_str(&quote(argument));
        }
        if config.alisim_log.is_some() {
            command.push_str(" --alisim-log \"$replay_dir/alisim.log\"");
        }
        let gamma_provenance: std::collections::BTreeMap<_, _> =
            pf_simulator::site_rates::EMPIRICAL_GAMMA
                .lines()
                .filter_map(|line| line.strip_prefix("# ")?.split_once('='))
                .collect();
        let msa_policy = json!({
            "length": config.alignment.length,
            "sequence_type": value(config.alignment.sequence_type),
            "base_model": config.alignment.model.as_deref().unwrap_or("LG"),
            "gamma": config.alignment.gamma.unwrap().to_string(),
            "empirical_gamma": {
                "source": "hogenom_alphas.txt",
                "active": !config.trees_only && config.alignment.gamma.unwrap() != pf_simulator::site_rates::GammaMode::None,
                "source_sha256": gamma_provenance["source_sha256"],
                "values_sha256": gamma_provenance["values_sha256"],
                "count": gamma_provenance["count"].parse::<u64>().expect("bundled count"),
                "law": "Uniform observation with replacement; Normal(a, a/10); max(draw, 0.05); fresh draw per attempt",
                "rng": "v1: rand 0.8.5 StdRng; rand_distr 0.4.3 Normal; seed=sample_seed(sample_seed(master_seed,candidate),0x47414d4d41000000 XOR attempt)",
                "replay": "candidate is numeric sample filename ID; attempt is zero-based; sample_alpha(master_seed,candidate,attempt) reconstructs alpha; optional AliSim log records each model and engine seed"
            },
            "engine_seed": "1 + ((sample_seed(sample_seed(master_seed,candidate),0) % 2147483646 + attempt) % 2147483646); zero-based candidate and attempt",
            "window_selection": if config.checks().window { "first_unique_tip_window" } else { "prefix; uniqueness not checked" },
            "window_endpoint_policy": value(config.checks().window_policy),
            "duplicate_check": config.checks().gapped,
            "tip_checks": { "distinct_gapped": config.checks().gapped, "distinct_ungapped": config.checks().ungapped, "unique_window": config.checks().window },
            "allow_duplicate_sequences": config.allow_duplicate_sequences,
            "materialize": config.materialize.iter().map(|&kind| value(kind)).collect::<Vec<_>>(),
            "max_attempts": if config.allow_duplicate_sequences { 1 } else { config.max_attempts },
            "replace_failed_candidates": config.ensure_samples,
            "indel": config.alignment.indel,
            "indel_size": config.alignment.indel_size,
            "no_unaligned": config.alignment.no_unaligned,
            "ancestral": config.alignment.ancestral,
        });
        let document = json!({
            "schema_version": 1,
            "status": "running",
            "seed": config.seed.expect("resolved seed"),
            "seed_source": if provided_seed { "provided" } else { "generated" },
            "requested_samples": config.samples,
            "ensure_samples": config.ensure_samples,
            "sample_count_mode": if config.ensure_samples { "successful_outputs" } else { "attempted_candidates" },
            "completed_samples": null,
            "written_samples": null,
            "failed_samples": null,
            "failed_duplicate_samples": null,
            "failed_window_samples": null,
            "failed_samples_by_last_rejection": null,
            "output_layout": { "full_msa": "msa_untrimmed", "sample_metadata": "metadata", "trimmed": "msa_trimmed (when materialized)", "unaligned": "unaligned (when materialized)" },
            "failed_duplicate_percentage": null,
            "failed_duplicate_percentage_denominator": if config.ensure_samples { "completed_samples" } else { "requested_samples" },
            "wall_time_seconds": null,
            "wall_time_scope": "Coordinator start before CLI parsing through simulation, output I/O, worker shutdown and log merging; excludes final summary serialization and terminal completion output.",
            "started_unix_seconds": started_at.duration_since(UNIX_EPOCH).map_err(io::Error::other)?.as_secs_f64(),
            "working_directory": directory.to_string_lossy(),
            "executable": executable.to_string_lossy(),
            "version": env!("CARGO_PKG_VERSION"),
            "build_profile": if cfg!(debug_assertions) { "debug" } else { "release" },
            "jobs": config.jobs,
            "effective_jobs": u64::from(config.jobs).min(config.samples),
            "trees_only": config.trees_only,
            "resolved_arguments": arguments.iter().map(|arg| arg.to_string_lossy()).collect::<Vec<_>>(),
            "reproduce_command": command,
            "reproduction_requirements": "Use the same executable/build and unchanged external inputs (such as custom model or empirical files), if supplied. Sample files are reproducible; timings, summary metadata and parallel log ordering vary.",
            "msa_policy": msa_policy,
            "error": null,
        });
        let summary = Self {
            path: config.output.join(FILE_NAME),
            started,
            document,
        };
        summary.save()?;
        Ok(summary)
    }

    pub fn finish(
        &mut self,
        counts: &Counts,
        error: Option<&dyn std::error::Error>,
    ) -> io::Result<()> {
        self.document["status"] = json!(if error.is_some() {
            "failed"
        } else {
            "completed"
        });
        self.document["written_samples"] = json!(counts.written);
        self.document["completed_samples"] = json!(counts.completed());
        self.document["failed_samples"] = json!(counts.failed);
        self.document["failed_duplicate_samples"] =
            json!(counts.failed_gapped + counts.failed_ungapped);
        self.document["failed_window_samples"] = json!(counts.failed_window);
        self.document["failed_samples_by_last_rejection"] = json!({ "gapped": counts.failed_gapped, "ungapped": counts.failed_ungapped, "window": counts.failed_window });
        let denominator = if self.document["ensure_samples"] == true {
            counts.completed()
        } else {
            self.document["requested_samples"]
                .as_u64()
                .expect("requested samples")
        };
        self.document["failed_duplicate_percentage"] = json!(if denominator == 0 {
            0.0
        } else {
            100.0 * (counts.failed_gapped + counts.failed_ungapped) as f64 / denominator as f64
        });
        let percentage = |count: u64| {
            if denominator == 0 {
                0.0
            } else {
                100.0 * count as f64 / denominator as f64
            }
        };
        self.document["failed_selection_percentage"] = json!(percentage(counts.failed));
        self.document["failed_window_percentage"] = json!(percentage(counts.failed_window));
        self.document["failed_samples_percentage_by_last_rejection"] = json!({
            "gapped": percentage(counts.failed_gapped),
            "ungapped": percentage(counts.failed_ungapped),
            "window": percentage(counts.failed_window),
        });
        self.document["wall_time_seconds"] = json!(self.started.elapsed().as_secs_f64());
        self.document["error"] = json!(error.map(ToString::to_string));
        self.save()
    }

    fn save(&self) -> io::Result<()> {
        let temporary = self.path.with_extension("json.tmp");
        let mut text = serde_json::to_vec_pretty(&self.document).map_err(io::Error::other)?;
        text.push(b'\n');
        std::fs::write(&temporary, text)?;
        std::fs::rename(temporary, &self.path)
    }
}
