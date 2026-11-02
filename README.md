# PF simulator


A Rust executable generating phylogenetic trees with `phylotree` and MSAs using
IQ-TREE/AliSim through a small C ABI bridge. Defaults: 10 tips, 500 amino acids, LG,
independent exponential branch lengths with mean 0.1 substitutions/site.

***Disclaimer:** This piece of software is being developed with the help of an LLM.*

## Build

Requires Rust/Cargo (edition 2024), CMake ≥3.20, a C++17 compiler, Eigen3,
Boost headers, and zlib. Tested on macOS ARM64; a Linux build path is included.
Cargo builds and statically links the native library automatically. The first
build takes a few minutes; Rust-only edits reuse the native build.
The local `phylotree` checkout must include the new `simulation` module; Cargo
uses `../../Development/rust/phylotree` (adjust the path for another layout).
IQ-TREE source must be in `dependencies/iqtree2` (tested commit
`a00094e03d1ae984e1497e16738f91514df8c366`). RevBayes is not required.

```sh
cargo build -j 8
cargo test
```

## Generate

```sh
cargo run -- --out example_data --samples 10 --tips 20 --length 500 --seed 42
cargo run -- --out parallel_data --samples 1000 --tips 200 --jobs 4
cargo run -- --out dna_data --seqtype DNA --model JC --tips 10 --length 1000
cargo run -- --out gamma_data --model 'LG+G4{0.5}' --samples 10
cargo run -- --out bd_data --preset python --samples 10 --tips 50
cargo run -- --out ete_data --tree-prior ete3 --empirical-diameter
cargo run -- --out gamma_trees --branch-prior gamma --gamma-shape 2 --gamma-scale 0.05
cargo run -- --out indel_data --samples 10 --indel 0.002,0.002 \
  --indel-size 'POW{1.7/100},POW{1.7/100}' --no-unaligned
```

Output must not exist. Each pair is `trees/sample_000000.nwk` and
`msa_untrimmed/sample_000000.fa`, with matching `T1`…`Tn` tip IDs. The first tree
has index zero. `--branch-mean` controls the exponential branch-length mean;
`-h` shows the common options; `--help` lists all options with details and examples.
Model and sequence type must agree.

`--seed` is optional. When omitted, the parent generates one random master seed
and passes it to every worker. Supply `--seed N` when you want a known,
reproducible dataset. The resolved seed is printed and saved before simulation.

Every output directory contains `simulation_summary.json`, with:

- The master `seed`, its source (`provided` or `generated`), resolved CLI options,
  executable/version/build profile, and a Bash `reproduce_command`.
- Requested, completed, and written sample counts; `failed_duplicate_samples`
  and `failed_duplicate_percentage` (percentage of requested candidates by default;
  percentage of completed candidates with `--ensure-samples`). The summary records
  `sample_count_mode`, `ensure_samples`, and the percentage denominator explicitly.
- `wall_time_seconds`, measured with a monotonic clock from coordinator startup
  before CLI parsing through generation, retries, output writes, worker shutdown,
  and log merging. Final summary serialization and completion printing are excluded.

The replay command spells out defaults, restores the working directory, and
creates a fresh output directory under `mktemp -d`; it also redirects any log to
that fresh directory. It requires the same executable/build and unchanged external
inputs (such as custom model or empirical files), if supplied. Tree and sequence files reproduce exactly; timing,
summary metadata, and parallel log ordering vary.

The summary starts with `status: running` and null result/timing fields, then is
atomically updated to `completed`, or `failed` for errors returned to the parent.
On fatal errors, counts reflect results acknowledged before stopping. An abrupt
termination (including upstream AliSim calling `exit()` in single-job mode) can
leave the initial `running` record, preserving the seed and replay command;
its null timing must not be treated as a completed measurement.

`-j N` / `--jobs N` runs up to N persistent worker processes. The default is 1,
which uses the direct sequential path. Workers generate trees, retry/check MSAs,
and write their own files; the parent collects progress and failure counts.
The same options and seed produce identical sample files for every job count.
Workers each have independent IQ-TREE state; larger job counts also use more RAM.
See the [parallelism benchmarks](benchmarks/parallelism/REPORT.md) for measured
whole-run speedups, including output I/O. The
[Python reference comparison](benchmarks/python/REPORT.md) measures complete
no-indel, duplicate-checked runs at 1, 4 and 8 worker/thread budgets, including
both the original Python threading strategy and a Python process pool. An
[indel follow-up](benchmarks/python_indels/REPORT.md) measures raw aligned and
unaligned output with duplicate checks disabled, including ancestral-output
and tip-only Python references.
The [full pipeline benchmark](benchmarks/python_pipeline/REPORT.md) additionally
matches empirical GC, ancestral simulation, tip-window selection and materialized
full/selected/unaligned output on both sides, with one attempt per tree.

`--preset default` (also used when omitted) selects uniform unrooted topologies
with exponential branch lengths of mean 0.1. `--preset python` selects the Python
birth-death pipeline: process-time branches, inherited rate changes, empirical
diameter scaling, then terminal repair below 0.001. It also defaults to LG with
empirical eight-category Gamma (G8) site rates. An explicit `--model` overrides this
alignment default; `--gamma` overrides the empirical rate mode. Numeric parameters such as `--branch-mean`,
`--birth-rate`, and `--minimum-terminal` remain configurable. The Python preset
requires birth-death topology and process-time branches; use the default preset
with individual tree options to build other pipelines.

Tree priors include uniform labeled rooted/unrooted topologies, ETE3, Yule and
DendroPy-compatible first-hit birth–death. `--python-bd` remains a compatibility
shorthand for `--preset python`. ETE3’s sampler is distinct from Yule.
See [tree priors](docs/tree-priors.md)
and the [distribution comparison](validation/tree_priors/REPORT.md).
`--trees-only` skips MSA generation.

`--gamma G8` or `--gamma GC` enables empirical site rates with any base model
(e.g. `--model WAG --gamma G8`). `--gamma none` disables empirical sampling;
fixed shapes remain available through `--model 'LG+G4{0.5}'`. On **each attempt**,
a uniformly selected bundled HOGENOM alpha is jittered by Normal(alpha, alpha/10)
and floored at 0.05, matching `alisim.py`. All 12,408 observations are included
without filtering or deduplication; no external alpha file or Python is needed.
This is separate from `--branch-prior gamma`.

**GC can be much slower than G8**, especially on large trees and alignments;
see the [Gamma benchmark](benchmarks/gamma_costs/REPORT.md). Use `--gamma GC`
explicitly to match the Python pipeline's continuous site-rate prior. G8 uses
eight discrete rate categories and changes that prior.

Gamma draws depend only on the master seed, numeric sample ID, and zero-based
attempt; they reproduce across worker counts. The summary records the empirical
source/payload SHA-256 hashes, resolved MSA and selection settings, and the exact
seed derivation. `--alisim-log` additionally records each attempted model and
engine seed. To reconstruct one alpha (candidate 12, attempt 2):

```bash
cargo run --example sample_gamma -- --seed 42 --candidate 12 --attempt 2 --count 1
```

Rust and Python use different RNGs: they sample the same prior, but equal numeric
seeds do not produce the same draws across languages. See the
[empirical Gamma comparison](validation/empirical_gamma/REPORT.md).

MSAs must have unique aligned tip sequences by default. `--max-attempts 20`
sets the maximum number of attempts per fixed tree, including the first attempt.
Retries use distinct deterministic seeds. If the limit is exhausted, that sample
is skipped and counted, and the run continues. No tree or sequence files are
written for skipped samples. The completion summary reports the number written
and the number that failed due to duplicates; the run exits successfully even
if all samples were skipped. `--samples` counts requested simulations, so the
output can contain fewer pairs and gaps in sample IDs. Later sample IDs and
seeds are unchanged. Other simulation and I/O errors still stop the run with a
nonzero exit status. `--allow-duplicate-sequences` accepts the first
attempt and overrides **all** tip checks. Checks are independent:

- `--require-distinct-gapped[=true|false]`: distinct full aligned tip rows;
  enabled by default. Use `--require-distinct-gapped=false` to disable it.
- `--require-distinct-ungapped[=true|false]`: distinct full tip rows after gap
  removal; disabled by default.
- `--require-unique-window[=true|false]`: a length-L window containing distinct
  aligned tip rows; disabled by default. Internal sequences never participate.

Each failed check retries the MSA on the fixed tree, using the same maximum-attempt
and optional replacement policy. Checks run in gapped, ungapped, window order;
logs identify the first failed check. The summary separates exhausted samples by
last rejection (`gapped`, `ungapped`, `window`) and records total failed samples.
Window failures are recoverable just like full-sequence duplicate failures.

`--window-policy inclusive` searches every valid start, including the last one.
`--window-policy python` preserves the legacy excluded final start for M > L.
The default is `python` under the Python preset, otherwise `inclusive`. Both
accept a unique full window at start zero when M == L, and reject M < L when
window checking is enabled. No columns are removed before checking.

```sh
cargo run --release -- --out selected_data --samples 100 --tips 50 \
  --indel 0.002,0.002 --require-distinct-ungapped --require-unique-window \
  --materialize trimmed,unaligned
```

Use `--ensure-samples` to make `--samples` the number of successful outputs instead:

```sh
cargo run --release -- --out complete_data --samples 100 --ensure-samples --jobs 4
```

After a tree exhausts `--max-attempts`, a fresh tree and seed are generated until
exactly the requested number of outputs is written. Only successful outputs
advance the progress bar. Failed candidates remain counted in the summary and
logs; `completed_samples` counts successes plus failures and can exceed the target.
Candidate IDs keep increasing, so output filenames can still have gaps. With a
fixed seed, outputs are identical across worker counts, without extra speculative
output files. Other simulation or I/O errors still stop the run.

**Warning:** `--ensure-samples` can significantly increase runtime and may never
finish if the settings cannot generate duplicate-free outputs. The CLI prints
this warning when the flag is used. With `--allow-duplicate-sequences` or
`--trees-only`, every successful candidate already counts toward the target.

AliSim returns numeric sequence rows to Rust through an in-memory adapter.
Duplicate checking compares those rows directly, including gap positions. Only
accepted attempts are converted to FASTA and written to disk. Indel reconstruction
also stays numeric, avoiding the intermediate text round trip. Rows remain in RAM.
The Rust tree is passed as numeric edges and branch lengths; the bridge builds
IQ-TREE nodes directly without Newick conversion or tree-file reads. Accepted
trees are serialized once to Newick for the output dataset.

`--indel INS,DEL` and `--indel-size 'INS_DIST,DEL_DIST'` pass their values
directly to AliSim, which parses and validates them. Rates are relative to
substitution rates. Omitting `--indel-size` uses AliSim's default size
distributions; it requires `--indel` when supplied. With indels, `--length`
sets the root sequence length; the final gapped alignment can be longer and is
kept in full. The default output stores one MSA plus per-sample metadata:

```text
simulation_summary.json
trees/sample_000000.nwk
msa_untrimmed/sample_000000.fa
metadata/sample_000000.json
```

Metadata records the actual model, tree/engine seeds, accepted attempt, earlier
rejection counts, tip identities, alignment length and selected window. Window
coordinates are **zero-based, half-open full-alignment column offsets**. If window
checking is disabled (including `--allow-duplicate-sequences`), the window is an
explicitly unchecked prefix `[0, min(L, M))`; it may be shorter than L.

`--materialize trimmed` additionally writes `msa_trimmed/sample_000000.fa`.
`--materialize unaligned` writes `unaligned/sample_000000.fa`, removing gaps from
**full** rows, not from the selected window. Use `--materialize trimmed,unaligned`
or `--materialize all` for both. These exports also work without indels and with
duplicates allowed. They neither change simulation RNGs nor acceptance. Neither
additional directory is created unless requested. The older `--no-unaligned`
flag remains accepted but conflicts with requesting an unaligned export.

A training loader can reproduce these views from the full MSA and metadata;
materializing them is optional. All requested sample files have the same stem.
For now metadata is a JSON sidecar per sample; a future packed format should embed
it in each sample record to avoid millions of extra filesystem entries. The
current `packed_data.py` tree/sequence format does **not** yet preserve these
sidecars, so retain them separately until that format is extended.

To pack raw, trimmed or unaligned sequences with the existing packer, select
`--sequence-dir msa_untrimmed`, `msa_trimmed` or `unaligned`, respectively.
See [AliSim's indel documentation](https://iqtree.github.io/doc/AliSim#insertion-and-deletion-models)
for distribution syntax.

A progress bar tracks completed pairs on the terminal. AliSim's informational
output is hidden by default; save it with `--alisim-log alisim.log` (the log file
must be new and its parent directory must exist). Native errors remain visible.
When stderr is not a terminal, the progress bar is hidden and the completion
summary is still printed. Logs include attempt seeds, acceptance/rejection, and
skipped sample IDs with their exhausted attempt counts.
With multiple jobs, workers use separate temporary log fragments next to the
requested log, then merge them into that log in worker order after stopping.
The merged log groups samples by worker rather than by sample ID. Fragments
are removed after a successful merge; a failed merge preserves diagnostics.

Pack with the existing PFceiver utility (Python ≥3.12 and tqdm):

```sh
python /Users/lucblassel/Work/PFceiver/PFceiver/packed_data.py \
  example_data packed_data --tree-dir trees --sequence-dir msa_untrimmed
```

The same seed and options reproduce the same sample files with the same
build; increasing `--samples` preserves the existing sample prefix. The Rust
tree generator now uses seeded `phylotree` samplers; its seed-to-tree mapping
differs from earlier PF-only and Python versions. Sample seeds still use
SplitMix64 derivation. Runs stopped by fatal errors
can leave partial output; use a new output directory when retrying.
In parallel mode, a fatal error stops dispatch and terminates/reaps all workers;
already completed outputs remain, and in-flight writes may be incomplete.

Tip-only output is the default. `--ancestral` additionally exports simulated
internal rows in full, selected and gap-free outputs. Their `I_<id>` labels match
the internal labels in the saved tree; metadata lists the ancestral IDs. All
acceptance checks remain tip-only. The native numeric transport includes ancestors
without serializing them through FASTA first. No tar/SQLite shards are written
directly. See [architecture](docs/architecture.md).

Application code is in `src/main.rs` (CLI/output), `src/tree.rs` (trees), and
`src/alignment.rs` (safe AliSim API). `src/workers.rs` coordinates process workers.
The C++ bridge is isolated in `native/`;
`build.rs` and CMake handle its compilation. The CLI uses Clap derive.
Use `cargo build --release -j 8` for an optimized binary in `target/release/pf-sim`.
`cargo test` runs Rust unit tests and the end-to-end checks (Python 3 required).

## Distribution validation

The [MSA comparison report](validation/msa_priors/REPORT.md) measures the linked
engine against independent standalone AliSim 2.4.0 and 2.2.2.6 executables. It also
compares the current Python preset with the actual legacy Python MSA pipeline,
including empirical gamma sampling, indels, and downstream window selection.
That report measured the earlier tree-only preset. The current Python preset
now defaults to empirical G8; use `--gamma GC` to reproduce the Python pipeline's
continuous site-rate prior. Indels, ancestral output and window selection require
explicit options.
See the report for measured differences, scope, raw metrics, and reproduction steps.
A [raw Python comparison](validation/python_raw/REPORT.md) separately matches
empirical Gamma/indel settings and compares `.untrimmed` and unaligned tip
outputs before duplicate rejection or window selection.
A [fixed-tree probe](validation/msa_priors/branch_probe/REPORT.md) diagnoses a
measured unrooted-tree discrepancy in the installed older binary.

## Bridge performance comparisons

[Measured results](benchmarks/REPORT.md) separately compare FASTA versus numeric
output and file-based versus direct tree input. Build the harness once, then run:

```sh
cargo build --release --features raw-benchmark --example benchmark_output
./target/release/examples/benchmark_output --comparison tree-input --rounds 9 --target-ms 200
./target/release/examples/benchmark_output --comparison output --rounds 9 --target-ms 200
python3 benchmarks/summarize.py
```

Numeric checking is enabled in normal builds; the feature only enables the
comparison harness. The original FASTA adapter remains as a reference.
The harness checks sequence equivalence and includes FASTA conversion of accepted
numeric output in its timings. It excludes final filesystem writes.

# TODOs
- [x] Compare output to previous python pipeline thoroughly
- [x] Implement topology choices  
    (need revbayes or MrBayes dep?, implement these in [`phylotree-rs`](https://github.com/lucblassel/phylotree-rs) and use that as a dependecy?)
    - [x] Birth death
    - [x] Uniform
    - [x] ete3's "uniform"
    - [x] Gamma branch lengths
    - [x] Rescaled birth-death branches
    - [x] rescaling to empirical diameters
- [x] Implement parallel simulation
- [ ] Simulate directly into new packed and sharded dataset format
