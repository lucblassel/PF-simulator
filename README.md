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
Cargo fetches `phylotree` from Git; `Cargo.lock` pins the revision.
IQ-TREE source must be in `dependencies/iqtree2` (tested commit
`a00094e03d1ae984e1497e16738f91514df8c366`).

```sh
cargo build -j 8
cargo test
```

### Build on Jean-Zay

```sh
# Clone repo
git clone --recursive https://github.com/lucblassel/PF-simulator.git && cd PF-simulator

# Load correct modules and deps
module load gcc/15.2.0  eigen/3.4.0  boost/1.86.0 rust/1.96.1  cmake/3.31.4

# Compile to fresh dir to avoid re-unsing stale CMAKE
CC=gcc CXX=g++ CARGO_TARGET_DIR=target-jeanzay cargo build --release --locked -j 4

# Use compiled binary
./target-jeanzay/release/pf-sim --help
```

Keep the matching GCC runtime loaded when **running** the executable, including
inside Slurm jobs. For the GCC 15.2 build above, put `module load gcc/15.2.0` in
the job script before launching `pf-sim`. Rust and CMake are build-time tools;
the linked C++ runtime is also needed at execution time. Check the actual loader
choice with `ldd /path/to/pf-sim`: `libstdc++.so.6` should resolve to the GCC module,
not an older system library.

An incompatible runtime can crash inside `std::ostream::sentry` with no stderr
message; GCC documents this stream-initialization incompatibility in its
[GCC 13 release notes](https://gcc.gnu.org/gcc-13/changes.html). A failed simulation
can leave a partial output directory; retries require a new output directory.
When a worker dies before replying, the coordinator reports its exit status or
signal, including SIGSEGV when available, before terminating the remaining workers.

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
Directory output stores metadata as JSON sidecars. Direct packed output embeds
it in SQLite instead. Post-packing directories with PFceiver's existing utility
still omits these sidecars; retain them separately when using that utility.

### Cluster throughput benchmark

```bash
sbatch --partition=YOUR_PARTITION scripts/benchmark_slurm.sh \
  /absolute/path/to/pf-sim /scratch/results
```

The 96-task array tests 50 tips × 500 sites / 200 tips × 5,000 sites,
1/4/8/16 workers, G8/GC, packed/files, and three paired seeds. Each task requests
16 CPUs; add your cluster's account/partition options and optionally `--exclusive`
for isolated CPU measurements. Use `--array=0-95%1` to limit the array to one
active task; the script's `%4` allows concurrent storage load
([Slurm arrays](https://slurm.schedmd.com/job_array.html)).
Defaults are 2,048/64 candidates, ancestors, POW indels, gapped/window checks and
one attempt per tree. Set `PF_BENCH_SAMPLES` to override candidate counts.

Each run saves `metrics.tsv`, its simulation summary and log. Metrics include
accepted/candidate throughput, failures, wall/CPU time, output bytes/files,
node and binary hash. Timing includes binary startup and output/index writing,
excluding Slurm launch, metric collection and cleanup. MiB/s measures logical
output over total runtime; writes are buffered and caches are not flushed.
Generated datasets are deleted after successful measurement; set
`PF_BENCH_KEEP_DATA=1` to retain them. Failed-run outputs are preserved.
Combine results after the array finishes:

```bash
awk 'FNR == 1 && NR != 1 {next} {print}' /scratch/results/bench-*/metrics.tsv > results.tsv
```

### Sustained buffer-size benchmark (8 Slurm jobs)

Rebuild the binary first; this benchmark requires `--shard-buffer-bytes`:

```bash
cargo build --release --locked
sbatch scripts/benchmark_buffers_slurm.sh /absolute/path/to/pf-sim \
  "$WORK/pf-buffer-bench" "$SCRATCH/pf-buffer-bench"
```

Add your cluster's partition/account options to `sbatch`. The array contains only
eight jobs and defaults to **one active job** (`0-7%1`) to reduce self-contention.
Each requests 16 CPUs, 32 GiB RAM, and one hour. Tasks 0/1 use 256 KiB per-worker
buffers, 2/3 use 1 MiB, 4/5 use 4 MiB, and 6/7 use 16 MiB. Even tasks write to
the supplied WORK directory; odd tasks use SCRATCH. All use the same seed and
simulation options: Python preset with G8, 200 tips, 5,000 sites, ancestors,
POW indels, gapped/window checks, 1 GiB shards, and one attempt per tree.

`--ensure-samples` requests **10,000 accepted outputs per job**; rejected candidates
are replaced. Allow roughly 25 GiB of temporary output per active job, potentially
more with indels. Successful datasets are deleted after measurement. Set
`PF_BENCH_KEEP_DATA=1` to retain them; failed output always remains. Small metadata,
logs and summaries remain under each filesystem's results directory.

Each point records `metrics.csv`, per-run commands/summaries/logs, CPU model and
affinity, timestamps, executable hash, worker/buffer settings, failures, wall/CPU
time, shard/file counts, and logical output MiB/s. This is **end-to-end simulator
throughput**, including final index merging, not a raw filesystem bandwidth test.
Data enumeration and cleanup are outside timing; OS caches are not flushed.

Optional environment overrides are `PF_BENCH_WORKERS` (within allocated CPUs),
`PF_BENCH_SAMPLES`, and `PF_BENCH_REPEATS`. Repeats run sequentially within each
of the same eight jobs; matching repeats use matching seeds across all points.
The default is one run per point, so it does not estimate run-to-run variability.
Combine results after the array finishes:

```bash
awk 'FNR == 1 && NR != 1 {next} {print}' \
  "$WORK/pf-buffer-bench"/buffer-*/metrics.csv \
  "$SCRATCH/pf-buffer-bench"/buffer-*/metrics.csv > buffer-metrics.csv
```

### Local output and SQLite benchmark

```bash
cargo build --release --locked --example benchmark_packed_io --bin pf-sim
python3 benchmarks/benchmark_local_sqlite.py \
  --writer current=target/release/examples/benchmark_packed_io \
  --simulator current=target/release/pf-sim \
  --out /path/on/test/filesystem/local-benchmark
```

This runs sequentially on macOS or Linux, randomizing the order of configurations.
It compares packed/files output at 1/4/8 workers with three repeats. The writer-only
cases use prebuilt 64 KiB / 2.6 MB payloads and report database setup, write,
finalization and index-merge timings. Their workers are independent threads; the
optional full simulator cases use the actual process pool and the cluster's G8
simulation settings. For before/after comparisons, save copies of both release
executables before editing, then supply additional `--writer old=/path/to/copy`
and `--simulator old=/path/to/copy` arguments.

Results include `metrics.csv`, `REPORT.md`, executable hashes, and a `runs.tar.gz`
archive of logs and timing JSON files. Output datasets are removed after each successful run unless
`--keep-data` is supplied; failed output is retained. Choose a new output directory
on the filesystem you want to measure. For sustained throughput, increase
`--large-samples` (for example to 512). `--shard-size-bytes` can exercise rotation.
Use `--workloads small --formats packed`, `--small-samples 8192` and
`--small-payload-bytes 256` for an index-heavy case that crosses the former
1,000-row commit interval. These payload size overrides affect writer-only runs.
SQLite durability settings remain unchanged; tar/FASTA writes are buffered, and
the benchmark does not flush the OS cache. Mac results do not predict Lustre latency.

### Direct packed output

```bash
pf-sim --out dataset/batches/job-000 --output-format packed \
  --preset python --samples 10000 --jobs 8 --seed 42 \
  --batch-id job-000 --file-prefix job-000
```

Workers write uncompressed tar files into one `shards/` directory, named
`worker-000-shard-000000.tar`, `worker-001-shard-000000.tar`, etc. There are no
per-worker directories. The coordinator creates the root `index.sqlite3` after
all workers finish and removes their temporary indices.
Tar members include their writer ID and have matching stems under `trees/`
and `sequences/`. Optional `--file-prefix job-000` adds `job-000-` to shard
and member filenames, for example:

```text
shards/job-000-worker-003-shard-000000.tar
  trees/job-000-worker-003-sample_000042.nwk
  sequences/job-000-worker-003-sample_000042.fa
```

Use a **unique file prefix per batch** to copy all archives into one directory
or extract them into a shared `trees/` and `sequences/` filetree without
collisions. `--batch-id` namespaces index records only; it does not prefix files.
Prefixes accept 1–40 ASCII letters, digits, dots, underscores or hyphens.
In directory mode the prefix applies to sample files and metadata sidecars,
e.g. `msa_untrimmed/job-000-sample_000042.fa`. Index and run-summary filenames
stay fixed. Naming does not affect sampling or RNGs; logical sample names remain
worker-independent, while `pairs.tree_name` and `sequence_name` record physical
archive paths. Saved replay commands preserve the prefix.

Only one full MSA and its tree are stored per example; `--ancestral` includes
internal sequences in that MSA. `--materialize` and `--trees-only` require file
output. The default `--output-format files` remains available.

`--shard-size-bytes` defaults to 1 GiB, including tar headers, padding and end
markers. Pairs stay together; an oversized pair gets its own shard. The index
uses PFceiver's version-3 `pairs` schema and is readable by `PackedPairs`.
Each worker keeps one SQLite transaction open for its current shard, including
initial schema creation. The transaction commits after closing that shard, or
at worker completion. Worker indices are merged only after all workers finish.
Each active shard uses a 1 MiB output buffer to coalesce the tar library's small
streaming writes. Set `--shard-buffer-bytes` to tune it independently of shard
size; the setting is retained in the summary and replay command.
Window coordinates and provenance are stored in
`sample_metadata(name, candidate_id, metadata_json)`; `dataset_metadata` records
the metadata extension version. The run summary remains at the dataset root.
No per-example files or extra FASTA views are written.

For cluster jobs, use separate output directories, distinct seeds and unique
batch IDs. Without explicit values, seeds are generated and batch IDs default
to `seed-<seed>`. Batch IDs name examples during merging; they do not change RNGs.
Merge completed batches without copying or repacking their tar files:

```bash
python3 scripts/merge_packed.py dataset/batches/job-* --out dataset
```

Inputs must be below `dataset`; keep the whole directory tree when relocating
it. The merger assigns dense IDs and names examples as `<batch-id>/<sample>`.
It preserves batch summaries in `batches` and maps examples in `sample_batches`.
Repeated identities, incomplete runs and conflicting policies are rejected;
`--allow-mixed-configs` permits an intentional mixture. Existing indices are
never overwritten. Sample contents are reproducible across worker counts;
worker assignment and shard bytes/order can vary.

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
without serializing them through FASTA first. Packed mode serializes each
accepted full MSA once into a tar member. See [architecture](docs/architecture.md).

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
- [x] Simulate directly into the packed and sharded dataset format
