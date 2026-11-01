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
`sequences/sample_000000.fa`, with matching `T1`…`Tn` tip IDs. The first tree
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
  and `failed_duplicate_percentage` (percentage of requested samples).
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
whole-run speedups, including output I/O.

`--preset default` (also used when omitted) selects uniform unrooted topologies
with exponential branch lengths of mean 0.1. `--preset python` selects the Python
birth-death pipeline: process-time branches, inherited rate changes, empirical
diameter scaling, then terminal repair below 0.001. Presets affect tree generation;
alignment settings remain independent. Numeric parameters such as `--branch-mean`,
`--birth-rate`, and `--minimum-terminal` remain configurable. The Python preset
requires birth-death topology and process-time branches; use the default preset
with individual tree options to build other pipelines.

Tree priors include uniform labeled rooted/unrooted topologies, ETE3, Yule and
DendroPy-compatible first-hit birth–death. `--python-bd` remains a compatibility
shorthand for `--preset python`. ETE3’s sampler is distinct from Yule.
See [tree priors](docs/tree-priors.md)
and the [distribution comparison](validation/tree_priors/REPORT.md).
`--trees-only` skips MSA generation.

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
attempt without checking. Duplicate comparisons include gap positions and use
the full alignment; gap-free outputs can still match when gap placement differs.

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
kept in full. Gap-free sequences are saved separately unless you pass
`--no-unaligned`:

```text
simulation_summary.json
trees/sample_000000.nwk
sequences/sample_000000.fa
unaligned_sequences/sample_000000.fa
```

`unaligned_sequences/` is created only when gap-free output is generated.
To pack gap-free sequences, use `--sequence-dir unaligned_sequences` instead
of `--sequence-dir sequences` in the packing command below.
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
  example_data packed_data --tree-dir trees --sequence-dir sequences
```

The same seed and options reproduce the same sample files with the same
build; increasing `--samples` preserves the existing sample prefix. The Rust
tree generator now uses seeded `phylotree` samplers; its seed-to-tree mapping
differs from earlier PF-only and Python versions. Sample seeds still use
SplitMix64 derivation. Runs stopped by fatal errors
can leave partial output; use a new output directory when retrying.
In parallel mode, a fatal error stops dispatch and terminates/reaps all workers;
already completed outputs remain, and in-flight writes may be incomplete.

This example emits tip sequences only and does not yet write tar/SQLite shards
directly. See [architecture](docs/architecture.md).

Application code is in `src/main.rs` (CLI/output), `src/tree.rs` (trees), and
`src/alignment.rs` (safe AliSim API). `src/workers.rs` coordinates process workers.
The C++ bridge is isolated in `native/`;
`build.rs` and CMake handle its compilation. The CLI uses Clap derive.
Use `cargo build --release -j 8` for an optimized binary in `target/release/pf-sim`.
`cargo test` runs Rust unit tests and the end-to-end checks (Python 3 required).

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
- [ ] Compare output to previous python pipeline thoroughly
- [ ] Implement topology choices  
    (need revbayes or MrBayes dep?, implement these in [`phylotree-rs`](https://github.com/lucblassel/phylotree-rs) and use that as a dependecy?)
    - [ ] Birth death
    - [x] Uniform
    - [ ] ete3's "uniform"
    - [ ] Gamma branch lengths
    - [ ] Rescaled birth-death branches
    - [ ] rescaling to empirical diameters
- [ ] Implement parallel simulation
- [ ] Simulate directly into new packed and sharded dataset format
