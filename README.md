# PF simulator

A Rust executable generating uniform unrooted binary trees and MSAs using
IQ-TREE/AliSim through a small C ABI bridge. Defaults: 10 tips, 500 amino acids, LG,
independent exponential branch lengths with mean 0.1 substitutions/site.

## Build

Requires Rust/Cargo (edition 2024), CMake ≥3.20, a C++17 compiler, Eigen3,
Boost headers, and zlib. Tested on macOS ARM64; a Linux build path is included.
Cargo builds and statically links the native library automatically. The first
build takes a few minutes; Rust-only edits reuse the native build.
IQ-TREE source must be in `dependencies/iqtree2` (tested commit
`a00094e03d1ae984e1497e16738f91514df8c366`). RevBayes is not required.

```sh
cargo build -j 8
cargo test
```

## Generate

```sh
cargo run -- --out example_data --samples 10 --tips 20 --length 500 --seed 42
cargo run -- --out dna_data --seqtype DNA --model JC --tips 10 --length 1000
cargo run -- --out gamma_data --model 'LG+G4{0.5}' --samples 10
cargo run -- --out indel_data --samples 10 --indel 0.002,0.002 \
  --indel-size 'POW{1.7/100},POW{1.7/100}' --no-unaligned
```

Output must not exist. Each pair is `trees/sample_000000.nwk` and
`sequences/sample_000000.fa`, with matching `T1`…`Tn` tip IDs. The first tree
has index zero. `--branch-mean` controls the exponential branch-length mean;
`--help` lists options. Model and sequence type must agree.

MSAs must have unique aligned tip sequences by default. `--max-attempts 20`
sets the maximum number of attempts per fixed tree, including the first attempt.
Retries use distinct deterministic seeds. If the limit is exhausted, the run
exits with an error and leaves no tree or sequence files for that failed sample;
earlier accepted samples remain. `--allow-duplicate-sequences` accepts the first
attempt without checking. Duplicate comparisons include gap positions and use
the full alignment; gap-free outputs can still match when gap placement differs.

AliSim returns numeric sequence rows to Rust through an in-memory adapter.
Duplicate checking compares those rows directly, including gap positions. Only
accepted attempts are converted to FASTA and written to disk. Indel reconstruction
also stays numeric, avoiding the intermediate text round trip. Rows remain in RAM;
tree files are still read by AliSim on each attempt.

`--indel INS,DEL` and `--indel-size 'INS_DIST,DEL_DIST'` pass their values
directly to AliSim, which parses and validates them. Rates are relative to
substitution rates. Omitting `--indel-size` uses AliSim's default size
distributions; it requires `--indel` when supplied. With indels, `--length`
sets the root sequence length; the final gapped alignment can be longer and is
kept in full. Gap-free sequences are saved separately unless you pass
`--no-unaligned`:

```text
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
summary is still printed. Logs include attempt seeds and acceptance/rejection.

Pack with the existing PFceiver utility (Python ≥3.12 and tqdm):

```sh
python /Users/lucblassel/Work/PFceiver/PFceiver/packed_data.py \
  example_data packed_data --tree-dir trees --sequence-dir sequences
```

Runs are sequential. Identical options reproduce the same files with the same
build; increasing `--samples` preserves the existing sample prefix. The Rust
tree generator uses xoshiro256** with SplitMix64 seeding, so its results differ
from the earlier C++ prototype for the same seed. Failed runs
can leave partial output; use a new output directory when retrying.

This example emits tip sequences only and has no
empirical diameter rescaling or birth–death prior. It does not yet
write tar/SQLite shards directly. See [architecture](docs/architecture.md).

Application code is in `src/main.rs` (CLI/output), `src/tree.rs` (trees), and
`src/alignment.rs` (safe AliSim API). The C++ bridge is isolated in `native/`;
`build.rs` and CMake handle its compilation. The CLI uses Clap derive.
Use `cargo build --release -j 8` for an optimized binary in `target/release/pf-sim`.
`cargo test` runs Rust unit tests and the end-to-end checks (Python 3 required).

## Output performance comparison

[Measured results](benchmarks/REPORT.md) compare the original FASTA transport
with the numeric transport now used by default. Run the comparison with:

```sh
cargo run --release --features raw-benchmark --example benchmark_output -- --rounds 9
```

Numeric checking is enabled in normal builds; the feature only enables the
comparison harness. The original FASTA adapter remains as a reference.
The harness checks sequence equivalence and includes FASTA conversion of accepted
numeric output in its timings. It excludes final filesystem writes.
