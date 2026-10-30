# AliSim bridge benchmarks

Recorded on 2026-10-30 on an Apple M2 Pro (10 physical cores, 32 GiB RAM),
macOS 27.0.1, Rust 1.98.1. Rust/native code use release optimization; AliSim
runs on one thread with the same portable likelihood kernel in both paths.

## What each comparison measures

- **Tree input:** numeric output on both sides; file/Newick input versus direct
  Rust edges converted to native IQ-TREE nodes. This isolates the latest bridge change.
- **Output transport:** FASTA versus numeric output, with file/Newick tree input
  on both sides. Rebuilding and rerunning the original harness does **not** measure
  direct tree input. These results describe the earlier output optimization.

The CLI now uses direct tree input and numeric output. Accepted trees are still
written once as Newick dataset artifacts. That final output cost is excluded here.

## Direct tree input: isolated effect

Plain LG paired gains: **+1.2%** at 10 tips × 500 sites (+58 µs/attempt); **+0.9%** at 50 tips × 500 sites (+53 µs/attempt); **+0.6%** at 200 tips × 5,000 sites (+133 µs/attempt).

At 200 tips, JC DNA gains **1.7%**; the duplicate-heavy LG case gains **2.1%**.

All continuous-gamma cases have gain IQRs crossing zero: there is no consistent
acceleration in those workloads. Direct input is a modest overhead reduction;
the large earlier LG-with-indels gain came from numeric output, not tree input.

Milliseconds are medians of round averages. Positive gain/saved time means direct
input is faster. Paired statistics use matching seeds within each round; the median
paired gain need not equal the ratio of the two displayed median times.

| Tips × root sites | Model | Branch mean | Numeric + file tree ms | Numeric + direct tree ms | Paired gain | Gain IQR | Saved µs | Rounds | Rejected |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 4 × 2 | JC DNA | 0.5 | 4.820 | 4.788 | +0.7% | +0.4…+1.9% | +34.9 | 7 | 93.4% |
| 10 × 500 | JC DNA | 0.1 | 4.893 | 4.848 | +0.7% | +0.5…+1.2% | +32.7 | 9 | 0.0% |
| 10 × 500 | LG | 0.1 | 4.978 | 4.935 | +1.2% | +0.8…+1.7% | +58.2 | 9 | 0.0% |
| 10 × 500 | LG + continuous gamma | 0.1 | 18.757 | 18.649 | +0.3% | -0.2…+0.7% | +55.7 | 9 | 0.0% |
| 10 × 500 | LG + continuous gamma + indels | 0.1 | 18.817 | 18.757 | +0.2% | -0.0…+0.3% | +38.5 | 9 | 0.0% |
| 10 × 500 | LG + indels | 0.1 | 5.031 | 4.987 | +0.9% | +0.6…+1.0% | +46.5 | 9 | 0.0% |
| 50 × 500 | JC DNA | 0.1 | 5.271 | 5.216 | +0.9% | +0.5…+1.5% | +47.6 | 9 | 0.0% |
| 50 × 500 | LG | 1e-08 | 5.006 | 4.922 | +1.7% | +0.9…+2.0% | +83.2 | 7 | 100.0% |
| 50 × 500 | LG | 0.1 | 5.548 | 5.480 | +0.9% | +0.8…+1.1% | +52.6 | 9 | 0.0% |
| 50 × 500 | LG + continuous gamma | 0.1 | 82.406 | 82.441 | -0.1% | -0.2…+0.2% | -44.1 | 9 | 0.0% |
| 50 × 500 | LG + continuous gamma + indels | 0.1 | 83.044 | 82.880 | -0.1% | -0.3…+0.4% | -88.2 | 9 | 0.0% |
| 50 × 500 | LG + indels | 0.1 | 5.674 | 5.588 | +1.1% | +0.2…+1.9% | +61.2 | 9 | 0.0% |
| 200 × 5,000 | JC DNA | 0.1 | 19.571 | 19.255 | +1.7% | +1.5…+2.2% | +326.9 | 9 | 0.0% |
| 200 × 5,000 | LG | 1e-08 | 8.778 | 8.615 | +2.1% | +1.3…+2.4% | +184.3 | 7 | 100.0% |
| 200 × 5,000 | LG | 0.1 | 21.454 | 21.299 | +0.6% | +0.4…+1.0% | +133.3 | 9 | 0.0% |
| 200 × 5,000 | LG + continuous gamma | 0.1 | 4053.221 | 4054.790 | +0.1% | -0.2…+0.7% | +5820.3 | 9 | 0.0% |
| 200 × 5,000 | LG + continuous gamma + indels | 0.1 | 4070.169 | 4067.293 | +0.1% | -0.0…+1.0% | +3291.0 | 9 | 0.0% |
| 200 × 5,000 | LG + indels | 0.1 | 28.527 | 28.558 | +0.8% | -0.8…+1.0% | +223.1 | 9 | 0.0% |

Sources: [tree-input-comparison-1793378508.csv](tree-input-comparison-1793378508.csv), [tree-input-comparison-1793378737.csv](tree-input-comparison-1793378737.csv), [tree-input-comparison-1793378739.csv](tree-input-comparison-1793378739.csv).

File input reads an already-created, usually cached Newick file. Tree generation
and the original tree-file write happen outside timing. Direct input includes
edge-array preparation, validation, native-node construction and cleanup. Both
sides include numeric row validation/checking and FASTA serialization of accepted
attempts. These are per-attempt costs, not complete dataset throughput.

## FASTA versus numeric output

LG with indels paired gains from numeric output: **+10.2%** at 10 × 500; **+31.5%** at 50 × 500; **+80.9%** at 200 × 5,000.

These gains include avoiding AliSim's indel numeric-text round trip; they are
not solely savings from duplicate checking. Both paths still read a tree file.

| Tips × root sites | Model | Branch mean | FASTA ms | Numeric ms | Paired gain | Gain IQR | FASTA check share | Rounds | Rejected |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 4 × 2 | JC DNA | 0.5 | 4.710 | 4.695 | -0.1% | -1.5…+1.5% | 0.01% | 9 | 94.6% |
| 10 × 500 | JC DNA | 0.1 | 4.877 | 4.861 | -0.1% | -0.7…+0.6% | 0.13% | 9 | 0.0% |
| 10 × 500 | LG | 0.1 | 4.955 | 4.941 | -0.3% | -0.8…+1.6% | 0.13% | 9 | 0.0% |
| 10 × 500 | LG + continuous gamma | 0.1 | 18.959 | 18.977 | -0.4% | -0.6…-0.0% | 0.04% | 9 | 0.0% |
| 10 × 500 | LG + continuous gamma + indels | 0.1 | 19.545 | 19.012 | +3.0% | +2.8…+3.0% | 0.04% | 9 | 0.0% |
| 10 × 500 | LG + indels | 0.1 | 5.588 | 5.062 | +10.2% | +8.8…+10.4% | 0.12% | 9 | 0.0% |
| 50 × 500 | JC DNA | 0.1 | 5.282 | 5.215 | +1.0% | +0.7…+1.4% | 0.59% | 9 | 0.0% |
| 50 × 500 | LG | 1e-08 | 4.813 | 4.802 | +0.1% | -0.3…+0.7% | 0.03% | 9 | 100.0% |
| 50 × 500 | LG | 0.1 | 5.548 | 5.523 | +0.3% | -0.4…+0.5% | 0.56% | 9 | 0.0% |
| 50 × 500 | LG + continuous gamma | 0.1 | 83.777 | 83.741 | +0.1% | -0.2…+0.2% | 0.04% | 9 | 0.0% |
| 50 × 500 | LG + continuous gamma + indels | 0.1 | 87.779 | 84.465 | +2.8% | +2.7…+3.3% | 0.04% | 9 | 0.0% |
| 50 × 500 | LG + indels | 0.1 | 8.343 | 5.719 | +31.5% | +30.9…+31.8% | 0.41% | 9 | 0.0% |
| 200 × 5,000 | JC DNA | 0.1 | 18.914 | 19.058 | -1.2% | -1.3…-0.7% | 5.79% | 9 | 0.0% |
| 200 × 5,000 | LG | 1e-08 | 8.034 | 8.133 | -0.8% | -1.4…-0.3% | 0.10% | 9 | 100.0% |
| 200 × 5,000 | LG | 0.1 | 21.075 | 21.198 | -0.7% | -0.9…-0.4% | 5.17% | 9 | 0.0% |
| 200 × 5,000 | LG + continuous gamma | 0.1 | 4130.770 | 4136.657 | -0.0% | -0.4…+0.0% | 0.03% | 9 | 0.0% |
| 200 × 5,000 | LG + continuous gamma + indels | 0.1 | 4297.968 | 4170.629 | +2.8% | +2.6…+2.9% | 0.04% | 9 | 0.0% |
| 200 × 5,000 | LG + indels | 0.1 | 145.072 | 27.833 | +80.9% | +80.7…+81.0% | 1.05% | 9 | 0.0% |

Sources: [output-comparison-indel-lg.csv](output-comparison-indel-lg.csv), [output-comparison-retries-dna.csv](output-comparison-retries-dna.csv), [output-comparison-retries-lg.csv](output-comparison-retries-lg.csv), [output-comparison.csv](output-comparison.csv).

## Method and limits

- One fixed uniform tree per workload (tree seed 42). Matching alignment seeds
  within each pair; backend order alternates each round. Round counts are shown above.
- Batches target 120–200 ms depending on invocation, bounded to 1–128 attempts.
  Slow simulations use one attempt per round.
- Timing includes native initialization, simulation, transport, duplicate checking,
  accepted FASTA/gap-free serialization and cleanup. It excludes final filesystem
  writes, tree generation, progress display and optional logging. Rejected numeric
  attempts are not formatted as FASTA. Disk tree reads are included only in file paths.
- Indels use rates 0.002,0.002 and POW{1.7/100},POW{1.7/100} sizes; continuous
  gamma shape is 0.5. Branch mean is 0.1 except where specified in the tables.
- The corrected tree-input harness compares aligned/gap-free warm-up sequences across
  FASTA, numeric/file and numeric/direct paths. Timed pairs assert identical duplicate
  decisions for every seed. The earlier output runs compare FASTA and numeric/file.
  Duplicate-heavy workloads time individual attempts, not full retry loops.
- Gain IQR is the middle half of paired gains, not a confidence interval. Positive
  values mean the candidate is faster; a range crossing zero indicates inconsistent
  gains across rounds. No CPU pinning or memory measurements; small differences
  should be treated cautiously. Cross-run timing changes do not isolate a code change.

## Reproduce

Build once before timing; do not run builds or competing benchmarks during trials.

```sh
cargo build --release --features raw-benchmark --example benchmark_output -j 8
# Isolate direct tree input (numeric output on both sides).
./target/release/examples/benchmark_output --comparison tree-input --rounds 9 --target-ms 200
./target/release/examples/benchmark_output --comparison tree-input --cases JC --workloads 4x2 --branch-mean 0.5
./target/release/examples/benchmark_output --comparison tree-input --cases LG --workloads 50x500,200x5000 --branch-mean 1e-8
# Original output-transport comparison (file input on both sides).
./target/release/examples/benchmark_output --comparison output --rounds 9 --target-ms 200
python3 benchmarks/summarize.py
```

The harness creates a timestamped CSV; `--output FILE` selects a new filename.
With no arguments the summary script selects the latest CSV for each workload and
comparison, using modification time; it never merges rounds from different files.
It writes this report plus `summary.csv` (output) and `tree-input-summary.csv`.
Pass explicit CSV paths to select a particular run; overlapping workloads are rejected.
The `indel` case is LG+GC{0.5} with indels; `indel_LG` is LG with indels.
