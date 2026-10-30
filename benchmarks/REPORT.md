# FASTA versus numeric output benchmark

Measured on 2026-10-30: Apple M2 Pro, 10 physical cores, 32 GiB RAM, macOS 27.0.1.
Rust 1.98.1; Rust and native code optimized in release mode. AliSim uses one thread
and the same portable likelihood kernel in both paths. No builds ran during the recorded trials.

## Conclusion

The current FASTA path is adequate for the existing LG + continuous-gamma workflow:
the numeric prototype saves approximately 3% with indels and offers no material gain
without them. Plain LG with indels is different: numeric output saves about 10% at
10 tips × 500 root sites, 31% at 50 × 500, and 81% at 200 × 5,000.
Rust FASTA duplicate parsing is a small component of those indel workloads; the
prototype also eliminates AliSim's native indel numeric-text storage and parsing.
A faster duplicate checker alone would not reproduce the large indel gains.

Decision after measurement: numeric checking is now the production default.
These timings were collected before promotion. Numeric output now has
regression comparisons against the FASTA
reference, including 200-tip trees. The benchmark feature enables the comparison
harness; ordinary builds use numeric checking.

## Method

- Nine paired rounds per workload. Backend order alternates each round.
- One fixed uniform tree per workload, tree seed 42; matching alignment seeds in each pair.
- Warm-up and semantic equivalence checks precede timing. Batch size is calibrated
  to about 200 ms per backend/round, bounded to 1–128 attempts; slow cases use one attempt.
- Times include native initialization, tree reading, simulation, transport, exact
  duplicate checking, buffer cleanup, and FASTA serialization of accepted numeric results.
  Indel cases also serialize accepted gap-free output. Rejected numeric results are not serialized.
- Tree generation, final file writes, terminal progress and optional logging are excluded.
  The numbers are compute time per attempt, not total CLI or durable-storage throughput.
- Indel rates: 0.002,0.002; size distributions: POW{1.7/100},POW{1.7/100}.
  Continuous gamma uses shape 0.5. Default mean branch length is 0.1 substitutions/site.
- Tip names and all aligned/gap-free residues match for the warm-up seed in every
  workload. Duplicate decisions match for every timed paired seed. Additional tests
  compare both transports across DNA/protein, gamma, indels and multiple seeds.
- IQR describes the middle half of the nine paired gains; it is not a confidence interval.
  Results are from one machine with no CPU pinning; sub-percent differences should
  not drive an architecture decision. Memory use was not measured.

## Normal branch lengths

Milliseconds per attempt are medians of round averages. Positive paired gain means numeric is faster.

| Tips × root sites | Model | FASTA ms | Numeric ms | Paired gain | Gain IQR | FASTA check share |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 10 × 500 | JC DNA | 4.877 | 4.861 | -0.1% | -0.7…+0.6% | 0.13% |
| 10 × 500 | LG | 4.955 | 4.941 | -0.3% | -0.8…+1.6% | 0.13% |
| 10 × 500 | LG + continuous gamma | 18.959 | 18.977 | -0.4% | -0.6…-0.0% | 0.04% |
| 10 × 500 | LG + continuous gamma + indels | 19.545 | 19.012 | +3.0% | +2.8…+3.0% | 0.04% |
| 10 × 500 | LG + indels | 5.588 | 5.062 | +10.2% | +8.8…+10.4% | 0.12% |
| 50 × 500 | JC DNA | 5.282 | 5.215 | +1.0% | +0.7…+1.4% | 0.59% |
| 50 × 500 | LG | 5.548 | 5.523 | +0.3% | -0.4…+0.5% | 0.56% |
| 50 × 500 | LG + continuous gamma | 83.777 | 83.741 | +0.1% | -0.2…+0.2% | 0.04% |
| 50 × 500 | LG + continuous gamma + indels | 87.779 | 84.465 | +2.8% | +2.7…+3.3% | 0.04% |
| 50 × 500 | LG + indels | 8.343 | 5.719 | +31.5% | +30.9…+31.8% | 0.41% |
| 200 × 5,000 | JC DNA | 18.914 | 19.058 | -1.2% | -1.3…-0.7% | 5.79% |
| 200 × 5,000 | LG | 21.075 | 21.198 | -0.7% | -0.9…-0.4% | 5.17% |
| 200 × 5,000 | LG + continuous gamma | 4130.770 | 4136.657 | -0.0% | -0.4…+0.0% | 0.03% |
| 200 × 5,000 | LG + continuous gamma + indels | 4297.968 | 4170.629 | +2.8% | +2.6…+2.9% | 0.04% |
| 200 × 5,000 | LG + indels | 145.072 | 27.833 | +80.9% | +80.7…+81.0% | 1.05% |

These ordinary workloads had no duplicate rejections among the timed seeds.

## Duplicate-heavy workloads

Each attempted MSA is timed once; this isolates per-attempt cost rather than running a retry loop.

| Tips × root sites | Model | Branch mean | Rejected | FASTA ms | Numeric ms | Paired gain |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 4 × 2 | JC DNA | 0.5 | 94.6% | 4.710 | 4.695 | -0.1% |
| 50 × 500 | LG | 1e-08 | 100.0% | 4.813 | 4.802 | +0.1% |
| 200 × 5,000 | LG | 1e-08 | 100.0% | 8.034 | 8.133 | -0.8% |

## Reproduce

```sh
cargo build --release --features raw-benchmark --example benchmark_output -j 8
./target/release/examples/benchmark_output --rounds 9 --target-ms 200
./target/release/examples/benchmark_output --cases JC --workloads 4x2 --branch-mean 0.5
./target/release/examples/benchmark_output --cases LG --workloads 50x500,200x5000 --branch-mean 1e-8
cargo test --release --features raw-benchmark
python3 benchmarks/summarize.py
```

The harness creates a fresh timestamped CSV by default. Pass `--output FILE` to select a new filename.
The summary script reads `benchmarks/output-comparison*.csv` by default; pass explicit CSV paths
to select one captured run. Its report metadata and conclusion describe the recorded M2 Pro run
and should be updated when collecting results on another machine or workload matrix.
Raw observations are in the `output-comparison*.csv` files; `summary.csv` is regenerated by the script.
The `indel` case means LG+GC{0.5} with indels; `indel_LG` means LG with indels.
