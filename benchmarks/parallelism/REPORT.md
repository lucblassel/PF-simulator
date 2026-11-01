# Process-worker benchmark

Recorded 2026-11-01T20:00:15.319709+00:00 on macOS-27.0.1-arm64-arm-64bit-Mach-O, 10 logical CPUs; release binary. 5 rounds per case/job count.

Each measurement runs the complete CLI: process startup, tree generation, AliSim, duplicate retries, FASTA conversion, final file writes, and worker shutdown. One job uses the direct sequential path; higher counts use persistent worker processes. Logs are disabled. AliSim uses one native thread per worker.

## Batch wall times

| Workload | Samples | Tips × root sites | 1 job (s) | 2 jobs (s / speedup) | 4 jobs (s / speedup) | 8 jobs (s / speedup) |
| --- | ---: | --- | ---: | ---: | ---: | ---: |
| tiny | 8 | 10 × 500 | 0.048 | 0.031 / 1.56× | 0.021 / 2.28× | 0.020 / 2.55× |
| default | 512 | 10 × 500 | 2.690 | 1.364 / 1.98× | 0.705 / 3.84× | 0.445 / 6.19× |
| large | 128 | 200 × 5000 | 2.862 | 1.439 / 1.99× | 0.749 / 3.82× | 0.503 / 5.78× |
| python | 128 | 200 × 1000 | 2.249 | 1.180 / 1.90× | 0.617 / 3.65× | 0.424 / 4.84× |
| indels | 64 | 200 × 2000 | 1.004 | 0.507 / 1.99× | 0.264 / 3.80× | 0.183 / 5.47× |
| continuous_gamma | 16 | 50 × 500 | 1.372 | 0.717 / 1.91× | 0.376 / 3.65× | 0.243 / 5.52× |
| retries | 64 | 4 × 2 | 1.542 | 0.792 / 1.97× | 0.435 / 3.66× | 0.288 / 5.37× |

Speedups are medians of same-round time ratios against one job; displayed times are independent medians. Samples/second count requested samples, including duplicate failures.

## Paired speedup variability

| Workload | 2 jobs speedup IQR | 4 jobs speedup IQR | 8 jobs speedup IQR |
| --- | ---: | ---: | ---: |
| tiny | 1.54–1.58× | 2.28–2.28× | 2.49–2.58× |
| default | 1.95–1.99× | 3.72–3.90× | 5.42–6.48× |
| large | 1.98–2.01× | 3.76–4.02× | 5.69–6.28× |
| python | 1.83–1.91× | 3.37–3.77× | 4.83–5.69× |
| indels | 1.99–2.02× | 3.75–3.89× | 3.75–6.18× |
| continuous_gamma | 1.91–1.93× | 3.64–3.72× | 5.27–5.64× |
| retries | 1.96–2.00× | 3.64–3.76× | 5.24–5.62× |

## Workloads and correctness

All cases use LG amino acids, the default uniform/exponential tree prior, and duplicate checking unless overridden below. All job counts use the same sample IDs and seed in each round.

- `tiny`: eight default samples, to expose worker startup overhead.
- `default`: a larger batch of the default 10-tip, 500-site simulation.
- `large`: larger trees and MSAs, including their full output I/O.
- `python`: the full `--preset python` tree pipeline.
- `indels`: insertion/deletion rates 0.002/0.002 and AliSim's default indel sizes; both aligned and unaligned FASTA are written.
- `continuous_gamma`: model `LG+GC{0.5}`; substitution-rate sampling is expensive.
- `retries`: 4 tips × 2 DNA sites under JC, branch mean 0.5; at most 20 attempts per tree. Exhausted samples are skipped and counted.

All 140 measured runs completed. Dataset hashes, written counts, and failure counts were checked against the same-round one-job reference for every worker count. Hashes cover every relative file path and every output byte; logs and timing are excluded.

## Method and limits

- Five rounds by default, with rotating job-count and workload order; one unmeasured two-sample warmup per workload. See `metadata.json` for exact options and run settings.
- The same executable is used for every job count; no per-sample subprocess startup. At most one in-flight sample per worker, with small status messages sent to the parent.
- Runs are not concurrent with one another. Each run writes to a fresh directory and includes startup and shutdown. File counting, hashing and deletion are excluded from timing.
- Writes use normal buffered filesystem I/O without fsync. This measures CLI completion on a warm local filesystem, not durable storage throughput or remote/shared storage.
- CPU seconds in `runs.csv` come from `RUSAGE_CHILDREN`, including waited-for descendants. Aggregate peak RAM was not measured; each worker holds its own native simulation state.
- Small batches are particularly sensitive to scheduling and launch overhead. These results describe this machine and these workload sizes; benchmark representative production runs.

Reproduce from the repository root:

```sh
cargo build --release -j 8
python3 benchmarks/benchmark_parallelism.py --output benchmarks/parallelism-new
```

Raw results: [runs.csv](runs.csv), [summary.csv](summary.csv), [metadata.json](metadata.json).
