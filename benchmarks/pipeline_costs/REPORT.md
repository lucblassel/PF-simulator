# Why the matched-pipeline speedup is smaller

The dominant extra cost is continuous Gamma inside AliSim. Window selection and ancestral export are much smaller in these measurements. Both backends pay the expensive native-engine cost, reducing the relative benefit of avoiding Python startup, tree-file input and FASTA parsing.

Single native thread, Python-prior trees, POW indels 0.002/0.002 throughout. Each cell is the median of three round means, in milliseconds per candidate. Each round uses 16 candidates at 50 tips × 500 sites, or four at 200 tips × 5,000 sites. Trees, native seeds and empirical alpha draws are paired across variants; Gamma changes the simulated data, as intended.

| Tips × root sites | Model | Ancestors | Linked engine | Gapped check | Gapped + window | Full FASTA formatting | Extra FASTA formatting | Write all three | Standalone AliSim |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 50 × 500 | LG + empirical GC | false | 72.422 | 0.021 | 0.051 | 0.022 | 0.051 | 0.390 | 84.041 |
| 50 × 500 | LG + empirical GC | true | 77.208 | 0.026 | 0.058 | 0.038 | 0.103 | 0.372 | 86.103 |
| 50 × 500 | LG | false | 9.720 | 0.029 | 0.051 | 0.035 | 0.067 | 0.611 | 37.491 |
| 50 × 500 | LG | true | 11.218 | 0.041 | 0.063 | 0.064 | 0.141 | 0.671 | 33.566 |
| 200 × 5000 | LG + empirical GC | false | 3579.644 | 0.621 | 1.828 | 0.697 | 2.000 | 1.486 | 3725.537 |
| 200 × 5000 | LG + empirical GC | true | 3679.011 | 0.867 | 2.092 | 1.474 | 4.000 | 2.464 | 3763.363 |
| 200 × 5000 | LG | false | 20.366 | 0.617 | 1.863 | 0.727 | 2.038 | 2.313 | 149.037 |
| 200 × 5000 | LG | true | 28.589 | 0.863 | 2.087 | 1.457 | 3.984 | 1.971 | 170.518 |

## Interpretation

- **50 tips:** continuous Gamma increases linked-engine time by 7.5×; standalone AliSim increases by 2.2×. With Gamma and ancestors, the combined gapped/window check costs 0.058 ms, versus 77.2 ms in the engine (0.076%).
- **200 tips:** continuous Gamma increases linked-engine time by 175.8×; standalone AliSim increases by 25.0×. With Gamma and ancestors, the combined gapped/window check costs 2.092 ms, versus 3679.0 ms in the engine (0.057%).

The earlier ~3× result used LG without Gamma and without indels. The later near-parity result also added continuous Gamma, ancestors, indels and materialized outputs. Those reports compare different workloads. This ablation holds indels fixed and isolates the large Gamma cost; it does not assign every historical timing difference to a single change.

AliSim already simulates ancestral states while propagating sequences down a tree. Exporting them adds retained rows, gap reconstruction, transport and output; it does not add an independent evolutionary simulation at each internal node.

The native implementation explicitly excludes continuous Gamma from its transition-matrix and substitution-rate caching paths (`dependencies/iqtree2/simulator/alisimulatorheterogeneity.cpp`). The timings demonstrate the engine cost; identifying which native subroutine accounts for it would require a native profiler.

## Method and limits

- Each candidate has exactly one attempt. No retries or replacement. All candidates are formatted and written even if selection fails, preventing rejection rates from changing the measured I/O workload. This deliberately differs from production, which skips rejected output.
- The gapped-only and gapped-plus-window checks run on the identical numeric rows. Each is repeated ten times and divided by ten for timing stability. The window column includes structural validation and the gapped check; it is not an additional cost to add to the gapped column. Rejected full duplicates exit before searching windows.
- The linked-engine column includes bridge setup and numeric transport, but excludes tree generation, selection and FASTA work. Full FASTA formatting is measured separately from formatting the trimmed and ungapped views. The write column covers all three byte buffers together; buffered overwrites of reused temporary files, no fsync. This isolates byte formatting/writing and does not reproduce production directory/metadata creation costs.
- Standalone uses the same generated tree (through Newick), exact model string, length and native seed. It includes native process startup, FASTA outputs and native auxiliary files. It is an engine/workflow control, not a new end-to-end timing of the Python script. Shared model/seed does not independently certify identical output files.
- Three independent master seeds, 20261102–20261104. Variant order and standalone/linked order alternate by round. No other benchmark runs concurrently. No warmup is subtracted. Small negative-looking ancestry differences reflect timing noise and should not be interpreted as ancestry making simulation faster.
- These short, serial phase measurements diagnose the bottleneck. The larger worker-scaling results remain in the [matched pipeline report](../python_pipeline/REPORT.md). No prior, acceptance rule or production behavior was changed.

## Reproduce

```bash
cargo build --release --example benchmark_pipeline_costs
target/release/examples/benchmark_pipeline_costs \
  --output benchmarks/pipeline_costs-new \
  --reference target/alisim-reference/build/iqtree/iqtree2
python3 benchmarks/summarize_pipeline_costs.py --input benchmarks/pipeline_costs-new
```

[Raw timings](runs.csv), [summary](summary.csv), [provenance](metadata.json).

## Discrete Gamma follow-up

A [paired G4/G8/GC benchmark](../gamma_costs/REPORT.md) confirms that G8 is much faster than GC for these workloads. At 200 tips × 5,000 sites with ancestors and indels, linked-engine times are **29 ms for LG, 34 ms for G4, 38 ms for G8 and 3,772 ms for GC**. G8 is about 100× faster than GC here. It approximates Gamma rates with eight categories and therefore changes the simulation prior relative to GC.
