# Discrete versus continuous Gamma

The dominant extra cost is continuous Gamma inside AliSim. Window selection and ancestral export are much smaller in these measurements. Both backends pay the expensive native-engine cost, reducing the relative benefit of avoiding Python startup, tree-file input and FASTA parsing.

Single native thread, Python-prior trees, POW indels 0.002/0.002 throughout. Each cell is the median of three round means, in milliseconds per candidate. Each round uses 16 candidates at 50 tips × 500 sites, or four at 200 tips × 5,000 sites. Trees, native seeds and empirical alpha draws are paired across variants; Gamma changes the simulated data, as intended.

| Tips × root sites | Model | Ancestors | Linked engine | Gapped check | Gapped + window | Full FASTA formatting | Extra FASTA formatting | Write all three | Standalone AliSim |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 50 × 500 | LG + empirical G4 | false | 6.313 | 0.020 | 0.043 | 0.022 | 0.051 | 0.214 | 16.487 |
| 50 × 500 | LG + empirical G4 | true | 6.350 | 0.025 | 0.050 | 0.037 | 0.095 | 0.200 | 16.484 |
| 50 × 500 | LG + empirical G8 | false | 7.025 | 0.020 | 0.046 | 0.020 | 0.049 | 0.372 | 18.464 |
| 50 × 500 | LG + empirical G8 | true | 7.169 | 0.026 | 0.052 | 0.041 | 0.099 | 0.293 | 18.560 |
| 50 × 500 | LG + empirical GC | false | 64.329 | 0.019 | 0.043 | 0.023 | 0.050 | 0.341 | 81.999 |
| 50 × 500 | LG + empirical GC | true | 68.006 | 0.027 | 0.050 | 0.039 | 0.099 | 0.337 | 88.670 |
| 50 × 500 | LG | false | 5.605 | 0.018 | 0.046 | 0.020 | 0.049 | 0.221 | 15.622 |
| 50 × 500 | LG | true | 5.800 | 0.026 | 0.052 | 0.040 | 0.099 | 0.241 | 15.893 |
| 200 × 5000 | LG + empirical G4 | false | 24.533 | 0.608 | 1.857 | 0.673 | 1.942 | 3.438 | 164.151 |
| 200 × 5000 | LG + empirical G4 | true | 33.743 | 0.842 | 2.086 | 1.427 | 4.057 | 2.299 | 177.273 |
| 200 × 5000 | LG + empirical G8 | false | 28.170 | 0.613 | 1.830 | 0.684 | 2.023 | 2.844 | 171.830 |
| 200 × 5000 | LG + empirical G8 | true | 37.546 | 0.853 | 2.102 | 1.529 | 3.908 | 2.332 | 170.846 |
| 200 × 5000 | LG + empirical GC | false | 3616.232 | 0.623 | 1.868 | 0.797 | 2.098 | 1.041 | 3762.047 |
| 200 × 5000 | LG + empirical GC | true | 3772.347 | 0.874 | 2.092 | 1.805 | 4.186 | 4.808 | 3936.032 |
| 200 × 5000 | LG | false | 19.974 | 0.630 | 1.834 | 0.702 | 2.066 | 2.717 | 150.483 |
| 200 × 5000 | LG | true | 28.553 | 0.875 | 2.113 | 1.563 | 4.015 | 4.548 | 165.500 |

## Interpretation

- **50 tips:** continuous Gamma increases linked-engine time by 11.5×; standalone AliSim increases by 5.2×. With Gamma and ancestors, the combined gapped/window check costs 0.050 ms, versus 68.0 ms in the engine (0.073%).
- **200 tips:** continuous Gamma increases linked-engine time by 181.1×; standalone AliSim increases by 25.0×. With Gamma and ancestors, the combined gapped/window check costs 2.092 ms, versus 3772.3 ms in the engine (0.055%).

## Discrete Gamma follow-up

Discrete Gamma uses a finite set of rate categories; continuous Gamma gives sites continuously varying rates. The shape draw is paired across G4/G8/GC, but these are different site-rate distributions, so changing GC to G8 changes the simulation prior.

- **50 tips, ancestors exported:** G4 takes 6.35 ms in the linked engine versus 68.01 ms for GC (10.7× faster). Standalone takes 16.48 versus 88.67 ms.
- **50 tips, ancestors exported:** G8 takes 7.17 ms in the linked engine versus 68.01 ms for GC (9.5× faster). Standalone takes 18.56 versus 88.67 ms.
- **200 tips, ancestors exported:** G4 takes 33.74 ms in the linked engine versus 3772.35 ms for GC (111.8× faster). Standalone takes 177.27 versus 3936.03 ms.
- **200 tips, ancestors exported:** G8 takes 37.55 ms in the linked engine versus 3772.35 ms for GC (100.5× faster). Standalone takes 170.85 versus 3936.03 ms.

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
  --reference target/alisim-reference/build/iqtree/iqtree2 --gammas none,G4,G8,GC
python3 benchmarks/summarize_pipeline_costs.py --input benchmarks/pipeline_costs-new
```

[Raw timings](runs.csv), [summary](summary.csv), [provenance](metadata.json).
