# Local SQLite and packed-output benchmark

Platform: macOS-27.0.1-arm64-arm-64bit-Mach-O.

Writer runs use prebuilt payloads and concurrent threads with independent worker connections. They isolate output handling and exclude AliSim/FASTA formatting. Simulation runs, when requested, use the simulator's actual process pool and the cluster benchmark's G8 options.

Run order is randomized with a fixed seed. All datasets use the same output filesystem. Data enumeration, integrity validation and cleanup are outside the internal timer. OS caching is not flushed. CPU times include validation and child startup. SQLite retains its normal durability settings; raw FASTA/tar writes are buffered. Worker-stage maxima are not additive wall-time components because workers overlap.

| Kind | Workload | Format | Workers | Variant | Median internal s | Min–max s | Accepted/s | Max worker init s* | Max worker finish s* | Merge s |
|---|---|---|---:|---|---:|---|---:|---:|---:|---:|
| simulation | large | files | 8 | buffer1m | 0.5095 | 0.4555–0.5230 | 125.6 | — | — | — |
| simulation | large | files | 8 | shard8k | 0.5601 | 0.4647–0.5671 | 114.3 | — | — | — |
| simulation | large | packed | 8 | buffer1m | 0.5122 | 0.4664–0.6531 | 125.0 | — | — | — |
| simulation | large | packed | 8 | shard8k | 0.5891 | 0.4979–0.7589 | 108.6 | — | — | — |
| simulation | small | files | 8 | buffer1m | 2.8631 | 2.6258–3.0874 | 476.8 | — | — | — |
| simulation | small | files | 8 | shard8k | 2.8043 | 2.7547–3.0270 | 496.7 | — | — | — |
| simulation | small | packed | 8 | buffer1m | 2.7033 | 2.5519–2.8393 | 520.8 | — | — | — |
| simulation | small | packed | 8 | shard8k | 2.5815 | 2.4830–2.8771 | 545.4 | — | — | — |
| writer | large | files | 8 | buffer1m | 0.0163 | 0.0160–0.0190 | 3923.1 | 0.0000 | 0.0000 | 0.0000 |
| writer | large | packed | 8 | buffer1m | 0.0307 | 0.0304–0.0324 | 2086.2 | 0.0037 | 0.0152 | 0.0031 |
| writer | large | packed | 8 | shard8k | 0.2121 | 0.2085–0.2889 | 301.7 | 0.0036 | 0.0034 | 0.0031 |
| writer | small | files | 8 | buffer1m | 0.3074 | 0.3032–0.3371 | 6662.9 | 0.0000 | 0.0000 | 0.0000 |
| writer | small | packed | 8 | buffer1m | 0.0727 | 0.0677–0.0873 | 28153.3 | 0.0042 | 0.0289 | 0.0346 |
| writer | small | packed | 8 | shard8k | 0.2975 | 0.2686–0.3091 | 6885.0 | 0.0046 | 0.0173 | 0.0353 |

*Median across runs of the slowest worker's stage duration. These figures diagnose setup/finalization overhead; they do not establish cluster filesystem performance.

Executable hashes and run settings are in `environment.json`; individual commands, logs, and timing/summary JSON files are retained in `runs.tar.gz`. `metrics.csv` contains every run.
