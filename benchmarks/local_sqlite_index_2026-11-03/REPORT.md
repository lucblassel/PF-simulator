# Local SQLite and packed-output benchmark

Platform: macOS-27.0.1-arm64-arm-64bit-Mach-O.

Writer runs use prebuilt payloads and concurrent threads with independent worker connections. They isolate output handling and exclude AliSim/FASTA formatting. Simulation runs, when requested, use the simulator's actual process pool and the cluster benchmark's G8 options.

Run order is randomized with a fixed seed. All datasets use the same output filesystem. Data enumeration, integrity validation and cleanup are outside the internal timer. OS caching is not flushed. CPU times include validation and child startup. SQLite retains its normal durability settings; raw FASTA/tar writes are buffered. Worker-stage maxima are not additive wall-time components because workers overlap.

| Kind | Workload | Format | Workers | Variant | Median internal s | Min–max s | Accepted/s | Max worker init s* | Max worker finish s* | Merge s |
|---|---|---|---:|---|---:|---|---:|---:|---:|---:|
| writer | small | packed | 1 | atomic | 0.1990 | 0.1940–0.2023 | 41166.7 | 0.0007 | 0.0023 | 0.1071 |
| writer | small | packed | 1 | legacy | 0.2009 | 0.1988–0.2028 | 40775.1 | 0.0023 | 0.0024 | 0.1084 |
| writer | small | packed | 1 | shard | 0.1901 | 0.1895–0.1992 | 43100.0 | 0.0007 | 0.0063 | 0.1081 |
| writer | small | packed | 8 | atomic | 0.2151 | 0.2091–0.2195 | 38092.5 | 0.0037 | 0.0041 | 0.1137 |
| writer | small | packed | 8 | legacy | 0.2237 | 0.2212–0.2260 | 36614.3 | 0.0098 | 0.0040 | 0.1150 |
| writer | small | packed | 8 | shard | 0.2098 | 0.2091–0.2214 | 39038.2 | 0.0038 | 0.0443 | 0.1112 |

*Median across runs of the slowest worker's stage duration. These figures diagnose setup/finalization overhead; they do not establish cluster filesystem performance.

Executable hashes and run settings are in `environment.json`; individual commands, logs, and timing/summary JSON files are retained in `runs.tar.gz`. `metrics.csv` contains every run.
