# Local SQLite and packed-output benchmark

Platform: macOS-27.0.1-arm64-arm-64bit-Mach-O.

Writer runs use prebuilt payloads and concurrent threads with independent worker connections. They isolate output handling and exclude AliSim/FASTA formatting. Simulation runs, when requested, use the simulator's actual process pool and the cluster benchmark's G8 options.

Run order is randomized with a fixed seed. All datasets use the same output filesystem. Data enumeration, integrity validation and cleanup are outside the internal timer. OS caching is not flushed. CPU times include validation and child startup. SQLite retains its normal durability settings; raw FASTA/tar writes are buffered. Worker-stage maxima are not additive wall-time components because workers overlap.

| Kind | Workload | Format | Workers | Variant | Median internal s | Min–max s | Accepted/s | Max worker init s* | Max worker finish s* | Merge s |
|---|---|---|---:|---|---:|---|---:|---:|---:|---:|
| simulation | large | files | 1 | legacy | 2.8176 | 2.5016–3.0204 | 22.7 | — | — | — |
| simulation | large | files | 1 | shard | 2.8195 | 2.8081–3.3940 | 22.7 | — | — | — |
| simulation | large | files | 4 | legacy | 0.7838 | 0.7226–0.8119 | 81.7 | — | — | — |
| simulation | large | files | 4 | shard | 0.7446 | 0.6617–0.7792 | 85.9 | — | — | — |
| simulation | large | files | 8 | legacy | 0.5505 | 0.4685–1.0452 | 116.3 | — | — | — |
| simulation | large | files | 8 | shard | 0.5375 | 0.5197–0.5813 | 119.1 | — | — | — |
| simulation | large | packed | 1 | legacy | 2.9214 | 2.7968–2.9574 | 21.9 | — | — | — |
| simulation | large | packed | 1 | shard | 2.7826 | 2.6526–2.9275 | 23.0 | — | — | — |
| simulation | large | packed | 4 | legacy | 0.7671 | 0.6994–0.8314 | 83.4 | — | — | — |
| simulation | large | packed | 4 | shard | 1.1143 | 0.6935–1.2524 | 57.4 | — | — | — |
| simulation | large | packed | 8 | legacy | 0.5806 | 0.4832–0.6328 | 110.2 | — | — | — |
| simulation | large | packed | 8 | shard | 0.5958 | 0.4980–0.8260 | 107.4 | — | — | — |
| simulation | small | files | 1 | legacy | 15.6027 | 15.4768–15.6059 | 90.0 | — | — | — |
| simulation | small | files | 1 | shard | 15.6192 | 15.4543–19.0357 | 88.3 | — | — | — |
| simulation | small | files | 4 | legacy | 3.9466 | 3.9263–6.1470 | 345.9 | — | — | — |
| simulation | small | files | 4 | shard | 4.0849 | 3.9076–4.1135 | 344.7 | — | — | — |
| simulation | small | files | 8 | legacy | 2.7887 | 2.7110–3.1841 | 504.9 | — | — | — |
| simulation | small | files | 8 | shard | 2.6278 | 2.5559–3.2134 | 534.1 | — | — | — |
| simulation | small | packed | 1 | legacy | 15.6250 | 15.4493–23.9609 | 87.4 | — | — | — |
| simulation | small | packed | 1 | shard | 15.2863 | 15.2237–15.6415 | 90.0 | — | — | — |
| simulation | small | packed | 4 | legacy | 4.0703 | 3.8686–4.9804 | 335.4 | — | — | — |
| simulation | small | packed | 4 | shard | 4.0554 | 3.9141–4.2292 | 343.5 | — | — | — |
| simulation | small | packed | 8 | legacy | 2.7655 | 2.5711–2.8859 | 493.6 | — | — | — |
| simulation | small | packed | 8 | shard | 2.5492 | 2.4063–2.6634 | 535.5 | — | — | — |
| writer | large | files | 1 | shard | 0.0329 | 0.0310–0.1268 | 1946.5 | 0.0000 | 0.0000 | 0.0000 |
| writer | large | files | 4 | shard | 0.0146 | 0.0143–0.0164 | 4397.7 | 0.0000 | 0.0000 | 0.0000 |
| writer | large | files | 8 | shard | 0.0226 | 0.0205–0.0245 | 2833.3 | 0.0000 | 0.0000 | 0.0000 |
| writer | large | packed | 1 | atomic | 0.1599 | 0.1445–0.1674 | 400.3 | 0.0008 | 0.0019 | 0.0017 |
| writer | large | packed | 1 | legacy | 0.1506 | 0.1356–0.1635 | 425.0 | 0.0027 | 0.0016 | 0.0036 |
| writer | large | packed | 1 | shard | 0.1956 | 0.1341–0.2117 | 327.2 | 0.0008 | 0.0017 | 0.0019 |
| writer | large | packed | 4 | atomic | 0.2040 | 0.1957–0.2045 | 313.7 | 0.0019 | 0.0033 | 0.0026 |
| writer | large | packed | 4 | legacy | 0.1932 | 0.1763–0.1942 | 331.3 | 0.0080 | 0.0034 | 0.0040 |
| writer | large | packed | 4 | shard | 0.1931 | 0.1867–0.2068 | 331.4 | 0.0051 | 0.0036 | 0.0024 |
| writer | large | packed | 8 | atomic | 0.3254 | 0.2639–0.3602 | 196.7 | 0.0044 | 0.0035 | 0.0029 |
| writer | large | packed | 8 | legacy | 0.2745 | 0.2450–0.3100 | 233.1 | 0.0082 | 0.0038 | 0.0053 |
| writer | large | packed | 8 | shard | 0.2260 | 0.2189–0.2379 | 283.2 | 0.0036 | 0.0042 | 0.0028 |
| writer | small | files | 1 | shard | 0.4595 | 0.4404–0.4736 | 4456.8 | 0.0000 | 0.0000 | 0.0000 |
| writer | small | files | 4 | shard | 0.2833 | 0.2685–0.3235 | 7229.8 | 0.0000 | 0.0000 | 0.0000 |
| writer | small | files | 8 | shard | 0.3030 | 0.2905–0.4288 | 6759.0 | 0.0000 | 0.0000 | 0.0000 |
| writer | small | packed | 1 | atomic | 0.2155 | 0.1923–0.3667 | 9502.6 | 0.0007 | 0.0021 | 0.0351 |
| writer | small | packed | 1 | legacy | 0.1621 | 0.1607–0.2400 | 12634.0 | 0.0021 | 0.0012 | 0.0296 |
| writer | small | packed | 1 | shard | 0.1872 | 0.1601–0.1928 | 10942.4 | 0.0007 | 0.0042 | 0.0296 |
| writer | small | packed | 4 | atomic | 0.2116 | 0.2059–0.2332 | 9679.6 | 0.0022 | 0.0124 | 0.0342 |
| writer | small | packed | 4 | legacy | 0.2013 | 0.1825–0.2189 | 10175.3 | 0.0048 | 0.0104 | 0.0370 |
| writer | small | packed | 4 | shard | 0.2043 | 0.2039–0.2394 | 10024.1 | 0.0017 | 0.0112 | 0.0336 |
| writer | small | packed | 8 | atomic | 0.3101 | 0.2261–0.3248 | 6604.7 | 0.0037 | 0.0177 | 0.0313 |
| writer | small | packed | 8 | legacy | 0.2811 | 0.2355–0.3100 | 7285.9 | 0.0150 | 0.0266 | 0.0366 |
| writer | small | packed | 8 | shard | 0.2479 | 0.2279–0.2576 | 8260.1 | 0.0037 | 0.0150 | 0.0389 |

*Median across runs of the slowest worker's stage duration. These figures diagnose setup/finalization overhead; they do not establish cluster filesystem performance.

Executable hashes and run settings are in `environment.json`; individual commands, logs, and timing/summary JSON files are retained in `runs.tar.gz`. `metrics.csv` contains every run.
