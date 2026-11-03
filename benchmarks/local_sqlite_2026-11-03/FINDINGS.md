# Local packed-output optimization findings

216 measured runs on this ARM64 Mac: 144 schema/commit comparisons, 42 buffer comparisons, and 30 index-heavy comparisons. All ran sequentially with randomized configuration order. Three repeats per configuration, except five in the index-heavy experiment. All measurements retained normal SQLite durability settings.

## Changes retained

1. Start `BEGIN IMMEDIATE` before creating tables and setting the schema version. Schema creation and initial inserts now share a transaction.
2. Commit worker index rows after closing each tar shard, or when the worker finishes. Remove the 1,000-example periodic commit.
3. Give each active tar shard a 1 MiB buffer. This adds roughly 1 MiB of output buffering per worker and coalesces the tar library's streaming writes.

The index schema/version and biological simulation behavior are unchanged. Root-index merging remains single-threaded. No queue or asynchronous SQLite writer was introduced.

## Why the buffer matters

The tar crate forwards data through `dyn Read`/`dyn Write` and `io::copy`, which streams a large byte slice in small chunks. A persistent file handle therefore did not imply large writes. A deterministic transport test recorded **319 underlying writer calls with an 8 KiB buffer versus 3 with a 1 MiB buffer** for a 2,600,000-byte entry, with equal total archive bytes. This counts calls at the underlying `Write` boundary; it is not a kernel syscall trace.

The following comparison holds schema/commit policy constant and changes only the output buffer. Payloads are prebuilt; this isolates writing and indexing rather than AliSim, selection or FASTA formatting.

| Writer workload, 8 workers | 8 KiB buffer, median s | 1 MiB buffer, median s | Ratio of median times |
|---|---:|---:|---:|
| 2,048 × 65,536-byte payloads | 0.2975 | 0.0727 | 4.09× |
| 64 × 2,600,000-byte payloads | 0.2121 | 0.0307 | 6.92× |

## SQLite policy experiment

The index-heavy case uses 8,192 examples with only 256 payload bytes each and approximately 1.5 KiB of metadata. It crosses the old 1,000-row interval even at eight workers. All three variants retain the old 8 KiB tar buffer, so buffer changes do not explain these differences. Five repeats per configuration.

| Workers | Legacy setup + 1,000-row commits, s | Atomic setup + 1,000-row commits, s | Atomic setup + shard commits, s | Time reduction, legacy to shard |
|---:|---:|---:|---:|---:|
| 1 | 0.2009 | 0.1990 | 0.1901 | 5.4% |
| 8 | 0.2237 | 0.2151 | 0.2098 | 6.2% |

For the original 2,048/64-example workloads, SQLite-only changes did not consistently improve total writing or end-to-end time. Atomic schema setup did reduce measured initialization time; at eight workers in the small writer case, the median slowest-worker initialization fell from about 15 ms to 3.7 ms. That is a small component of total local runtime. Most worker configurations in that grid never reached 1,000 rows per worker, so reducing periodic commit frequency cannot help them.

## Full simulator check

The buffer comparison also ran the real process pool with G8, empirical gamma shapes, ancestors, indels and duplicate/window checks, using the same three seeds. Accepted counts matched across builds, formats and worker counts.

| Workload, 8 workers | 8 KiB packed, median s | 1 MiB packed, median s | 8 KiB build files control, s | 1 MiB build files control, s |
|---|---:|---:|---:|---:|
| 2,048 candidates, 50 tips × 500 sites | 2.5815 | 2.7033 | 2.8043 | 2.8631 |
| 64 candidates, 200 tips × 5,000 sites | 0.5891 | 0.5122 | 0.5601 | 0.5095 |

The large packed case improved from about 0.59 to 0.51 s, but the unchanged file-mode control also varied (0.56 to 0.51 s). The small packed case did not improve. These short, three-repeat measurements do not establish a reliable end-to-end speedup. The isolated writer improvement is much clearer and should be rechecked on the cluster filesystem.

## Why no SQLite thread/queue yet

A queue would not by itself reduce database setup, transaction count or index finalization. A central coordinator index could eliminate per-worker index files and their merge, but would require sending index records across worker-process IPC and coordinating failures and publication. The measured local SQLite costs do not currently justify that complexity. The simple transaction changes and larger output buffer preserve worker ownership of shards.

## Method and limits

- Writer-only measurements use threads with independent production `PackedWriter` instances; full simulator measurements use its process pool. They answer different questions.
- Worker initialization/write/finalization durations and coordinator merge time are recorded separately. Maxima across workers are not additive wall-time components.
- Writes are OS-buffered. Integrity checks, output enumeration and cleanup happen after the internal timer. SQLite durability settings were not relaxed; normal FASTA/tar writes do not force durability. No OS cache flush or exclusive-machine guarantee.
- The toy payloads do not simulate phylogenetic data or reproduce its FASTA-formatting cost. Metadata/payload sizes and run settings are recorded.
- One transaction per shard can commit more often than every 1,000 rows if very small shards are requested. A crash can lose the current uncommitted shard index; incomplete runs still do not publish the completed root index. Neither the old nor new pipeline promises power-loss durability for raw tar data.
- Timings on local macOS storage do not quantify SQLite synchronization or metadata latency on Lustre. Three repeats use min–max ranges, not confidence intervals.

## Reproduce

See `benchmarks/benchmark_local_sqlite.py --help` and the README's local benchmark section. Pass saved before/after executables with repeated `--writer NAME=PATH` and optional `--simulator NAME=PATH` arguments. The script records executable hashes and full commands; binary paths for historical variants may point to temporary directories. Historical packed-writer/schema source snapshots are in `sources/`.

Detailed results:

- [Schema/commit grid](REPORT.md), [raw metrics](metrics.csv).
- [Buffer comparison](../local_sqlite_buffer_2026-11-03/REPORT.md), [raw metrics](../local_sqlite_buffer_2026-11-03/metrics.csv).
- [Index-heavy case](../local_sqlite_index_2026-11-03/REPORT.md), [raw metrics](../local_sqlite_index_2026-11-03/metrics.csv).

Validation: all debug tests/doctests passed for transaction changes; release tests across all targets passed after the buffer change, including the Python executable contract and packed-offset/rotation/merge tests. Clippy with `--all-targets -- -D warnings`, Ruff and formatting checks passed.
