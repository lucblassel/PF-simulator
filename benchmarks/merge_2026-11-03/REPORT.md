# Packed index merge benchmark

Measured locally on macOS ARM64 with Python 3.15.1, using twenty synthetic
100,000-sample batches (2,000,000 examples). Each sample has approximately
1.5 KiB of JSON metadata. There are 200 sparse placeholder shards; these are
not usable training archives. Source generation and logical-content checks are
excluded from merge timings. No tree/MSA payloads are read or copied.

The previous script is from commit `99a5981239ced24129c12a74f4ef402a7897ff60`.
The updated script caches each validated shard path and size, and inserts
1,000 samples at a time with `executemany`. Both versions use a single SQLite
transaction and keep the same publication, policy, and bounds checks.

| Implementation | Wall time | Examples/s | Index size |
| --- | ---: | ---: | ---: |
| Previous | 120.55 s | 16,591 | 4.41 GiB |
| Updated | 35.87 s | 55,753 | 4.41 GiB |

The updated merger was **3.36× faster** (70.2% less wall time). The two indices
had identical logical contents and the same size (4,732,235,776 bytes).

Full timings are in `metrics.json`; environment and merger hashes are in
`environment.json`. One full-scale measurement per implementation, old first;
OS caches were not flushed. This is a local comparison, not a measurement of
Lustre performance or a predicted cluster speedup.

The benchmark compares SHA-256 digests of every row in all five tables in a
fixed order, including metadata and batch provenance. Fixture and merged index
files are deleted after verification.

Reproduce after checking out the old script and its schema into a separate
checkout:

```bash
python3 benchmarks/benchmark_merge.py --out benchmarks/merge-repeat \
  --merger before=/path/to/old-checkout/scripts/merge_packed.py \
  --merger after=scripts/merge_packed.py --rounds 3
```

The merger previously resolved each shard path for every sample, despite
caching its size. For this fixture, shard-path resolution falls from 2,000,000
to 200 calls. On shared storage, repeated path/symlink lookups can add metadata
latency; the exact contribution to the reported ten-minute cluster run has
not been profiled. Merging still reads and writes several GiB of index data,
updates SQLite indices, and duplicates the per-sample metadata in the combined
index. It is not just concatenating twenty small headers.

Validation: all 35 executable/merger tests pass. Added regression coverage for
one path resolution per shard, batch-boundary content preservation, rollback
on invalid offsets or missing metadata, unsafe paths, and external symlinks.
