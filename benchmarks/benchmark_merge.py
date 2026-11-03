"""Measure index merging with synthetic rows and sparse placeholder shards.

This isolates the merger, without simulating data or creating usable training
archives. Inputs and generated indices are temporary. Example:

    python3 benchmarks/benchmark_merge.py --out benchmarks/merge-results \
        --merger before=/path/to/old/scripts/merge_packed.py \
        --merger after=scripts/merge_packed.py

Each script needs its adjacent ../src/packed_schema.sql. Defaults match twenty
100,000-sample batches. Allow several GB of temporary disk space. All versions
must produce identical logical SQLite contents; digest checks are outside timing.
"""

import argparse
import hashlib
import json
import platform
import sqlite3
import subprocess
import sys
import tempfile
from contextlib import closing
from pathlib import Path
from time import perf_counter

ROOT = Path(__file__).resolve().parents[1]


def fixture(root, batches, samples, metadata_bytes):
    sources = []
    metadata = json.dumps({"synthetic": True, "padding": "x" * metadata_bytes})
    for batch in range(batches):
        source = root / "batches" / f"batch-{batch:03}"
        (source / "shards").mkdir(parents=True)
        for shard in range((samples + 9999) // 10000):
            with (source / "shards" / f"{shard:04}.tar").open("wb") as output:
                output.truncate(10000 * 2048)
        summary = {
            "status": "completed",
            "batch_id": source.name,
            "output_layout": {"format": "packed"},
            "written_samples": samples,
            "resolved_arguments": ["--model=LG", "--gamma=G8"],
        }
        (source / "simulation_summary.json").write_text(json.dumps(summary))
        with closing(sqlite3.connect(source / "index.sqlite3")) as db:
            db.executescript((ROOT / "src/packed_schema.sql").read_text())
            db.executemany(
                "INSERT INTO pairs VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (
                    (
                        i,
                        f"sample_{i:06}",
                        f"shards/{i // 10000:04}.tar",
                        f"trees/worker-000-sample_{i:06}.nwk",
                        (i % 10000) * 2048 + 512,
                        100,
                        f"msa_untrimmed/worker-000-sample_{i:06}.fa",
                        (i % 10000) * 2048 + 1536,
                        500,
                    )
                    for i in range(samples)
                ),
            )
            db.executemany(
                "INSERT INTO sample_metadata VALUES (?, ?, ?)",
                ((f"sample_{i:06}", i, metadata) for i in range(samples)),
            )
            db.commit()
        sources.append(source)
    return sources


def digest_index(path):
    digest = hashlib.sha256()
    with closing(sqlite3.connect(path)) as db:
        for table, order in (
            ("pairs", "id"),
            ("sample_metadata", "name"),
            ("sample_batches", "name"),
            ("batches", "batch_id"),
            ("dataset_metadata", "key"),
        ):
            digest.update(table.encode())
            for row in db.execute(f"SELECT * FROM {table} ORDER BY {order}"):
                digest.update(json.dumps(row, separators=(",", ":")).encode())
                digest.update(b"\n")
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument(
        "--merger", action="append", required=True, metavar="LABEL=SCRIPT"
    )
    parser.add_argument("--batches", type=int, default=20)
    parser.add_argument("--samples-per-batch", type=int, default=100000)
    parser.add_argument("--metadata-bytes", type=int, default=1500)
    parser.add_argument("--rounds", type=int, default=1)
    args = parser.parse_args()
    if (
        min(args.batches, args.samples_per_batch, args.rounds) < 1
        or args.metadata_bytes < 0
    ):
        parser.error("Counts must be positive; metadata bytes must be nonnegative")
    args.out.mkdir(parents=True, exist_ok=False)
    mergers = [
        (label, Path(script).resolve())
        for label, script in (item.split("=", 1) for item in args.merger)
    ]
    environment = {
        "platform": platform.platform(),
        "python": sys.version,
        "batches": args.batches,
        "samples_per_batch": args.samples_per_batch,
        "metadata_padding_bytes": args.metadata_bytes,
        "note": "Synthetic indices; sparse placeholders are not training archives",
        "mergers": {
            label: {
                "path": str(script),
                "sha256": hashlib.sha256(script.read_bytes()).hexdigest(),
            }
            for label, script in mergers
        },
    }
    (args.out / "environment.json").write_text(json.dumps(environment, indent=2) + "\n")
    results = []
    expected = None
    with tempfile.TemporaryDirectory(prefix="pf-merge-benchmark-") as directory:
        root = Path(directory)
        print("Creating synthetic source indices…", flush=True)
        sources = fixture(
            root, args.batches, args.samples_per_batch, args.metadata_bytes
        )
        for repeat in range(args.rounds):
            for label, script in mergers[:: 1 if repeat % 2 == 0 else -1]:
                start = perf_counter()
                result = subprocess.run(
                    [
                        sys.executable,
                        str(script),
                        *map(str, sources),
                        "--out",
                        str(root),
                    ],
                    capture_output=True,
                    text=True,
                )
                elapsed = perf_counter() - start
                (args.out / f"{label}-{repeat}.log").write_text(
                    result.stdout + result.stderr
                )
                result.check_returncode()
                index = root / "index.sqlite3"
                record = {
                    "label": label,
                    "repeat": repeat,
                    "wall_seconds": elapsed,
                    "samples_per_second": args.batches
                    * args.samples_per_batch
                    / elapsed,
                    "index_bytes": index.stat().st_size,
                }
                print(json.dumps(record), flush=True)
                record["logical_sha256"] = digest_index(index)
                if expected is None:
                    expected = record["logical_sha256"]
                if record["logical_sha256"] != expected:
                    raise RuntimeError("Merged index contents differ")
                results.append(record)
                (args.out / "metrics.json").write_text(
                    json.dumps(results, indent=2) + "\n"
                )
                index.unlink()


if __name__ == "__main__":
    main()
