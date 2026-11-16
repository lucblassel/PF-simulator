"""Merge completed simulations into one index and a flat shards/ directory.

Shards are hard linked, so sources and destination must share a filesystem.
Source datasets are retained. Requires Python 3.10+; no extra packages.
"""

from __future__ import annotations

import argparse
import errno
import json
import os
import sqlite3
import sys
import tempfile
from contextlib import closing
from io import BytesIO
from pathlib import Path, PurePosixPath
from time import monotonic
from typing import Any

IGNORED_ARGUMENTS = {
    "--samples",
    "--seed",
    "--jobs",
    "--output-format",
    "--shard-size-bytes",
    "--shard-buffer-bytes",
    "--batch-id",
    "--file-prefix",
}

INSERT_BATCH_SIZE = 1_000
SUMMARY_NAME = "simulation_summary.json"


def insert_samples(target: sqlite3.Connection, samples: list[tuple]) -> None:
    """Insert a bounded batch while retaining the enclosing merge transaction."""
    target.executemany(
        "INSERT INTO pairs VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
        (sample[:9] for sample in samples),
    )
    target.executemany(
        "INSERT INTO sample_metadata VALUES (?, ?, ?)",
        ((sample[1], sample[9], sample[10]) for sample in samples),
    )
    target.executemany(
        "INSERT INTO sample_batches VALUES (?, ?)",
        ((sample[1], sample[11]) for sample in samples),
    )
    target.executemany(
        "INSERT INTO sequence_lengths VALUES (?, ?)",
        ((sample[0], sample[12]) for sample in samples),
    )
    samples.clear()


def policy(summary: dict[str, Any]) -> tuple[str, ...]:
    return tuple(
        argument
        for argument in summary["resolved_arguments"]
        if argument.split("=", 1)[0] not in IGNORED_ARGUMENTS
    )


def legacy_length(metadata: str, shard: Path, offset: int, size: int, name: str) -> int:
    # Older simulator outputs recorded the full gapped MSA width in JSON.
    columns = json.loads(metadata).get("alignment_columns")
    if type(columns) is int and 0 <= columns <= size:
        return columns

    with shard.open("rb") as handle:
        handle.seek(offset)
        payload = handle.read(size)
    if len(payload) != size:
        raise ValueError(f"Truncated sequence payload: {name}")
    maximum = length = 0
    found_header = False
    for line in BytesIO(payload):
        line = line.strip()
        if line.startswith(b">"):
            maximum = max(maximum, length)
            length = 0
            found_header = True
        elif found_header:
            length += len(line)
        else:
            raise ValueError(f"Sequence data before first FASTA header in {name}")
    if not found_header:
        raise ValueError(f"No FASTA header in {name}")
    return max(maximum, length)


def merge_indices(
    sources: list[Path],
    destination: Path,
    *,
    allow_mixed_configs: bool = False,
    show_progress: bool = False,
) -> int:
    """Publish a new dense index; preserve source batch identity and provenance."""
    destination = destination.resolve()
    destination.mkdir(parents=True, exist_ok=True)
    index = destination / "index.sqlite3"
    summary_path = destination / SUMMARY_NAME
    for output in (index, summary_path):
        if output.exists() or output.is_symlink():
            raise FileExistsError(output)
    resolved_sources = sorted(source.resolve() for source in sources)
    if len(set(resolved_sources)) != len(sources):
        raise ValueError("Repeated source directory")
    expected_policy: tuple[str, ...] | None = None
    identities: set[str] = set()
    summaries: list[dict[str, Any]] = []
    total = 0
    for source in resolved_sources:
        summary = json.loads((source / SUMMARY_NAME).read_text(encoding="utf-8"))
        if (
            summary.get("status") != "completed"
            or summary.get("output_layout", {}).get("format") != "packed"
        ):
            raise ValueError(f"Not a completed packed simulation: {source}")
        identity = summary["batch_id"]
        if identity in identities:
            raise ValueError(f"Repeated batch identity: {identity}")
        identities.add(identity)
        current_policy = policy(summary)
        if (
            expected_policy is not None
            and current_policy != expected_policy
            and not allow_mixed_configs
        ):
            raise ValueError(
                "Simulation policies differ; use --allow-mixed-configs for an intentional mixture"
            )
        expected_policy = current_policy
        written = summary.get("written_samples")
        if type(written) is not int or written < 0:
            raise ValueError(f"Invalid summary sample count: {source}")
        total += written
        summaries.append(summary)

    shard_dir = destination / "shards"
    shard_dir.mkdir()  # Exclusive ownership; never remove preexisting shards.
    staging: Path | None = None
    published = False
    summary_published = False
    next_shard = 0
    count = 0
    interactive = show_progress and sys.stderr.isatty()
    last_update = 0.0
    last_bucket = -1
    line_width = 0

    def report(phase: str | None = None) -> None:
        nonlocal last_update, last_bucket, line_width
        if not show_progress:
            return
        percent = min(99, count * 100 // total) if total else 0
        bucket = percent // 10
        now = monotonic()
        if phase is None:
            if interactive and now - last_update < 0.2:
                return
            if not interactive and bucket <= last_bucket:
                return
        last_update, last_bucket = now, bucket
        if phase == "finalizing":
            label = "Finalizing index"
        else:
            label = f"Merging: {100 if phase == 'done' else percent:3d}%"
        line = f"{label} ({count:,}/{total:,} samples)"
        if interactive:
            print(
                "\r" + line.ljust(line_width),
                end="\n" if phase == "done" else "",
                file=sys.stderr,
                flush=True,
            )
            line_width = len(line) if phase != "done" else 0
        else:
            print(line, file=sys.stderr, flush=True)

    try:
        report("start")
        descriptor, temporary = tempfile.mkstemp(
            prefix=".index-", suffix=".sqlite3", dir=destination
        )
        staging = Path(temporary)
        os.close(descriptor)
        with closing(sqlite3.connect(staging)) as target:
            schema = Path(__file__).resolve().parents[1] / "src" / "packed_schema.sql"
            target.executescript(schema.read_text())
            target.execute(
                "CREATE TABLE batches (batch_id TEXT PRIMARY KEY, path TEXT NOT NULL, summary_json TEXT NOT NULL)"
            )
            target.execute(
                "CREATE TABLE sample_batches (name TEXT PRIMARY KEY REFERENCES pairs(name), batch_id TEXT NOT NULL REFERENCES batches(batch_id))"
            )
            for source, summary in zip(resolved_sources, summaries, strict=True):
                identity = summary["batch_id"]
                target.execute(
                    "INSERT INTO batches VALUES (?, ?, ?)",
                    (
                        identity,
                        Path(os.path.relpath(source, destination)).as_posix(),
                        json.dumps(summary),
                    ),
                )
                uri = (source / "index.sqlite3").as_uri() + "?mode=ro&immutable=1"
                with closing(sqlite3.connect(uri, uri=True)) as db:
                    if db.execute("PRAGMA user_version").fetchone()[0] != 3:
                        raise ValueError(f"Unsupported packed index version: {source}")
                    source_count = db.execute("SELECT COUNT(*) FROM pairs").fetchone()[
                        0
                    ]
                    if source_count != summary["written_samples"]:
                        raise ValueError(
                            f"Index/summary sample count mismatch: {source}"
                        )
                    seen = 0
                    # Validate each shard once; keep memory proportional to shard count.
                    shards: dict[str, tuple[str, int]] = {}
                    samples: list[tuple] = []
                    has_lengths = (
                        db.execute(
                            "SELECT 1 FROM sqlite_master WHERE type='table' AND name='sequence_lengths'"
                        ).fetchone()
                        is not None
                    )
                    length_column = "l.max_sequence_length" if has_lengths else "NULL"
                    length_join = (
                        "LEFT JOIN sequence_lengths l ON l.pair_id = p.id"
                        if has_lengths
                        else ""
                    )
                    rows = db.execute(
                        f"SELECT p.*, m.candidate_id, m.metadata_json, {length_column} "
                        f"FROM pairs p JOIN sample_metadata m USING(name) {length_join} ORDER BY p.id"
                    )
                    for row in rows:
                        (
                            _,
                            name,
                            shard,
                            tree_name,
                            tree_offset,
                            tree_size,
                            sequence_name,
                            sequence_offset,
                            sequence_size,
                            candidate,
                            metadata,
                            length,
                        ) = row
                        if shard not in shards:
                            relative = PurePosixPath(shard)
                            if relative.is_absolute() or ".." in relative.parts:
                                raise ValueError(f"Unsafe shard path: {shard}")
                            path = (source / shard).resolve()
                            if not path.is_relative_to(source):
                                raise ValueError(
                                    f"Shard path escapes source directory: {path}"
                                )
                            shard_name = f"shard-{next_shard:06d}.tar"
                            linked = shard_dir / shard_name
                            try:
                                os.link(path, linked)
                            except OSError as error:
                                if error.errno == errno.EXDEV:
                                    raise OSError(
                                        f"Cannot link {path}: source and output must share a filesystem"
                                    ) from error
                                raise
                            shards[shard] = (
                                f"shards/{shard_name}",
                                linked.stat().st_size,
                            )
                            next_shard += 1
                        shard_path, shard_size = shards[shard]
                        for offset, size in [
                            (tree_offset, tree_size),
                            (sequence_offset, sequence_size),
                        ]:
                            if offset < 512 or size < 0 or offset + size > shard_size:
                                raise ValueError(
                                    f"Invalid payload bounds in {source}: {name}"
                                )
                        qualified = f"{identity}/{name}"
                        if not has_lengths:
                            length = legacy_length(
                                metadata,
                                destination / shard_path,
                                sequence_offset,
                                sequence_size,
                                qualified,
                            )
                        if type(length) is not int or not 0 <= length <= sequence_size:
                            raise ValueError(
                                f"Missing or invalid sequence length: {qualified}"
                            )
                        samples.append(
                            (
                                count,
                                qualified,
                                shard_path,
                                tree_name,
                                tree_offset,
                                tree_size,
                                sequence_name,
                                sequence_offset,
                                sequence_size,
                                candidate,
                                metadata,
                                identity,
                                length,
                            ),
                        )
                        count += 1
                        seen += 1
                        if len(samples) >= INSERT_BATCH_SIZE:
                            insert_samples(target, samples)
                            if count < total:
                                report()
                    if seen != source_count:
                        raise ValueError(f"Missing sample metadata: {source}")
                    insert_samples(target, samples)
                    if count < total:
                        report()
            report("finalizing")
            target.execute(
                "INSERT INTO sequence_length_stats "
                "SELECT COALESCE(MAX(max_sequence_length), 0) FROM sequence_lengths"
            )
            target.commit()
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            prefix=".summary-",
            suffix=".json",
            dir=destination,
        ) as output:
            json.dump(
                {"written_samples": count, "batches": summaries}, output, indent=2
            )
            output.write("\n")
            output.flush()
            os.link(output.name, summary_path)
            summary_published = True
        # Exclusive publication prevents accidentally replacing an existing index.
        os.link(staging, index)
        published = True
        report("done")
    finally:
        if interactive and line_width:
            print(file=sys.stderr, flush=True)
        if staging is not None:
            staging.unlink(missing_ok=True)
        if not published:
            if summary_published:
                summary_path.unlink()
            for shard in shard_dir.iterdir():
                shard.unlink()
            shard_dir.rmdir()
    return count


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sources", nargs="+", type=Path)
    parser.add_argument(
        "--out",
        required=True,
        type=Path,
        help="Output directory for the merged index, summary and shards (same filesystem as sources)",
    )
    parser.add_argument(
        "--allow-mixed-configs",
        action="store_true",
        help="Allow an intentional mixture of simulation policies",
    )
    args = parser.parse_args()
    start = monotonic()
    count = merge_indices(
        args.sources,
        args.out,
        allow_mixed_configs=args.allow_mixed_configs,
        show_progress=True,
    )
    elapsed = monotonic() - start
    print(
        f"Indexed {count} examples in {args.out / 'index.sqlite3'} in {elapsed}sec; tar shards were not modified."
    )


if __name__ == "__main__":
    main()
