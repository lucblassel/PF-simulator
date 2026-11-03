"""Merge completed simulation indices without copying or rewriting tar shards.

Sources must live below the destination directory, making all shard references
relative and the resulting dataset relocatable. The complete directory tree must
be retained when moving the dataset. Requires Python 3.10+; no extra packages.
"""

from __future__ import annotations

import argparse
from contextlib import closing
import json
import os
from pathlib import Path, PurePosixPath
import sqlite3
import tempfile
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


def policy(summary: dict[str, Any]) -> tuple[str, ...]:
    return tuple(
        argument
        for argument in summary["resolved_arguments"]
        if argument.split("=", 1)[0] not in IGNORED_ARGUMENTS
    )


def merge_indices(
    sources: list[Path], destination: Path, *, allow_mixed_configs: bool = False
) -> int:
    """Publish a new dense index; preserve source batch identity and provenance."""
    destination = destination.resolve()
    destination.mkdir(parents=True, exist_ok=True)
    index = destination / "index.sqlite3"
    if index.exists():
        raise FileExistsError(index)
    if len(set(source.resolve() for source in sources)) != len(sources):
        raise ValueError("Repeated source directory")
    descriptor, temporary = tempfile.mkstemp(
        prefix=".index-", suffix=".sqlite3", dir=destination
    )
    os.close(descriptor)
    staging = Path(temporary)
    count = 0
    expected_policy: tuple[str, ...] | None = None
    identities: set[str] = set()
    try:
        with closing(sqlite3.connect(staging)) as target:
            schema = Path(__file__).resolve().parents[1] / "src" / "packed_schema.sql"
            target.executescript(schema.read_text())
            target.execute(
                "CREATE TABLE batches (batch_id TEXT PRIMARY KEY, path TEXT NOT NULL, summary_json TEXT NOT NULL)"
            )
            target.execute(
                "CREATE TABLE sample_batches (name TEXT PRIMARY KEY REFERENCES pairs(name), batch_id TEXT NOT NULL REFERENCES batches(batch_id))"
            )
            for source in sorted(
                (p.resolve() for p in sources), key=lambda p: p.as_posix()
            ):
                source.relative_to(
                    destination
                )  # Reject external/absolute shard dependencies.
                summary = json.loads((source / "simulation_summary.json").read_text())
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
                target.execute(
                    "INSERT INTO batches VALUES (?, ?, ?)",
                    (
                        identity,
                        source.relative_to(destination).as_posix(),
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
                    sizes: dict[str, int] = {}
                    rows = db.execute(
                        "SELECT p.*, m.candidate_id, m.metadata_json FROM pairs p JOIN sample_metadata m USING(name) ORDER BY p.id"
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
                        ) = row
                        relative = PurePosixPath(shard)
                        if relative.is_absolute() or ".." in relative.parts:
                            raise ValueError(f"Unsafe shard path: {shard}")
                        path = (source / shard).resolve()
                        shard_path = path.relative_to(destination).as_posix()
                        if shard not in sizes:
                            sizes[shard] = path.stat().st_size
                        for offset, size in [
                            (tree_offset, tree_size),
                            (sequence_offset, sequence_size),
                        ]:
                            if offset < 512 or size < 0 or offset + size > sizes[shard]:
                                raise ValueError(
                                    f"Invalid payload bounds in {source}: {name}"
                                )
                        qualified = f"{identity}/{name}"
                        target.execute(
                            "INSERT INTO pairs VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
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
                            ),
                        )
                        target.execute(
                            "INSERT INTO sample_metadata VALUES (?, ?, ?)",
                            (qualified, candidate, metadata),
                        )
                        target.execute(
                            "INSERT INTO sample_batches VALUES (?, ?)",
                            (qualified, identity),
                        )
                        count += 1
                        seen += 1
                    if seen != source_count:
                        raise ValueError(f"Missing sample metadata: {source}")
            target.commit()
        # Exclusive publication prevents accidentally replacing an existing index.
        os.link(staging, index)
    finally:
        staging.unlink(missing_ok=True)
    return count


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sources", nargs="+", type=Path)
    parser.add_argument(
        "--out",
        required=True,
        type=Path,
        help="Common parent directory; receives index.sqlite3",
    )
    parser.add_argument(
        "--allow-mixed-configs",
        action="store_true",
        help="Allow an intentional mixture of simulation policies",
    )
    args = parser.parse_args()
    count = merge_indices(
        args.sources, args.out, allow_mixed_configs=args.allow_mixed_configs
    )
    print(
        f"Indexed {count} examples in {args.out / 'index.sqlite3'}; tar shards were not modified."
    )


if __name__ == "__main__":
    main()
