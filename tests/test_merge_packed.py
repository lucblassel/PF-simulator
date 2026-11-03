"""Index-only merger regression tests; placeholder shards need no simulator."""

import importlib.util
import json
import sqlite3
import tempfile
import unittest
from contextlib import closing
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "merge_packed", ROOT / "scripts/merge_packed.py"
)
merge_packed = importlib.util.module_from_spec(spec)
spec.loader.exec_module(merge_packed)


class MergeTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name).resolve()
        self.source = self.root / "batches" / "a"
        (self.source / "shards").mkdir(parents=True)
        (self.source / "shards/a.tar").write_bytes(bytes(4096))
        self.samples = merge_packed.INSERT_BATCH_SIZE + 3
        summary = {
            "status": "completed",
            "output_layout": {"format": "packed"},
            "batch_id": "a",
            "resolved_arguments": [],
            "written_samples": self.samples,
        }
        (self.source / "simulation_summary.json").write_text(json.dumps(summary))
        with closing(sqlite3.connect(self.source / "index.sqlite3")) as db:
            db.executescript((ROOT / "src/packed_schema.sql").read_text())
            db.executemany(
                "INSERT INTO pairs VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (
                    (
                        i,
                        f"sample_{i}",
                        "shards/a.tar",
                        f"trees/{i}.nwk",
                        512,
                        10,
                        f"msa/{i}.fa",
                        1536,
                        100,
                    )
                    for i in range(self.samples)
                ),
            )
            db.executemany(
                "INSERT INTO sample_metadata VALUES (?, ?, ?)",
                (
                    (f"sample_{i}", i, json.dumps({"candidate": i}))
                    for i in range(self.samples)
                ),
            )
            db.commit()

    def change(self, sql, parameters=()):
        with closing(sqlite3.connect(self.source / "index.sqlite3")) as db:
            db.execute(sql, parameters)
            db.commit()

    def assert_rejected(self, message):
        with self.assertRaisesRegex(ValueError, message):
            merge_packed.merge_indices([self.source], self.root)
        self.assertFalse((self.root / "index.sqlite3").exists())
        self.assertEqual(list(self.root.glob(".index-*.sqlite3")), [])

    def test_batched_rows_preserve_content_and_resolve_each_shard_once(self):
        original = Path.resolve
        resolved = []

        def resolve(path, *args, **kwargs):
            resolved.append(path)
            return original(path, *args, **kwargs)

        with patch.object(Path, "resolve", resolve):
            self.assertEqual(
                merge_packed.merge_indices([self.source], self.root), self.samples
            )
        self.assertEqual(resolved.count(self.source / "shards/a.tar"), 1)
        with closing(sqlite3.connect(self.root / "index.sqlite3")) as db:
            rows = list(
                db.execute(
                    "SELECT p.*, m.candidate_id, m.metadata_json, b.batch_id "
                    "FROM pairs p JOIN sample_metadata m USING(name) "
                    "JOIN sample_batches b USING(name) ORDER BY p.id"
                )
            )
            self.assertEqual(len(rows), self.samples)
            for i, row in enumerate(rows):
                self.assertEqual(
                    row,
                    (
                        i,
                        f"a/sample_{i}",
                        "batches/a/shards/a.tar",
                        f"trees/{i}.nwk",
                        512,
                        10,
                        f"msa/{i}.fa",
                        1536,
                        100,
                        i,
                        json.dumps({"candidate": i}),
                        "a",
                    ),
                )
            self.assertEqual(db.execute("PRAGMA integrity_check").fetchone(), ("ok",))
            self.assertEqual(list(db.execute("PRAGMA foreign_key_check")), [])
        with self.assertRaises(FileExistsError):
            merge_packed.merge_indices([self.source], self.root)

    def test_payload_bounds_checked_after_first_batch(self):
        self.change(
            "UPDATE pairs SET sequence_size=999999 WHERE id=?", (self.samples - 1,)
        )
        self.assert_rejected("Invalid payload bounds")

    def test_missing_metadata_does_not_publish_partial_index(self):
        self.change(
            "DELETE FROM sample_metadata WHERE candidate_id=?", (self.samples - 1,)
        )
        self.assert_rejected("Missing sample metadata")

    def test_unsafe_paths_rejected(self):
        for path in ("../outside.tar", "/absolute.tar"):
            with self.subTest(path=path):
                self.change("UPDATE pairs SET shard=?", (path,))
                self.assert_rejected("Unsafe shard path")

    def test_symlink_outside_dataset_rejected(self):
        with tempfile.TemporaryDirectory() as outside:
            shard = self.source / "shards/a.tar"
            shard.unlink()
            shard.symlink_to(Path(outside) / "outside.tar")
            self.assert_rejected("not in the subpath|is not a subpath")


if __name__ == "__main__":
    unittest.main()
