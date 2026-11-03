"""Exercise the executable and check tree/sequence consistency without dependencies."""

from __future__ import annotations

import math
from contextlib import closing
import sqlite3
import tarfile
import json
import os
import time
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest
from test_merge_packed import MergeTests  # noqa: F401 -- included by unittest.main


EXECUTABLE = str(Path(sys.argv.pop(1)).resolve())


def run(output: Path, *options: str, success: bool = True) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(
        [EXECUTABLE, "--out", str(output), *options], capture_output=True, text=True, timeout=30
    )
    if success and result.returncode != 0:
        raise AssertionError(f"Exit code {result.returncode}\n" + result.stdout + result.stderr)
    if not success and result.returncode == 0:
        raise AssertionError("Expected simulation to fail")
    return result


def fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    name = ""
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].strip()
            if name in sequences:
                raise AssertionError(f"Duplicate FASTA ID: {name}")
            sequences[name] = ""
        else:
            sequences[name] += line.strip()
    return sequences


def tree_paths(path: Path) -> dict[str, dict[int, float]]:
    """Parse our Newick subset into root paths keyed by unique edge IDs."""
    tokens = iter(re.findall(r"T\d+|[(),:;]|[0-9.eE+\-]+", path.read_text()))
    current = next(tokens)
    paths: dict[str, dict[int, float]] = {}
    edge_id = 0

    def advance() -> None:
        nonlocal current
        current = next(tokens, "")

    def subtree() -> list[str]:
        nonlocal edge_id
        tips: list[str] = []
        if current == "(":
            advance()
            tips.extend(subtree())
            while current == ",":
                advance()
                tips.extend(subtree())
            assert current == ")"
            advance()
        else:
            assert current.startswith("T") and current not in paths
            tips.append(current)
            paths[current] = {}
            advance()

        if current == ":":
            advance()
            length = float(current)
            assert math.isfinite(length) and length > 0
            for tip in tips:
                paths[tip][edge_id] = length
            edge_id += 1
            advance()
        return tips

    subtree()
    assert current == ";" and next(tokens, None) is None
    assert edge_id == 2 * len(paths) - 3
    return paths


def packed_snapshot(root: Path) -> dict[str, tuple[bytes, bytes, dict]]:
    """Independently check tar headers and direct-offset payloads, including metadata."""
    snapshot = {}
    members = {}
    for path in root.rglob("*.tar"):
        with tarfile.open(path) as archive:
            assert path.stat().st_size % 512 == 0
            with path.open("rb") as stream:
                stream.seek(-1024, 2)
                assert stream.read() == b"\0" * 1024
            entries = {}
            for member in archive:
                assert member.isfile()
                entries[member.name] = (member.offset_data, member.size, archive.extractfile(member).read())
            members[path.relative_to(root).as_posix()] = entries
    assert not list(root.rglob("*.partial"))
    with closing(sqlite3.connect(root / "index.sqlite3")) as db:
        assert db.execute("PRAGMA user_version").fetchone()[0] == 3
        rows = list(db.execute("SELECT * FROM pairs ORDER BY id"))
        assert [row[0] for row in rows] == list(range(len(rows)))
        for _, name, shard, tree_name, tree_offset, tree_size, sequence_name, sequence_offset, sequence_size in rows:
            with (root / shard).open("rb") as handle:
                handle.seek(tree_offset)
                tree = handle.read(tree_size)
                handle.seek(sequence_offset)
                sequence = handle.read(sequence_size)
            assert members[shard].pop(tree_name) == (tree_offset, tree_size, tree)
            assert members[shard].pop(sequence_name) == (sequence_offset, sequence_size, sequence)
            candidate, metadata = db.execute("SELECT candidate_id, metadata_json FROM sample_metadata WHERE name=?", (name,)).fetchone()
            metadata = json.loads(metadata)
            worker = re.search(r"worker-(\d+)-shard-\d+\.tar$", Path(shard).name).group(1)
            prefix = metadata["sample"].removesuffix(f"sample_{candidate:06}")
            assert tree_name == f"trees/{prefix}worker-{worker}-sample_{candidate:06}.nwk"
            assert sequence_name == f"sequences/{prefix}worker-{worker}-sample_{candidate:06}.fa"
            assert candidate == metadata["candidate_id"]
            snapshot[name] = (tree, sequence, metadata)
        assert all(not entries for entries in members.values())
    return snapshot


class SimulationTests(unittest.TestCase):
    def test_packed_outputs_offsets_rotation_and_worker_determinism(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            options = ("--preset", "python", "--ancestral", "--samples", "6", "--tips", "8",
                       "--length", "500", "--seed", "42", "--indel", "0.002,0.002",
                       "--indel-size", "POW{1.7/100},POW{1.7/100}", "--require-unique-window",
                       "--max-attempts", "2", "--file-prefix", "batch.001")
            files = root / "files"
            run(files, *options)
            expected = {path.stem: ((files / "trees" / f"{path.stem}.nwk").read_bytes(), path.read_bytes(),
                                   json.loads((files / "metadata" / f"{path.stem}.json").read_text()))
                        for path in (files / "msa_untrimmed").glob("*.fa")}
            self.assertTrue(expected)
            for jobs, target, buffer in ((1, 1024, 262144), (3, 1024, 1048576), (3, 32768, 16777216)):
                output = root / f"{jobs}-{target}"
                run(output, *options, "--jobs", str(jobs), "--output-format", "packed", "--shard-size-bytes", str(target),
                    "--shard-buffer-bytes", str(buffer))
                self.assertEqual(packed_snapshot(output), expected)
                for name in ("trees", "msa_untrimmed", "metadata", "unaligned", "msa_trimmed"):
                    self.assertFalse((output / name).exists())
                shards = list(output.rglob("*.tar"))
                if target == 1024:
                    self.assertEqual(len(shards), len(expected))  # Oversize pairs each get their own shard.
                else:
                    self.assertTrue(all(path.stat().st_size <= target for path in shards))
                with closing(sqlite3.connect(output / "index.sqlite3")) as db:
                    paths = [row[0] for row in db.execute("SELECT DISTINCT shard FROM pairs")]
                    matches = [re.fullmatch(r"shards/batch\.001-worker-(\d{3,})-shard-\d{6,}\.tar", path) for path in paths]
                    self.assertTrue(all(matches))
                    owners = {match.group(1) for match in matches}
                    self.assertLessEqual(len(owners), jobs)
                self.assertEqual({path.name for path in output.iterdir()}, {"shards", "index.sqlite3", "simulation_summary.json"})
                self.assertTrue(all(path.parent == output / "shards" for path in shards))
                self.assertTrue(all(path.is_file() and path.suffix == ".tar" for path in (output / "shards").iterdir()))
                self.assertFalse((output / "workers").exists())
                summary = json.loads((output / "simulation_summary.json").read_text())
                self.assertEqual(summary["output_layout"]["shard_size_bytes"], target)
                self.assertEqual(summary["output_layout"]["shard_buffer_bytes"], buffer)
                self.assertIn(f"--shard-buffer-bytes={buffer}", summary["resolved_arguments"])
                self.assertEqual(summary["file_prefix"], "batch.001")
                self.assertIn("--file-prefix=batch.001", summary["resolved_arguments"])
                command = summary["reproduce_command"] + ' && printf "\nREPLAY=%s\n" "$replay_dir"'
                result = subprocess.run(["bash", "-c", command], capture_output=True, text=True, timeout=30)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                replay = Path(re.search(r"^REPLAY=(.+)$", result.stdout, re.MULTILINE).group(1))
                try:
                    self.assertEqual(packed_snapshot(replay / "data"), expected)
                finally:
                    import shutil
                    shutil.rmtree(replay)

    def test_packed_prefixes_extract_and_merge_without_collisions(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            extracted = root / "extracted"
            seen_members = set()
            seen_shards = set()
            expected = {}
            sources = []
            for batch, prefix, seed in (("a", "batch-a", "1"), ("b", "b" * 40, "2")):
                output = root / "batches" / batch
                sources.append(output)
                run(output, "--output-format", "packed", "--samples", "4", "--jobs", "3", "--tips", "6",
                    "--length", "100", "--allow-duplicate-sequences", "--ancestral", "--preset", "python",
                    "--seed", seed, "--batch-id", batch, "--file-prefix", prefix, "--shard-size-bytes", "1024")
                snapshot = packed_snapshot(output)
                self.assertEqual(len(snapshot), 4)
                expected.update({f"{batch}/{name}": row for name, row in snapshot.items()})
                for shard in sorted((output / "shards").glob("*.tar")):
                    self.assertNotIn(shard.name, seen_shards)
                    self.assertTrue(shard.name.startswith(f"{prefix}-worker-"))
                    seen_shards.add(shard.name)
                    with tarfile.open(shard) as archive:
                        for member in archive.getmembers():
                            self.assertNotIn(member.name, seen_members)
                            seen_members.add(member.name)
                        if hasattr(tarfile, "data_filter"):
                            archive.extractall(extracted, filter="data")
                        else:
                            archive.extractall(extracted)
                with closing(sqlite3.connect(output / "index.sqlite3")) as db:
                    for name, tree_name, sequence_name in db.execute("SELECT name, tree_name, sequence_name FROM pairs"):
                        tree, sequence, _ = snapshot[name]
                        self.assertEqual((extracted / tree_name).read_bytes(), tree)
                        self.assertEqual((extracted / sequence_name).read_bytes(), sequence)
                        self.assertEqual(Path(tree_name).stem, Path(sequence_name).stem)
            self.assertEqual(len(seen_members), 16)
            self.assertEqual({p.relative_to(extracted).as_posix() for p in extracted.rglob("*") if p.is_file()}, seen_members)
            script = Path(__file__).resolve().parents[1] / "scripts" / "merge_packed.py"
            result = subprocess.run([sys.executable, str(script), *map(str, sources), "--out", str(root)],
                                    capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertEqual(packed_snapshot(root), expected)

    def test_file_prefix_directory_views_trees_only_and_validation(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            options = ("--samples", "2", "--tips", "4", "--seqtype", "DNA", "--model", "JC", "--length", "50",
                       "--indel", "0.1,0.05", "--indel-size", "GEO{2},GEO{3}", "--allow-duplicate-sequences", "--seed", "42", "--materialize", "all")
            baseline, prefixed = root / "baseline", root / "prefixed"
            run(baseline, *options)
            run(prefixed, *options, "--file-prefix", "batch.x")
            for kind in ("trees", "msa_untrimmed", "msa_trimmed", "unaligned"):
                for path in (baseline / kind).iterdir():
                    self.assertEqual(path.read_bytes(), (prefixed / kind / f"batch.x-{path.name}").read_bytes())
            for path in (prefixed / "metadata").glob("*.json"):
                metadata = json.loads(path.read_text())
                self.assertEqual(path.stem, metadata["sample"])
                self.assertTrue(path.name.startswith("batch.x-sample_"))
            trees = root / "trees"
            run(trees, "--trees-only", "--samples", "1", "--file-prefix", "trees", "--seed", "1")
            self.assertTrue((trees / "trees" / "trees-sample_000000.nwk").is_file())
            for i, invalid in enumerate(("", "../escape", "a/b", "a\\b", "a b", "a" * 41)):
                output = root / f"invalid{i}"
                result = run(output, "--file-prefix", invalid, success=False)
                self.assertIn("--file-prefix", result.stderr)
                self.assertFalse(output.exists())

    def test_packed_quota_empty_failures_and_validation(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            options = ("--output-format", "packed", "--samples", "4", "--tips", "4",
                       "--length", "2", "--seqtype", "DNA", "--model", "JC", "--branch-mean", "0.5",
                       "--seed", "1", "--max-attempts", "1")
            snapshots = []
            for jobs in (1, 3):
                output = root / f"quota{jobs}"
                run(output, *options, "--ensure-samples", "--jobs", str(jobs))
                snapshots.append(packed_snapshot(output))
                self.assertEqual(len(snapshots[-1]), 4)
                self.assertNotIn("sample_000000", snapshots[-1])
            self.assertEqual(*snapshots)
            empty = root / "empty"
            run(empty, "--output-format", "packed", "--samples", "3", "--jobs", "2", "--tips", "5",
                "--length", "1", "--seqtype", "DNA", "--model", "JC", "--max-attempts", "1", "--seed", "1")
            self.assertEqual(packed_snapshot(empty), {})
            fatal = root / "fatal"
            run(fatal, "--output-format", "packed", "--jobs", "2", "--samples", "3", "--model", "INVALID_MODEL", success=False)
            self.assertFalse((fatal / "index.sqlite3").exists())
            for flag in ("--trees-only", "--materialize=all"):
                out = root / flag.replace("=", "_")
                result = run(out, "--output-format", "packed", flag, success=False)
                self.assertIn("one full MSA/tree pair", result.stderr)
                self.assertFalse(out.exists())
            invalid = root / "invalid-buffer"
            run(invalid, "--shard-buffer-bytes", "0", success=False)
            self.assertFalse(invalid.exists())

    def test_packed_batch_merge_without_repacking(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            options = ("--output-format", "packed", "--samples", "3", "--jobs", "2", "--tips", "5",
                       "--length", "100", "--allow-duplicate-sequences")
            batches = []
            expected = {}
            for batch, seed in (("a", "1"), ("b", "2")):
                output = root / "batches" / batch
                run(output, *options, "--seed", seed, "--batch-id", batch,
                    "--shard-buffer-bytes", str(262144 if batch == "a" else 16777216))
                batches.append(output)
                expected.update({f"{batch}/{name}": row for name, row in packed_snapshot(output).items()})
            original = {path: path.read_bytes() for path in root.rglob("*.tar")}
            script = Path(__file__).resolve().parents[1] / "scripts" / "merge_packed.py"
            result = subprocess.run([sys.executable, str(script), *map(str, batches), "--out", str(root)],
                                    capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertEqual(packed_snapshot(root), expected)
            self.assertEqual({path: path.read_bytes() for path in root.rglob("*.tar")}, original)
            with closing(sqlite3.connect(root / "index.sqlite3")) as db:
                self.assertEqual(db.execute("SELECT COUNT(*) FROM batches").fetchone()[0], 2)
            duplicate = root / "batches" / "duplicate"
            run(duplicate, *options, "--seed", "3", "--batch-id", "a")
            result = subprocess.run([sys.executable, str(script), str(batches[0]), str(duplicate), "--out", str(root / "batches")],
                                    capture_output=True, text=True, timeout=30)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Repeated batch identity", result.stderr)
            self.assertFalse((root / "batches" / "index.sqlite3").exists())
            summary_path = batches[1] / "simulation_summary.json"
            original_summary = summary_path.read_text()
            incomplete = json.loads(original_summary)
            incomplete["status"] = "running"
            summary_path.write_text(json.dumps(incomplete))
            result = subprocess.run([sys.executable, str(script), *map(str, batches), "--out", str(root / "batches")],
                                    capture_output=True, text=True, timeout=30)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Not a completed packed simulation", result.stderr)
            summary_path.write_text(original_summary)
            changed = root / "batches" / "changed"
            run(changed, *options, "--seed", "4", "--batch-id", "changed", "--gamma", "G8")
            command = [sys.executable, str(script), str(batches[0]), str(changed), "--out", str(root / "batches")]
            result = subprocess.run(command, capture_output=True, text=True, timeout=30)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Simulation policies differ", result.stderr)
            self.assertFalse((root / "batches" / "index.sqlite3").exists())
            result = subprocess.run([*command, "--allow-mixed-configs"], capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            with closing(sqlite3.connect(root / "batches" / "index.sqlite3")) as db:
                self.assertEqual(db.execute("SELECT COUNT(*) FROM pairs").fetchone()[0], 6)


    def test_ancestral_exports_map_to_tree_and_preserve_worker_determinism(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            options = ("--preset", "python", "--ancestral", "--samples", "6", "--tips", "8",
                       "--length", "500", "--seed", "42", "--indel", "0.002,0.002",
                       "--indel-size", "POW{1.7/100},POW{1.7/100}", "--materialize", "all",
                       "--require-unique-window", "--max-attempts", "1")
            outputs = []
            for jobs in [1, 3]:
                output = root / str(jobs)
                run(output, *options, "--jobs", str(jobs))
                files = sorted((output / "msa_untrimmed").glob("*.fa"))
                self.assertGreater(len(files), 0)
                for path in files:
                    rows = fasta(path)
                    self.assertEqual(len(rows), 15)
                    internal = set(rows) - {f"T{i}" for i in range(1, 9)}
                    tree = (output / "trees" / f"{path.stem}.nwk").read_text()
                    self.assertEqual(set(re.findall(r"I_\d+", tree)), internal)
                    self.assertEqual(fasta(output / "unaligned" / path.name), {name: row.replace("-", "") for name, row in rows.items()})
                    meta = json.loads((output / "metadata" / f"{path.stem}.json").read_text())
                    self.assertEqual(set(meta["ancestral_names"]), internal)
                    self.assertTrue(meta["ancestral"])
                    window = meta["window"]
                    self.assertEqual(fasta(output / "msa_trimmed" / path.name), {name: row[window["start"]:window["end"]] for name, row in rows.items()})
                outputs.append({p.relative_to(output): p.read_bytes() for p in output.rglob("*") if p.is_file() and p.name != "simulation_summary.json"})
            self.assertEqual(*outputs)

    def test_default_output_and_optional_materialized_views(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            options = ("--tips", "8", "--samples", "3", "--length", "100", "--seed", "42",
                       "--branch-mean", "0.3", "--indel", "0.1,0.05", "--indel-size", "GEO{2},GEO{3}")
            plain, full = root / "plain", root / "full"
            run(plain, *options)
            run(full, *options, "--materialize", "trimmed,unaligned")
            self.assertFalse((plain / "msa_trimmed").exists())
            self.assertFalse((plain / "unaligned").exists())
            self.assertFalse((plain / "sequences").exists())
            for source in (plain / "msa_untrimmed").glob("*.fa"):
                raw = fasta(source)
                self.assertEqual(source.read_bytes(), (full / "msa_untrimmed" / source.name).read_bytes())
                self.assertEqual(fasta(full / "msa_trimmed" / source.name), {name: row[:100] for name, row in raw.items()})
                self.assertEqual(fasta(full / "unaligned" / source.name), {name: row.replace("-", "") for name, row in raw.items()})
                metadata = json.loads((plain / "metadata" / f"{source.stem}.json").read_text())
                self.assertEqual(metadata["window"]["start"], 0)
                self.assertEqual(metadata["window"]["end"], 100)
                self.assertFalse(metadata["window"]["uniqueness_checked"])
                self.assertEqual(metadata["checks"], {"distinct_gapped": True, "distinct_ungapped": False, "unique_window": False})
                self.assertEqual(metadata["alignment_columns"], len(next(iter(raw.values()))))
                self.assertEqual(metadata["tip_names"], [f"T{i}" for i in range(1, 9)])
            conflict = run(root / "conflict", *options, "--materialize", "unaligned", "--no-unaligned", success=False)
            self.assertIn("conflicts", conflict.stderr)

    def test_independent_checks_allow_override_and_window_failures(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            impossible = ("--tips", "5", "--length", "1", "--seqtype", "DNA", "--model", "JC",
                          "--max-attempts", "2", "--samples", "2", "--seed", "42")
            for name, checks, reason in [
                ("gapped", (), "gapped"),
                ("ungapped", ("--require-distinct-gapped=false", "--require-distinct-ungapped"), "ungapped"),
                ("window", ("--require-distinct-gapped=false", "--require-unique-window"), "window"),
            ]:
                output = root / name
                run(output, *impossible, *checks, "--jobs", "2")
                summary = json.loads((output / "simulation_summary.json").read_text())
                self.assertEqual(summary["failed_samples"], 2)
                self.assertEqual(summary["failed_samples_by_last_rejection"][reason], 2)
                self.assertEqual(summary["failed_samples_percentage_by_last_rejection"][reason], 100.0)
                self.assertEqual(summary["written_samples"], 0)
                self.assertEqual(list((output / "metadata").iterdir()), [])
                self.assertEqual(summary["failed_duplicate_samples"], 0 if reason == "window" else 2)
            all_checks = ("--require-distinct-gapped", "--require-distinct-ungapped", "--require-unique-window")
            allowed = root / "allowed"
            run(allowed, *impossible, *all_checks, "--allow-duplicate-sequences", "--materialize", "all", "--jobs", "2")
            for path in (allowed / "metadata").glob("*.json"):
                metadata = json.loads(path.read_text())
                self.assertFalse(any(metadata["checks"].values()))
                self.assertEqual(metadata["attempt_number"], 1)
                self.assertEqual(metadata["window"]["selection"], "prefix")
                self.assertFalse(metadata["window"]["uniqueness_checked"])
                msa = fasta(allowed / "msa_untrimmed" / f"{path.stem}.fa")
                self.assertLess(len(set(msa.values())), len(msa))
                self.assertEqual(fasta(allowed / "unaligned" / f"{path.stem}.fa"), msa)
            unchecked = root / "unchecked"
            run(unchecked, *impossible, "--require-distinct-gapped=false")
            self.assertEqual(len(list((unchecked / "msa_untrimmed").glob("*.fa"))), 2)

    def test_window_selection_materialization_replay_and_workers(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            options = ("--tips", "8", "--samples", "4", "--length", "60", "--seed", "17",
                       "--branch-mean", "0.3", "--indel", "0.1,0.05", "--indel-size", "GEO{2},GEO{3}",
                       "--require-distinct-ungapped", "--require-unique-window", "--window-policy", "inclusive",
                       "--materialize", "all")
            serial = root / "serial"
            run(serial, *options)
            self.assertEqual(len(list((serial / "metadata").glob("*.json"))), 4)
            for path in (serial / "metadata").glob("*.json"):
                meta = json.loads(path.read_text())
                msa = fasta(serial / "msa_untrimmed" / f"{path.stem}.fa")
                rows = list(msa.values())
                width = len(rows[0])
                start = next(i for i in range(width - 60 + 1) if len({row[i:i+60] for row in rows}) == len(rows))
                self.assertEqual(meta["window"]["start"], start)
                self.assertEqual(meta["window"]["end"], start + 60)
                self.assertTrue(meta["window"]["uniqueness_checked"])
                self.assertEqual(fasta(serial / "msa_trimmed" / f"{path.stem}.fa"), {name: row[start:start+60] for name, row in msa.items()})
                self.assertEqual(len({row.replace("-", "") for row in rows}), len(rows))
            parallel = root / "parallel"
            run(parallel, *options, "--jobs", "3")
            def snapshot(output):
                return {p.relative_to(output): p.read_bytes() for p in output.rglob("*") if p.is_file() and p.name != "simulation_summary.json"}
            self.assertEqual(snapshot(serial), snapshot(parallel))
            summary = json.loads((serial / "simulation_summary.json").read_text())
            self.assertIn("--require-unique-window=true", summary["resolved_arguments"])
            self.assertIn("--materialize=all", summary["resolved_arguments"])
            result = subprocess.run(["bash", "-c", summary["reproduce_command"] + '\nprintf "REPLAY=%s\\n" "$replay_dir"'], capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stderr)
            replay = Path(re.search(r"REPLAY=(.+)", result.stdout).group(1))
            try:
                self.assertEqual(snapshot(serial), snapshot(replay / "data"))
            finally:
                import shutil
                shutil.rmtree(replay)
            # Failed windows must also trigger replacement candidates under --ensure-samples.
            quota_options = ("--tips", "4", "--length", "2", "--seqtype", "DNA", "--model", "JC",
                             "--samples", "5", "--seed", "1", "--branch-mean", "0.5", "--max-attempts", "1",
                             "--require-distinct-gapped=false", "--require-unique-window", "--ensure-samples")
            quotas = []
            for jobs in [1, 3]:
                out = root / f"quota_{jobs}"
                run(out, *quota_options, "--jobs", str(jobs))
                info = json.loads((out / "simulation_summary.json").read_text())
                self.assertEqual(info["written_samples"], 5)
                self.assertGreater(info["failed_window_samples"], 0)
                quotas.append(snapshot(out))
            self.assertEqual(*quotas)

    def test_empirical_gamma_presets_replay_and_workers(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            options = ("--preset", "python", "--samples", "5", "--tips", "8",
                       "--length", "80", "--seed", "20261102", "--allow-duplicate-sequences")
            baseline = root / "baseline"
            run(baseline, *options, "--alisim-log", str(root / "baseline.log"))
            summary = json.loads((baseline / "simulation_summary.json").read_text())
            self.assertEqual(summary["msa_policy"]["gamma"], "G8")
            self.assertIn("--gamma=G8", summary["resolved_arguments"])
            self.assertEqual(summary["msa_policy"]["empirical_gamma"]["source_sha256"],
                             "1527b12722751a395ee9483ee15af2f92715a25dc2560011e47d269090355ccb")
            for name, extra in [("parallel", ("--jobs", "3")),
                                ("explicit", ("--model", "LG", "--gamma", "G8"))]:
                output = root / name
                run(output, *options, *extra)
                for path in baseline.rglob("*"):
                    if path.is_file() and path.name != "simulation_summary.json":
                        self.assertEqual(path.read_bytes(), (output / path.relative_to(baseline)).read_bytes())
            # Exercise the recorded command: its explicit --model must not cancel G8.
            command = summary["reproduce_command"] + '\nprintf "REPLAY=%s\\n" "$replay_dir"'
            result = subprocess.run(["bash", "-c", command], capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stderr)
            replay = Path(re.search(r"REPLAY=(.+)", result.stdout).group(1))
            try:
                for path in baseline.rglob("*.fa"):
                    self.assertEqual(path.read_bytes(), (replay / "data" / path.relative_to(baseline)).read_bytes())
            finally:
                import shutil
                shutil.rmtree(replay)
            continuous = root / "continuous"
            run(continuous, *options, "--gamma", "GC")
            continuous_summary = json.loads((continuous / "simulation_summary.json").read_text())
            self.assertEqual(continuous_summary["msa_policy"]["gamma"], "GC")
            for path in (continuous / "metadata").glob("*.json"):
                self.assertTrue(json.loads(path.read_text())["model"].startswith("LG+GC{"))
            for name, extra in [("off", ("--gamma", "none")),
                                ("override", ("--model", "LG+G4{0.5}")),
                                ("bare", ("--model", "LG"))]:
                output = root / name
                run(output, *options, *extra)
                metadata = json.loads((output / "simulation_summary.json").read_text())
                self.assertEqual(metadata["msa_policy"]["gamma"], "none")

    def test_empirical_gamma_retries_and_validation(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            log = root / "attempts.log"
            run(root / "retries", "--gamma", "G8", "--seqtype", "DNA", "--model", "JC",
                "--tips", "5", "--length", "1", "--samples", "2", "--seed", "42",
                "--max-attempts", "3", "--alisim-log", str(log))
            # Five DNA tips and one site guarantee duplicates on every attempt.
            models = re.findall(r"Sample sample_\d+, attempt \d+/3, seed \d+, model (JC\+G8\{[^}]+\})", log.read_text())
            self.assertEqual(len(models), 6)
            self.assertEqual(len(set(models)), 6)
            for gamma in ["G0", "G", "invalid"]:
                run(root / gamma, "--gamma", gamma, success=False)
            result = run(root / "conflict", "--gamma", "GC", "--model", "LG+G4{0.5}", success=False)
            self.assertIn("existing +G/+R", result.stderr)
            self.assertFalse((root / "conflict").exists())

    def test_ensure_samples_replaces_failures_without_extra_outputs(self) -> None:
        with tempfile.TemporaryDirectory(prefix="pf ensure samples ") as directory:
            root = Path(directory)
            # Candidate zero fails with this seed; later trees can succeed.
            options = ("--samples", "8", "--tips", "4", "--length", "2",
                       "--branch-mean", "0.5", "--seqtype", "DNA", "--model", "JC",
                       "--seed", "1", "--max-attempts", "1")
            expected = None
            completed = None
            for jobs in (1, 2, 4):
                output, log = root / f"jobs_{jobs}", root / f"jobs_{jobs}.log"
                result = run(output, *options, "--jobs", str(jobs), "--ensure-samples",
                             "--alisim-log", str(log))
                self.assertEqual(result.stderr.count("significantly increase runtime"), 1)
                summary = json.loads((output / "simulation_summary.json").read_text())
                self.assertEqual(summary["status"], "completed")
                self.assertTrue(summary["ensure_samples"])
                self.assertEqual(summary["sample_count_mode"], "successful_outputs")
                self.assertEqual(summary["requested_samples"], 8)
                self.assertEqual(summary["written_samples"], 8)
                self.assertGreater(summary["failed_duplicate_samples"], 0)
                total = summary["completed_samples"]
                self.assertEqual(total, 8 + summary["failed_duplicate_samples"])
                self.assertEqual(summary["failed_duplicate_percentage_denominator"], "completed_samples")
                self.assertAlmostEqual(summary["failed_duplicate_percentage"],
                                       100 * summary["failed_duplicate_samples"] / total)
                self.assertIn("--ensure-samples", summary["reproduce_command"])
                self.assertEqual(log.read_text().count("[Alignment Simulator] Executing"), total)
                self.assertEqual(log.read_text().count("Skipped sample:"), total - 8)
                self.assertEqual(list(root.glob(f"{log.name}.worker-*")), [])
                trees = sorted((output / "trees").glob("*.nwk"))
                msas = sorted((output / "msa_untrimmed").glob("*.fa"))
                self.assertEqual(len(trees), 8)
                self.assertEqual(len(msas), 8)
                self.assertEqual({p.stem for p in trees}, {p.stem for p in msas})
                self.assertFalse((output / "trees" / "sample_000000.nwk").exists())
                self.assertEqual(msas[-1].stem, f"sample_{total - 1:06}")
                for path in msas:
                    rows = fasta(path)
                    self.assertEqual(len(set(rows.values())), 4)
                files = {p.relative_to(output): p.read_bytes() for p in output.rglob("*")
                         if p.is_file() and p.name != "simulation_summary.json"}
                if expected is None:
                    expected, completed = files, total
                self.assertEqual(files, expected)
                self.assertEqual(total, completed)

            # The enforced dataset is exactly the accepted prefix of the same
            # candidate stream, including fresh trees beyond the original quota.
            reference = root / "bounded_reference"
            run(reference, *options[2:], "--samples", str(completed))
            self.assertEqual({p.relative_to(reference): p.read_bytes() for p in reference.rglob("*")
                              if p.is_file() and p.name != "simulation_summary.json"}, expected)

            # Increasing the output target preserves all previous accepted files.
            larger = root / "larger"
            run(larger, *options[2:], "--samples", "10", "--jobs", "4", "--ensure-samples")
            self.assertEqual(len(list((larger / "msa_untrimmed").glob("*.fa"))), 10)
            for path, content in expected.items():
                self.assertEqual((larger / path).read_bytes(), content)

    def test_ensure_samples_bypass_and_fatal_errors(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            # Impossible uniqueness is harmless when duplicate checking is bypassed.
            output = root / "bypass"
            run(output, "--samples", "3", "--jobs", "4", "--ensure-samples",
                "--allow-duplicate-sequences", "--tips", "5", "--length", "1",
                "--seqtype", "DNA", "--model", "JC", "--seed", "1")
            summary = json.loads((output / "simulation_summary.json").read_text())
            self.assertEqual(summary["written_samples"], 3)
            self.assertEqual(summary["completed_samples"], 3)
            self.assertEqual(summary["failed_duplicate_samples"], 0)

            output = root / "fatal"
            run(output, "--samples", "3", "--jobs", "2", "--ensure-samples",
                "--model", "NOT_A_REAL_MODEL", success=False)
            summary = json.loads((output / "simulation_summary.json").read_text())
            self.assertEqual(summary["status"], "failed")
            self.assertEqual(summary["completed_samples"], 0)
            self.assertEqual(summary["failed_duplicate_percentage"], 0)
            self.assertIn("Worker", summary["error"])

    def test_generated_seed_summary_and_executable_replay_command(self) -> None:
        with tempfile.TemporaryDirectory(prefix="pf summary ' ") as directory:
            root = Path(directory)
            prior = root / "empirical's prior $sample.tsv"
            prior.write_text("1\n2 3\n")
            seeds = set()
            for jobs in (1, 3):
                output, log = root / f"random_{jobs}", root / f"random_{jobs}.log"
                start = time.perf_counter()
                run(output, "--jobs", str(jobs), "--samples", "5", "--tips", "8",
                    "--length", "100", "--model", "LG+G4{0.5}", "--indel", "0.01,0.01",
                    "--indel-size", "GEO{2},GEO{3}", "--empirical-diameters", str(prior),
                    "--minimum-terminal", "0.001", "--alisim-log", str(log))
                elapsed = time.perf_counter() - start
                summary = json.loads((output / "simulation_summary.json").read_text())
                self.assertEqual(summary["status"], "completed")
                self.assertEqual(summary["seed_source"], "generated")
                self.assertGreaterEqual(summary["seed"], 0)
                self.assertLess(summary["seed"], 2**64)
                self.assertNotIn(summary["seed"], seeds)
                seeds.add(summary["seed"])
                self.assertEqual(summary["requested_samples"], 5)
                self.assertEqual(summary["written_samples"] + summary["failed_duplicate_samples"], 5)
                self.assertAlmostEqual(summary["failed_duplicate_percentage"],
                                       100 * summary["failed_duplicate_samples"] / 5)
                self.assertGreater(summary["wall_time_seconds"], 0)
                self.assertLessEqual(summary["wall_time_seconds"], elapsed)
                self.assertIn("--seed=", summary["reproduce_command"])
                self.assertFalse((output / "simulation_summary.json.tmp").exists())

                # Actually execute the saved Bash, from another cwd, with quoted
                # filenames and a redirected log. Replays get fresh directories.
                # Read the directory from the shell: mktemp implementations and
                # shell environments need not place it under our supplied TMPDIR.
                # Use && so printing the path cannot hide a failed replay.
                command = summary["reproduce_command"] + ' && printf "\nREPLAY=%s\n" "$replay_dir"'
                replay = subprocess.run(["bash", "-c", command],
                                        capture_output=True, text=True, timeout=30, cwd="/",
                                        env={**os.environ, "TMPDIR": str(root)})
                self.assertEqual(replay.returncode, 0, replay.stdout + replay.stderr)
                marker = re.search(r"^REPLAY=(.+)$", replay.stdout, re.MULTILINE)
                self.assertIsNotNone(marker, replay.stdout + replay.stderr)
                replay_parent = Path(marker.group(1))
                try:
                    reproduced = replay_parent / "data"
                    replay_summary = json.loads((reproduced / "simulation_summary.json").read_text())
                    self.assertEqual(replay_summary["seed"], summary["seed"])
                    self.assertEqual(replay_summary["seed_source"], "provided")
                    self.assertEqual(replay_summary["failed_duplicate_samples"], summary["failed_duplicate_samples"])
                    original_files = {p.relative_to(output): p.read_bytes() for p in output.rglob("*")
                                      if p.is_file() and p.name != "simulation_summary.json"}
                    replay_files = {p.relative_to(reproduced): p.read_bytes() for p in reproduced.rglob("*")
                                    if p.is_file() and p.name != "simulation_summary.json"}
                    self.assertEqual(original_files, replay_files)
                finally:
                    import shutil
                    shutil.rmtree(replay_parent)

    def test_explicit_seed_boundaries_and_zero_failure_summary(self) -> None:
        with tempfile.TemporaryDirectory(prefix="pf explicit seed ") as directory:
            for seed in (0, 2**64 - 1):
                output = Path(directory) / str(seed)
                run(output, f"--seed={seed}", "--samples", "3", "--jobs", "2", "--trees-only")
                summary = json.loads((output / "simulation_summary.json").read_text())
                self.assertEqual(summary["seed"], seed)
                self.assertEqual(summary["seed_source"], "provided")
                self.assertEqual(summary["written_samples"], 3)
                self.assertEqual(summary["failed_duplicate_samples"], 0)
                self.assertEqual(summary["failed_duplicate_percentage"], 0)
                self.assertGreater(summary["wall_time_seconds"], 0)

    def test_summary_names_cannot_be_used_as_alisim_logs(self) -> None:
        with tempfile.TemporaryDirectory(prefix="pf reserved summary ") as directory:
            for index, name in enumerate(("simulation_summary.json", "simulation_summary.json.tmp")):
                output = Path(directory) / str(index)
                result = run(output, "--alisim-log", str(output / name), success=False)
                self.assertIn("summary filename", result.stderr)
                summary = json.loads((output / "simulation_summary.json").read_text())
                self.assertEqual(summary["status"], "failed")
                self.assertEqual(summary["completed_samples"], 0)
                self.assertFalse((output / "simulation_summary.json.tmp").exists())

    def test_process_workers_preserve_outputs_failures_and_logs(self) -> None:
        cases = [
            ("protein", ("--tips", "20", "--length", "100", "--model", "LG+G4{0.5}")),
            ("python", ("--tips", "20", "--length", "100", "--preset", "python")),
            ("empirical_indels", ("--tips", "20", "--length", "100", "--model", "LG",
                                  "--gamma", "G8", "--indel", "0.01,0.01",
                                  "--indel-size", "POW{1.7/100},POW{1.7/100}")),
            ("indels", ("--tips", "20", "--length", "100", "--indel", "0.02,0.01",
                        "--indel-size", "GEO{2},GEO{3}")),
            ("mixed", ("--tips", "4", "--length", "2", "--branch-mean", "0.5",
                       "--seqtype", "DNA", "--model", "JC", "--max-attempts", "1")),
            ("all_fail", ("--tips", "5", "--length", "1", "--seqtype", "DNA",
                          "--model", "JC", "--max-attempts", "2")),
            ("trees", ("--tips", "20", "--preset", "python", "--trees-only")),
        ]
        with tempfile.TemporaryDirectory(prefix="pf parallel ") as directory:
            root = Path(directory)
            for name, flags in cases:
                baseline = root / f"{name}_serial"
                serial_log = root / f"{name}_serial.log"
                options = ("--samples", "7", "--seed", "1", *flags)
                serial = run(baseline, *options, "--alisim-log", str(serial_log))
                expected = {p.relative_to(baseline): p.read_bytes()
                            for p in baseline.rglob("*") if p.is_file() and p.name != "simulation_summary.json"}
                for jobs in (2, 4):
                    output = root / f"{name}_{jobs}"
                    log = root / f"{name}_{jobs}.log"
                    # Exercise both forms when forwarding options to workers.
                    log_options = ("--alisim-log", str(log)) if jobs == 2 else (f"--alisim-log={log}",)
                    parallel = run(output, *options, "--jobs", str(jobs), *log_options)
                    self.assertEqual(parallel.stderr, "")
                    self.assertEqual(parallel.stdout.replace(str(output), "OUT"),
                                     serial.stdout.replace(str(baseline), "OUT"))
                    actual = {p.relative_to(output): p.read_bytes()
                              for p in output.rglob("*") if p.is_file() and p.name != "simulation_summary.json"}
                    self.assertEqual(actual, expected)
                    for marker in ("[Alignment Simulator] Executing", "Rejected alignment:",
                                   "Accepted alignment", "Skipped sample:"):
                        self.assertEqual(log.read_text().count(marker), serial_log.read_text().count(marker))
                    self.assertEqual(list(root.glob(f"{log.name}.worker-*")), [])
            # A single sample falls back to serial without spawning excess jobs.
            run(root / "one", "--samples", "1", "--jobs", "100", "--trees-only")
            invalid = root / "invalid"
            run(invalid, "--jobs", "0", success=False)
            self.assertFalse(invalid.exists())

    def test_process_worker_fatal_errors_and_startup_cleanup(self) -> None:
        with tempfile.TemporaryDirectory(prefix="pf worker failure ") as directory:
            root = Path(directory)
            output, log = root / "fatal", root / "fatal.log"
            result = run(output, "--jobs", "2", "--samples", "20", "--model", "NOT_A_REAL_MODEL",
                         "--alisim-log", str(log), success=False)
            self.assertIn("Worker", result.stderr)
            self.assertNotIn("Wrote", result.stdout)
            self.assertNotIn("Sample sample_000002", log.read_text())
            self.assertEqual(list(root.glob(f"{log.name}.worker-*")), [])
            self.assertEqual(list((output / "msa_untrimmed").iterdir()), [])
            summary = json.loads((output / "simulation_summary.json").read_text())
            self.assertEqual(summary["status"], "failed")
            self.assertIn("Worker", summary["error"])
            self.assertGreater(summary["wall_time_seconds"], 0)

            # A late worker-startup error must stop already spawned workers,
            # and must never overwrite a pre-existing log fragment.
            log = root / "collision.log"
            collision = root / "collision.log.worker-001.log"
            collision.write_text("keep")
            run(root / "collision", "--jobs", "2", "--samples", "100",
                "--alisim-log", str(log), success=False)
            self.assertEqual(collision.read_text(), "keep")
            self.assertFalse((root / "collision.log.worker-000.log").exists())

    def test_short_and_long_help(self) -> None:
        def help_text(flag: str) -> str:
            result = subprocess.run([EXECUTABLE, flag], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(result.stderr, "")
            return result.stdout

        short, detailed = help_text("-h"), help_text("--help")
        self.assertIn("--preset", short)
        self.assertIn("--samples", short)
        self.assertIn("--help", short)
        self.assertNotIn("--small-event-mean", short)
        self.assertNotIn("--gamma-shape", short)
        self.assertIn("--small-event-mean", detailed)
        self.assertIn("--output-format", detailed)
        self.assertIn("--shard-size-bytes", detailed)
        self.assertIn("--gamma-shape", detailed)
        self.assertIn("--python-bd", detailed)
        self.assertIn("--preset python", detailed)
        self.assertIn("--ensure-samples", detailed)
        self.assertIn("significantly increase runtime", detailed)
        self.assertIn("Tree priors:", detailed)
        self.assertIn("Default: G8", detailed)
        self.assertIn("GC (continuous Gamma) can be much slower than G8", detailed)
        self.assertLess(len(short.splitlines()), 25)

    def test_presets_preserve_tree_and_alignment_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name, original, preset in [
                ("default", (), ("--preset", "default")),
                ("python", ("--python-bd",), ("--preset", "python")),
            ]:
                first, second = root / f"{name}_old", root / f"{name}_preset"
                options = ("--tips", "6", "--length", "100", "--seed", "42")
                run(first, *original, *options)
                run(second, *preset, *options)
                for path in first.rglob("*"):
                    if path.is_file() and path.name != "simulation_summary.json":
                        self.assertEqual(path.read_bytes(), (second / path.relative_to(first)).read_bytes())
            invalid = root / "invalid"
            run(invalid, "--preset", "unknown", success=False)
            self.assertFalse(invalid.exists())

    def test_protein_pairs_and_reproducibility(self) -> None:
        with tempfile.TemporaryDirectory(prefix="pf sim ") as directory:
            first, second = Path(directory) / "a", Path(directory) / "b"
            options = ("--samples", "3", "--tips", "8", "--length", "120", "--seed", "42",
                       "--model", "LG+G4{0.5}")
            result = run(first, *options)
            self.assertNotIn("[Alignment Simulator]", result.stdout + result.stderr)
            self.assertNotIn("SUBSTITUTION PROCESS", result.stdout + result.stderr)
            self.assertEqual(result.stderr, "")  # Progress is hidden outside a terminal.
            run(second, *options)
            trees = sorted((first / "trees").glob("*.nwk"))
            self.assertEqual(len(trees), 3)
            self.assertEqual(len(list((first / "msa_untrimmed").iterdir())), 3)
            self.assertFalse((first / "unaligned").exists())

            for tree in trees:
                paths = tree_paths(tree)
                sequences = fasta(first / "msa_untrimmed" / f"{tree.stem}.fa")
                self.assertEqual(set(paths), {f"T{i}" for i in range(1, 9)})
                self.assertEqual(set(paths), set(sequences))
                self.assertEqual({len(seq) for seq in sequences.values()}, {120})
                self.assertLessEqual(set("".join(sequences.values())), set("ACDEFGHIKLMNPQRSTVWY"))

            for path in first.rglob("*"):
                if path.is_file() and path.name != "simulation_summary.json":
                    self.assertEqual(path.read_bytes(), (second / path.relative_to(first)).read_bytes())

            # Increasing sample count must preserve the previous prefix.
            prefix = Path(directory) / "prefix"
            run(prefix, "--samples", "1", *options[2:])
            for folder, suffix in (("trees", ".nwk"), ("msa_untrimmed", ".fa")):
                filename = Path(folder) / f"sample_000000{suffix}"
                self.assertEqual((prefix / filename).read_bytes(), (first / filename).read_bytes())

    def test_jc_distances_follow_tree(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "dna"
            run(output, "--tips", "4", "--length", "50000", "--model", "JC",
                "--seqtype", "DNA", "--seed", "99")
            paths = tree_paths(output / "trees" / "sample_000000.nwk")
            sequences = fasta(output / "msa_untrimmed" / "sample_000000.fa")
            names = sorted(paths)
            for i, a in enumerate(names):
                for b in names[i + 1:]:
                    distance = sum(paths[a].get(edge, paths[b].get(edge, 0.0))
                                   for edge in paths[a].keys() ^ paths[b].keys())
                    expected = 0.75 * (1 - math.exp(-4 * distance / 3))
                    observed = sum(x != y for x, y in zip(sequences[a], sequences[b])) / 50000
                    self.assertAlmostEqual(observed, expected, delta=0.015)

    def test_validation_and_no_overwrite(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "existing"
            output.mkdir()
            marker = output / "keep.txt"
            marker.write_text("keep")
            run(output, success=False)
            self.assertEqual(marker.read_text(), "keep")
            for index, options in enumerate([
                ("--tips", "2"), ("--samples", "0"), ("--seed", "-1"),
                ("--length", "2147483648"), ("--branch-mean", "nan"),
                ("--branch-mean", "0.1junk"), ("--unknown", "1"),
                ("--indel-size", "GEO{2},GEO{3}"),
                ("--max-attempts", "0"),
            ]):
                target = Path(directory) / str(index)
                run(target, *options, success=False)
                self.assertFalse(target.exists())

    def test_optional_alisim_log(self) -> None:
        with tempfile.TemporaryDirectory(prefix="pf log ") as directory:
            root = Path(directory)
            quiet, logged = root / "quiet", root / "logged"
            log = root / "alisim output.log"
            options = ("--samples", "2", "--tips", "4", "--length", "50", "--seed", "42")
            run(quiet, *options)
            result = run(logged, *options, "--alisim-log", str(log))
            content = log.read_text()
            self.assertEqual(content.count("[Alignment Simulator] Executing"), 2)
            self.assertIn("Model of substitution: LG", content)
            self.assertIn("sample_000000", content)
            self.assertIn("sample_000001", content)
            self.assertNotIn("[Alignment Simulator]", result.stdout + result.stderr)
            self.assertNotIn("\x1b", content)
            self.assertNotIn("Wrote 2", content)
            for path in quiet.rglob("*"):
                if path.is_file() and path.name != "simulation_summary.json":
                    self.assertEqual(path.read_bytes(), (logged / path.relative_to(quiet)).read_bytes())

            original = log.read_bytes()
            run(root / "existing_log", *options, "--alisim-log", str(log), success=False)
            self.assertEqual(log.read_bytes(), original)

            # Logging must not create a file where a tree or MSA will be written.
            conflict = root / "conflict"
            run(conflict, *options, "--alisim-log",
                str(conflict / "msa_untrimmed" / "sample_000000.fa"), success=False)
            self.assertFalse((conflict / "msa_untrimmed" / "sample_000000.fa").exists())

            # Native errors still reach stderr when informational output is hidden.
            result = run(root / "invalid_model", "--model", "INVALID_MODEL", success=False)
            self.assertIn("INVALID_MODEL", result.stderr)

    def test_indel_alignments_and_gap_free_output(self) -> None:
        with tempfile.TemporaryDirectory(prefix="pf indels ") as directory:
            root = Path(directory)
            for seqtype, model, alphabet in [
                ("AA", "LG+G4{0.5}", "ACDEFGHIKLMNPQRSTVWY"),
                ("DNA", "JC", "ACGT"),
            ]:
                with self.subTest(seqtype=seqtype):
                    first, repeat, aligned_only = (root / f"{seqtype}_{name}"
                                                    for name in ("first", "repeat", "aligned"))
                    options = ("--samples", "2", "--tips", "6", "--length", "500",
                               "--branch-mean", "0.3", "--seed", "42",
                               "--seqtype", seqtype, "--model", model,
                               "--indel", "0.1,0.05", "--indel-size", "GEO{2},GEO{3}")
                    result = run(first, *options, "--materialize", "unaligned")
                    self.assertNotIn("[Alignment Simulator]", result.stdout + result.stderr)
                    run(repeat, *options, "--materialize", "unaligned")
                    run(aligned_only, *options, "--no-unaligned")

                    self.assertEqual(len(list((first / "msa_untrimmed").iterdir())), 2)
                    self.assertEqual({path.name for path in (first / "msa_untrimmed").iterdir()},
                                     {path.name for path in (first / "unaligned").iterdir()})
                    self.assertEqual(len(list((aligned_only / "msa_untrimmed").iterdir())), 2)
                    self.assertFalse((aligned_only / "unaligned").exists())
                    for tree in (first / "trees").glob("*.nwk"):
                        paths = tree_paths(tree)
                        aligned_path = Path("msa_untrimmed") / f"{tree.stem}.fa"
                        aligned = fasta(first / aligned_path)
                        unaligned = fasta(first / "unaligned" / f"{tree.stem}.fa")
                        self.assertEqual(set(paths), set(aligned))
                        self.assertEqual(set(paths), set(unaligned))
                        self.assertEqual(len({len(seq) for seq in aligned.values()}), 1)
                        self.assertGreater(len(next(iter(aligned.values()))), 500)
                        self.assertIn("-", "".join(aligned.values()))
                        self.assertLessEqual(set("".join(aligned.values())), set(alphabet + "-"))
                        self.assertEqual({name: seq.replace("-", "") for name, seq in aligned.items()},
                                         unaligned)
                        self.assertEqual((first / aligned_path).read_bytes(),
                                         (aligned_only / aligned_path).read_bytes())
                    for path in first.rglob("*"):
                        if path.is_file() and path.name != "simulation_summary.json":
                            self.assertEqual(path.read_bytes(), (repeat / path.relative_to(first)).read_bytes())

            # Exercise the original workflow's power-law syntax, including '/'.
            run(root / "powerlaw", "--samples", "2", "--length", "200", "--seed", "42",
                "--indel", "0.002,0.002", "--indel-size", "POW{1.7/100},POW{1.7/100}",
                "--no-unaligned")
            result = run(root / "bad_rates", "--indel", "invalid", success=False)
            self.assertIn("--indel", result.stderr)

    def test_tree_priors_and_python_bd_pipeline(self) -> None:
        with tempfile.TemporaryDirectory(prefix="pf tree priors ") as directory:
            root = Path(directory)
            variants = {
                "uniform": ("--tree-prior", "uniform"),
                "uniform_rooted": ("--tree-prior", "uniform-rooted"),
                "ete3": ("--tree-prior", "ete3"),
                "yule": ("--tree-prior", "yule"),
                "bd_raw": ("--tree-prior", "birth-death"),
                "bd_python": ("--python-bd",),
                "gamma": ("--branch-prior", "gamma", "--gamma-shape", "2", "--gamma-scale", "0.05"),
            }
            for name, flags in variants.items():
                out = root / name
                run(out, "--tips", "6", "--length", "100", "--allow-duplicate-sequences", *flags)
                sequences = fasta(out / "msa_untrimmed" / "sample_000000.fa")
                self.assertEqual(len(sequences), 6)
                self.assertEqual(set(sequences), {f"T{i}" for i in range(1, 7)})
            run(root / "trees_only", "--samples", "3", "--tips", "20", "--python-bd", "--trees-only")
            self.assertEqual(len(list((root / "trees_only" / "trees").glob("*.nwk"))), 3)
            self.assertEqual(list((root / "trees_only" / "msa_untrimmed").iterdir()), [])
            run(root / "conflict", "--python-bd", "--tree-prior", "ete3", success=False)
            run(root / "bad_gamma", "--branch-prior", "gamma", "--gamma-shape", "-1", success=False)

    def test_duplicate_retry_limit_and_bypass(self) -> None:
        with tempfile.TemporaryDirectory(prefix="pf retries ") as directory:
            root = Path(directory)
            # Five one-site DNA sequences cannot all be unique (four states).
            impossible = ("--tips", "5", "--length", "1", "--seqtype", "DNA", "--model", "JC")
            failed, log = root / "failed", root / "failed.log"
            result = run(failed, *impossible, "--samples", "3", "--max-attempts", "3",
                         "--alisim-log", str(log))
            self.assertEqual(result.stderr, "")
            self.assertIn("Wrote 0 tree/FASTA pairs", result.stdout)
            self.assertIn("3 failed due to duplicate sequences; 3 requested", result.stdout)
            self.assertEqual(log.read_text().count("[Alignment Simulator] Executing"), 9)
            self.assertEqual(log.read_text().count("Rejected alignment:"), 9)
            self.assertEqual(log.read_text().count("Skipped sample:"), 3)
            self.assertIn("No duplicate-free MSA for sample_000002 after 3 attempts", log.read_text())
            self.assertEqual(list((failed / "trees").iterdir()), [])
            self.assertEqual(list((failed / "msa_untrimmed").iterdir()), [])
            self.assertFalse((failed / "unaligned").exists())
            summary = json.loads((failed / "simulation_summary.json").read_text())
            self.assertEqual(summary["status"], "completed")
            self.assertEqual(summary["written_samples"], 0)
            self.assertEqual(summary["failed_duplicate_samples"], 3)
            self.assertEqual(summary["failed_duplicate_percentage"], 100)

            allowed, allowed_log = root / "allowed", root / "allowed.log"
            run(allowed, *impossible, "--max-attempts", "3", "--allow-duplicate-sequences",
                "--alisim-log", str(allowed_log))
            self.assertEqual(allowed_log.read_text().count("[Alignment Simulator] Executing"), 1)
            sequences = fasta(allowed / "msa_untrimmed" / "sample_000000.fa")
            self.assertLess(len(set(sequences.values())), len(sequences))

    def test_retries_resimulate_a_fixed_tree_and_are_reproducible(self) -> None:
        with tempfile.TemporaryDirectory(prefix="pf retry success ") as directory:
            root = Path(directory)
            # This fixed seed needs retries for this small alignment.
            options = ("--tips", "4", "--length", "2", "--branch-mean", "0.5",
                       "--seqtype", "DNA", "--model", "JC", "--seed", "1")
            one_attempt = run(root / "one_attempt", *options, "--max-attempts", "1")
            self.assertIn("1 failed due to duplicate sequences; 1 requested", one_attempt.stdout)
            self.assertEqual(list((root / "one_attempt" / "msa_untrimmed").iterdir()), [])
            accepted, log = root / "accepted", root / "attempts.log"
            run(accepted, *options, "--max-attempts", "20", "--alisim-log", str(log))
            content = log.read_text()
            count = content.count("[Alignment Simulator] Executing")
            self.assertGreater(count, 1)
            self.assertLessEqual(count, 20)
            self.assertEqual(content.count("Rejected alignment:"), count - 1)
            seeds = re.findall(r"attempt \d+/\d+, seed (\d+)", content)
            self.assertEqual(len(set(seeds)), count)
            sequences = fasta(accepted / "msa_untrimmed" / "sample_000000.fa")
            self.assertEqual(len(set(sequences.values())), len(sequences))
            self.assertEqual(len(list((accepted / "msa_untrimmed").iterdir())), 1)

            # Exactly enough attempts succeeds; one fewer skips this sample.
            short_limit = run(root / "short_limit", *options, "--max-attempts", str(count - 1))
            self.assertIn("1 failed due to duplicate sequences; 1 requested", short_limit.stdout)
            repeat = root / "repeat"
            run(repeat, *options, "--max-attempts", str(count))
            first = root / "first"
            run(first, *options, "--allow-duplicate-sequences")
            for folder, extension in (("trees", "nwk"), ("msa_untrimmed", "fa")):
                relative = Path(folder) / f"sample_000000.{extension}"
                self.assertEqual((accepted / relative).read_bytes(), (repeat / relative).read_bytes())
            self.assertEqual((accepted / "trees" / "sample_000000.nwk").read_bytes(),
                             (first / "trees" / "sample_000000.nwk").read_bytes())

    def test_mixed_failures_keep_sample_ids_and_seeds(self) -> None:
        with tempfile.TemporaryDirectory(prefix="pf mixed failures ") as directory:
            root = Path(directory)
            options = ("--samples", "12", "--tips", "4", "--length", "2",
                       "--branch-mean", "0.5", "--seqtype", "DNA", "--model", "JC", "--seed", "1")
            reference, mixed = root / "reference", root / "mixed"
            run(reference, *options, "--allow-duplicate-sequences")
            result = run(mixed, *options, "--max-attempts", "1")
            self.assertEqual(result.stderr, "")
            accepted, failed = 0, 0
            for sample in range(12):
                name = f"sample_{sample:06}"
                sequences = fasta(reference / "msa_untrimmed" / f"{name}.fa")
                if len(set(sequences.values())) < len(sequences):
                    failed += 1
                    self.assertFalse((mixed / "trees" / f"{name}.nwk").exists())
                    self.assertFalse((mixed / "msa_untrimmed" / f"{name}.fa").exists())
                else:
                    accepted += 1
                    for folder, suffix in (("trees", "nwk"), ("msa_untrimmed", "fa")):
                        relative = Path(folder) / f"{name}.{suffix}"
                        self.assertEqual((mixed / relative).read_bytes(), (reference / relative).read_bytes())
            self.assertGreater(accepted, 0)
            self.assertGreater(failed, 0)
            self.assertIn(f"Wrote {accepted} tree/FASTA pairs", result.stdout)
            self.assertIn(f"{failed} failed due to duplicate sequences; 12 requested", result.stdout)
            summary = json.loads((mixed / "simulation_summary.json").read_text())
            self.assertEqual(summary["written_samples"], accepted)
            self.assertEqual(summary["failed_duplicate_samples"], failed)
            self.assertAlmostEqual(summary["failed_duplicate_percentage"], 100 * failed / 12)
            self.assertEqual(len(list((mixed / "trees").iterdir())), accepted)
            self.assertEqual(len(list((mixed / "msa_untrimmed").iterdir())), accepted)

    def test_non_duplicate_simulation_errors_still_stop_the_run(self) -> None:
        with tempfile.TemporaryDirectory(prefix="pf fatal error ") as directory:
            root = Path(directory)
            output, log = root / "output", root / "attempts.log"
            result = run(output, "--samples", "3", "--model", "NOT_A_REAL_MODEL",
                         "--alisim-log", str(log), success=False)
            self.assertNotIn("Wrote", result.stdout)
            self.assertIn("Sample sample_000000, attempt 1/20", log.read_text())
            self.assertNotIn("Sample sample_000001", log.read_text())
            self.assertNotIn("Skipped sample:", log.read_text())
            self.assertEqual(list((output / "trees").iterdir()), [])
            self.assertEqual(list((output / "msa_untrimmed").iterdir()), [])
            # IQ-TREE exits the process here, so the initial record survives
            # without inventing a completed wall-time measurement.
            summary = json.loads((output / "simulation_summary.json").read_text())
            self.assertEqual(summary["status"], "running")
            self.assertIsNone(summary["wall_time_seconds"])
            self.assertIsNone(summary["written_samples"])


if __name__ == "__main__":
    unittest.main()
