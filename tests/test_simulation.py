"""Exercise the executable and check tree/sequence consistency without dependencies."""

from __future__ import annotations

import math
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest


EXECUTABLE = str(Path(sys.argv.pop(1)).resolve())


def run(output: Path, *options: str, success: bool = True) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(
        [EXECUTABLE, "--out", str(output), *options], capture_output=True, text=True
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


class SimulationTests(unittest.TestCase):
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
            self.assertEqual(len(list((first / "sequences").iterdir())), 3)
            self.assertFalse((first / "unaligned_sequences").exists())

            for tree in trees:
                paths = tree_paths(tree)
                sequences = fasta(first / "sequences" / f"{tree.stem}.fa")
                self.assertEqual(set(paths), {f"T{i}" for i in range(1, 9)})
                self.assertEqual(set(paths), set(sequences))
                self.assertEqual({len(seq) for seq in sequences.values()}, {120})
                self.assertLessEqual(set("".join(sequences.values())), set("ACDEFGHIKLMNPQRSTVWY"))

            for path in first.rglob("*"):
                if path.is_file():
                    self.assertEqual(path.read_bytes(), (second / path.relative_to(first)).read_bytes())

            # Increasing sample count must preserve the previous prefix.
            prefix = Path(directory) / "prefix"
            run(prefix, "--samples", "1", *options[2:])
            for folder, suffix in (("trees", ".nwk"), ("sequences", ".fa")):
                filename = Path(folder) / f"sample_000000{suffix}"
                self.assertEqual((prefix / filename).read_bytes(), (first / filename).read_bytes())

    def test_jc_distances_follow_tree(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "dna"
            run(output, "--tips", "4", "--length", "50000", "--model", "JC",
                "--seqtype", "DNA", "--seed", "99")
            paths = tree_paths(output / "trees" / "sample_000000.nwk")
            sequences = fasta(output / "sequences" / "sample_000000.fa")
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
                if path.is_file():
                    self.assertEqual(path.read_bytes(), (logged / path.relative_to(quiet)).read_bytes())

            original = log.read_bytes()
            run(root / "existing_log", *options, "--alisim-log", str(log), success=False)
            self.assertEqual(log.read_bytes(), original)

            # Logging must not create a file where a tree or MSA will be written.
            conflict = root / "conflict"
            run(conflict, *options, "--alisim-log",
                str(conflict / "sequences" / "sample_000000.fa"), success=False)
            self.assertFalse((conflict / "sequences" / "sample_000000.fa").exists())

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
                    result = run(first, *options)
                    self.assertNotIn("[Alignment Simulator]", result.stdout + result.stderr)
                    run(repeat, *options)
                    run(aligned_only, *options, "--no-unaligned")

                    self.assertEqual(len(list((first / "sequences").iterdir())), 2)
                    self.assertEqual({path.name for path in (first / "sequences").iterdir()},
                                     {path.name for path in (first / "unaligned_sequences").iterdir()})
                    self.assertEqual(len(list((aligned_only / "sequences").iterdir())), 2)
                    self.assertFalse((aligned_only / "unaligned_sequences").exists())
                    for tree in (first / "trees").glob("*.nwk"):
                        paths = tree_paths(tree)
                        aligned_path = Path("sequences") / f"{tree.stem}.fa"
                        aligned = fasta(first / aligned_path)
                        unaligned = fasta(first / "unaligned_sequences" / f"{tree.stem}.fa")
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
                        if path.is_file():
                            self.assertEqual(path.read_bytes(), (repeat / path.relative_to(first)).read_bytes())

            # Exercise the original workflow's power-law syntax, including '/'.
            run(root / "powerlaw", "--samples", "2", "--length", "200", "--seed", "42",
                "--indel", "0.002,0.002", "--indel-size", "POW{1.7/100},POW{1.7/100}",
                "--no-unaligned")
            result = run(root / "bad_rates", "--indel", "invalid", success=False)
            self.assertIn("--indel", result.stderr)

    def test_duplicate_retry_limit_and_bypass(self) -> None:
        with tempfile.TemporaryDirectory(prefix="pf retries ") as directory:
            root = Path(directory)
            # Five one-site DNA sequences cannot all be unique (four states).
            impossible = ("--tips", "5", "--length", "1", "--seqtype", "DNA", "--model", "JC")
            failed, log = root / "failed", root / "failed.log"
            result = run(failed, *impossible, "--max-attempts", "3",
                         "--alisim-log", str(log), success=False)
            self.assertIn("after 3 attempts", result.stderr)
            self.assertEqual(log.read_text().count("[Alignment Simulator] Executing"), 3)
            self.assertEqual(log.read_text().count("Rejected alignment:"), 3)
            self.assertEqual(list((failed / "trees").iterdir()), [])
            self.assertEqual(list((failed / "sequences").iterdir()), [])
            self.assertFalse((failed / "unaligned_sequences").exists())

            allowed, allowed_log = root / "allowed", root / "allowed.log"
            run(allowed, *impossible, "--max-attempts", "3", "--allow-duplicate-sequences",
                "--alisim-log", str(allowed_log))
            self.assertEqual(allowed_log.read_text().count("[Alignment Simulator] Executing"), 1)
            sequences = fasta(allowed / "sequences" / "sample_000000.fa")
            self.assertLess(len(set(sequences.values())), len(sequences))

    def test_retries_resimulate_a_fixed_tree_and_are_reproducible(self) -> None:
        with tempfile.TemporaryDirectory(prefix="pf retry success ") as directory:
            root = Path(directory)
            # This fixed seed needs retries for this small alignment.
            options = ("--tips", "4", "--length", "2", "--branch-mean", "0.5",
                       "--seqtype", "DNA", "--model", "JC", "--seed", "1")
            run(root / "one_attempt", *options, "--max-attempts", "1", success=False)
            accepted, log = root / "accepted", root / "attempts.log"
            run(accepted, *options, "--max-attempts", "20", "--alisim-log", str(log))
            content = log.read_text()
            count = content.count("[Alignment Simulator] Executing")
            self.assertGreater(count, 1)
            self.assertLessEqual(count, 20)
            self.assertEqual(content.count("Rejected alignment:"), count - 1)
            seeds = re.findall(r"attempt \d+/\d+, seed (\d+)", content)
            self.assertEqual(len(set(seeds)), count)
            sequences = fasta(accepted / "sequences" / "sample_000000.fa")
            self.assertEqual(len(set(sequences.values())), len(sequences))
            self.assertEqual(len(list((accepted / "sequences").iterdir())), 1)

            # Exactly enough attempts succeeds; one fewer must fail.
            run(root / "short_limit", *options, "--max-attempts", str(count - 1), success=False)
            repeat = root / "repeat"
            run(repeat, *options, "--max-attempts", str(count))
            first = root / "first"
            run(first, *options, "--allow-duplicate-sequences")
            for folder, extension in (("trees", "nwk"), ("sequences", "fa")):
                relative = Path(folder) / f"sample_000000.{extension}"
                self.assertEqual((accepted / relative).read_bytes(), (repeat / relative).read_bytes())
            self.assertEqual((accepted / "trees" / "sample_000000.nwk").read_bytes(),
                             (first / "trees" / "sample_000000.nwk").read_bytes())


if __name__ == "__main__":
    unittest.main()
