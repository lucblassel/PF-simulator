"""Exercise the executable and check tree/sequence consistency without dependencies."""

from __future__ import annotations

import math
import json
import os
import time
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest


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


class SimulationTests(unittest.TestCase):
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
                before = set(root.iterdir())
                replay = subprocess.run(["bash", "-c", summary["reproduce_command"]],
                                        capture_output=True, text=True, timeout=30, cwd="/",
                                        env={**os.environ, "TMPDIR": str(root)})
                self.assertEqual(replay.returncode, 0, replay.stdout + replay.stderr)
                created = set(root.iterdir()) - before
                self.assertEqual(len(created), 1)
                reproduced = created.pop() / "data"
                replay_summary = json.loads((reproduced / "simulation_summary.json").read_text())
                self.assertEqual(replay_summary["seed"], summary["seed"])
                self.assertEqual(replay_summary["seed_source"], "provided")
                self.assertEqual(replay_summary["failed_duplicate_samples"], summary["failed_duplicate_samples"])
                original_files = {p.relative_to(output): p.read_bytes() for p in output.rglob("*")
                                  if p.is_file() and p.name != "simulation_summary.json"}
                replay_files = {p.relative_to(reproduced): p.read_bytes() for p in reproduced.rglob("*")
                                if p.is_file() and p.name != "simulation_summary.json"}
                self.assertEqual(original_files, replay_files)

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
            self.assertEqual(list((output / "sequences").iterdir()), [])
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
        self.assertIn("--gamma-shape", detailed)
        self.assertIn("--python-bd", detailed)
        self.assertIn("--preset python", detailed)
        self.assertIn("Tree priors:", detailed)
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
                if path.is_file() and path.name != "simulation_summary.json":
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
                if path.is_file() and path.name != "simulation_summary.json":
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
                sequences = fasta(out / "sequences" / "sample_000000.fa")
                self.assertEqual(len(sequences), 6)
                self.assertEqual(set(sequences), {f"T{i}" for i in range(1, 7)})
            run(root / "trees_only", "--samples", "3", "--tips", "20", "--python-bd", "--trees-only")
            self.assertEqual(len(list((root / "trees_only" / "trees").glob("*.nwk"))), 3)
            self.assertEqual(list((root / "trees_only" / "sequences").iterdir()), [])
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
            self.assertEqual(list((failed / "sequences").iterdir()), [])
            self.assertFalse((failed / "unaligned_sequences").exists())
            summary = json.loads((failed / "simulation_summary.json").read_text())
            self.assertEqual(summary["status"], "completed")
            self.assertEqual(summary["written_samples"], 0)
            self.assertEqual(summary["failed_duplicate_samples"], 3)
            self.assertEqual(summary["failed_duplicate_percentage"], 100)

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
            one_attempt = run(root / "one_attempt", *options, "--max-attempts", "1")
            self.assertIn("1 failed due to duplicate sequences; 1 requested", one_attempt.stdout)
            self.assertEqual(list((root / "one_attempt" / "sequences").iterdir()), [])
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

            # Exactly enough attempts succeeds; one fewer skips this sample.
            short_limit = run(root / "short_limit", *options, "--max-attempts", str(count - 1))
            self.assertIn("1 failed due to duplicate sequences; 1 requested", short_limit.stdout)
            repeat = root / "repeat"
            run(repeat, *options, "--max-attempts", str(count))
            first = root / "first"
            run(first, *options, "--allow-duplicate-sequences")
            for folder, extension in (("trees", "nwk"), ("sequences", "fa")):
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
                sequences = fasta(reference / "sequences" / f"{name}.fa")
                if len(set(sequences.values())) < len(sequences):
                    failed += 1
                    self.assertFalse((mixed / "trees" / f"{name}.nwk").exists())
                    self.assertFalse((mixed / "sequences" / f"{name}.fa").exists())
                else:
                    accepted += 1
                    for folder, suffix in (("trees", "nwk"), ("sequences", "fa")):
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
            self.assertEqual(len(list((mixed / "sequences").iterdir())), accepted)

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
            self.assertEqual(list((output / "sequences").iterdir()), [])
            # IQ-TREE exits the process here, so the initial record survives
            # without inventing a completed wall-time measurement.
            summary = json.loads((output / "simulation_summary.json").read_text())
            self.assertEqual(summary["status"], "running")
            self.assertIsNone(summary["wall_time_seconds"])
            self.assertIsNone(summary["written_samples"])


if __name__ == "__main__":
    unittest.main()
