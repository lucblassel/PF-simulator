"""Measure complete dataset throughput and verify worker-count reproducibility.

Uses only the standard library. Run after `cargo build --release`.
All jobs run sequentially across measurements; no competing benchmark batches.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import statistics
import subprocess
import tempfile
import time

CASES = [
    ("tiny", 8, 10, 500, []),
    ("default", 512, 10, 500, []),
    ("large", 128, 200, 5000, []),
    ("python", 128, 200, 1000, ["--preset", "python"]),
    ("indels", 64, 200, 2000, ["--indel", "0.002,0.002"]),
    ("continuous_gamma", 16, 50, 500, ["--model", "LG+GC{0.5}"]),
    ("retries", 64, 4, 2, ["--seqtype", "DNA", "--model", "JC", "--branch-mean", "0.5"]),
]
JOBS = [1, 2, 4, 8]
FIELDS = ["case", "round", "jobs", "seed", "requested", "written", "failed", "wall_seconds",
          "cpu_seconds", "output_bytes", "dataset_sha256"]


def digest_dataset(directory: Path) -> tuple[str, int]:
    digest = hashlib.sha256()
    size = 0
    for path in sorted(p for p in directory.rglob("*") if p.is_file() and p.name != "simulation_summary.json"):
        digest.update(str(path.relative_to(directory)).encode() + b"\0")
        digest.update(path.stat().st_size.to_bytes(8, "little"))
        with path.open("rb") as stream:
            while chunk := stream.read(1024 * 1024):
                size += len(chunk)
                digest.update(chunk)
    return digest.hexdigest(), size


def measure(binary: Path, case: tuple, jobs: int, seed: int, temporary_root: Path) -> dict:
    name, samples, tips, length, flags = case
    with tempfile.TemporaryDirectory(prefix="pf-parallel-", dir=temporary_root) as temp:
        output = Path(temp) / "data"
        command = [str(binary), "--out", str(output), "--samples", str(samples),
                   "--tips", str(tips), "--length", str(length), "--seed", str(seed),
                   "--jobs", str(jobs), *flags]
        before = resource.getrusage(resource.RUSAGE_CHILDREN)
        start = time.perf_counter()
        result = subprocess.run(command, capture_output=True, text=True, timeout=180, check=True)
        elapsed = time.perf_counter() - start
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        written = len(list((output / "sequences").glob("*.fa")))
        failed = samples - written
        expected_summary = f"{failed} failed due to duplicate sequences; {samples} requested"
        if expected_summary not in result.stdout or result.stderr:
            raise RuntimeError(f"Unexpected run output: {result.stdout}\n{result.stderr}")
        # Hashing, counting, and deleting files are outside the timed region.
        digest, size = digest_dataset(output)
        return dict(case=name, jobs=jobs, seed=seed, requested=samples, written=written,
                    failed=failed, wall_seconds=elapsed,
                    cpu_seconds=(after.ru_utime + after.ru_stime - before.ru_utime - before.ru_stime),
                    output_bytes=size, dataset_sha256=digest)


def quartiles(values: list[float]) -> tuple[float, float]:
    q = statistics.quantiles(values, n=4, method="inclusive")
    return q[0], q[2]


def summarize(output: Path) -> None:
    with (output / "runs.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    metadata = json.loads((output / "metadata.json").read_text())
    groups = {(name, jobs): [row for row in rows if row["case"] == name and int(row["jobs"]) == jobs]
              for name, *_ in CASES for jobs in JOBS}
    summary = []
    for name, *_ in CASES:
        baseline = {row["round"]: float(row["wall_seconds"]) for row in groups[name, 1]}
        for jobs in JOBS:
            group = groups[name, jobs]
            if len(group) != metadata["repeats"]:
                raise RuntimeError(f"Incomplete benchmark: {name}, {jobs} jobs")
            times = [float(row["wall_seconds"]) for row in group]
            paired = [baseline[row["round"]] / float(row["wall_seconds"]) for row in group]
            low, high = quartiles(paired)
            summary.append(dict(case=name, jobs=jobs, median_seconds=statistics.median(times),
                                paired_speedup=statistics.median(paired), speedup_q1=low, speedup_q3=high,
                                requested_per_second=statistics.median([int(row["requested"]) / float(row["wall_seconds"]) for row in group]),
                                median_written=statistics.median([int(row["written"]) for row in group]),
                                median_failed=statistics.median([int(row["failed"]) for row in group]),
                                mean_cpu_cores=statistics.median([float(row["cpu_seconds"]) / float(row["wall_seconds"]) for row in group])))
    with (output / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)
    lines = ["# Process-worker benchmark", "", f"Recorded {metadata['started_utc']} on {metadata['platform']}, "
             f"{metadata['logical_cpus']} logical CPUs; release binary. {metadata['repeats']} rounds per case/job count.", "",
             "Each measurement runs the complete CLI: process startup, tree generation, AliSim, "
             "duplicate retries, FASTA conversion, final file writes, and worker shutdown. "
             "One job uses the direct sequential path; higher counts use persistent worker processes. "
             "Logs are disabled. AliSim uses one native thread per worker.", "",
             "## Batch wall times", "",
             "| Workload | Samples | Tips × root sites | 1 job (s) | 2 jobs (s / speedup) | 4 jobs (s / speedup) | 8 jobs (s / speedup) |",
             "| --- | ---: | --- | ---: | ---: | ---: | ---: |"]
    for name, samples, tips, length, _ in CASES:
        group = [next(row for row in summary if row["case"] == name and row["jobs"] == jobs) for jobs in JOBS]
        times = [f"{row['median_seconds']:.3f}" + (f" / {row['paired_speedup']:.2f}×" if row['jobs'] > 1 else "") for row in group]
        lines.append(f"| {name} | {samples} | {tips} × {length} | " + " | ".join(times) + " |")
    lines += ["", "Speedups are medians of same-round time ratios against one job; displayed times are "
              "independent medians. Samples/second count requested samples, including duplicate failures.", "",
              "## Paired speedup variability", "",
              "| Workload | 2 jobs speedup IQR | 4 jobs speedup IQR | 8 jobs speedup IQR |",
              "| --- | ---: | ---: | ---: |"]
    for name, *_ in CASES:
        cells = []
        for jobs in JOBS[1:]:
            row = next(row for row in summary if row['case'] == name and row['jobs'] == jobs)
            cells.append(f"{row['speedup_q1']:.2f}–{row['speedup_q3']:.2f}×")
        lines.append(f"| {name} | " + " | ".join(cells) + " |")
    lines += ["", "## Workloads and correctness", "",
              "All cases use LG amino acids, the default uniform/exponential tree prior, and duplicate "
              "checking unless overridden below. All job counts use the same sample IDs and seed in each round.", "",
              "- `tiny`: eight default samples, to expose worker startup overhead.",
              "- `default`: a larger batch of the default 10-tip, 500-site simulation.",
              "- `large`: larger trees and MSAs, including their full output I/O.",
              "- `python`: the full `--preset python` tree pipeline.",
              "- `indels`: insertion/deletion rates 0.002/0.002 and AliSim's default indel sizes; "
              "both aligned and unaligned FASTA are written.",
              "- `continuous_gamma`: model `LG+GC{0.5}`; substitution-rate sampling is expensive.",
              "- `retries`: 4 tips × 2 DNA sites under JC, branch mean 0.5; at most 20 attempts per tree. "
              "Exhausted samples are skipped and counted.", "",
              f"All {len(rows)} measured runs completed. Dataset hashes, written counts, and failure counts "
              "were checked against the same-round one-job reference for every worker count. "
              "Hashes cover every relative sample file path and every sample output byte; logs and run summaries are excluded.", "",
              "## Method and limits", "",
              "- Five rounds by default, with rotating job-count and workload order; one unmeasured "
              "two-sample warmup per workload. See `metadata.json` for exact options and run settings.",
              "- The same executable is used for every job count; no per-sample subprocess startup. "
              "At most one in-flight sample per worker, with small status messages sent to the parent.",
              "- Runs are not concurrent with one another. Each run writes to a fresh directory and "
              "includes startup and shutdown. File counting, hashing and deletion are excluded from timing.",
              "- Writes use normal buffered filesystem I/O without fsync. This measures CLI completion "
              "on a warm local filesystem, not durable storage throughput or remote/shared storage.",
              "- CPU seconds in `runs.csv` come from `RUSAGE_CHILDREN`, including waited-for descendants. "
              "Aggregate peak RAM was not measured; each worker holds its own native simulation state.",
              "- Small batches are particularly sensitive to scheduling and launch overhead. "
              "These results describe this machine and these workload sizes; benchmark representative production runs.", "",
              "Reproduce from the repository root:", "", "```sh", "cargo build --release -j 8",
              "python3 benchmarks/benchmark_parallelism.py --output benchmarks/parallelism-new", "```", "",
              "Raw results: [runs.csv](runs.csv), [summary.csv](summary.csv), [metadata.json](metadata.json).", ""]
    (output / "REPORT.md").write_text("\n".join(lines))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=Path("target/release/pf-sim"))
    parser.add_argument("--output", type=Path, default=Path("benchmarks/parallelism"))
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--temporary-root", type=Path, default=Path(tempfile.gettempdir()))
    parser.add_argument("--summarize-only", action="store_true")
    args = parser.parse_args()
    if args.summarize_only:
        summarize(args.output)
        return
    if args.repeats < 2:
        parser.error("--repeats must be at least 2")
    binary = args.binary.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    metadata = dict(started_utc=datetime.now(timezone.utc).isoformat(), platform=platform.platform(),
                    logical_cpus=os.cpu_count(), python=platform.python_version(),
                    rustc=subprocess.check_output(["rustc", "--version"], text=True).strip(),
                    binary=str(binary), binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                    harness_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    repeats=args.repeats, jobs=JOBS, cases=CASES,
                    temporary_root=str(args.temporary_root.resolve()), seed_base=1729,
                    order="rotate cases by round, rotate jobs by round + case index")
    (args.output / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    for case in CASES:
        warmup = (case[0], 2, *case[2:])
        measure(binary, warmup, 1, 123, args.temporary_root)
    with (args.output / "runs.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        for repeat in range(args.repeats):
            order = CASES[repeat % len(CASES):] + CASES[:repeat % len(CASES)]
            for case in order:
                index = CASES.index(case)
                offset = (repeat + index) % len(JOBS)
                jobs_order = JOBS[offset:] + JOBS[:offset]
                reference = None
                for jobs in jobs_order:
                    row = measure(binary, case, jobs, 1729 + repeat, args.temporary_root)
                    row["round"] = repeat
                    signature = (row["dataset_sha256"], row["written"], row["failed"])
                    if reference is None:
                        reference = signature
                    elif signature != reference:
                        raise RuntimeError(f"Worker-count-dependent result: {case[0]}, round {repeat}")
                    writer.writerow(row)
                    stream.flush()
                    print(f"round {repeat + 1}/{args.repeats} {case[0]:>16} jobs={jobs}: "
                          f"{row['wall_seconds']:.3f}s, {row['written']} written/{row['failed']} failed", flush=True)
    metadata["completed_utc"] = datetime.now(timezone.utc).isoformat()
    (args.output / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    summarize(args.output)
    print(f"Report: {args.output / 'REPORT.md'}")


if __name__ == "__main__":
    main()
