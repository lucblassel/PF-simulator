#!/usr/bin/env python3
"""Compare packed writer builds locally, with optional complete simulator runs.

Build `cargo build --release --example benchmark_packed_io --bin pf-sim` first.
Pass --writer NAME=PATH repeatedly to compare saved before/after executables.
The benchmark writes temporary datasets beside its results, then removes them.
"""

import argparse
import csv
import hashlib
import itertools
import json
import os
import platform
import random
import resource
import shutil
import statistics
import subprocess
import tarfile
import time
from pathlib import Path


def executable(value):
    name, separator, path = value.partition("=")
    if not separator or not name:
        raise argparse.ArgumentTypeError("Use NAME=/path/to/executable")
    path = Path(path).resolve()
    if not path.is_file() or not os.access(path, os.X_OK):
        raise argparse.ArgumentTypeError(f"Not executable: {path}")
    return name, path


def archive_runs(output):
    paths = sorted([*output.glob("run-*.log"), *output.glob("run-*.json")])
    with tarfile.open(output / "runs.tar.gz", "x:gz") as archive:
        for path in paths:
            archive.add(path, arcname=path.name)
    for path in paths:
        path.unlink()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--writer", action="append", type=executable, required=True)
    parser.add_argument("--simulator", action="append", type=executable, default=[])
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--workers", default="1,4,8")
    parser.add_argument("--rounds", type=int, default=3)
    parser.add_argument("--small-samples", type=int, default=2048)
    parser.add_argument("--large-samples", type=int, default=64)
    parser.add_argument("--small-payload-bytes", type=int, default=65_536)
    parser.add_argument("--large-payload-bytes", type=int, default=2_600_000)
    parser.add_argument("--workloads", default="small,large")
    parser.add_argument("--formats", default="packed,files")
    parser.add_argument("--shard-size-bytes", type=int, default=1 << 30)
    parser.add_argument("--keep-data", action="store_true")
    args = parser.parse_args()
    workers = [int(value) for value in args.workers.split(",")]
    selected_workloads = args.workloads.split(",")
    formats = args.formats.split(",")
    if not set(selected_workloads) <= {"small", "large"} or not set(formats) <= {"packed", "files"}:
        parser.error("Select small/large workloads and packed/files formats")
    if min(args.small_payload_bytes, args.large_payload_bytes) < 16:
        parser.error("Payloads must be at least 16 bytes")
    if (args.rounds < 1 or min(workers) < 1 or
            max(workers) > min(args.small_samples, args.large_samples)):
        parser.error("Positive rounds/workers required; each workload must cover every worker")
    if len({name for name, _ in args.writer}) != len(args.writer):
        parser.error("Writer variant names must be unique")
    if len({name for name, _ in args.simulator}) != len(args.simulator):
        parser.error("Simulator variant names must be unique")
    args.out.mkdir(parents=True, exist_ok=False)
    data_root = args.out / "data"
    data_root.mkdir()
    workloads = [("small", 50, 500, args.small_payload_bytes, args.small_samples),
                 ("large", 200, 5000, args.large_payload_bytes, args.large_samples)]
    workloads = [workload for workload in workloads if workload[0] in selected_workloads]
    variants = {("writer", name): path for name, path in args.writer}
    variants.update({("simulation", name): path for name, path in args.simulator})
    metadata = {"platform": platform.platform(), "machine": platform.machine(),
                "logical_cpus": os.cpu_count(), "arguments": vars(args) | {
                    "writer": [(n, str(p)) for n, p in args.writer],
                    "simulator": [(n, str(p)) for n, p in args.simulator], "out": str(args.out)},
                "executables": {f"{kind}/{name}": {"path": str(path),
                    "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
                    for (kind, name), path in variants.items()}}
    (args.out / "environment.json").write_text(json.dumps(metadata, indent=2))
    tasks = []
    for kind, name in variants:
        for workload, count, repeat, fmt in itertools.product(workloads, workers, range(args.rounds), formats):
            # File-mode source is identical across writer revisions: one control.
            if kind == "writer" and fmt == "files" and name != args.writer[-1][0]:
                continue
            tasks.append((kind, name, workload, count, repeat, fmt))
    random.Random(20261103).shuffle(tasks)
    rows = []
    accepted = {}
    for i, (kind, name, workload, count, repeat, fmt) in enumerate(tasks):
        case, tips, length, payload, samples = workload
        output = data_root / f"run-{i:04}"
        command = [str(variants[(kind, name)]), "--out", str(output), "--samples", str(samples)]
        if kind == "writer":
            command += ["--workers", str(count), "--format", fmt, "--payload-bytes", str(payload)]
        else:
            command += ["--jobs", str(count), "--output-format", fmt, "--preset", "python",
                        "--gamma", "G8", "--seed", str(20261102 + repeat),
                        "--tips", str(tips), "--length", str(length), "--ancestral",
                        "--indel", "0.002,0.002", "--indel-size", "POW{1.7/100},POW{1.7/100}",
                        "--require-distinct-gapped", "--require-distinct-ungapped=false",
                        "--require-unique-window", "--window-policy", "python", "--max-attempts", "1"]
        command += ["--shard-size-bytes", str(args.shard_size_bytes)]
        environment = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
        usage_before = resource.getrusage(resource.RUSAGE_CHILDREN)
        started = time.perf_counter()
        result = subprocess.run(command, text=True, capture_output=True, env=environment)
        wall = time.perf_counter() - started
        usage_after = resource.getrusage(resource.RUSAGE_CHILDREN)
        (args.out / f"run-{i:04}.log").write_text(" ".join(command) + "\n" + result.stdout + result.stderr)
        if result.returncode:
            raise RuntimeError(f"Run {i} failed; see run-{i:04}.log. Data retained.")
        cpu = usage_after.ru_utime + usage_after.ru_stime - usage_before.ru_utime - usage_before.ru_stime
        row = dict(kind=kind, variant=name, workload=case, workers=count, repeat=repeat,
                   format=fmt, samples=samples, wall_s=wall, cpu_s=cpu)
        if kind == "writer":
            detail = json.loads(result.stdout)
            row.update(timed_s=detail["total_s"], accepted=samples, merge_s=detail["merge_s"])
            for stage in ["init_s", "write_s", "finish_s"]:
                row[f"worker_{stage}_max"] = max(t[stage] for t in detail["worker_timings"])
            (args.out / f"run-{i:04}.json").write_text(json.dumps(detail, indent=2))
        else:
            detail = json.loads((output / "simulation_summary.json").read_text())
            row.update(timed_s=detail["wall_time_seconds"], accepted=detail["written_samples"],
                       merge_s="", worker_init_s_max="", worker_write_s_max="", worker_finish_s_max="")
            identity = case, repeat
            expected = accepted.setdefault(identity, detail["written_samples"])
            assert detail["written_samples"] == expected, "Accepted counts changed across builds/formats/workers"
            (args.out / f"run-{i:04}.json").write_text(json.dumps(detail, indent=2))
        row["accepted_per_s"] = row["accepted"] / row["timed_s"]
        paths = [path for path in output.rglob("*") if path.is_file()]
        row.update(output_bytes=sum(path.stat().st_size for path in paths), output_files=len(paths))
        rows.append(row)
        with (args.out / "metrics.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        if not args.keep_data:
            shutil.rmtree(output)
        print(f"{i + 1}/{len(tasks)} {kind} {name} {case} {fmt} {count} workers: "
              f"{row['timed_s']:.4f}s", flush=True)

    groups = {}
    for row in rows:
        key = tuple(row[c] for c in ["kind", "workload", "format", "workers", "variant"])
        groups.setdefault(key, []).append(row)
    report = ["# Local SQLite and packed-output benchmark", "", f"Platform: {platform.platform()}.", "",
              "Writer runs use prebuilt payloads and concurrent threads with independent worker connections. "
              "They isolate output handling and exclude AliSim/FASTA formatting. Simulation runs, when "
              "requested, use the simulator's actual process pool and the cluster benchmark's G8 options.", "",
              "Run order is randomized with a fixed seed. All datasets use the same output filesystem. "
              "Data enumeration, integrity validation and cleanup are outside the internal timer. "
              "OS caching is not flushed. CPU times include validation and child startup. "
              "SQLite retains its normal durability settings; raw FASTA/tar writes are buffered. "
              "Worker-stage maxima are not additive wall-time components because workers overlap.", "",
              "| Kind | Workload | Format | Workers | Variant | Median internal s | Min–max s | Accepted/s | Max worker init s* | Max worker finish s* | Merge s |",
              "|---|---|---|---:|---|---:|---|---:|---:|---:|---:|"]
    for key, group in sorted(groups.items()):
        durations = [r["timed_s"] for r in group]
        stages = [f"{statistics.median(r[c] for r in group):.4f}" if group[0][c] != "" else "—"
                  for c in ["worker_init_s_max", "worker_finish_s_max", "merge_s"]]
        report.append("| " + " | ".join(map(str, key)) +
                      f" | {statistics.median(durations):.4f} | {min(durations):.4f}–{max(durations):.4f} | "
                      f"{statistics.median(r['accepted_per_s'] for r in group):.1f} | " + " | ".join(stages) + " |")
    report += ["", "*Median across runs of the slowest worker's stage duration. "
               "These figures diagnose setup/finalization overhead; they do not establish cluster filesystem performance.", "",
               "Executable hashes and run settings are in `environment.json`; individual commands, logs, "
               "and timing/summary JSON files are retained in `runs.tar.gz`. `metrics.csv` contains every run.", ""]
    (args.out / "REPORT.md").write_text("\n".join(report))
    archive_runs(args.out)
    print(f"Results: {args.out / 'REPORT.md'}")


if __name__ == "__main__":
    main()
