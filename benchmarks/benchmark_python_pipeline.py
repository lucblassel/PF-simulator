"""Time matched empirical-GC, ancestral, indel and window-selected pipelines.

The Python reference calls the actual legacy tree/alignment/window functions.
One attempt per tree, tip gapped/window checks, no ungapped check, no replacement.
Both materialize full, selected and ungapped rows including ancestors. Python's
additional prefix-trimmed/raw-rejected/native-tree intermediates remain timed.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
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
import sys
import tempfile
import time

import benchmark_python as base

CASES = [("pipeline_50", 128, 50, 500), ("pipeline_200", 32, 200, 5000)]
SETTINGS = [(backend, jobs) for jobs in [1, 4, 8] for backend in ["rust", "python"]]
POW = "POW{1.7/100},POW{1.7/100}"
WINDOW = None
LAST_ALPHA = None


def initialize(legacy_root):
    global WINDOW
    base.initialize(legacy_root)
    root = Path(legacy_root)
    WINDOW = base.load_module("actual_window_selection", root / "find_unduplicated.py")
    base.MSA_MODULE.args.gamma = "GC"
    base.MSA_MODULE.args.allow_duplicate_sequences = True
    base.MSA_MODULE.alphas = base.MSA_MODULE.load_list(root / "hogenom_alphas.txt")
    original_sample = base.MSA_MODULE.sample_scale
    original_simulate = base.MSA_MODULE.simulate_alignment

    def sample_alpha(values):
        global LAST_ALPHA
        LAST_ALPHA = float(original_sample(values))
        return LAST_ALPHA

    def simulate(*args):
        args = list(args)
        args[2], args[3], args[6], args[9] = "GC", True, True, 1
        return original_simulate(*args)

    base.MSA_MODULE.sample_scale = sample_alpha
    base.MSA_MODULE.simulate_alignment = simulate


def python_sample(task):
    from Bio import SeqIO

    row = base.python_sample(task)
    assert row["attempts"] == row["written"] == 1
    out = Path(task["output"])
    name = f"sample_{task['sample']:06}"
    prefix = out / "sequences" / name
    Path(str(prefix) + ".fa").rename(out / "prefix_trimmed" / f"{name}.fa")
    Path(str(prefix) + ".fa.untrimmed").rename(out / "msa_untrimmed" / f"{name}.fa")
    Path(str(prefix) + ".unaligned.fa").rename(out / "unaligned" / f"{name}.fa")
    Path(str(prefix) + ".full.treefile").rename(out / "labelled_trees" / f"{name}.nwk")
    start = time.perf_counter()
    records = list(SeqIO.parse(out / "msa_untrimmed" / f"{name}.fa", "fasta"))
    tips = list(filter(WINDOW.seq_is_tip, records))
    if WINDOW.has_duplicates(tips):
        reason, offset = "gapped", None
    else:
        offset = WINDOW.find_unduplicated(tips, task["length"])
        reason = "window" if offset is None else None
    if reason is None:
        offset = int(offset)  # The legacy M == L case returns False for offset 0.
        SeqIO.write(
            (record[offset : offset + task["length"]] for record in records),
            out / "msa_trimmed" / f"{name}.fa",
            "fasta",
        )
    metadata = dict(
        candidate_id=task["sample"],
        model=f"LG+GC{{{LAST_ALPHA}}}",
        selected=reason is None,
        rejection=reason,
        window=None
        if offset is None
        else dict(start=offset, end=offset + task["length"]),
    )
    (out / "metadata" / f"{name}.json").write_text(json.dumps(metadata) + "\n")
    row.update(
        written=int(reason is None),
        rejection=reason,
        selection_seconds=time.perf_counter() - start,
    )
    return row


def reference_run(args):
    args.data.mkdir()
    for name in [
        "trees",
        "sequences",
        "msa_untrimmed",
        "msa_trimmed",
        "unaligned",
        "metadata",
        "labelled_trees",
        "prefix_trimmed",
    ]:
        (args.data / name).mkdir()
    tasks = [
        dict(
            sample=i,
            seed=args.seed,
            tips=args.tips,
            length=args.length,
            output=str(args.data),
            reference=str(args.reference),
            native_threads=1,
        )
        for i in range(args.samples)
    ]
    if args.jobs == 1:
        initialize(str(args.legacy_root))
        rows = [python_sample(task) for task in tasks]
    else:
        with ProcessPoolExecutor(
            max_workers=min(args.jobs, args.samples),
            initializer=initialize,
            initargs=(str(args.legacy_root),),
        ) as pool:
            rows = list(pool.map(python_sample, tasks, chunksize=1))
    rejected = {
        name: sum(row["rejection"] == name for row in rows)
        for name in ["gapped", "window"]
    }
    summary = dict(
        requested_samples=args.samples,
        written_samples=sum(row["written"] for row in rows),
        failed_samples=sum(rejected.values()),
        failed_samples_by_last_rejection=rejected,
        attempts=args.samples,
        tree_seconds_sum=sum(row["tree_seconds"] for row in rows),
        msa_seconds_sum=sum(row["msa_seconds"] for row in rows),
        selection_seconds_sum=sum(row["selection_seconds"] for row in rows),
    )
    (args.data / "simulation_summary.json").write_text(json.dumps(summary) + "\n")


def read(path):
    from Bio import SeqIO

    records = list(SeqIO.parse(path, "fasta"))
    rows = {record.id: str(record.seq) for record in records}
    assert len(rows) == len(records)
    return rows


def validate(output, samples, tips, length, backend):
    """Excluded from timing: independent exact checks on retained data and metadata."""
    full_files = sorted((output / "msa_untrimmed").glob("*.fa"))
    selected_files = sorted((output / "msa_trimmed").glob("*.fa"))
    selected_ids = {path.stem for path in selected_files}
    summary = json.loads((output / "simulation_summary.json").read_text())
    assert len(selected_ids) == summary["written_samples"]
    assert len(full_files) == (samples if backend == "python" else len(selected_ids))
    assert len(list((output / "trees").glob("*.nwk"))) == len(full_files)
    assert len(list((output / "unaligned").glob("*.fa"))) == len(full_files)
    expected_tips = {f"T{i}" for i in range(1, tips + 1)}
    reconstruction_checks = retained_checks = raw_duplicate_count = window_failures = 0
    retained_sig = hashlib.sha256()
    for path in full_files:
        full = read(path)
        ungapped = read(output / "unaligned" / path.name)
        assert set(full) == set(ungapped) and expected_tips < full.keys()
        assert len(full) == 2 * tips - 1, (backend, path, len(full), 2 * tips - 1)
        assert len({len(row) for row in full.values()}) == 1
        for name, row in full.items():
            assert row.replace("-", "") == ungapped[name]
            reconstruction_checks += 1
        rows = [full[name] for name in sorted(expected_tips)]
        duplicate = len(set(rows)) != tips
        raw_duplicate_count += duplicate
        width = len(rows[0])
        starts = (
            range(width - length)
            if width > length
            else range(1 if width == length else 0)
        )
        offset = (
            None
            if duplicate
            else next(
                (
                    i
                    for i in starts
                    if len({row[i : i + length] for row in rows}) == tips
                ),
                None,
            )
        )
        if offset is None:
            window_failures += not duplicate
        assert (offset is not None) == (path.stem in selected_ids), (
            backend,
            path,
            offset,
        )
        meta = json.loads((output / "metadata" / f"{path.stem}.json").read_text())
        if offset is not None:
            selected = read(output / "msa_trimmed" / path.name)
            assert selected == {
                name: row[offset : offset + length] for name, row in full.items()
            }
            assert (
                meta["window"]["start"] == offset
                and meta["window"]["end"] == offset + length
            )
            retained_checks += len(selected)
            for directory in ["trees", "msa_untrimmed", "msa_trimmed", "unaligned"]:
                filename = f"{path.stem}.nwk" if directory == "trees" else path.name
                retained_sig.update((directory + "/" + filename).encode())
                retained_sig.update((output / directory / filename).read_bytes())
        if backend == "python":
            assert read(output / "prefix_trimmed" / path.name) == {
                name: row[:length] for name, row in full.items()
            }
    all_sig = hashlib.sha256()
    size = 0
    for path in sorted(
        p
        for p in output.rglob("*")
        if p.is_file() and p.name != "simulation_summary.json"
    ):
        data = path.read_bytes()
        all_sig.update(str(path.relative_to(output)).encode() + b"\0" + data)
        size += len(data)
    return dict(
        written=len(selected_files),
        reconstruction_checks=reconstruction_checks,
        selected_row_checks=retained_checks,
        raw_files=len(full_files),
        observed_raw_gapped_duplicates=raw_duplicate_count,
        observed_raw_window_failures=window_failures,
        dataset_sha256=all_sig.hexdigest(),
        retained_sha256=retained_sig.hexdigest(),
        output_bytes=size,
    )


def measure(args, case, backend, jobs, seed):
    name, samples, tips, length = case
    with tempfile.TemporaryDirectory(prefix="pf-pipeline-bench-") as temp:
        out = Path(temp) / "data"
        if backend == "rust":
            command = [
                str(args.binary),
                "--out",
                str(out),
                "--preset",
                "python",
                "--model",
                "LG",
                "--gamma",
                "GC",
                "--samples",
                str(samples),
                "--tips",
                str(tips),
                "--length",
                str(length),
                "--seed",
                str(seed),
                "--jobs",
                str(jobs),
                "--ancestral",
                "--indel",
                "0.002,0.002",
                "--indel-size",
                POW,
                "--require-distinct-gapped",
                "--require-distinct-ungapped=false",
                "--require-unique-window",
                "--window-policy",
                "python",
                "--max-attempts",
                "1",
                "--materialize",
                "all",
            ]
        else:
            command = [
                sys.executable,
                str(Path(__file__).resolve()),
                "--reference-run",
                "--data",
                str(out),
                "--samples",
                str(samples),
                "--tips",
                str(tips),
                "--length",
                str(length),
                "--seed",
                str(seed),
                "--jobs",
                str(jobs),
                "--reference",
                str(args.reference),
                "--legacy-root",
                str(args.legacy_root),
            ]
        before = resource.getrusage(resource.RUSAGE_CHILDREN)
        start = time.perf_counter()
        process = subprocess.run(command, capture_output=True, text=True, timeout=900)
        elapsed = time.perf_counter() - start
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        if process.returncode:
            raise RuntimeError(f"{command}\n{process.stdout}\n{process.stderr}")
        summary = json.loads((out / "simulation_summary.json").read_text())
        measurements = validate(out, samples, tips, length, backend)
        return dict(
            case=name,
            backend=backend,
            jobs=jobs,
            seed=seed,
            requested=samples,
            **measurements,
            failed=samples - measurements["written"],
            failed_gapped=summary["failed_samples_by_last_rejection"].get("gapped", 0),
            failed_window=summary["failed_samples_by_last_rejection"].get("window", 0),
            wall_seconds=elapsed,
            cpu_seconds=after.ru_utime
            + after.ru_stime
            - before.ru_utime
            - before.ru_stime,
            requested_per_second=samples / elapsed,
            accepted_per_second=measurements["written"] / elapsed,
            command=json.dumps(command),
        )


def summarize(output):
    meta = json.loads((output / "metadata.json").read_text())
    rows = list(csv.DictReader((output / "runs.csv").open()))
    result = []
    lines = [
        "# Python pipeline timing: empirical GC, ancestors, indels and windows",
        "",
        "Both generate Python-prior trees and LG + empirical GC with POW indels 0.002/0.002, "
        "export ancestors, check full gapped tips and the first unique length-L tip window, "
        "and materialize full/selected/ungapped FASTAs. One attempt per tree; no ungapped "
        "duplicate check or replacement.",
        "",
        "| Case (pairs × tips × root sites) | Workers | Rust seconds | Python seconds | Python / Rust | Rust retained | Python retained |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for name, samples, tips, length in meta["cases"]:
        for jobs in meta["jobs"]:
            groups = {
                b: [
                    r
                    for r in rows
                    if r["case"] == name
                    and int(r["jobs"]) == jobs
                    and r["backend"] == b
                ]
                for b in ["rust", "python"]
            }
            assert all(len(group) == meta["repeats"] for group in groups.values())
            times = {
                b: statistics.median(float(r["wall_seconds"]) for r in group)
                for b, group in groups.items()
            }
            written = {
                b: statistics.median(int(r["written"]) for r in group)
                for b, group in groups.items()
            }
            per_round_py = {
                r["round"]: float(r["wall_seconds"]) for r in groups["python"]
            }
            ratios = [
                per_round_py[r["round"]] / float(r["wall_seconds"])
                for r in groups["rust"]
            ]
            ratio = statistics.median(ratios)
            lines.append(
                f"| {samples} × {tips} × {length} | {jobs} | {times['rust']:.3f} | {times['python']:.3f} | {ratio:.2f}× | {written['rust']:g}/{samples} | {written['python']:g}/{samples} |"
            )
            for backend, group in groups.items():
                serial = {
                    r["round"]: float(r["wall_seconds"])
                    for r in rows
                    if r["case"] == name
                    and r["backend"] == backend
                    and int(r["jobs"]) == 1
                }
                result.append(
                    dict(
                        case=name,
                        samples=samples,
                        tips=tips,
                        length=length,
                        backend=backend,
                        jobs=jobs,
                        median_seconds=times[backend],
                        min_seconds=min(float(r["wall_seconds"]) for r in group),
                        max_seconds=max(float(r["wall_seconds"]) for r in group),
                        median_written=written[backend],
                        median_accepted_per_second=statistics.median(
                            float(r["accepted_per_second"]) for r in group
                        ),
                        speedup_vs_own_serial=statistics.median(
                            serial[r["round"]] / float(r["wall_seconds"]) for r in group
                        ),
                        python_over_rust=ratio,
                    )
                )
    with (output / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(result[0]))
        writer.writeheader()
        writer.writerows(result)
    lines += [
        "",
        "## Scope and interpretation",
        "",
        f"- Each timing is the median of {meta['repeats']} interleaved rounds. Timing includes process startup/imports, "
        "tree generation, simulation, selection, output I/O and worker shutdown. Independent validation "
        "and temporary-directory cleanup are outside timing. No Snakemake scheduling overhead is charged.",
        "- Python calls the actual legacy functions, with only invocation parameters and native seeds "
        "provided by the harness. It uses the unadapted pinned AliSim 2.4.0 executable. One worker "
        "is the script’s serial tree loop; 4/8 workers are an optimized process-pool reference, each "
        "with one native thread. Indels force one native thread in the original script.",
        "- Python writes raw and ungapped files for rejected candidates, plus its initial first-L "
        "truncation and native labelled tree. Rust skips serialization of rejected candidates and "
        "has no redundant prefix file; Rust additionally writes provenance/window metadata. These "
        "are production workflow differences, not identical disk workloads. All are included in timing.",
        "- Full/selected/ungapped rows include ancestors for BOTH backends. Window and duplicate "
        "checks use tips only, with the legacy excluded-final-window endpoint. The M == L unique "
        "case selects offset zero, matching the actual script after its full-tip duplicate check.",
        "- Within each backend, dataset hashes must match across worker counts for the same seed. "
        "Across backends, trees, alpha draws and sequence RNG streams differ: equal numeric seeds "
        "do not imply identical datasets. Retention and accepted throughput are reported explicitly. "
        "These timing batches do not replace the larger distribution validations.",
        "- Filesystem page caches are warm; there is no forced fsync or cache eviction. "
        "No other benchmark is run concurrently by the harness.",
        "",
        "## Verification",
        "",
        f"All {len(rows)} timed runs passed independent validation. "
        f"{sum(int(r['reconstruction_checks']) for r in rows):,} exact aligned-to-ungapped checks and "
        f"{sum(int(r['selected_row_checks']) for r in rows):,} selected-row checks passed, including ancestors. "
        "Retained offsets were independently checked by exhaustive exact tip-row comparisons. "
        "All within-backend worker-count hashes agreed.",
        "",
        "## Reproduce",
        "",
        "```bash",
        "cargo build --release",
        "/path/to/legacy/.venv/bin/python benchmarks/benchmark_python_pipeline.py \\\n  --binary target/release/pf-sim \\\n  --reference target/alisim-reference/build/iqtree/iqtree2 \\\n  --legacy-root /path/to/bayesnj_simulation/scripts \\\n  --output benchmarks/python_pipeline",
        "```",
        "",
        "See [runs.csv](runs.csv) for all timings, counts, output sizes, hashes and commands; "
        "[metadata.json](metadata.json) records input/source/binary identities, platform and policy; "
        "[summary.csv](summary.csv) contains grouped statistics.",
        "",
    ]
    (output / "REPORT.md").write_text("\n".join(lines))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--binary", type=Path, default=Path("target/release/pf-sim").resolve()
    )
    parser.add_argument(
        "--reference",
        type=Path,
        default=Path("target/alisim-reference/build/iqtree/iqtree2").resolve(),
    )
    parser.add_argument(
        "--legacy-root",
        type=Path,
        default=Path("/Users/lucblassel/Work/bayesnj_simulation/scripts"),
    )
    parser.add_argument(
        "--output", type=Path, default=Path("benchmarks/python_pipeline")
    )
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--cases", nargs="*", choices=[case[0] for case in CASES])
    parser.add_argument("--summarize-only", action="store_true")
    parser.add_argument("--reference-run", action="store_true", help=argparse.SUPPRESS)
    for flag, typ in [
        ("data", Path),
        ("samples", int),
        ("tips", int),
        ("length", int),
        ("seed", int),
        ("jobs", int),
        ("backend", str),
    ]:
        parser.add_argument("--" + flag, type=typ, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.reference_run:
        reference_run(args)
        return
    if args.summarize_only:
        summarize(args.output)
        return
    if args.repeats < 2:
        parser.error("at least two repeats required")
    args.binary = args.binary.resolve()
    args.reference = args.reference.resolve()
    cases = [case for case in CASES if not args.cases or case[0] in args.cases]
    args.output.mkdir(parents=True, exist_ok=False)
    sources = [
        Path(__file__),
        Path(base.__file__),
        args.legacy_root / "alisim.py",
        args.legacy_root / "find_unduplicated.py",
        args.legacy_root / "SimulateAlterRescale.py",
        args.legacy_root / "hogenom_alphas.txt",
        args.legacy_root / "hogenom_diams.txt",
        args.legacy_root / "raxml_diams.txt",
        *sorted(Path("src").glob("*.rs")),
        *sorted(Path("native").glob("*")),
        Path("Cargo.lock"),
        Path("Cargo.toml"),
        Path("build.rs"),
        Path("CMakeLists.txt"),
        Path("data/empirical_gamma.tsv"),
        Path("data/empirical_diameters.tsv"),
        *sorted(
            (base.ROOT / "../../Development/rust/phylotree/src").resolve().rglob("*.rs")
        ),
    ]
    meta = dict(
        started_utc=datetime.now(timezone.utc).isoformat(),
        platform=platform.platform(),
        logical_cpus=os.cpu_count(),
        python=platform.python_version(),
        repeats=args.repeats,
        cases=cases,
        jobs=[1, 4, 8],
        settings=SETTINGS,
        seed_base=20261102,
        policy=dict(
            model="LG",
            gamma="empirical GC",
            indel="0.002,0.002",
            indel_size=POW,
            ancestral=True,
            gapped_check=True,
            ungapped_check=False,
            window_check=True,
            window_policy="python",
            max_attempts=1,
            ensure_samples=False,
            materialize="all",
        ),
        binaries={str(p): base.sha(p) for p in [args.binary, args.reference]},
        sources={str(p): base.sha(p) for p in sources},
        reference_version=subprocess.check_output(
            [str(args.reference), "--version"], text=True
        ).strip(),
    )
    (args.output / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    for backend, jobs in SETTINGS:
        measure(args, ("warmup", 2, 50, 500), backend, jobs, 19937)
    hashes = {}
    with (args.output / "runs.csv").open("w", newline="") as stream:
        writer = None
        for repeat in range(args.repeats):
            for ci, case in enumerate(
                cases[repeat % len(cases) :] + cases[: repeat % len(cases)]
            ):
                shift = (repeat + ci) % len(SETTINGS)
                for backend, jobs in SETTINGS[shift:] + SETTINGS[:shift]:
                    row = measure(args, case, backend, jobs, 20261102 + repeat)
                    row["round"] = repeat
                    key = (case[0], backend, repeat)
                    assert key not in hashes or hashes[key] == row["dataset_sha256"], (
                        "Worker-dependent data",
                        key,
                    )
                    hashes[key] = row["dataset_sha256"]
                    if writer is None:
                        writer = csv.DictWriter(stream, fieldnames=list(row))
                        writer.writeheader()
                    writer.writerow(row)
                    stream.flush()
                    print(
                        f"round {repeat + 1}/{args.repeats} {case[0]} {backend}/{jobs}: {row['wall_seconds']:.3f}s; {row['written']}/{row['requested']} retained",
                        flush=True,
                    )
    meta["completed_utc"] = datetime.now(timezone.utc).isoformat()
    (args.output / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    summarize(args.output)


if __name__ == "__main__":
    main()
