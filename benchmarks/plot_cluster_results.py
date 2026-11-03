#!/usr/bin/env python3
"""Summarize the Slurm benchmark CSV and plot paired worker/format comparisons.

Requires matplotlib and numpy. Example:
    python benchmarks/plot_cluster_results.py all-metrics.csv --out benchmarks/cluster_2026-11-03
"""

import argparse
import csv
import hashlib
import itertools
from collections import defaultdict
from pathlib import Path
from statistics import median

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--array-job", help="Default: the only complete 96-task array")
    args = parser.parse_args()
    raw = list(csv.DictReader(args.input.open()))
    jobs = defaultdict(list)
    for row in raw:
        jobs[row["array_job"]].append(row)
    complete = [job for job, values in jobs.items()
                if len(values) == 96 and {int(r["task"]) for r in values} == set(range(96))]
    if args.array_job is None and len(complete) != 1:
        parser.error("Select an array explicitly with --array-job; no unique complete array")
    job = args.array_job or complete[0]
    rows = jobs[job]
    if len(rows) != 96 or {int(r["task"]) for r in rows} != set(range(96)):
        parser.error("Expected one complete 96-task array without duplicate tasks")

    integers = ["tips", "length", "workers", "effective_workers", "seed", "repeat",
                "requested", "completed", "written", "failed", "failed_duplicates",
                "failed_windows", "output_files", "output_bytes", "exit_code"]
    floats = ["wall_s", "simulator_wall_s", "cpu_s", "avg_cpu_cores", "accepted_per_s",
              "candidates_per_s", "output_mib_per_s"]
    for row in rows:
        for name in integers:
            row[name] = int(row[name])
        for name in floats:
            row[name] = float(row[name])
        assert row["exit_code"] == 0 and row["status"] == "completed", row
        assert row["completed"] == row["requested"] == row["written"] + row["failed"]
        assert row["failed"] == row["failed_duplicates"] + row["failed_windows"]
        assert row["effective_workers"] == row["workers"]
        assert np.isclose(row["accepted_per_s"], row["written"] / row["wall_s"])
        assert np.isclose(row["candidates_per_s"], row["completed"] / row["wall_s"])
        assert np.isclose(row["avg_cpu_cores"], row["cpu_s"] / row["wall_s"])
    assert len({r["binary_sha256"] for r in rows}) == 1
    workloads = sorted({(r["tips"], r["length"], r["requested"]) for r in rows})
    workers = sorted({r["workers"] for r in rows})
    gammas, formats = ["G8", "GC"], ["packed", "files"]
    assert len(workloads) == 2 and workers == [1, 4, 8, 16]
    groups = defaultdict(list)
    index = {}
    for row in rows:
        key = (row["tips"], row["length"], row["gamma"], row["format"], row["workers"])
        groups[key].append(row)
        identity = (*key, row["seed"])
        assert identity not in index
        index[identity] = row
    seeds = sorted({r["seed"] for r in rows})
    assert len(seeds) == 3
    for tips, length, _ in workloads:
        for gamma, fmt, count in itertools.product(gammas, formats, workers):
            group = groups[(tips, length, gamma, fmt, count)]
            assert sorted(r["seed"] for r in group) == seeds
            for row in group:
                reference = index[(tips, length, gamma, "packed", 1, row["seed"])]
                for name in ["completed", "written", "failed_duplicates", "failed_windows"]:
                    assert row[name] == reference[name], (name, row)

    def values(tips, length, gamma, fmt, count, metric):
        return [r[metric] for r in groups[(tips, length, gamma, fmt, count)]]

    def speeds(tips, length, gamma, fmt, count):
        return [index[(tips, length, gamma, fmt, 1, seed)]["wall_s"] /
                index[(tips, length, gamma, fmt, count, seed)]["wall_s"] for seed in seeds]

    def gains(tips, length, gamma, count):
        return [index[(tips, length, gamma, "files", count, seed)]["wall_s"] /
                index[(tips, length, gamma, "packed", count, seed)]["wall_s"] for seed in seeds]

    args.out.mkdir(parents=True, exist_ok=True)
    summary = []
    for key, group in sorted(groups.items()):
        tips, length, gamma, fmt, count = key
        record = dict(tips=tips, length=length, gamma=gamma, format=fmt, workers=count, repeats=len(group))
        for metric in [*floats, "output_files", "output_bytes", "written"]:
            numbers = [r[metric] for r in group]
            record.update({f"{metric}_median": median(numbers), f"{metric}_min": min(numbers),
                           f"{metric}_max": max(numbers)})
        paired = speeds(tips, length, gamma, fmt, count)
        record.update(speedup_median=median(paired), speedup_min=min(paired), speedup_max=max(paired))
        summary.append(record)
    with (args.out / "summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)

    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.titleweight": "bold", "pdf.fonttype": 42})
    colors = {"G8": "#0072B2", "GC": "#D55E00"}

    def line(axis, series, label, color, style="-", marker="o"):
        center = np.array([median(v) for v in series])
        lower, upper = np.array([min(v) for v in series]), np.array([max(v) for v in series])
        axis.errorbar(workers, center, yerr=[center - lower, upper - center],
                      label=label, color=color, linestyle=style, marker=marker,
                      capsize=3, linewidth=1.8, markersize=5, elinewidth=1)

    def axes_style(axis):
        axis.set_xticks(workers)
        axis.set_xlabel("Simulation workers")
        axis.grid(axis="y", alpha=0.2)
        axis.set_axisbelow(True)

    figure, axes = plt.subplots(2, 2, figsize=(12, 8))
    for col, (tips, length, requested) in enumerate(workloads):
        for gamma, fmt in itertools.product(gammas, formats):
            label = f"{gamma}, {fmt}"
            style, marker = ("-", "o") if fmt == "packed" else ("--", "s")
            line(axes[0, col], [values(tips, length, gamma, fmt, w, "accepted_per_s") for w in workers],
                 label, colors[gamma], style, marker)
            line(axes[1, col], [speeds(tips, length, gamma, fmt, w) for w in workers],
                 label, colors[gamma], style, marker)
        axes[0, col].set_title(f"{tips} tips × {length:,} sites\n{requested:,} candidates per run", pad=10)
        axes[0, col].set_ylabel("Accepted examples / second")
        axes[0, col].set_ylim(bottom=0)
        axes[1, col].plot(workers, workers, ":", color="#777777", linewidth=1.2, label="Linear scaling")
        axes[1, col].set_ylabel("Wall-time speedup vs 1 worker")
        axes[1, col].set_ylim(bottom=0)
        for axis in axes[:, col]:
            axes_style(axis)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.suptitle("PF-simulator throughput and worker scaling", fontsize=17, y=0.99)
    figure.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.945), ncol=4, frameon=False)
    figure.text(0.5, 0.015, "Median of 3 seeds; whiskers show min–max. Speedups pair runs by seed. Dotted line: linear scaling.",
                ha="center", fontsize=9, color="#444444")
    figure.tight_layout(rect=(0, 0.045, 1, 0.885), h_pad=2)
    for extension in ["png", "pdf"]:
        figure.savefig(args.out / f"throughput_scaling.{extension}", dpi=180)
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    work_colors = ["#0072B2", "#8E5EA2"]
    for (tips, length, _), color in zip(workloads, work_colors):
        for gamma in gammas:
            line(axes[0], [gains(tips, length, gamma, w) for w in workers],
                 f"{tips} tips × {length:,}, {gamma}", color,
                 "-" if gamma == "G8" else "--", "o" if gamma == "G8" else "s")
        for fmt in formats:
            line(axes[1], [values(tips, length, "G8", fmt, w, "output_files") for w in workers],
                 f"{tips} tips × {length:,}, {fmt}", color,
                 "-" if fmt == "packed" else "--", "o" if fmt == "packed" else "s")
    axes[0].axhline(1, color="#777777", linestyle=":", linewidth=1)
    axes[0].set_ylabel("Files wall time / packed wall time")
    axes[0].set_title("Packed output speed advantage", pad=10)
    axes[0].set_ylim(bottom=0)
    axes[1].set_title("Files per dataset (G8)", pad=10)
    axes[1].set_ylabel("Regular file count, log scale")
    axes[1].set_yscale("log")
    for axis in axes:
        axes_style(axis)
        axis.legend(frameon=False, fontsize=8, loc="upper right")
    axes[1].legend(frameon=False, fontsize=8, loc="center", bbox_to_anchor=(0.5, 0.73))
    figure.suptitle("Output format tradeoffs", fontsize=17, y=0.99)
    figure.text(0.5, 0.015, "3 seeds per configuration. Left: paired median and min–max; values above 1 favor packed output.",
                ha="center", fontsize=9, color="#444444")
    figure.tight_layout(rect=(0, 0.045, 1, 0.92))
    for extension in ["png", "pdf"]:
        figure.savefig(args.out / f"output_formats.{extension}", dpi=180)
    plt.close(figure)

    md = ["# Jean Zay benchmark results", "",
          f"Source: `{args.input.name}` (SHA-256 `{hashlib.sha256(args.input.read_bytes()).hexdigest()}`).",
          f"Selected array **{job}**: 96 completed runs, 32 configurations × 3 seeds. "
          f"Excluded {len(raw) - len(rows)} earlier run(s) to keep the grid balanced.", "",
          f"Node(s): {', '.join(sorted({r['node'] for r in rows}))}. All runs used the same binary "
          f"(`{rows[0]['binary_sha256']}`) and reported exit code 0.", "",
          "## Throughput and scaling", "",
          "Rates and times below are medians across seeds. Speedups are medians of "
          "the three seed-paired wall-time ratios, so they need not equal ratios of the displayed medians.", "",
          "| Tips × sites | Gamma | Format | Wall s, 1 worker | Wall s, 16 workers | Accepted/s, 1 worker | Accepted/s, 16 workers | Paired speedup, 16 workers |",
          "|---|---|---|---:|---:|---:|---:|---:|"]
    for tips, length, _ in workloads:
        for gamma, fmt in itertools.product(gammas, formats):
            metrics = [median(values(tips, length, gamma, fmt, w, metric))
                       for metric in ["wall_s", "accepted_per_s"] for w in [1, 16]]
            md.append(f"| {tips} × {length:,} | {gamma} | {fmt} | {metrics[0]:.2f} | {metrics[1]:.2f} | "
                      f"{metrics[2]:.2f} | {metrics[3]:.2f} | {median(speeds(tips, length, gamma, fmt, 16)):.2f}× |")
    md += ["", "![Throughput and scaling](throughput_scaling.png)", "", "## Acceptance", "",
           "Each seed is counted once per workload/model, using the 1-worker packed run. "
           "Accepted counts and rejection categories match exactly across worker counts and formats for every seed.", "",
           "| Tips × sites | Gamma | Candidates | Accepted | Acceptance | Duplicate failures | Window failures |",
           "|---|---|---:|---:|---:|---:|---:|"]
    for tips, length, _ in workloads:
        for gamma in gammas:
            group = groups[(tips, length, gamma, "packed", 1)]
            totals = {name: sum(r[name] for r in group)
                      for name in ["completed", "written", "failed_duplicates", "failed_windows"]}
            md.append(f"| {tips} × {length:,} | {gamma} | {totals['completed']:,} | {totals['written']:,} | "
                      f"{100 * totals['written'] / totals['completed']:.2f}% | "
                      f"{totals['failed_duplicates']:,} | {totals['failed_windows']:,} |")
    md += ["", "## Output format and CPU use", "",
           "| Tips × sites | Workers | Packed/files wall-time gain, G8 | Packed/files wall-time gain, GC | Packed file count, G8 | Separate file count, G8 |",
           "|---|---:|---:|---:|---:|---:|"]
    for tips, length, _ in workloads:
        for w in workers:
            md.append(f"| {tips} × {length:,} | {w} | {median(gains(tips, length, 'G8', w)):.2f}× | "
                      f"{median(gains(tips, length, 'GC', w)):.2f}× | "
                      f"{median(values(tips, length, 'G8', 'packed', w, 'output_files')):,.0f} | "
                      f"{median(values(tips, length, 'G8', 'files', w, 'output_files')):,.0f} |")
    md += ["", "![Output formats](output_formats.png)", "",
           "## Interpretation", "",
           "- For 50-tip, 500-site examples, packed output is 3.61× faster than separate files at one "
           "worker with G8, falling to 1.41× at 16 workers. Its file count falls from about 4,168 to 3–18. "
           "At 16 workers it reaches about 491 accepted examples/s with G8.",
           "- For 200-tip, 5,000-site examples, packed output is not consistently faster. At 16 workers, "
           "G8 reaches about 35 accepted examples/s packed versus 43 with separate files. Both use about "
           "159 MiB of logical output; packed output reduces 193 files to 18. This grid has only 64 "
           "candidates per large run, leaving little work per worker at 16 workers.",
           "- G8 is faster than continuous gamma in every matched configuration in this grid. "
           "For large examples the median seed-paired G8 throughput advantage ranges from "
           f"{min(median([index[(200, 5000, 'G8', fmt, w, s)]['accepted_per_s'] / index[(200, 5000, 'GC', fmt, w, s)]['accepted_per_s'] for s in seeds]) for fmt, w in itertools.product(formats, workers)):.2f}× to "
           f"{max(median([index[(200, 5000, 'G8', fmt, w, s)]['accepted_per_s'] / index[(200, 5000, 'GC', fmt, w, s)]['accepted_per_s'] for s in seeds]) for fmt, w in itertools.product(formats, workers)):.2f}×. "
           "These are different rate models, not interchangeable simulation priors.",
           "- The low serial CPU use for small separate-file runs (median 0.28 cores for G8, "
           "0.38 for GC), compared with roughly one core for packed runs, is consistent with substantial "
           "I/O or scheduling waits. The CSV alone cannot isolate the cause. The apparent superlinear "
           "scaling of separate files is relative to this slow serial baseline, not evidence of superlinear CPU computation.",
           "- For large G8 packed runs, median wall time is effectively flat from 8 to 16 workers "
           "(1.83 to 1.84 s). More workers increase total CPU cost without improving median throughput "
           "for these short batches. Larger batches are needed to assess sustained scaling.", "",
           "## Scope and limitations", "",
           "- This compares this simulator's worker counts, rate models and output formats. There are "
           "no Python reference measurements in this CSV.",
           "- The supplied benchmark script uses indels, ancestral sequences, gapped and window "
           "duplicate checks, and one attempt per candidate. It requests candidate counts rather than "
           "enforcing accepted output counts. The CSV itself does not retain the full command, so exact "
           "settings should be confirmed from the per-run simulation summaries.",
           "- All rows identify `jean-zay-pp3`. The CSV does not establish exclusive CPU allocation "
           "or whether other array jobs overlapped. The current script permits four concurrent tasks "
           "(`%4`), so shared-resource contention is possible. Validate on the intended compute partition "
           "with isolated runs before treating these as production capacity estimates.",
           "- Three different seeds vary both simulation cost and the number of accepted samples. "
           "Whiskers are observed min–max ranges, not confidence intervals. Format comparisons and "
           "speedups pair the same seeds; G8 and GC may retain different numbers of small examples.",
           "- Timing includes executable startup and buffered output writes. Data enumeration and "
           "cleanup happen after timing. Logical output bytes are not physical filesystem usage or a "
           "measurement of sustained disk bandwidth. Tar rotation at the 1 GiB target is not exercised "
           "by these datasets.", "",
           "All 32 configuration summaries, including min/max timings and CPU use, are in `summary.csv`. "
           "Figures are also available as vector PDFs.", ""]
    (args.out / "REPORT.md").write_text("\n".join(md))
    print(f"Analyzed {len(rows)} runs from array {job}; excluded {len(raw) - len(rows)} earlier run(s).")
    print(f"Report, summary table and plots saved to {args.out}")


if __name__ == "__main__":
    main()
