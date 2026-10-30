"""Summarize paired output-transport trials; Python standard library only."""

import csv
import sys
from collections import defaultdict
from pathlib import Path
from statistics import median, quantiles


root = Path(__file__).resolve().parent
groups = defaultdict(lambda: defaultdict(dict))
paths = [Path(value) for value in sys.argv[1:]] if len(sys.argv) > 1 else sorted(root.glob("output-comparison*.csv"))
for path in paths:
    with path.open() as stream:
        for row in csv.DictReader(stream):
            key = (int(row["tips"]), int(row["root_sites"]), row["case"],
                   float(row.get("branch_mean", "0.1")))
            trials = groups[key][row["mode"]]
            if int(row["round"]) in trials:
                raise ValueError("Repeated workload/round: pass only one run's CSV files to summarize.py")
            trials[int(row["round"])] = row

labels = {"LG": "LG", "JC": "JC DNA", "LG+GC{0.5}": "LG + continuous gamma",
          "indel": "LG + continuous gamma + indels", "indel_LG": "LG + indels"}
summary = []
for (tips, sites, case, branch_mean), modes in sorted(groups.items()):
    assert modes["fasta"].keys() == modes["numeric"].keys()
    fasta = [modes["fasta"][index] for index in sorted(modes["fasta"])]
    numeric = [modes["numeric"][index] for index in sorted(modes["numeric"])]
    gains = [100 * (1 - float(b["pipeline_ms"]) / float(a["pipeline_ms"]))
             for a, b in zip(fasta, numeric)]
    quartiles = quantiles(gains, n=4, method="inclusive")
    summary.append({
        "tips": tips, "root_sites": sites, "case": case, "branch_mean": branch_mean,
        "fasta_ms": median(float(row["pipeline_ms"]) for row in fasta),
        "numeric_ms": median(float(row["pipeline_ms"]) for row in numeric),
        "paired_gain_pct": median(gains), "gain_q25": quartiles[0], "gain_q75": quartiles[2],
        "fasta_check_pct": median(100 * float(row["check_ms"]) / float(row["pipeline_ms"]) for row in fasta),
        "rejection_pct": 100 * sum(int(row["rejections"]) for row in fasta) / sum(int(row["iterations"]) for row in fasta),
        "attempts_per_backend": sum(int(row["iterations"]) for row in fasta),
        "rounds": len(fasta),
    })

with (root / "summary.csv").open("w") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(summary[0]))
    writer.writeheader()
    writer.writerows(summary)

report = [
    "# FASTA versus numeric output benchmark", "",
    "Measured on 2026-10-30: Apple M2 Pro, 10 physical cores, 32 GiB RAM, macOS 27.0.1.",
    "Rust 1.98.1; Rust and native code optimized in release mode. AliSim uses one thread",
    "and the same portable likelihood kernel in both paths. No builds ran during the recorded trials.", "",
    "## Conclusion", "",
    "The current FASTA path is adequate for the existing LG + continuous-gamma workflow:",
    "the numeric prototype saves approximately 3% with indels and offers no material gain",
    "without them. Plain LG with indels is different: numeric output saves about 10% at",
    "10 tips × 500 root sites, 31% at 50 × 500, and 81% at 200 × 5,000.",
    "Rust FASTA duplicate parsing is a small component of those indel workloads; the",
    "prototype also eliminates AliSim's native indel numeric-text storage and parsing.",
    "A faster duplicate checker alone would not reproduce the large indel gains.", "",
    "Decision after measurement: numeric checking is now the production default.",
    "These timings were collected before promotion. Numeric output now has",
    "regression comparisons against the FASTA",
    "reference, including 200-tip trees. The benchmark feature enables the comparison",
    "harness; ordinary builds use numeric checking.", "",
    "## Method", "",
    "- Nine paired rounds per workload. Backend order alternates each round.",
    "- One fixed uniform tree per workload, tree seed 42; matching alignment seeds in each pair.",
    "- Warm-up and semantic equivalence checks precede timing. Batch size is calibrated",
    "  to about 200 ms per backend/round, bounded to 1–128 attempts; slow cases use one attempt.",
    "- Times include native initialization, tree reading, simulation, transport, exact",
    "  duplicate checking, buffer cleanup, and FASTA serialization of accepted numeric results.",
    "  Indel cases also serialize accepted gap-free output. Rejected numeric results are not serialized.",
    "- Tree generation, final file writes, terminal progress and optional logging are excluded.",
    "  The numbers are compute time per attempt, not total CLI or durable-storage throughput.",
    "- Indel rates: 0.002,0.002; size distributions: POW{1.7/100},POW{1.7/100}.",
    "  Continuous gamma uses shape 0.5. Default mean branch length is 0.1 substitutions/site.",
    "- Tip names and all aligned/gap-free residues match for the warm-up seed in every",
    "  workload. Duplicate decisions match for every timed paired seed. Additional tests",
    "  compare both transports across DNA/protein, gamma, indels and multiple seeds.",
    "- IQR describes the middle half of the nine paired gains; it is not a confidence interval.",
    "  Results are from one machine with no CPU pinning; sub-percent differences should",
    "  not drive an architecture decision. Memory use was not measured.", "",
    "## Normal branch lengths", "",
    "Milliseconds per attempt are medians of round averages. Positive paired gain means numeric is faster.", "",
    "| Tips × root sites | Model | FASTA ms | Numeric ms | Paired gain | Gain IQR | FASTA check share |",
    "| --- | --- | ---: | ---: | ---: | ---: | ---: |",
]
for row in summary:
    if row["branch_mean"] != 0.1:
        continue
    report.append(f"| {row['tips']} × {row['root_sites']:,} | {labels[row['case']]} | "
                  f"{row['fasta_ms']:.3f} | {row['numeric_ms']:.3f} | {row['paired_gain_pct']:+.1f}% | "
                  f"{row['gain_q25']:+.1f}…{row['gain_q75']:+.1f}% | {row['fasta_check_pct']:.2f}% |")
report.extend(["", "These ordinary workloads had no duplicate rejections among the timed seeds.", "",
               "## Duplicate-heavy workloads", "",
               "Each attempted MSA is timed once; this isolates per-attempt cost rather than running a retry loop.", "",
               "| Tips × root sites | Model | Branch mean | Rejected | FASTA ms | Numeric ms | Paired gain |",
               "| --- | --- | ---: | ---: | ---: | ---: | ---: |"])
for row in summary:
    if row["branch_mean"] == 0.1:
        continue
    report.append(f"| {row['tips']} × {row['root_sites']:,} | {labels[row['case']]} | {row['branch_mean']:g} | "
                  f"{row['rejection_pct']:.1f}% | {row['fasta_ms']:.3f} | {row['numeric_ms']:.3f} | {row['paired_gain_pct']:+.1f}% |")
report.extend(["", "## Reproduce", "", "```sh",
               "cargo build --release --features raw-benchmark --example benchmark_output -j 8",
               "./target/release/examples/benchmark_output --rounds 9 --target-ms 200",
               "./target/release/examples/benchmark_output --cases JC --workloads 4x2 --branch-mean 0.5",
               "./target/release/examples/benchmark_output --cases LG --workloads 50x500,200x5000 --branch-mean 1e-8",
               "cargo test --release --features raw-benchmark",
               "python3 benchmarks/summarize.py", "```", "",
               "The harness creates a fresh timestamped CSV by default. Pass `--output FILE` to select a new filename.",
               "The summary script reads `benchmarks/output-comparison*.csv` by default; pass explicit CSV paths",
               "to select one captured run. Its report metadata and conclusion describe the recorded M2 Pro run",
               "and should be updated when collecting results on another machine or workload matrix.",
               "Raw observations are in the `output-comparison*.csv` files; `summary.csv` is regenerated by the script.",
               "The `indel` case means LG+GC{0.5} with indels; `indel_LG` means LG with indels.", ""])
(root / "REPORT.md").write_text("\n".join(report))
print(f"Wrote {len(summary)} workload summaries to {root / 'REPORT.md'}")
