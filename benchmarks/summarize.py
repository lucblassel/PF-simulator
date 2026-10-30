"""Summarize output and tree-input comparisons; Python standard library only.

With no arguments, use the latest CSV for each workload/comparison.
Explicit paths select a run and reject overlapping workload/round data.
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path
from statistics import median, quantiles

ROOT = Path(__file__).resolve().parent
LABELS = {
    "LG": "LG",
    "JC": "JC DNA",
    "LG+GC{0.5}": "LG + continuous gamma",
    "indel": "LG + continuous gamma + indels",
    "indel_LG": "LG + indels",
}


def summarize(paths, strict=False):
    selected = {}
    for path in paths:
        workloads = defaultdict(lambda: defaultdict(dict))
        with path.open() as stream:
            for row in csv.DictReader(stream):
                key = (
                    int(row["tips"]), int(row["root_sites"]), row["case"],
                    float(row.get("branch_mean", "0.1")),
                )
                trials = workloads[key][row["mode"]]
                index = int(row["round"])
                if index in trials:
                    raise ValueError(f"Repeated round in {path}: {key}, {index}")
                trials[index] = row
        if not workloads:
            raise ValueError(f"Empty benchmark: {path}")
        for key, modes in workloads.items():
            if set(modes) == {"fasta", "numeric"}:
                comparison, baseline, candidate = "output", "fasta", "numeric"
            elif set(modes) == {"numeric", "numeric_tree"}:
                comparison, baseline, candidate = "tree-input", "numeric", "numeric_tree"
            else:
                raise ValueError(f"Unknown or incomplete comparison in {path}: {set(modes)}")
            if modes[baseline].keys() != modes[candidate].keys():
                raise ValueError(f"Unpaired rounds in {path}: {key}")
            if len(modes[baseline]) < 3:
                raise ValueError(f"Need at least three complete rounds in {path}: {key}")
            selection = (comparison, key)
            if strict and selection in selected:
                raise ValueError(f"Repeated workload across selected files: {key}")
            # Select the entire workload, never merge rounds from different runs.
            selected[selection] = (path, modes, baseline, candidate)

    results = defaultdict(list)
    for (comparison, (tips, sites, case, branch_mean)), (path, modes, baseline, candidate) in sorted(selected.items()):
        indices = sorted(modes[baseline])
        a = [modes[baseline][index] for index in indices]
        b = [modes[candidate][index] for index in indices]
        if any((x["iterations"], x["rejections"]) != (y["iterations"], y["rejections"]) for x, y in zip(a, b)):
            raise ValueError(f"Unmatched iteration/rejection counts: {path}, {case}")
        gains = [100 * (1 - float(y["pipeline_ms"]) / float(x["pipeline_ms"])) for x, y in zip(a, b)]
        q25, _, q75 = quantiles(gains, n=4, method="inclusive")
        results[comparison].append({
            "tips": tips, "root_sites": sites, "case": case, "branch_mean": branch_mean,
            "baseline_ms": median(float(row["pipeline_ms"]) for row in a),
            "candidate_ms": median(float(row["pipeline_ms"]) for row in b),
            "paired_gain_pct": median(gains), "gain_q25": q25, "gain_q75": q75,
            "paired_saved_us": median(1000 * (float(x["pipeline_ms"]) - float(y["pipeline_ms"])) for x, y in zip(a, b)),
            "baseline_check_pct": median(100 * float(row["check_ms"]) / float(row["pipeline_ms"]) for row in a),
            "rejection_pct": 100 * sum(int(row["rejections"]) for row in a) / sum(int(row["iterations"]) for row in a),
            "attempts_per_backend": sum(int(row["iterations"]) for row in a),
            "rounds": len(a), "source": path.name,
        })
    return results


def write_csv(path, rows):
    if rows:
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)


def table(rows, baseline, candidate, tree_input=False):
    extra = "Saved µs" if tree_input else "FASTA check share"
    lines = [
        f"| Tips × root sites | Model | Branch mean | {baseline} ms | {candidate} ms | Paired gain | Gain IQR | {extra} | Rounds | Rejected |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in rows:
        extra_value = f"{row['paired_saved_us']:+.1f}" if tree_input else f"{row['baseline_check_pct']:.2f}%"
        lines.append(
            f"| {row['tips']} × {row['root_sites']:,} | {LABELS[row['case']]} | {row['branch_mean']:g} | "
            f"{row['baseline_ms']:.3f} | {row['candidate_ms']:.3f} | {row['paired_gain_pct']:+.1f}% | "
            f"{row['gain_q25']:+.1f}…{row['gain_q75']:+.1f}% | {extra_value} | {row['rounds']} | {row['rejection_pct']:.1f}% |"
        )
    return lines


def sources(rows):
    return ", ".join(f"[{name}]({name})" for name in sorted({row["source"] for row in rows}))


def main():
    explicit = bool(sys.argv[1:])
    paths = [Path(value) for value in sys.argv[1:]] if explicit else sorted(
        [*ROOT.glob("output-comparison*.csv"), *ROOT.glob("tree-input-comparison*.csv")],
        key=lambda path: (path.stat().st_mtime_ns, path.name),
    )
    if not paths:
        raise ValueError("No benchmark CSV files found")
    results = summarize(paths, strict=explicit)
    output = results["output"]
    tree = results["tree-input"]
    write_csv(ROOT / "summary.csv", output)
    write_csv(ROOT / "tree-input-summary.csv", tree)
    report = [
        "# AliSim bridge benchmarks", "",
        "Recorded on 2026-10-30 on an Apple M2 Pro (10 physical cores, 32 GiB RAM),",
        "macOS 27.0.1, Rust 1.98.1. Rust/native code use release optimization; AliSim",
        "runs on one thread with the same portable likelihood kernel in both paths.", "",
        "## What each comparison measures", "",
        "- **Tree input:** numeric output on both sides; file/Newick input versus direct",
        "  Rust edges converted to native IQ-TREE nodes. This isolates the latest bridge change.",
        "- **Output transport:** FASTA versus numeric output, with file/Newick tree input",
        "  on both sides. Rebuilding and rerunning the original harness does **not** measure",
        "  direct tree input. These results describe the earlier output optimization.", "",
        "The CLI now uses direct tree input and numeric output. Accepted trees are still",
        "written once as Newick dataset artifacts. That final output cost is excluded here.", "",
        "## Direct tree input: isolated effect", "",
    ]
    if tree:
        normal_lg = [row for row in tree if row["case"] == "LG" and row["branch_mean"] == 0.1]
        if normal_lg:
            report.append("Plain LG paired gains: " + "; ".join(
                f"**{row['paired_gain_pct']:+.1f}%** at {row['tips']} tips × {row['root_sites']:,} sites ({row['paired_saved_us']:+.0f} µs/attempt)"
                for row in normal_lg) + ".")
            report.append("")
        dna = [row for row in tree if row["case"] == "JC" and row["branch_mean"] == 0.1]
        rejected = [row for row in tree if row["case"] == "LG" and row["branch_mean"] == 1e-8]
        if dna and rejected:
            large_dna = max(dna, key=lambda row: row["tips"])
            large_rejected = max(rejected, key=lambda row: row["tips"])
            report.extend([
                f"At {large_dna['tips']} tips, JC DNA gains **{large_dna['paired_gain_pct']:.1f}%**; "
                f"the duplicate-heavy LG case gains **{large_rejected['paired_gain_pct']:.1f}%**.", "",
            ])
        gamma = [row for row in tree if row["case"] in {"indel", "LG+GC{0.5}"}]
        if gamma and all(row["gain_q25"] <= 0 <= row["gain_q75"] for row in gamma):
            report.extend([
                "All continuous-gamma cases have gain IQRs crossing zero: there is no consistent",
                "acceleration in those workloads. Direct input is a modest overhead reduction;",
                "the large earlier LG-with-indels gain came from numeric output, not tree input.", "",
            ])
        report.extend([
            "Milliseconds are medians of round averages. Positive gain/saved time means direct",
            "input is faster. Paired statistics use matching seeds within each round; the median",
            "paired gain need not equal the ratio of the two displayed median times.", "",
            *table(tree, "Numeric + file tree", "Numeric + direct tree", True), "",
            "Sources: " + sources(tree) + ".", "",
            "File input reads an already-created, usually cached Newick file. Tree generation",
            "and the original tree-file write happen outside timing. Direct input includes",
            "edge-array preparation, validation, native-node construction and cleanup. Both",
            "sides include numeric row validation/checking and FASTA serialization of accepted",
            "attempts. These are per-attempt costs, not complete dataset throughput.", "",
        ])
    else:
        report.extend(["No tree-input measurements selected yet; see the reproduction commands below.", ""])
    report.extend(["## FASTA versus numeric output", ""])
    if output:
        indels = [row for row in output if row["case"] == "indel_LG" and row["branch_mean"] == 0.1]
        if indels:
            report.extend(["LG with indels paired gains from numeric output: " + "; ".join(
                f"**{row['paired_gain_pct']:+.1f}%** at {row['tips']} × {row['root_sites']:,}"
                for row in indels) + ".", ""])
        report.extend([
            "These gains include avoiding AliSim's indel numeric-text round trip; they are",
            "not solely savings from duplicate checking. Both paths still read a tree file.", "",
            *table(output, "FASTA", "Numeric"), "",
            "Sources: " + sources(output) + ".", "",
        ])
    report.extend([
        "## Method and limits", "",
        "- One fixed uniform tree per workload (tree seed 42). Matching alignment seeds",
        "  within each pair; backend order alternates each round. Round counts are shown above.",
        "- Batches target 120–200 ms depending on invocation, bounded to 1–128 attempts.",
        "  Slow simulations use one attempt per round.",
        "- Timing includes native initialization, simulation, transport, duplicate checking,",
        "  accepted FASTA/gap-free serialization and cleanup. It excludes final filesystem",
        "  writes, tree generation, progress display and optional logging. Rejected numeric",
        "  attempts are not formatted as FASTA. Disk tree reads are included only in file paths.",
        "- Indels use rates 0.002,0.002 and POW{1.7/100},POW{1.7/100} sizes; continuous",
        "  gamma shape is 0.5. Branch mean is 0.1 except where specified in the tables.",
        "- The corrected tree-input harness compares aligned/gap-free warm-up sequences across",
        "  FASTA, numeric/file and numeric/direct paths. Timed pairs assert identical duplicate",
        "  decisions for every seed. The earlier output runs compare FASTA and numeric/file.",
        "  Duplicate-heavy workloads time individual attempts, not full retry loops.",
        "- Gain IQR is the middle half of paired gains, not a confidence interval. Positive",
        "  values mean the candidate is faster; a range crossing zero indicates inconsistent",
        "  gains across rounds. No CPU pinning or memory measurements; small differences",
        "  should be treated cautiously. Cross-run timing changes do not isolate a code change.", "",
        "## Reproduce", "",
        "Build once before timing; do not run builds or competing benchmarks during trials.", "",
        "```sh",
        "cargo build --release --features raw-benchmark --example benchmark_output -j 8",
        "# Isolate direct tree input (numeric output on both sides).",
        "./target/release/examples/benchmark_output --comparison tree-input --rounds 9 --target-ms 200",
        "./target/release/examples/benchmark_output --comparison tree-input --cases JC --workloads 4x2 --branch-mean 0.5",
        "./target/release/examples/benchmark_output --comparison tree-input --cases LG --workloads 50x500,200x5000 --branch-mean 1e-8",
        "# Original output-transport comparison (file input on both sides).",
        "./target/release/examples/benchmark_output --comparison output --rounds 9 --target-ms 200",
        "python3 benchmarks/summarize.py", "```", "",
        "The harness creates a timestamped CSV; `--output FILE` selects a new filename.",
        "With no arguments the summary script selects the latest CSV for each workload and",
        "comparison, using modification time; it never merges rounds from different files.",
        "It writes this report plus `summary.csv` (output) and `tree-input-summary.csv`.",
        "Pass explicit CSV paths to select a particular run; overlapping workloads are rejected.",
        "The `indel` case is LG+GC{0.5} with indels; `indel_LG` is LG with indels.", "",
    ])
    (ROOT / "REPORT.md").write_text("\n".join(report))
    print(f"Wrote {len(output)} output and {len(tree)} tree-input summaries to {ROOT / 'REPORT.md'}")


if __name__ == "__main__":
    main()
