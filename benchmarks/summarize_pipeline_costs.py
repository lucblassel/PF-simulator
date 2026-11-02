"""Summarize the fixed-input pipeline-cost ablation (no simulation changes)."""

import argparse
import collections
import csv
from pathlib import Path
import statistics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("benchmarks/pipeline_costs"))
    args = parser.parse_args()
    rows = list(csv.DictReader((args.input / "runs.csv").open()))
    grouped = collections.defaultdict(list)
    for row in rows:
        grouped[
            (int(row["tips"]), row["gamma"], row["ancestral"], int(row["round"]))
        ].append(row)
    phases = [
        "tree_s",
        "engine_s",
        "gapped_s",
        "window_s",
        "full_format_s",
        "extra_format_s",
        "write_all_s",
        "standalone_s",
    ]
    rounds = {}
    for key, batch in grouped.items():
        rounds[key] = {
            phase: statistics.mean(float(r[phase]) for r in batch) for phase in phases
        }
    keys = sorted({key[:3] for key in grouped})
    summary = []
    for tips, gamma, ancestral in keys:
        batches = [v for k, v in rounds.items() if k[:3] == (tips, gamma, ancestral)]
        assert len(batches) == 3, (tips, gamma, ancestral, len(batches))
        summary.append(
            dict(
                tips=tips,
                gamma=gamma,
                ancestral=ancestral,
                **{
                    phase: statistics.median(b[phase] for b in batches)
                    for phase in phases
                },
            )
        )
    with (args.input / "summary.csv").open("w") as file:
        writer = csv.DictWriter(file, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)
    lines = [
        "# Discrete versus continuous Gamma"
        if any(r["gamma"] not in {"none", "GC"} for r in summary)
        else "# Why the matched-pipeline speedup is smaller",
        "",
        "The dominant extra cost is continuous Gamma inside AliSim. Window selection and ancestral export are much smaller in these measurements. Both backends pay the expensive native-engine cost, reducing the relative benefit of avoiding Python startup, tree-file input and FASTA parsing.",
        "",
        "Single native thread, Python-prior trees, POW indels 0.002/0.002 throughout. Each cell is the median of three round means, in milliseconds per candidate. Each round uses 16 candidates at 50 tips × 500 sites, or four at 200 tips × 5,000 sites. Trees, native seeds and empirical alpha draws are paired across variants; Gamma changes the simulated data, as intended.",
        "",
        "| Tips × root sites | Model | Ancestors | Linked engine | Gapped check | Gapped + window | Full FASTA formatting | Extra FASTA formatting | Write all three | Standalone AliSim |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in summary:
        values = [f"{row[p] * 1000:.3f}" for p in phases[1:]]
        lines.append(
            f"| {row['tips']} × {500 if row['tips'] == 50 else 5000} | {'LG' if row['gamma'] == 'none' else 'LG + empirical ' + row['gamma']} | {row['ancestral']} | "
            + " | ".join(values)
            + " |"
        )
    lines += ["", "## Interpretation", ""]
    for tips in [50, 200]:
        fixed = next(
            r
            for r in summary
            if r["tips"] == tips and r["gamma"] == "none" and r["ancestral"] == "false"
        )
        gc = next(
            r
            for r in summary
            if r["tips"] == tips and r["gamma"] == "GC" and r["ancestral"] == "false"
        )
        anc = next(
            r
            for r in summary
            if r["tips"] == tips and r["gamma"] == "GC" and r["ancestral"] == "true"
        )
        lines += [
            f"- **{tips} tips:** continuous Gamma increases linked-engine time by {gc['engine_s'] / fixed['engine_s']:.1f}×; standalone AliSim increases by {gc['standalone_s'] / fixed['standalone_s']:.1f}×. With Gamma and ancestors, the combined gapped/window check costs {anc['window_s'] * 1000:.3f} ms, versus {anc['engine_s'] * 1000:.1f} ms in the engine ({100 * anc['window_s'] / anc['engine_s']:.3f}%)."
        ]
    discrete = sorted({r["gamma"] for r in summary if r["gamma"] not in {"none", "GC"}})
    if discrete:
        lines += [
            "",
            "## Discrete Gamma follow-up",
            "",
            "Discrete Gamma uses a finite set of rate categories; continuous Gamma gives sites continuously varying rates. The shape draw is paired across G4/G8/GC, but these are different site-rate distributions, so changing GC to G8 changes the simulation prior.",
            "",
        ]
        for tips in [50, 200]:
            gc = next(
                r
                for r in summary
                if r["tips"] == tips and r["gamma"] == "GC" and r["ancestral"] == "true"
            )
            for mode in discrete:
                row = next(
                    r
                    for r in summary
                    if r["tips"] == tips
                    and r["gamma"] == mode
                    and r["ancestral"] == "true"
                )
                lines.append(
                    f"- **{tips} tips, ancestors exported:** {mode} takes {row['engine_s'] * 1000:.2f} ms in the linked engine versus {gc['engine_s'] * 1000:.2f} ms for GC ({gc['engine_s'] / row['engine_s']:.1f}× faster). Standalone takes {row['standalone_s'] * 1000:.2f} versus {gc['standalone_s'] * 1000:.2f} ms."
                )
    lines += [
        "",
        "The earlier ~3× result used LG without Gamma and without indels. The later near-parity result also added continuous Gamma, ancestors, indels and materialized outputs. Those reports compare different workloads. This ablation holds indels fixed and isolates the large Gamma cost; it does not assign every historical timing difference to a single change.",
        "",
        "AliSim already simulates ancestral states while propagating sequences down a tree. Exporting them adds retained rows, gap reconstruction, transport and output; it does not add an independent evolutionary simulation at each internal node.",
        "",
        "The native implementation explicitly excludes continuous Gamma from its transition-matrix and substitution-rate caching paths (`dependencies/iqtree2/simulator/alisimulatorheterogeneity.cpp`). The timings demonstrate the engine cost; identifying which native subroutine accounts for it would require a native profiler.",
        "",
        "## Method and limits",
        "",
        "- Each candidate has exactly one attempt. No retries or replacement. All candidates are formatted and written even if selection fails, preventing rejection rates from changing the measured I/O workload. This deliberately differs from production, which skips rejected output.",
        "- The gapped-only and gapped-plus-window checks run on the identical numeric rows. Each is repeated ten times and divided by ten for timing stability. The window column includes structural validation and the gapped check; it is not an additional cost to add to the gapped column. Rejected full duplicates exit before searching windows.",
        "- The linked-engine column includes bridge setup and numeric transport, but excludes tree generation, selection and FASTA work. Full FASTA formatting is measured separately from formatting the trimmed and ungapped views. The write column covers all three byte buffers together; buffered overwrites of reused temporary files, no fsync. This isolates byte formatting/writing and does not reproduce production directory/metadata creation costs.",
        "- Standalone uses the same generated tree (through Newick), exact model string, length and native seed. It includes native process startup, FASTA outputs and native auxiliary files. It is an engine/workflow control, not a new end-to-end timing of the Python script. Shared model/seed does not independently certify identical output files.",
        "- Three independent master seeds, 20261102–20261104. Variant order and standalone/linked order alternate by round. No other benchmark runs concurrently. No warmup is subtracted. Small negative-looking ancestry differences reflect timing noise and should not be interpreted as ancestry making simulation faster.",
        "- These short, serial phase measurements diagnose the bottleneck. The larger worker-scaling results remain in the [matched pipeline report](../python_pipeline/REPORT.md). No prior, acceptance rule or production behavior was changed.",
        "",
        "## Reproduce",
        "",
        "```bash",
        "cargo build --release --example benchmark_pipeline_costs",
        "target/release/examples/benchmark_pipeline_costs \\",
        "  --output benchmarks/pipeline_costs-new \\",
        "  --reference target/alisim-reference/build/iqtree/iqtree2 --gammas "
        + ",".join(dict.fromkeys(r["gamma"] for r in rows)),
        "python3 benchmarks/summarize_pipeline_costs.py --input benchmarks/pipeline_costs-new",
        "```",
        "",
        "[Raw timings](runs.csv), [summary](summary.csv), [provenance](metadata.json).",
    ]
    (args.input / "REPORT.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
