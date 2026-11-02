"""Plot matched-pipeline wall time and accepted-pair throughput."""

import argparse
import csv
import os
from pathlib import Path
import tempfile

os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="pf-pipeline-mpl-"))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "directory", type=Path, nargs="?", default=Path("benchmarks/python_pipeline")
    )
    args = parser.parse_args()
    rows = list(csv.DictReader((args.directory / "summary.csv").open()))
    cases = list(dict.fromkeys(row["case"] for row in rows))
    figure, axes = plt.subplots(
        2, len(cases), figsize=(10, 7), squeeze=False, layout="constrained"
    )
    colors = {"rust": "#117a65", "python": "#b05a25"}
    for column, case in enumerate(cases):
        group = [row for row in rows if row["case"] == case]
        for backend in ["python", "rust"]:
            points = sorted(
                (row for row in group if row["backend"] == backend),
                key=lambda row: int(row["jobs"]),
            )
            x = [int(row["jobs"]) for row in points]
            times = [float(row["median_seconds"]) for row in points]
            errors = [
                [
                    float(row["median_seconds"]) - float(row["min_seconds"])
                    for row in points
                ],
                [
                    float(row["max_seconds"]) - float(row["median_seconds"])
                    for row in points
                ],
            ]
            label = "Rust" if backend == "rust" else "Python reference"
            axes[0, column].errorbar(
                x,
                times,
                yerr=errors,
                marker="o",
                capsize=3,
                color=colors[backend],
                label=label,
            )
            axes[1, column].plot(
                x,
                [float(row["median_accepted_per_second"]) for row in points],
                marker="o",
                color=colors[backend],
                label=label,
            )
        row = group[0]
        axes[0, column].set_title(
            f"{row['samples']} pairs · {row['tips']} tips · {row['length']} root sites"
        )
        axes[0, column].set_ylabel("Wall time (seconds)")
        axes[1, column].set_ylabel("Accepted pairs / second")
        for axis in axes[:, column]:
            axis.set_xticks([1, 4, 8])
            axis.set_xlabel("Worker processes")
            axis.grid(alpha=0.2)
            axis.legend(frameon=False)
            axis.set_ylim(bottom=0)
    figure.suptitle(
        "Empirical GC + indels + ancestral rows + tip-window selection\nAll FASTA views materialized; median and range of three rounds"
    )
    figure.savefig(args.directory / "scaling.png", dpi=180)
    figure.savefig(args.directory / "scaling.pdf")


if __name__ == "__main__":
    main()
