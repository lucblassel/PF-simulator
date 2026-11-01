"""Compare Rust tree-level distributions with the actual legacy Python simulator.

Requires numpy, ete3, dendropy, matplotlib (the bayesnj_simulation environment).
Uses equal-size two-sample DKWM bounds (n >= 458), Bonferroni family alpha=.01.
Raw tree-level samples are retained as gzip CSV; no pooled-edge pseudoreplication.
"""

import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
import math
import os
import pickle
import random
import subprocess
import tempfile
from pathlib import Path
import numpy as np
from ete3 import Tree
from dendropy.utility import GLOBAL_RNG

ROOT = Path(__file__).resolve().parents[1]
FIELDS = [
    "sample",
    "diameter",
    "total_length",
    "crown_height",
    "height_span",
    "terminal_mean",
    "terminal_cv",
    "internal_mean",
    "min_terminal",
    "zero_terminals",
    "cherries",
    "sackin",
    "colless",
    "root_balance",
    "parent_child_correlation",
    "diameter_target",
    "diameter_increase",
]


def metrics(tree, target=None):
    root = tree.get_tree_root()
    order = list(tree.traverse("preorder"))
    heights = {root: 0.0}
    terminal = []
    internal = []
    pairs = []
    leaf_heights = []
    sackin = cherries = 0
    for node in order:
        if not node.is_root():
            heights[node] = heights[node.up] + node.dist
            if node.up is not root:
                pairs.append((node.dist, node.up.dist))
        if node.is_leaf():
            terminal.append(node.dist)
            leaf_heights.append(heights[node])
            sackin += len(node.get_ancestors())
        else:
            if node is not root:
                internal.append(node.dist)
            if len(node.children) == 2 and all(c.is_leaf() for c in node.children):
                cherries += 1
    sizes = {}
    longest = {}
    diameter = 0.0
    colless = 0
    for node in reversed(order):
        counts = [sizes[c] for c in node.children]
        sizes[node] = sum(counts) if counts else 1
        if len(counts) == 2:
            colless += abs(counts[0] - counts[1])
        ds = sorted(
            [longest[c] + c.dist for c in node.children] + [0.0, 0.0], reverse=True
        )
        longest[node] = ds[0]
        diameter = max(diameter, sum(ds[:2]))
    high = max(leaf_heights)
    span = high - min(leaf_heights)
    if span < 1e-10 * max(high, 1.0):
        span = 0.0
    x = np.asarray([p[0] for p in pairs])
    y = np.asarray([p[1] for p in pairs])
    x = x - x.mean()
    y = y - y.mean()
    v = float((x * x).sum() * (y * y).sum())
    corr = float((x * y).sum() / math.sqrt(v)) if v > 0 else 0.0
    increase = diameter - target if target is not None else math.nan
    if target is not None and abs(increase) < 1e-10 * max(1.0, diameter):
        increase = 0.0
    return dict(
        zip(
            FIELDS[1:],
            [
                diameter,
                sum(terminal) + sum(internal),
                high,
                span,
                float(np.mean(terminal)),
                float(np.std(terminal) / np.mean(terminal)),
                float(np.mean(internal)),
                min(terminal),
                sum(x == 0 for x in terminal),
                cherries,
                sackin,
                colless,
                min(sizes[c] for c in root.children) / len(terminal),
                corr,
                target if target is not None else math.nan,
                increase,
            ],
        )
    )


def ks(a, b):
    a = np.sort(a)
    b = np.sort(b)
    grid = np.unique(np.concatenate((a, b)))
    return float(
        np.max(
            np.abs(
                np.searchsorted(a, grid, side="right") / len(a)
                - np.searchsorted(b, grid, side="right") / len(b)
            )
        )
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--python-script",
        type=Path,
        default=Path(
            "/Users/lucblassel/Work/bayesnj_simulation/scripts/SimulateAlterRescale.py"
        ),
    )
    parser.add_argument(
        "--rust-binary",
        type=Path,
        default=ROOT / "target/release/examples/tree_prior_samples",
    )
    parser.add_argument("--samples", type=int, default=3000)
    parser.add_argument("--tips", type=int, nargs="+", default=[10, 50, 200])
    parser.add_argument("--out", type=Path, default=ROOT / "validation/tree_priors")
    args = parser.parse_args()
    if args.samples < 458:
        parser.error("DKWM comparison requires at least 458 samples per cohort")
    spec = importlib.util.spec_from_file_location(
        "legacy_simulation", args.python_script
    )
    legacy = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(legacy)
    parent = args.python_script.parent
    h = pickle.load((parent / "hogenom_diams.txt").open("rb"))
    r = pickle.load((parent / "raxml_diams.txt").open("rb"))
    diameters = [float(v) for v in h + r * 10 if 0.02 < float(v) < 15]
    args.out.mkdir(parents=True, exist_ok=True)
    # Positive comparisons and intentionally wrong-prior controls are separate families.
    cases = [
        ("bd_raw", ["--tree-prior", "birth-death"]),
        (
            "bd_altered",
            [
                "--tree-prior",
                "birth-death",
                "--python-rescaling",
                "--minimum-terminal",
                "0",
            ],
        ),
        ("bd_python", ["--python-bd"]),
        ("bd_python_bl", ["--python-bd", "--use-branch-times"]),
        ("ete3", ["--tree-prior", "ete3", "--empirical-diameter"]),
    ]
    all_results = []
    cohorts = {}
    invariants = []
    files = []
    serialization_errors = []
    with tempfile.TemporaryDirectory(prefix="pf-prior-comparison-") as tmp:
        tmp = Path(tmp)
        for case_index, (case, flags) in enumerate(
            cases
            + [
                ("control_bd_wrong_death", ["--tree-prior", "birth-death"]),
                (
                    "control_ete_vs_yule",
                    ["--tree-prior", "yule", "--empirical-diameter"],
                ),
            ]
        ):
            tips_list = (
                [50]
                if case == "control_bd_wrong_death"
                else [4]
                if case == "control_ete_vs_yule"
                else args.tips
            )
            for tips in tips_list:
                print(
                    f"{case}: {tips} tips, {args.samples} trees per implementation",
                    flush=True,
                )
                rust_path = tmp / "rust.csv"
                subprocess.run(
                    [
                        str(args.rust_binary),
                        "--samples",
                        str(args.samples),
                        "--tips",
                        str(tips),
                        "--seed",
                        str(271828 + case_index * 100003 + tips),
                        "--output",
                        str(rust_path),
                        *flags,
                    ],
                    check=True,
                )
                with rust_path.open() as stream:
                    rust = list(csv.DictReader(stream))
                python = []
                target = [None]
                captured = [None]
                original_write = Tree.write

                def capture_write(tree, *a, **kw):
                    captured[0] = tree
                    return original_write(tree, *a, **kw)

                Tree.write = capture_write
                original_scale = legacy.sample_scale
                original_rescale = legacy.rescale_tree

                def record_scale(diams):
                    target[0] = original_scale(diams)
                    return target[0]

                legacy.sample_scale = record_scale
                if case == "bd_altered":
                    legacy.rescale_tree = lambda tree, scale: tree
                seed = 314159 + case_index * 100003 + tips
                random.seed(seed)
                np.random.seed(seed % (2**32))
                GLOBAL_RNG.seed(seed + 1)
                try:
                    for i in range(args.samples):
                        if case in {"bd_raw", "control_bd_wrong_death"}:
                            t = legacy.treesim.birth_death_tree(
                                birth_rate=1.0,
                                death_rate=0.0 if case.startswith("control") else 0.5,
                                num_extant_tips=tips,
                            )
                            tree = Tree(
                                t.as_string(schema="newick", suppress_rooting=True)
                            )
                            draw = None
                        else:
                            legacy.simulateTree(
                                i=i,
                                numtips=tips,
                                treeType="uniform"
                                if case in {"ete3", "control_ete_vs_yule"}
                                else "birth-death",
                                outdir=str(tmp),
                                verbose=False,
                                diams=diameters,
                                use_bl=case == "bd_python_bl",
                                rate_small=0.03,
                                multiplier_small=0.1,
                                rate_big=1.0,
                                multiplier_big=0.2,
                                minimum_value=0.0 if case == "bd_altered" else 0.001,
                            )
                            p = tmp / f"{i}_{tips}_tips.nwk"
                            parsed = Tree(str(p), format=1)
                            p.unlink()
                            tree = captured[0]
                            draw = None if case == "bd_altered" else target[0]
                            # Compare the stochastic prior before ETE's 6-significant-digit writer.
                            error = max(
                                abs(a.dist - b.dist) / max(abs(a.dist), 1e-300)
                                for a, b in zip(
                                    list(tree.traverse("preorder"))[1:],
                                    list(parsed.traverse("preorder"))[1:],
                                )
                            )
                            if error > 5.1e-6:
                                raise AssertionError(
                                    "Unexpected Newick rounding discrepancy"
                                )
                            serialization_errors.append(error)
                        row = metrics(tree, draw)
                        row["sample"] = i
                        python.append(row)
                finally:
                    legacy.sample_scale = original_scale
                    legacy.rescale_tree = original_rescale
                    Tree.write = original_write
                for implementation, rows in [("rust", rust), ("python", python)]:
                    path = args.out / f"{case}-{tips}-{implementation}.csv.gz"
                    files.append(path.name)
                    with gzip.open(path, "wt", newline="") as stream:
                        writer = csv.DictWriter(
                            stream, fieldnames=FIELDS, lineterminator="\n"
                        )
                        writer.writeheader()
                        writer.writerows(rows)
                cohorts[case, tips] = (rust, python)
                if case == "bd_raw":
                    for implementation, rows in [("rust", rust), ("python", python)]:
                        invariants.append(
                            {
                                "case": case,
                                "tips": tips,
                                "implementation": implementation,
                                "two_newborn_zero_edges": all(
                                    float(row["zero_terminals"]) == 2 for row in rows
                                ),
                                "ultrametric": all(
                                    float(row["height_span"]) == 0 for row in rows
                                ),
                            }
                        )
                if case in {"bd_python", "bd_python_bl"}:
                    for implementation, rows in [("rust", rust), ("python", python)]:
                        invariants.append(
                            {
                                "case": case,
                                "tips": tips,
                                "implementation": implementation,
                                "terminal_floor": all(
                                    float(row["min_terminal"]) >= 0.001 - 1e-9
                                    for row in rows
                                ),
                                "diameter_not_shrunk_by_floor": all(
                                    float(row["diameter_increase"]) >= -1e-4
                                    for row in rows
                                ),
                            }
                        )
                for field in FIELDS[1:]:
                    a = np.array([float(row[field]) for row in rust])
                    b = np.array([float(row[field]) for row in python])
                    if not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
                        continue
                    d = ks(a, b)
                    p_bound = min(1.0, 2 * math.exp(-args.samples * d * d))
                    all_results.append(
                        dict(
                            case=case,
                            tips=tips,
                            metric=field,
                            n=args.samples,
                            rust_mean=float(a.mean()),
                            python_mean=float(b.mean()),
                            rust_median=float(np.median(a)),
                            python_median=float(np.median(b)),
                            ks_d=d,
                            p_upper_bound=p_bound,
                        )
                    )
    main_tests = [r for r in all_results if not r["case"].startswith("control")]
    m = len(main_tests)
    cutoff = 0.01 / m
    threshold = math.sqrt(math.log(2 * m / 0.01) / args.samples)
    for row in all_results:
        row["family_alpha"] = 0.01
        row["bonferroni_tests"] = m
        row["detected_difference"] = row["p_upper_bound"] < cutoff
    with (args.out / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=list(all_results[0]), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(all_results)
    versions = {}
    import ete3
    import dendropy

    versions.update(
        ete3=ete3.__version__, dendropy=dendropy.__version__, numpy=np.__version__
    )
    metadata = dict(
        rust_binary_sha256=hashlib.sha256(args.rust_binary.read_bytes()).hexdigest(),
        samples_per_case=args.samples,
        tips=args.tips,
        python_script=str(args.python_script),
        python_script_sha256=hashlib.sha256(
            args.python_script.read_bytes()
        ).hexdigest(),
        versions=versions,
        number_of_tests=m,
        family_alpha=0.01,
        ks_threshold=threshold,
        invariants=invariants,
        raw_files=files,
        rust_seed_base=271828,
        python_seed_base=314159,
        max_python_serialization_relative_error=max(serialization_errors, default=0.0),
    )
    (args.out / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    failures = [r for r in main_tests if r["detected_difference"]]
    report = [
        "# Tree-prior distribution comparison",
        "",
        f"{args.samples:,} independent trees per implementation per case/tip count; {m} tree-level comparisons.",
        f"Family-wise alpha 0.01; two-sample KS/DKWM Bonferroni critical distance **{threshold:.4f}**.",
        "",
        f"**{len(failures)} detected differences** in the intended matching-prior comparisons.",
        "",
        "This is evidence of agreement at the measured resolution, not proof of equality of all tree distributions.",
        "Prior agreement also relies on the explicit construction/conditioning rules and structural tests below.",
        "",
        "| Pipeline | Tips | Largest KS D | Metric | Detected differences |",
        "| --- | ---: | ---: | --- | ---: |",
    ]
    for case, _ in cases:
        for tips in args.tips:
            rows = [r for r in main_tests if r["case"] == case and r["tips"] == tips]
            worst = max(rows, key=lambda r: r["ks_d"])
            report.append(
                f"| {case} | {tips} | {worst['ks_d']:.4f} | {worst['metric']} | {sum(r['detected_difference'] for r in rows)}/{len(rows)} |"
            )
    report.extend(
        [
            "",
            "## Negative controls",
            "",
            "| Deliberate mismatch | Metric | KS D | Detected |",
            "| --- | --- | ---: | --- |",
        ]
    )
    for case, metric in [
        ("control_bd_wrong_death", "crown_height"),
        ("control_ete_vs_yule", "root_balance"),
    ]:
        row = next(
            r for r in all_results if r["case"] == case and r["metric"] == metric
        )
        report.append(
            f"| {case} | {metric} | {row['ks_d']:.4f} | {row['detected_difference']} |"
        )
    report.extend(
        [
            "",
            "## Priors and transformation order",
            "",
            "- Uniform unrooted: labeled PDA, probability 1/(2n-5)!!; uniform rooted additionally supported.",
            "- ETE3 3.1.3: coin flip between deque ends, not uniformly choosing a leaf. Yule is a separate sampler.",
            "- BD: start with one lineage, birth=1, death=0.5; restart extinct realizations, stop at the first n extant tips, prune extinct lineages and suppress unary nodes. No fixed crown age, GSA or extra terminal waiting time.",
            "- Raw BD leaves are contemporaneous and exactly two newborn terminal edges are zero. Both implementations satisfy these invariants.",
            "- bd_altered: the actual Python rate-change loop with use_bl=False, before diameter scaling and without terminal repair.",
            "- bd_python: actual SimulateAlterRescale.py default BD transformation. bd_python_bl also tests use_bl=True.",
            "- Rate events use exponential means 0.03 and 1.0; multipliers are LogNormal(0, sigma) with sigma=0.1 and 0.2. Descendants inherit the final parent rate. With use_bl=False, events run over the mean original branch time and integration is scaled back by original/mean.",
            "- Empirical prior: filtered HOGENOM values plus ten replicated copies of filtered RAxML values, not tenfold diameter values; strict 0.02 < mean < 15. Draw N(mean, mean/10), clamp below at 0.02.",
            "- Rescale altered edges to this target, then replace shorter BD terminal edges by N(0.001,0.005), rejecting draws below 0.001. This can increase the final diameter; do not rescale again. ETE3 does not apply this terminal repair.",
            "",
            "## What was compared",
            "",
            "All metrics are per tree: diameter, total length, crown height, height span, terminal mean/CV/minimum/zero count, internal mean, cherries, Sackin, Colless, root split balance and parent-child edge correlation. Empirical cases additionally compare sampled diameter targets and final-minus-target increments. Edges within a tree are not treated as independent samples.",
            "",
            f"Python versions: {versions}. The actual legacy script is imported; full-pipeline trees are captured immediately before its real Newick writer. The written files are also parsed and checked for the expected 6-significant-digit rounding (maximum relative error recorded in metadata). Rust uses the same TreeSimulation generator called by pf-sim, via examples/tree_prior_samples.rs.",
            "",
            "Near-zero height spans and diameter increments are coalesced at relative/absolute tolerance 1e-10 to avoid treating floating-point roundoff as a prior difference. Rust emits full-precision lengths; ETE3 writes only 6 significant digits, so the stochastic-prior comparison uses pre-serialization lengths on both sides.",
            "",
            "The two-sample DKWM bound p <= 2 exp(-n D²) applies to equal cohort sizes n >= 458; Bonferroni controls the declared family regardless of correlations among metrics. Ties make the distribution-free comparison conservative. Negative controls demonstrate sensitivity to distinct birth/death parameters and to substituting Yule for ETE3.",
            "",
            "[Two-sample DKWM reference](https://arxiv.org/abs/1107.5356), [ETE3 populate source](https://github.com/etetoolkit/ete/blob/3.1.3/ete3/coretype/tree.py), [DendroPy birth-death source](https://jeetsukumaran.github.io/DendroPy/_modules/dendropy/model/birthdeath.html).",
            "",
            "The PRNGs differ, so matching seeds do not produce identical trees. This comparison uses independent seeded ensembles. Finite marginal checks do not establish equality of every joint feature or exclude smaller effects.",
            "",
            "![Tree-prior empirical CDFs](comparison.png)",
            "",
            "## Reproduce",
            "",
            "```sh",
            "cargo build --release --example tree_prior_samples",
            f"/Users/lucblassel/Work/bayesnj_simulation/.venv/bin/python scripts/compare_tree_priors.py --samples {args.samples}",
            "```",
            "",
            "Raw tree-level samples: `*-rust.csv.gz` and `*-python.csv.gz`; detailed results: `summary.csv`; source hash, versions and structural checks: `metadata.json`.",
        ]
    )
    (args.out / "REPORT.md").write_text("\n".join(report) + "\n")
    os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "target/matplotlib"))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 3, figsize=(12, 7))
    for ax, (case, field) in zip(
        axes.flat,
        [
            ("bd_raw", "crown_height"),
            ("bd_raw", "root_balance"),
            ("bd_altered", "parent_child_correlation"),
            ("bd_python", "diameter"),
            ("bd_python", "terminal_cv"),
            ("bd_python", "colless"),
        ],
    ):
        for label, rows in zip(
            ["Rust", "Python"], cohorts[case, 50 if 50 in args.tips else args.tips[0]]
        ):
            values = np.sort([float(r[field]) for r in rows])
            ax.step(
                values,
                np.arange(1, len(values) + 1) / len(values),
                where="post",
                label=label,
                alpha=0.8,
            )
        ax.set_title(f"{case}: {field}", fontsize=10)
        ax.set_ylabel("Empirical CDF")
        ax.grid(alpha=0.2)
    axes.flat[0].legend()
    fig.suptitle(f"Tree priors: {args.samples:,} trees per implementation, 50 tips")
    fig.tight_layout()
    fig.savefig(args.out / "comparison.png", dpi=180)
    plt.close(fig)
    print(
        f"Report: {args.out / 'REPORT.md'}; {len(failures)} detected differences among {m} tests",
        flush=True,
    )
    if failures:
        raise SystemExit(1)
    if any(
        not all(
            v for k, v in inv.items() if k not in {"case", "tips", "implementation"}
        )
        for inv in invariants
    ):
        raise SystemExit("Structural invariant failed")


if __name__ == "__main__":
    main()
