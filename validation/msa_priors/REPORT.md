# Standalone AliSim and Python MSA-prior validation

3,000 independent MSAs per engine/case, using independent tree banks and sequence seeds. No production simulation behavior was changed.

## Findings

- **Pinned standalone 2.4.0:** no detected differences in 290 variable-feature
  comparisons across 10 cases (30,000 MSAs per backend). Maximum KS distance
  0.0433; the corrected detection threshold is 0.0636 at n=3,000.
- **Installed standalone 2.2.2.6:** 10 detected feature differences, concentrated
  in the 10-tip unrooted LG and JC cases. The [fixed-quartet follow-up](branch_probe/REPORT.md)
  is consistent with this build omitting the first terminal edge on that input.
  For T1–T2, the JC expectation is 0.24726 differences/site; observed means are
  0.24761 (bridge), 0.24702 (2.4.0), and 0.17574 (2.2.2.6).
- **Full emitted datasets do not currently match the Python pipeline.** The preset
  reproduces tree settings; empirical gamma shapes, indels, cropping and
  acceptance policies differ. The actual-function audit detects 14 feature
  differences for GC and 14 for G8. This is a measured mismatch, not merely a
  list of different defaults. No simulation behavior was changed.
- The deliberately changed Gamma prior is detected in 13 features, showing that
  this measurement can distinguish a substantive MSA-prior change. Failure to
  detect a difference elsewhere does not prove complete distributional equality.

## Follow-up: raw output with matched empirical priors

See the [implementation checklist](../python_raw/REPORT.md#implementation-checklist-for-python-parity)
for the remaining work toward Python-pipeline parity.

The [raw Python comparison](../python_raw/REPORT.md) subsequently isolates the
stage before cropping or duplicate filtering. With empirical Gamma sampling and
POW indels matched, it finds **0 / 552** detected variable-feature differences
in `.untrimmed` and unaligned tip outputs across 50/200-tip GC/G8 cases, using
both standalone versions. All 4.5 million gap-removal/reconstruction checks pass.

The full-policy audit below compared default LG without indels to the configured
Python pipeline. It mixed model-setting differences with selection differences;
its results do not identify window deduplication as the sole cause. The raw
follow-up also compares those mismatched defaults without filtering and still
finds differences. Production behavior remains unchanged.

## Engine comparisons

| Case | Standalone | Variable metric checks | Detected differences | Largest KS D |
| --- | --- | ---: | ---: | ---: |
| lg_u10 | standalone_240 | 31 | 0 | 0.0337 |
| lg_u10 | standalone_2226 | 31 | 5 | 0.1230 |
| lg_bd50 | standalone_240 | 31 | 0 | 0.0393 |
| lg_bd50 | standalone_2226 | 31 | 0 | 0.0490 |
| lg_u200 | standalone_240 | 29 | 0 | 0.0343 |
| lg_u200 | standalone_2226 | 29 | 0 | 0.0323 |
| g8_bd50 | standalone_240 | 31 | 0 | 0.0433 |
| g8_bd50 | standalone_2226 | 31 | 0 | 0.0573 |
| gc_bd50 | standalone_240 | 31 | 0 | 0.0343 |
| gc_bd50 | standalone_2226 | 31 | 0 | 0.0350 |
| invariant_bd50 | standalone_240 | 31 | 0 | 0.0353 |
| invariant_bd50 | standalone_2226 | 31 | 0 | 0.0290 |
| pow_bd50 | standalone_240 | 38 | 0 | 0.0373 |
| pow_bd50 | standalone_2226 | 38 | 0 | 0.0457 |
| gc_pow_bd50 | standalone_240 | 38 | 0 | 0.0413 |
| gc_pow_bd50 | standalone_2226 | 38 | 0 | 0.0417 |
| jc_u10 | standalone_240 | 15 | 0 | 0.0390 |
| jc_u10 | standalone_2226 | 15 | 5 | 0.1093 |
| duplicates_u4 | standalone_240 | 15 | 0 | 0.0267 |
| duplicates_u4 | standalone_2226 | 15 | 0 | 0.0557 |

## Actual Python-pipeline output audit

PF `--preset python`, default LG, 200 root sites: **2,154 / 3,000 retained
(71.80%)**, with 846 exhausted duplicate retries (28.20%). The actual Python
GC and G8 policies retain **23.37%** and **23.10%**, respectively. These Python
runs use standalone 2.4.0 to isolate policy differences from executable differences;
they do not assume which IQ-TREE version generated historical datasets.

For example, in the GC comparison's equal-sized retained ensembles, mean gap
fraction changes from **0% to 16.62%**, invariant-column fraction from **0.87% to
19.62%**, and within-MSA identity variance from **0.0228 to 0.0847**. Both policies
emit 200-column MSAs containing tips only; ancestral rows occur in Python's
intermediate files and are filtered out before the comparison.

The Python gamma prior samples one of 12,408 empirical alpha values uniformly,
adds Normal noise with standard deviation alpha/10, then floors at 0.05. Its
indel rates are 0.002/0.002 with POW{1.7/100} sizes. It simulates once with
duplicates allowed, drops full-alignment duplicates, then selects the first
unique-tip window using the actual legacy function (including its endpoint
behavior). PF instead retries whole alignments up to 20 times on each fixed tree.

| Python policy | Retained / requested | Detected metric differences | Worst KS D |
| --- | ---: | ---: | ---: |
| GC + POW indels + actual window selection | 701 / 3000 | 14 / 34 | 0.9957 |
| G8 + POW indels + actual window selection | 693 / 3000 | 14 / 34 | 0.9942 |

Negative control: changing Gamma shape 0.5 → 5 detected 13 / 31 metric differences.

## What is and is not established

- Engine tests hold the model/length/indel settings and tree prior fixed across implementations, while using independent draws. No duplicate rejection or alignment cropping is applied to these unconditional MSAs.
- Cases use LG, LG+G8{0.5}, LG+GC{0.5}, LG+I{0.2}, JC DNA, and the Python POW{1.7/100} indel law. Trees cover 4, 10, 50, and 200 tips, rooted Python BD and unrooted uniform trees.
- The pinned standalone 2.4.0 build has the original source writers/tree loading. It uses OpenMP with one runtime thread and the portable kernel. The bridge uses a sequential build with explicit per-call RNG initialization; identical numeric seeds need not produce identical sequences.
- The separately installed 2.2.2.6 binary is also tested, so version/build differences are not hidden.
- Each data point is a whole MSA. Pair distances are summarized within each MSA over fixed tip pairs; sites, taxa, and pairs are not treated as independent replicate samples.
- Measured features include lengths, gap fraction/columns/run lengths, invariant sites, entropy, between-site identity variation, tip distances, duplicate incidence/unique-tip fraction, and 20 state frequencies. Constant structural values are checked separately instead of inflating the statistical test count.
- There are 679 variable-metric comparisons and 241 constant-metric checks; all 920 count in the Bonferroni denominator.
- Statistical comparisons use two-sample KS with the equal-size DKWM bound `2 exp(-n D²)` (n >= 458); for smaller retained ensembles, the conservative union of one-sample DKW bounds `4 exp(-n D²/2)` is used. Bonferroni family alpha is 0.01 over all listed comparisons. Retained policy ensembles are truncated to equal sample counts. No detected difference is evidence at this finite resolution, not proof of equality of every aspect of the joint distribution.
- The Python-policy audit runs the actual imported `alisim.py` and `find_unduplicated.py` functions. Only process launching is instrumented to pin executable, explicit independent seeds, one thread, portable kernel and log prefix. Alpha sampling, POW indels, ancestral simulation, trimming and the downstream unique-window selection are the actual legacy code.
- The current preset only selects the Python **tree** pipeline. Its default MSA is LG without rate heterogeneity or indels. It retries duplicates on a fixed tree and retains full alignments. The Snakefile uses empirical gamma shapes (GC or configured G8), POW indels, and downstream unique-window selection/cropping to a fixed length, then filters ancestral records to tips.
- Therefore agreeing raw tree priors and conditional engines does not establish full emitted-dataset prior equality. Different retry/filter/window policies can also change which trees are retained.
- Optional user-supplied root sequences and custom model-definition files are not covered by these cases; they are not currently exposed by the bridge CLI. No universal claim over every AliSim option is made.

## Reproduce

```sh
python3 scripts/build_alisim_reference.py
cargo build --release --example validate_alisim
cargo build --release
/path/to/legacy/.venv/bin/python scripts/validate_msa_priors.py --output validation/msa_priors-new
/path/to/legacy/.venv/bin/python scripts/probe_alisim_branches.py --output validation/msa_priors-new/branch_probe
```

The validation environment needs NumPy and Biopython, plus the legacy scripts’
dependencies. Measurement helpers have six analytic checks, run with
`python -m unittest discover -s scripts -p test_validate_msa_priors.py`.

Raw per-MSA metrics/hashes are stored as gzip CSVs. See [comparisons.csv](comparisons.csv), [metadata.json](metadata.json), and [structural_constants.json](structural_constants.json).

References: [Wei & Dudley, two-sample DKWM](https://arxiv.org/abs/1107.5356), [AliSim documentation](https://iqtree.github.io/doc/AliSim).
