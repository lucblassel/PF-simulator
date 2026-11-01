# Tree simulation priors

PF-simulator now delegates topology and independent branch-length sampling to
`phylotree::simulation`. Its `Tree` is retained in Rust; a small adapter maps
named tips and the root to IQ-TREE's numeric edge interface. Neither direction
uses Newick as an intermediate. Both rooted and unrooted binary trees are
supported. Newick remains a final output artifact.

## CLI examples

```sh
# Existing default: equal-probability unrooted labeled topologies, Exp(mean=0.1).
pf-sim --out uniform --preset default --samples 100 --tips 50
# ETE3's actual populate algorithm, not Yule; independent exponential edges.
pf-sim --out ete --tree-prior ete3 --branch-mean 0.1
# Gamma mean = shape * scale = 0.1.
pf-sim --out gamma --tree-prior uniform --branch-prior gamma --gamma-shape 2 --gamma-scale 0.05
# Full default SimulateAlterRescale.py birth-death transformation.
pf-sim --out python-prior --preset python --samples 100 --tips 50
# Its alternate use_bl=True variant.
pf-sim --out python-bl --preset python --use-branch-times
# Draw only trees, with no alignment simulation.
pf-sim --out bd-times --tree-prior birth-death --trees-only
# Independently enable empirical diameter scaling for an ETE3 tree.
pf-sim --out ete-scaled --tree-prior ete3 --empirical-diameter
```

`-h` shows common settings, including `--preset`; `--help` includes all tree,
rate-process, diameter, indel, and output options. `--preset default` is used
when omitted and permits individual topology/branch-prior choices.

`--tree-prior` accepts `uniform` (unrooted PDA), `uniform-rooted` (rooted PDA),
`ete3`, `yule`, and `birth-death`. The uniform priors are over **labeled** binary
topologies, not equally weighted unlabeled shapes. ETE3's deque-end coin flip is
not the uniform-extant-leaf split used by Yule. Their four-tip balanced-tree
probabilities differ: rooted PDA 1/5, ETE3 1/2, Yule 1/3.

`--branch-prior` defaults to `process-time` for birth-death and `exponential`
otherwise. `--branch-mean` is an exponential mean, while `--gamma-shape` and
`--gamma-scale` use the shape/scale convention. Explicit independent branch
priors replace the time lengths after birth-death topology generation; that
joint branch prior is different from retaining the process's times.

Birth-death starts with one lineage, has default birth/death rates 1.0/0.5,
restarts complete extinctions, and stops immediately on first reaching the
requested number of living tips. Extinct lineages are removed and unary nodes
suppressed. It does **not** use GSA, condition on a crown age, or add another
waiting interval at the requested count. These conventions match the actual
DendroPy call in the Python script. `--birth-rate` and `--death-rate` override
rates explicitly. The root stem is excluded from branch transformations.

Raw first-hit BD has exactly two zero-length newborn terminal edges, so those
two tips have identical sequences. To retain that raw prior when simulating
MSAs, use `--allow-duplicate-sequences`. The Python-compatible pipeline repairs
short tips as described below; do not silently repair raw process times.

## Exact Python transformation order

`--preset python` (or the compatibility flag `--python-bd`) selects birth-death
process times, `--python-rescaling`, empirical
diameter scaling and the terminal repair. Incompatible topology/branch overrides
are rejected. The corresponding explicit pipeline is:

1. Draw an empirical mean, then `Normal(mean, mean/10)` and clamp the target
   below at 0.02. The candidate means are strictly between 0.02 and 15; the
   jittered target may exceed 15, just as in Python.
2. Generate the first-hit BD tree at fixed tip count.
3. Starting with root rate 1, inherit the final parent rate down the tree.
   Competing exponential waiting times have **means** 0.03 and 1.0; multiply
   the current rate by `LogNormal(0, sigma)` at events, with log standard
   deviations 0.1 and 0.2. Integrate rate over each edge.
4. With default `use_bl=False`, run events over the mean original edge time,
   then scale integration by original edge time / mean time. With
   `--use-branch-times`, use each original edge time instead.
5. Rescale all altered edges to the sampled target diameter.
6. For terminal edges below 0.001, draw `Normal(0.001, 0.005)` repeatedly until
   the draw is at least 0.001. This is a truncated normal, not a simple clamp.
   **Do not rescale again.** The final diameter can exceed the target.

The four process parameters are configurable through `--small-event-mean`,
`--large-event-mean`, `--small-log-sd` and `--large-log-sd`. `--minimum-terminal`
overrides the repair floor; 0 disables it. `--python-rescaling` enables inherited
rate changes and default terminal repair without implicitly enabling empirical
scaling. Other priors can use these transforms explicitly, but only the BD
preset is claimed to reproduce the Python script's BD pipeline.

## Empirical data

`data/empirical_diameters.tsv` contains the exact retained float values with
integer weights: 59,100 HOGENOM entries of weight 1 and 4,402 RAxML entries of
weight 10, for total weight 103,120. Python's `raxml_list * 10` repeats a list;
it does **not** multiply its diameters. Each source pickle's SHA256 is recorded
in the text header. Data are embedded in the executable, so simulation needs
neither Python nor external pickle files.

`--empirical-diameters FILE` enables an alternate empirical list with one
`diameter [integer_weight]` per line; blank lines/comments are ignored. Values
must already satisfy 0.02 < value < 15. Regenerate the bundled prior with:

```sh
python scripts/export_empirical_diameters.py /path/to/bayesnj_simulation/scripts
```

## Validation and reproducibility

[Distribution comparison](../validation/tree_priors/REPORT.md) compares independent
ensembles generated by this binary's tree generator with the actual Python
script. It covers raw BD, altered BD, both full BD use_bl settings and ETE3,
plus deliberately mismatched controls. Counts, topology balance, branch-time
features, inherited-rate correlation and empirical diameters are tested
separately. Tree-level samples are retained, avoiding pseudo-replication from
correlated branches within one tree.

PRNGs differ from Python and from the earlier PF-only xoshiro sampler. Seeds
reproduce files with the same implementation and dependency versions; old
seed-to-tree outputs are not preserved. Sample seed derivation and MSA retry
seed derivation remain independent across samples. The crate's old convenience
sampling APIs retain their original behavior.

The local `phylotree` checkout is currently a Cargo path dependency at
`../../Development/rust/phylotree`. It must contain the new simulation module;
use an appropriate path or released revision when moving this project.
