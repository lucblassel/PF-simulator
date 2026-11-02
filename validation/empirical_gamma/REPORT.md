# Empirical Gamma sampling

On 2026-11-02, the new Rust sampler and the actual `alisim.py::sample_scale`
function each generated **20,000 independent draws**. No distribution difference
was detected: two-sample KS D = **0.0070**, below the 1% DKWM threshold
**0.01628** (upper bound on p = 0.7506). This is evidence of consistency, not a
proof of distributional equality.

Both implementations uniformly select an observation from the empirical list,
draw Normal(alpha, alpha/10), and floor the result at 0.05. The observed fraction
at the floor was 0.175% for Python and 0.215% for Rust. The complete empirical CDF,
including this atom, is included in the KS comparison. No observations are
filtered or deduplicated, and the floor is applied after the Normal draw.

The bundled 12,408 values match the source list **exactly**, including order and
multiplicity. Source pickle SHA-256:
`1527b12722751a395ee9483ee15af2f92715a25dc2560011e47d269090355ccb`.
Canonical exported values SHA-256:
`6683cac6b11364056886f3f3ea243fa85148ea1db8559e86ccd76fa816eb4a77`.

Seeds were Rust 42, Python `random` 20261102, NumPy 20261103. RNG algorithms
are different; the same numeric seed is not expected to match Python draws.
The Rust stream is deterministic for each master seed, candidate ID and attempt,
independent of worker scheduling and of the tree/AliSim RNG streams. Its reference
vector is covered by a unit test. Integration checks verify identical serial and
parallel sample files (GC and G8), fresh shapes on retries, explicit overrides,
and execution of the summary's saved reproduction command.

## Reproduce

```bash
cargo build --example sample_gamma
/path/to/legacy/.venv/bin/python scripts/validate_empirical_gamma.py \
  --sampler target/debug/examples/sample_gamma \
  --python-script /path/to/bayesnj_simulation/scripts/alisim.py
```

The legacy environment supplies NumPy and the dependencies needed to import the
actual script. Runtime simulation does not require Python. To regenerate the
bundle from the trusted legacy pickle:

```bash
python3 scripts/export_empirical_gamma.py /path/to/hogenom_alphas.txt
```

[results.json](results.json) records input, script and sampler hashes, seeds and
measurements. [python.csv.gz](python.csv.gz) and [rust.csv.gz](rust.csv.gz) retain
all draws. This tests the implemented parameter prior; the earlier
[raw MSA comparison](../python_raw/REPORT.md) independently tested AliSim outputs
under matched GC/G8 priors and is not rerun or relabeled here.

## Build verification scope

At the time of the sampling measurement, iterator work in the external
`phylotree` checkout did not compile. Initial simulator tests, Clippy and the
sampling comparison therefore used an isolated build with `phylotree` commit
`f372bd004d4e082c6c66475333ffa651897b12f3` and the existing native AliSim
library. The sampler hash above identifies that measured build.

Follow-up: the traversal iterators have since been implemented in the live
checkout. The full simulator test suite (including 21 executable integration
cases and two new doctests) now also passes against that checkout; the earlier
build blocker is resolved.
