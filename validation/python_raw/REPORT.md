# Raw Python / bridge alignment comparison

**This comparison precedes duplicate rejection and window selection.** Each backend uses the actual empirical Gamma prior, LG, POW indels (rates 0.002/0.002), and 200 root sites. All attempts are retained. The measurements below used the production behavior at the time of this comparison; see the implementation update below for subsequent changes.

| Case | Python executable | Output | Variable checks | Detected differences | Largest KS D |
| --- | --- | --- | ---: | ---: | ---: |
| GC_bd50 | python_240 | aligned | 40 | 0 | 0.0397 |
| GC_bd50 | python_240 | unaligned | 29 | 0 | 0.0397 |
| GC_bd50 | python_2226 | aligned | 40 | 0 | 0.0460 |
| GC_bd50 | python_2226 | unaligned | 29 | 0 | 0.0460 |
| G8_bd50 | python_240 | aligned | 40 | 0 | 0.0350 |
| G8_bd50 | python_240 | unaligned | 29 | 0 | 0.0387 |
| G8_bd50 | python_2226 | aligned | 40 | 0 | 0.0457 |
| G8_bd50 | python_2226 | unaligned | 29 | 0 | 0.0503 |
| GC_bd200 | python_240 | aligned | 40 | 0 | 0.0367 |
| GC_bd200 | python_240 | unaligned | 29 | 0 | 0.0367 |
| GC_bd200 | python_2226 | aligned | 40 | 0 | 0.0490 |
| GC_bd200 | python_2226 | unaligned | 29 | 0 | 0.0360 |
| G8_bd200 | python_240 | aligned | 40 | 0 | 0.0390 |
| G8_bd200 | python_240 | unaligned | 29 | 0 | 0.0390 |
| G8_bd200 | python_2226 | aligned | 40 | 0 | 0.0523 |
| G8_bd200 | python_2226 | unaligned | 29 | 0 | 0.0343 |

## Findings

With matched empirical Gamma and indel priors, **none of the 552 variable-feature
comparisons detected a difference** between the bridge and the actual Python
script’s raw aligned/unaligned tip output. This covers **36,000 simulations**:
3,000 per GC/G8 × 50/200-tip case × three backends. The largest KS distance was
0.0523, below the corrected threshold of 0.0630.

All **4,500,000 tip-sequence checks** passed: removing gaps from each raw aligned
row reproduced its actual unaligned sequence exactly. Python’s additional
ancestral rows are a format difference; they are excluded by matching tip IDs.
Both standalone versions agree in these rooted birth-death cases; the earlier
2.2.2.6 unrooted-tree discrepancy is a separate, unchanged finding.

The earlier full-pipeline comparison changed both model settings and selection
policies. Its retention/statistic differences therefore cannot be assigned to
window deduplication alone. The raw default-LG baseline below quantifies
model-setting differences before any filtering. Once these settings are matched,
these experiments find no evidence of a raw-engine discrepancy; window/cropping
and retry policies remain relevant to the final emitted-dataset prior.

## What was compared

- Python: the actual `alisim.py.simulate_alignment` function with empirical alpha sampling, `--write-all`, `-mwopt`, POW indels and duplicates allowed. Read `.fa.untrimmed`, never the trimmed `.fa`, and the actual `.unaligned.fa`.
- Bridge: the production direct-tree/numeric-output functions through `validate_alisim`, with the same Gamma sampling law and indel/model settings, with no rejection or trimming. The CLI does not currently sample empirical Gamma shapes automatically; this harness supplies one sampled model string per MSA.
- There are 50 or 200 tips (`bd50`/`bd200`), using the Python tree preset. Independent tree banks, alpha draws and sequence seeds are used per backend. The raw-tree prior was validated separately; this is a comparison of the induced raw sequence prior.
- Python emits ancestors as well as tips; our output contains tips. Tip IDs are matched before comparing. Raw tip-projected alignments retain all columns, including all-gap tip columns. Source row counts and all-gap columns are recorded explicitly.
- Every unaligned tip sequence is checked, residue for residue, against the matching raw aligned row after removing gaps, for every simulation in every backend. Unaligned distributions are compared separately using lengths, composition and duplicate/unique-sequence fractions.
- Python process launch only substitutes the executable and adds an explicit seed/log prefix. It uses one thread and the standalone default CPU kernel; no portable-kernel override is added. All other arguments and the alpha sampling/trim implementation are the actual Python code.
- The bridge and standalone use different RNG initialization; matching seed values do not guarantee byte identity across executables. Distributional comparisons therefore use independent ensembles.

## Method and scope

3,000 independent whole MSAs per backend/case. KS comparisons use the equal-size DKWM upper bound `2 exp(-n D²)`, Bonferroni family alpha 0.01. Constant checks remain in the correction denominator. Sites and tips are not counted as independent replicates.

704 variable comparisons + 32 constant checks; family size 736; detection distance 0.0630.

No detected difference is evidence at this resolution, not proof of equality of all joint distributions or every possible AliSim option. Custom matrices and supplied root sequences are not covered.

## Raw defaults versus Python, without any filtering

- Default LG (no indels) versus raw Python GC, python_240, 50 tips: 19 / 38 variable features differ **before either deduplication or cropping**.
- Default LG (no indels) versus raw Python G8, python_240, 50 tips: 19 / 38 variable features differ **before either deduplication or cropping**.
- Default LG (no indels) versus raw Python GC, python_2226, 50 tips: 21 / 38 variable features differ **before either deduplication or cropping**.
- Default LG (no indels) versus raw Python G8, python_2226, 50 tips: 18 / 38 variable features differ **before either deduplication or cropping**.

This baseline reuses the earlier unconditional `lg_bd50-bridge.csv.gz` ensemble (never duplicate-filtered). It is separate from the fair, matched-prior comparisons above. Thus the earlier full-pipeline differences cannot be attributed to window deduplication alone.

## Implementation checklist for Python parity

Target the default **Snakefile pipeline** first (LG + empirical GC + POW indels),
then its configurable variants. Standalone `alisim.py` defaults differ from the
Snakefile: Gamma, indels and ancestral output are opt-in in that script.
Checked items already exist; unchecked items remain to implement or validate.

### Already available

- [x] Python-compatible birth-death trees, inherited branch-rate changes,
  empirical diameter scaling and terminal repair (validated separately).
- [x] Explicit GC/G8 model strings and POW indel parameters passed to AliSim;
  numeric duplicate checks and raw aligned/unaligned tip output.
- [x] First-attempt output without duplicate rejection via
  `--allow-duplicate-sequences`.
- [x] Raw tip-output distribution comparison under matched empirical priors,
  including exact aligned-to-unaligned reconstruction checks (this report).

### Match raw `alisim.py` generation automatically

- [x] Add empirical Gamma-shape sampling to the simulator: uniformly select an
  alpha from the same empirical values, draw Normal(alpha, alpha/10), and floor
  at 0.05. Sample once per attempt; resample alpha if an optional retry occurs.
  Support continuous GC and configurable discrete G[n], including G8.
- [x] Make `--preset python` (and `--python-bd`) select LG + empirical GC.
  An explicit `--model` overrides the preset MSA default; `--gamma GC`, `G8`
  (or another G[n]) and `none` select empirical rates or disable sampling.
- [ ] Extend the profile with default indels 0.002/0.002 and POW{1.7/100}
  insertion/deletion sizes. These remain explicit flags for now.
- [x] Make sampled model parameters reproducible across worker counts and save
  the empirical-input identity/hash and resolved MSA/selection policy in the run
  summary. Retain enough per-sample seed/parameter information to replay a sample.

Implementation update (2026-11-02): the binary embeds all 12,408 empirical
observations, preserving order and multiplicity. It samples a fresh alpha per
attempt using a dedicated deterministic stream. The summary saves source/payload
SHA-256 hashes, the resolved model/rate and selection policy, and versioned seed
recipes keyed by sample ID and attempt. Optional logs save actual attempted model
strings and engine seeds; the `sample_gamma` example reconstructs individual
alphas. Saved commands explicitly include the resolved `--gamma` mode, so replay
is independent of preset defaults. Serial/parallel identity, retry resampling,
fixed-model overrides and command replay are covered by integration checks.
The [Gamma sampling comparison](../empirical_gamma/REPORT.md) measures 20,000
Rust and 20,000 actual Python `sample_scale` draws: KS D = 0.0070, below the
1% threshold 0.01628; no difference detected. This validates the newly implemented
parameter sampler; the historical MSA measurements above are unchanged.

### Match final Snakefile dataset selection

- [ ] Add a Python selection mode: simulate **once per tree**, allowing raw
  duplicates, then drop rejected samples without retrying or replacing trees.
  Keep this policy distinct from the current maximum-attempts policy.
- [x] Reject alignments with duplicate full-length **tip** rows. Otherwise select
  the first contiguous L-column window whose aligned tip rows are unique, keeping
  gap positions in the comparison. Drop samples with no acceptable window.
- [x] Match the legacy window edge cases explicitly: for M > L the script tests
  starts in `range(M - L)`, excluding the last possible window; for M = L a
  full-unique alignment is retained at offset zero; M < L is rejected. Cover
  these cases with fixtures. Correcting the excluded endpoint would be a
  separately documented policy change, since it can change the retained prior.
- [x] Emit the selected L-column tip MSA while keeping raw full aligned and
  unaligned outputs clearly distinguishable. Unaligned output from `alisim.py`
  corresponds to the **full raw simulation**, not to the selected window.
  Preserve the raw/selected outputs requested by the user with matching sample IDs.
- [x] Report exhausted samples by their last gapped/ungapped/window rejection,
  including counts and percentages. Accepted sample metadata and optional logs
  retain the first failing check for earlier rejected attempts.
- [ ] Aggregate all rejected-attempt counts across accepted and exhausted samples
  in the run summary, in addition to its current failed-sample categories.
- [ ] Validate the completed Rust profile against the actual Python pipeline
  with matched settings: raw distributions, retained fractions, selected-window
  statistics and the distribution of **retained trees**. The raw checks here
  do not establish equivalence after selection.

Selection/output update: gapped uniqueness remains on by default; ungapped
uniqueness and the unique tip-window search are separate opt-in checks. All are
bypassed by `--allow-duplicate-sequences`. The default stores `msa_untrimmed/`
and `metadata/` with a window interval; `--materialize trimmed,unaligned` writes
both derived FASTA versions. Without a window check, trimming uses the first L
columns (or the full alignment if shorter), explicitly marked unchecked. Window
checking uses exact numeric row comparisons, with rolling hashes only as an
accelerator. `--window-policy python` selects legacy endpoints; `inclusive`
includes the final valid start. Metadata records actual accepted model/seeds and
prior attempt rejections. Failed-sample summary categories describe the last
attempt's first failing check; optional logs retain every attempt's reason.
These options do not yet implement the separate once-per-tree Python selection
profile or establish post-selection distributional parity.

### Optional script features and intermediate-file parity

These are needed for the corresponding Python options or exact intermediate
files, rather than for the default final tip-MSA distribution alone.

- [x] Add ancestral sequence export (`--ancestral`, native `--write-all`) for
  aligned/unaligned records. Exact numeric/reference-FASTA comparisons cover all
  rows for rooted/unrooted trees, with/without indels. Direct-input ancestor IDs
  match internal labels in the saved tree; selection remains tip-only. The
  [matched pipeline benchmark](../../benchmarks/python_pipeline/REPORT.md) exercises
  larger ancestral outputs and validates each selected/ungapped row.
- [ ] Offer the script's first-L-column `.fa` export alongside full
  `.fa.untrimmed` output when requested. The downstream unique-window search must
  still read the full alignment, not this initial truncation.
- [ ] Support supplied root FASTA sequences: choose a record uniformly per tree,
  keep that root record on retries, and let it determine root length as the
  Python script does. Validate this separate prior.
- [ ] Expose custom NEXUS model-definition files (`-mdef`) and the script's
  associated model-name handling; validate custom/mixture-model configurations.

Distributional parity is the target. Identical cross-executable sequences for
identical numeric seeds are not required; the RNG initialization differs.

## Files and reproduction

Raw aligned/unaligned per-MSA measurements are gzip CSVs. `*-structure.csv.gz` records original row counts. `fixtures/` preserves the first three tree/raw-alignment/unaligned files per backend/case, plus exact model strings and seeds; backend 0 is the bridge, 1 is Python + 2.4.0, 2 is Python + 2.2.2.6. Python fixture commands contain temporary paths; replace them with the retained tree and a new output prefix.

```sh
cargo build --release --example validate_alisim
/path/to/legacy/.venv/bin/python scripts/validate_raw_python.py --output validation/python_raw-new
```

See [comparisons.csv](comparisons.csv), [metadata.json](metadata.json), and [constants.json](constants.json). Reference executables are built as described in [the earlier report](../msa_priors/REPORT.md).
