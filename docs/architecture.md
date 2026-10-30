# Architecture

Use a Rust application with explicit boundaries between topology/branch
sampling, sequence simulation, and dataset storage:

```text
CLI / sample configuration / deterministic per-sample seeds
    -> tree prior -> Tree (edges, labels, lengths)
    -> AliSim adapter -> tip alignment
    -> dataset writer -> matched Newick + FASTA files
```

The implementation keeps the tree representation independent of vendor types.
`src/tree.rs` samples trees, `src/alignment.rs` exposes the safe AliSim API,
and `src/main.rs` controls sample IDs and file output. C++ is isolated in
`native/alisim_bridge.cpp`; vendor headers never enter the Rust application. The AliSim adapter reads the serialized
Newick and returns numeric tip rows to Rust before any sequence files are written.

## Rust and the FFI boundary

The bridge exports one `extern "C"` function. It accepts borrowed NUL-terminated
strings and fixed-width integer settings, returns an integer status, and writes
errors into a buffer allocated by Rust. A synchronous callback copies numeric
rows, names, the state alphabet and gap-free export metadata into Rust-owned vectors. C++ objects and
allocations remain on the C++ side; no pointers outlive a callback. The bridge catches C++ exceptions; none may cross this ABI.
The safe Rust wrapper owns the input strings through the call, validates scalar
settings, and locks a process-wide mutex before entering IQ-TREE. This protects
safe Rust callers from concurrent access to the upstream singleton and RNGs.

An indicatif progress bar counts completed pairs on stderr. The native bridge
temporarily redirects C++ `cout`/`clog` to an optional append-only log or
`/dev/null`, restoring the streams on return or exception. Rust's stderr is not
redirected, so the progress bar remains live during simulation. Native errors
on `cerr` remain visible. The CLI creates a new log once for the run and refuses
logs inside the tree/sequence directories to prevent collisions with data files.

Optional `--indel` and `--indel-size` strings pass unchanged through borrowed
C strings to the upstream parser; `--no-unaligned` forwards its output switch.
Rust does not reimplement the indel model or distribution grammar. The bridge
assembles these argument strings before creating its `argv` pointers, keeping
them valid through simulation. Each call resets AliSim's defaults through its
parser, including indel settings. Gapped alignments are retained without trimming;
the requested length refers to the root sequence, rather than final MSA width.
After accepting an MSA, the Rust writer saves aligned bytes to
`sequences/sample_000000.fa` and optional gap-free bytes to
`unaligned_sequences/sample_000000.fa`. The aligned and gap-free directories
therefore have identical filenames, paired by stem with `trees/`. The gap-free
directory is created only when AliSim actually generates a gap-free file.

`cargo build` invokes CMake through `build.rs`, discovers IQ-TREE's transitive
static archives, and bundles them with the bridge for static linking into the
Rust executable. The binary uses system C++/zlib libraries, and requires no
separate IQ-TREE executable or bridge shared library. No bindgen or C++ type
bindings are needed for this small API. Native build outputs live under Cargo's
`OUT_DIR`; changing Rust code alone does not rebuild IQ-TREE. The native build
currently supports host builds on macOS/Linux; macOS ARM64 has been tested.

See the [Cargo build script reference](https://doc.rust-lang.org/cargo/reference/build-scripts.html)
and [Rust FFI guide](https://doc.rust-lang.org/nomicon/ffi.html).

## Why link AliSim but avoid RevBayes for the first prior?

AliSim supplies the substitution machinery already used in the original workflow,
including LG and gamma models. The local IQ-TREE source supports `BUILD_LIB=ON`;
its `runAliSim(Params&, Checkpoint*)` entry point is usable through a thin adapter.
It is an internal C++ interface, so pin and test the vendor revision. Merely being
written in C++ does not make two projects independently reusable libraries.

The upstream `main` library includes its CLI entry point and supporting globals
in one translation unit. Our build renames that entry point to avoid colliding
with our own `main`, while retaining the required symbols. The vendor checkout
is unchanged; the build compiles generated copies with the stream adaptations
described below. A forced-include compatibility header supplies sequential OpenMP stubs
for unguarded calls in this revision's simulator and tree libraries. The adapter uses the
upstream argument parser to initialize defaults,
then calls `runAliSim` in process. It never invokes another executable.
It also initializes AliSim's per-thread RNG arrays explicitly, since upstream
only initializes those arrays in its OpenMP build, and selects the portable
scalar likelihood kernel instead of relying on the CLI's CPU detection.

A uniform unrooted labeled binary topology needs little code: start with the
unique three-tip topology, choose an edge uniformly, subdivide it, and attach
the next labeled tip. Repeating this samples all `(2n-5)!!` topologies uniformly:
every larger labeled tree has a unique predecessor obtained by deleting its
last tip and suppressing the resulting degree-two node. Branch lengths are drawn
independently after topology generation. The degree-three Newick serialization
root has no biological interpretation. This is distinct from a uniform rooted
prior, a Yule prior, and choosing arbitrary pairs of lineages to merge.

RevBayes becomes useful if its particular birth–death or constrained priors are
required. Add a separate adapter around its core distribution classes and convert
the resulting tree to the application representation. Keep its Rev interpreter,
DAG objects, RNG, and tree ownership inside that adapter. Match conditioning
explicitly: fixed tip count alone does not specify a birth–death distribution.
Reproducing the old Python workflow also requires its inherited branch-rate
process, empirical diameter draws, and terminal-branch adjustments.

## In-memory output and duplicate retries

Upstream `runAliSim` orchestrates the simulation. Without indels, output flows
through `AliSimulator::openOutputStream` and the sequence-writing functions in
`simulator/alisimulator.cpp`. With indels, `writeSequencesToFile` in
`main/alisim.cpp` writes the final aligned and gap-free FASTA; its
`writeSeqsFromTmpDataAndGenomeTreesIndels` helper reads temporary numeric states.
There is no public return-an-MSA API at this entry point.

`native/memory_io.cmake` generates adapted copies of just these two translation
units. It replaces uncompressed output streams in the relevant writer functions
with `MemoryOutputFile`, and the indel temporary-data reader with
`MemoryInputFile`; temporary-data cleanup erases the corresponding memory entry.
The reference FASTA mode preserves upstream serialization. Numeric mode instead
captures tip states and reconstructs indel alignments without the text round trip.
Input tree files and
optional logs still use real files. The adapter is limited to the current
single-threaded, uncompressed FASTA configuration; new compression, partition,
or output-format options would require auditing the additional I/O paths.

Source anchors must match during CMake configuration, so vendor changes that
break this adaptation fail the build. Keep the IQ-TREE revision pinned and verify
this adapter after upgrades. A `MemorySession` clears buffers on return or C++
exception. Callback inputs are borrowed; Rust reserves storage and copies the
bytes without unwinding across the ABI. Each attempt retains output and indel
working buffers in RAM, so memory use grows with alignment size. This removes
alignment disk I/O and rejected-attempt FASTA serialization; memory copies remain.

Rust's numeric `has_duplicates` compares complete aligned state rows, including
gaps, with exact hash-set equality (hash collisions cannot cause a false
duplicate). Codes exporting the same character are canonicalized in the bridge.
The FASTA checker remains available for reference comparisons. The CLI does
not trim the MSA to the requested root length, unlike the old Python indel
workflow. Only tip sequences are requested from AliSim.

`--max-attempts` defaults to 20 and includes the initial simulation. Trees and
model settings stay fixed during retries. The first attempt keeps the original
seed; subsequent attempts increment it modulo the supported seed range. Seeds
are distinct within a sample and independent of attempts consumed by other
samples. `--allow-duplicate-sequences` bypasses the check and runs once.

Only accepted MSAs are written. Exhaustion aborts the run with a nonzero exit,
removes the failed sample's input tree, and retains earlier accepted pairs.
This differs from the old Python batch script, which skipped exhausted samples
and reported them after the batch. AliSim errors that call `exit()` remain a
separate, unrecoverable upstream behavior.

Numeric transport is enabled in ordinary builds. `native/numeric_output.hpp`
captures tip states before FASTA conversion and retains numeric indel originals.
For indels it reuses `GenomeTree::exportNewGenome` with the upstream
insertion-history update order. Rust validates row dimensions and alphabet
indices before checking duplicates. Only accepted attempts are converted to
aligned FASTA and, when requested by AliSim, gap-free FASTA.

The `raw-benchmark` feature enables the comparison harness. The original FASTA
adapter remains available through `simulate_alignment_in_memory` as a reference;
`simulate_alignment_numeric` is the API used by the CLI. Regression comparisons
cover LG, discrete/continuous gamma, invariant sites and JC, with and without
indels, plus 200-tip trees and zero-rate/no-unaligned settings.
See `benchmarks/REPORT.md` for measured costs, methodology and limitations.

## Scaling and packed output

Use one isolated worker process per independent batch initially. IQ-TREE uses
singleton parameters, global RNGs, and error paths that can terminate the process;
the adapter is not a reentrant or recoverable service API. The example disables
OpenMP/MPI and calls it sequentially. Avoid concurrent simulations from threads without the Rust wrapper
until the integration has been refactored and tested for that purpose.

Each sample gets a deterministic seed derived from `(master_seed, sample_id)`
and a separate AliSim seed. Future workers should receive explicit, disjoint sample
ID ranges so scheduling does not affect the data. Pin compiler/library versions
for byte reproducibility; AliSim's random distributions and floating-point
results need not match across platforms. The Rust tree RNG uses xoshiro256** with SplitMix64 initialization, rejection
sampling for unbiased edge choices, and inverse-CDF exponential lengths. Reference
vectors and quartet/branch-length tests cover this implementation. Tree samples
therefore differ from the earlier C++/MT19937 prototype for the same seed.
A production run manifest should record
all parameters, sample seeds, vendor revisions, and output checksums.

Keep a future packed writer behind the same sample boundary. `packed_data.py`
format version 3 stores uncompressed tar shards and a SQLite index containing
pair IDs, names, shard names, original filenames, and each payload's byte offset
and size. Index IDs must be contiguous from zero; offsets point to payloads,
not tar headers. Tree/sequence pairs share a stem. The current two directories
can already be packed without changing the simulator.

For a direct writer, rotate shards at a configured pair count, write the Newick
and FASTA bytes, then record their payload offsets and sizes. Use a single writer
or independent worker shards followed by an index merge. Stage an entire dataset
and publish it only after successful completion. This avoids partially published
datasets and contention on a shared SQLite connection. The current
`simulate_alignment_in_memory` API already supplies the FASTA payload bytes
needed by a direct packed writer.

## Validation

The tests check approximately equal frequencies of the three quartet topologies,
branch lengths and degrees, matching Newick/FASTA labels and lengths, LG with
gamma rates, repeatability and sample-prefix stability, invalid arguments, and
refusal to overwrite existing output. A long JC DNA simulation compares observed
pairwise differences with `3/4 * (1 - exp(-4d/3))` for distances on the emitted tree.
Additional tests cover duplicate retry success/exhaustion, exact attempt limits,
the duplicate bypass, indels, and returning data with a nonexistent output
directory. The memory adapter has also been compared byte for byte against the
original disk backend for protein, DNA, and continuous-gamma indel runs.
These are sanity checks, not a complete scientific validation of all AliSim models.

References: [AliSim manual](https://iqtree.github.io/doc/AliSim),
[RevBayes uniform topology example](https://revbayes.github.io/tutorials/partition/).
Before distributing a linked application, check the licenses supplied with the
pinned dependencies; this example retains the upstream source and license files.
