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
`phylotree::simulation` samples trees; `src/tree_simulation.rs` applies the Python
rate/diameter pipeline and `src/tree.rs` adapts the resulting tree.
`src/alignment.rs` exposes the safe AliSim API,
and `src/main.rs` controls sample IDs and file output. C++ is isolated in
`native/alisim_bridge.cpp`; vendor headers never enter the Rust application. The
AliSim adapter accepts numeric edges and lengths, builds native IQ-TREE nodes, and
returns numeric tip rows. Newick is serialized only for accepted output trees.

## Rust and the FFI boundary

The bridge exports one `extern "C"` function. It accepts borrowed NUL-terminated
strings, fixed-width integer settings and an array of C-compatible tree edges,
returns an integer status, and writes
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

## AliSim and Rust tree priors

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

The `phylotree` crate now also implements the ETE3 deque process, actual Yule
topologies and DendroPy-compatible first-hit birth–death times. PF-simulator
applies the inherited-rate, empirical-diameter and terminal-repair steps in
the Python order. See [tree-prior definitions](tree-priors.md) and the
[distribution comparison](../validation/tree_priors/REPORT.md). Fixed tip count
alone does not specify a birth–death prior; the stopping rule, extinction
restarts, pruning and crown/stem convention are explicit.

## Direct tree input

`simulate_alignment_from_tree(&Tree, ...)` passes a numeric edge array through
FFI. Rust retains the `phylotree` tree behind a small application adapter. `native/tree_input.hpp` validates the edge
list and uses IQ-TREE's virtual `newNode` factory and `addNeighbor` methods to
construct the native tree AliSim expects. No Newick formatting, parsing or tree
file reads occur on this path. Native nodes belong to IQ-TREE and do not escape
the call. A scoped session clears input state on return or exception.

The conversion preserves the Newick reader's leaf IDs, neighbor order and root
selection, since traversal order changes RNG consumption. Topology, lengths and
sequence results are compared against file-based input with identical seeds.
Each retry builds a fresh native tree because AliSim modifies its nodes. The
accepted Rust tree is written once as a Newick output artifact. Direct input
represents rooted and unrooted binary trees with one length per
edge (rooted inputs add IQ-TREE's zero-length virtual root leaf), without partition trees, node annotations or class-specific branch lengths.

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
The CLI supplies a Rust tree directly; optional logs still use real files.
The file-based input APIs remain available for reference comparisons. The adapter is limited to the current
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

Only accepted trees and MSAs are written. On exhaustion,
`SimulationError::DuplicateSequences` identifies the sample and attempt count.
The CLI counts the failure, logs the skip when a log is enabled, advances the
progress bar, and continues. `--samples` counts requested simulations; skipped
IDs remain gaps, preserving later sample seeds. The final summary reports
written and failed counts, and duplicate exhaustion alone does not cause a
nonzero exit. This matches the old Python batch script's skip-and-report behavior.

`SimulationError::Other` wraps errors returned by AliSim or attempt logging;
these remain fatal. Tree-generation and output I/O errors also stop the run.
Failed attempts need no tree or alignment file cleanup. AliSim errors that call
`exit()` remain a separate, unrecoverable upstream behavior.

Numeric transport is enabled in ordinary builds. `native/numeric_output.hpp`
captures tip states before FASTA conversion and retains numeric indel originals.
For indels it reuses `GenomeTree::exportNewGenome` with the upstream
insertion-history update order. Rust validates row dimensions and alphabet
indices before checking duplicates. Only accepted attempts are converted to
aligned FASTA and, when requested by AliSim, gap-free FASTA.

The `raw-benchmark` feature enables the comparison harness. The original FASTA
adapter remains available through `simulate_alignment_in_memory` as a reference;
`simulate_alignment_from_tree` is the API used by the CLI;
`simulate_alignment_numeric` retains file-based tree input for comparisons. Regression comparisons
cover LG, discrete/continuous gamma, invariant sites and JC, with and without
indels, plus 200-tip trees and zero-rate/no-unaligned settings. Direct-tree comparisons
also exercise the heterotachy tree-reload path with LG+H2.
See `benchmarks/REPORT.md` for measured costs, methodology and limitations.

## Scaling and packed output

`--jobs N` uses N persistent child processes of the same executable, capped at
the requested sample count. One job keeps the original direct sequential path.
IQ-TREE uses singleton parameters, global RNGs, and error paths that can terminate
the process, so the adapter remains serialized inside each worker. OpenMP/MPI
remain disabled. The process boundary provides independent native state without
making the C++ integration reentrant.

The parent creates output directories once. One lightweight coordinator thread
per worker sends a sample ID over stdin, then reads its written/skipped result
from stdout. A shared atomic counter assigns the next ID when a worker finishes,
balancing samples with different retry costs. Workers generate trees, simulate,
check duplicates, serialize accepted results, and write distinct output files.
Only small completion messages pass through the parent's status channel; trees
and MSAs stay inside the worker. No dedicated writer or bulk-data IPC is needed.

The parent alone updates the progress bar and counts outcomes. Duplicate exhaustion
is reported as skipped and does not stop dispatch. Fatal worker errors, unexpected
exits, protocol errors, and startup failures cancel dispatch and kill/reap remaining
workers. Completed samples remain; fatal errors may leave partial in-flight writes.
With `--alisim-log`, each worker writes a separate fragment created with
`create_new`; the parent merges fragments after workers stop, including on fatal
errors. Successful merges remove fragments, and pre-existing files are preserved.
Log order depends on worker assignment; sample-file contents do not.

Each sample gets a deterministic seed derived from `(master_seed, sample_id)`
and a separate AliSim seed. `--seed` is optional: the coordinator resolves a
random seed once, then forwards it explicitly to every worker (replacing any
original seed argument). Workers cannot start without an explicit master seed.
Workers receive explicit, disjoint sample IDs so
scheduling does not affect the data. Pin compiler/library versions
for byte reproducibility; AliSim's random distributions and floating-point
results need not match across platforms. Tree sampling uses `phylotree::simulation` with a caller-seeded rand 0.8 StdRng.
SplitMix64 still derives independent sample/MSA seeds. Tree samples differ from
the previous PF-only xoshiro sampler and from Python for the same seed. Topology
frequencies, branch moments, stopping rules and cross-language ensemble
comparisons cover the new implementation.
`simulation_summary.json` records the master seed and its source, resolved CLI
arguments, executable/version/profile, replay command, requested/completed/written
counts, duplicate failure count/percentage, and monotonic wall time. It is written
before simulations and replaced atomically at completion or handled failure.
Only the parent writes this file. The timer includes setup, simulation, final
output I/O, worker startup/shutdown and log merging; final summary serialization
and the terminal completion message are excluded. Duplicate percentage uses
requested samples as its denominator; fatal-run counts cover acknowledged results.

Replay commands quote Bash arguments, spell out numeric defaults and effective
topology/branch priors, restore the working directory, and use a fresh temporary
output/log destination. Reproduction requires the same build and unchanged external
inputs. Abrupt exits can leave `status: running` with null result/timing
fields; the seed/command are still recoverable. A future production manifest can
also record vendor revisions, executable and external-input hashes, and sample
output checksums. Sample-file checks exclude summary metadata and timing.

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
