"""Compare complete no-indel, duplicate-checked Python and Rust datasets.

Run in the legacy Python environment. The reference imports the actual tree and
MSA functions. python_threads reproduces the script's sequential tree loop with
N AliSim threads; python_pool is a benchmark-only scheduler with N processes and
one AliSim thread each. Rust uses its normal CLI. No simulator code is changed.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import random
import resource
import shlex
import statistics
import subprocess
import sys
import tempfile
import time
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
CASES = [
    ('small_50', 8, 50, 500),
    ('medium_50', 128, 50, 500),
    ('batch_50', 512, 50, 500),
    ('small_200', 8, 200, 5000),
    ('batch_200', 128, 200, 5000),
]
SETTINGS = [('rust', 1), ('python_threads', 1), ('rust', 4), ('python_threads', 4),
            ('python_pool', 4), ('rust', 8), ('python_threads', 8), ('python_pool', 8)]
TREE_MODULE = MSA_MODULE = NP = GLOBAL_RNG = None
DIAMETERS = None


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def splitmix(seed, sample):
    mask = (1 << 64) - 1
    value = (seed + 0x9e3779b97f4a7c15 * (sample + 1)) & mask
    value = ((value ^ (value >> 30)) * 0xbf58476d1ce4e5b9) & mask
    value = ((value ^ (value >> 27)) * 0x94d049bb133111eb) & mask
    return value ^ (value >> 31)


def initialize(legacy_root):
    global TREE_MODULE, MSA_MODULE, NP, GLOBAL_RNG, DIAMETERS
    import numpy
    from dendropy.utility import GLOBAL_RNG as rng
    NP, GLOBAL_RNG = numpy, rng
    root = Path(legacy_root)
    TREE_MODULE = load_module('legacy_tree_benchmark', root / 'SimulateAlterRescale.py')
    MSA_MODULE = load_module('legacy_msa_benchmark', root / 'alisim.py')
    MSA_MODULE.args = SimpleNamespace(gamma=None, allow_duplicate_sequences=False, keep_logfiles=False)
    DIAMETERS = TREE_MODULE.read_list(root / 'hogenom_diams.txt') + TREE_MODULE.read_list(root / 'raxml_diams.txt') * 10
    DIAMETERS = [x for x in DIAMETERS if 0.02 < x < 15]


def python_sample(task):
    sample_seed = splitmix(task['seed'], task['sample'])
    random.seed(sample_seed)
    NP.random.seed(sample_seed & 0xffffffff)
    GLOBAL_RNG.seed(sample_seed + 1)
    tree_dir = Path(task['output']) / 'trees'
    name = f"sample_{task['sample']:06}"
    start = time.perf_counter()
    TREE_MODULE.simulateTree(task['sample'], task['tips'], 'birth-death', str(tree_dir), False,
                             DIAMETERS, False, 0.03, 0.1, 1.0, 0.2, 0.001)
    original = tree_dir / f"{task['sample']}_{task['tips']}_tips.nwk"
    tree = tree_dir / f'{name}.nwk'
    original.rename(tree)
    tree_seconds = time.perf_counter() - start
    attempts = 0
    first_seed = splitmix(sample_seed, 0) % 2147483646
    launch_original = subprocess.Popen

    def launch(command, **kwargs):
        nonlocal attempts
        words = shlex.split(command)
        words[0] = task['reference']
        words += ['--seed', str(1 + (first_seed + attempts) % 2147483646), '-redo']
        attempts += 1
        # Preserve the actual function's shell launch and default CPU kernel.
        return launch_original(shlex.join(words), **kwargs)

    MSA_MODULE.subprocess = SimpleNamespace(Popen=launch, PIPE=subprocess.PIPE)
    start = time.perf_counter()
    failure, _ = MSA_MODULE.simulate_alignment(str(tree), 'LG', None, False, '', '', False,
                                              str(Path(task['output']) / 'sequences'), task['length'],
                                              20, task['native_threads'], '', None)
    msa_seconds = time.perf_counter() - start
    path = Path(task['output']) / 'sequences' / f'{name}.fa'
    if failure is not None:
        if failure != str(tree) or attempts != 20:
            raise RuntimeError(f'Unexpected Python simulation failure: {failure}')
        assert not path.exists()
        tree.unlink()  # Keep complete pairs only, as the Rust CLI does.
    else:
        assert path.exists() and attempts >= 1
    return dict(sample=task['sample'], written=int(failure is None), attempts=attempts,
                tree_seconds=tree_seconds, msa_seconds=msa_seconds)


def reference_run(args):
    """Fresh reference process, including imports/worker startup in outer timing."""
    args.data.mkdir()
    (args.data / 'trees').mkdir()
    (args.data / 'sequences').mkdir()
    tasks = [dict(sample=i, seed=args.seed, tips=args.tips, length=args.length,
                  output=str(args.data), reference=str(args.reference),
                  native_threads=args.jobs if args.backend == 'python_threads' else 1)
             for i in range(args.samples)]
    if args.backend == 'python_pool':
        with ProcessPoolExecutor(max_workers=min(args.jobs, args.samples), initializer=initialize,
                                 initargs=(str(args.legacy_root),)) as pool:
            rows = list(pool.map(python_sample, tasks, chunksize=1))
    else:
        initialize(str(args.legacy_root))
        rows = [python_sample(task) for task in tasks]
    written = sum(row['written'] for row in rows)
    summary = dict(requested_samples=args.samples, written_samples=written,
                   failed_duplicate_samples=args.samples - written,
                   attempts=sum(row['attempts'] for row in rows),
                   tree_seconds_sum=sum(row['tree_seconds'] for row in rows),
                   msa_seconds_sum=sum(row['msa_seconds'] for row in rows))
    (args.data / 'simulation_summary.json').write_text(json.dumps(summary) + '\n')


def validate_dataset(output, samples, tips, length):
    """Excluded from timing: validate actual files and all accepted tip rows."""
    from Bio import SeqIO
    digest = hashlib.sha256()
    tree_digest = hashlib.sha256()
    written = size = 0
    sequence_dir = output / ('msa_untrimmed' if (output / 'msa_untrimmed').exists() else 'sequences')
    seqs = sorted(sequence_dir.glob('*.fa'))
    trees = sorted((output / 'trees').glob('*.nwk'))
    assert {p.stem for p in seqs} == {p.stem for p in trees}
    for p in seqs:
        records = list(SeqIO.parse(p, 'fasta'))
        assert len(records) == tips
        strings = [str(record.seq) for record in records]
        assert all(len(row) == length for row in strings)
        assert len(set(strings)) == tips, f'Duplicate accepted output: {p}'
        assert {record.id for record in records} == {f'T{i}' for i in range(1, tips + 1)}
        written += 1
    for p in trees + seqs:
        content = p.read_bytes()
        digest.update(str(p.relative_to(output)).encode() + b'\0' + content)
        size += len(content)
        if p.suffix == '.nwk':
            tree_digest.update(p.name.encode() + b'\0' + content)
    assert written <= samples
    return written, digest.hexdigest(), tree_digest.hexdigest(), size


def measure(args, case, backend, jobs, seed):
    name, samples, tips, length = case
    with tempfile.TemporaryDirectory(prefix='python-pf-bench-') as temp:
        output = Path(temp) / 'data'
        if backend == 'rust':
            command = [str(args.binary), '--out', str(output), '--samples', str(samples),
                       '--tips', str(tips), '--length', str(length), '--seed', str(seed),
                       '--preset', 'python', '--jobs', str(jobs), '--model', 'LG', '--max-attempts', '20']
        else:
            command = [sys.executable, str(Path(__file__).resolve()), '--reference-run',
                       '--data', str(output), '--samples', str(samples), '--tips', str(tips),
                       '--length', str(length), '--seed', str(seed), '--backend', backend,
                       '--jobs', str(jobs), '--reference', str(args.reference),
                       '--legacy-root', str(args.legacy_root)]
        before = resource.getrusage(resource.RUSAGE_CHILDREN)
        start = time.perf_counter()
        result = subprocess.run(command, capture_output=True, text=True, timeout=600)
        elapsed = time.perf_counter() - start
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        if result.returncode:
            raise RuntimeError(f'{command}\n{result.stdout}\n{result.stderr}')
        summary = json.loads((output / 'simulation_summary.json').read_text())
        written, digest, tree_digest, size = validate_dataset(output, samples, tips, length)
        assert written == summary['written_samples']
        row = dict(case=name, backend=backend, jobs=jobs, seed=seed, requested=samples,
                   written=written, failed=samples-written, wall_seconds=elapsed,
                   cpu_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,
                   output_bytes=size, dataset_sha256=digest, trees_sha256=tree_digest,
                   attempts=summary.get('attempts', ''),
                   tree_seconds_sum=summary.get('tree_seconds_sum', ''),
                   msa_seconds_sum=summary.get('msa_seconds_sum', ''))
        return row


def quartiles(values):
    q = statistics.quantiles(values, n=4, method='inclusive')
    return q[0], q[2]


def summarize(output):
    meta = json.loads((output / 'metadata.json').read_text())
    rows = list(csv.DictReader((output / 'runs.csv').open()))
    result = []
    for case, samples, tips, length in meta['cases']:
        for backend, jobs in meta['settings']:
            group = [r for r in rows if r['case'] == case and r['backend'] == backend and int(r['jobs']) == jobs]
            assert len(group) == meta['repeats'], (case, backend, jobs)
            same_round_py = {r['round']: float(r['wall_seconds']) for r in rows
                             if r['case'] == case and r['backend'] == 'python_threads' and int(r['jobs']) == jobs}
            serial = {r['round']: float(r['wall_seconds']) for r in rows
                      if r['case'] == case and r['backend'] == backend and int(r['jobs']) == 1}
            if backend == 'python_pool':
                serial = {r['round']: float(r['wall_seconds']) for r in rows
                          if r['case'] == case and r['backend'] == 'python_threads' and int(r['jobs']) == 1}
            times = [float(r['wall_seconds']) for r in group]
            speedups = [serial[r['round']] / float(r['wall_seconds']) for r in group]
            relative = [same_round_py[r['round']] / float(r['wall_seconds']) for r in group]
            low, high = quartiles(relative)
            result.append(dict(case=case, samples=samples, tips=tips, length=length, backend=backend, jobs=jobs,
                               median_seconds=statistics.median(times),
                               speedup_vs_own_serial=statistics.median(speedups),
                               speedup_vs_python_threads=statistics.median(relative),
                               relative_speedup_q1=low, relative_speedup_q3=high,
                               accepted_per_second=statistics.median([int(r['written'])/float(r['wall_seconds']) for r in group]),
                               median_written=statistics.median([int(r['written']) for r in group]),
                               median_failed=statistics.median([int(r['failed']) for r in group])))
    with (output / 'summary.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(result[0]))
        writer.writeheader()
        writer.writerows(result)
    lines = ['# Python reference versus PF simulator: no indels, duplicate checks enabled', '',
             f"Recorded {meta['started_utc']} on {meta['platform']}; {meta['logical_cpus']} logical CPUs. "
             f"{meta['repeats']} rounds, release Rust binary, Python {meta['python']}.", '',
             '## Matched single-threaded and parallel wall times', '',
             '| Case | Pairs requested | Tips × sites | Worker/thread budget | Rust (s) | Python native threads (s) | Python pool (s) | Python native / Rust | Python pool / Rust |',
             '| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |']
    for case, samples, tips, length in meta['cases']:
        for jobs in meta['jobs']:
            def get(backend, j=jobs):
                return next(r for r in result if r['case'] == case and r['backend'] == backend and r['jobs'] == j)
            rust, py = get('rust'), get('python_threads')
            pool = py if jobs == 1 else get('python_pool')
            same_round_rust = {r['round']: float(r['wall_seconds']) for r in rows
                               if r['case'] == case and r['backend'] == 'rust' and int(r['jobs']) == jobs}
            pool_ratios = [float(r['wall_seconds']) / same_round_rust[r['round']] for r in rows
                           if r['case'] == case and r['backend'] == ('python_threads' if jobs == 1 else 'python_pool') and int(r['jobs']) == jobs]
            lines.append(f"| {case} | {samples} | {tips} × {length} | {jobs} | {rust['median_seconds']:.3f} | "
                         f"{py['median_seconds']:.3f} | {pool['median_seconds']:.3f} | "
                         f"{rust['speedup_vs_python_threads']:.2f}× | {statistics.median(pool_ratios):.2f}× |")
    lines += ['', 'Times are independent medians; ratios are medians of same-round ratios. '
              'A value greater than 1 in the last columns favors Rust. The one-thread Python result '
              'is shared as both Python baselines, not measured twice.', '',
              '## Retention and throughput', '',
              '| Case | Backend | Budget | Median written / requested | Accepted MSAs/s | Speedup vs own serial |',
              '| --- | --- | ---: | ---: | ---: | ---: |']
    for r in result:
        lines.append(f"| {r['case']} | {r['backend']} | {r['jobs']} | {r['median_written']:g} / {r['samples']} | "
                     f"{r['accepted_per_second']:.1f} | {r['speedup_vs_own_serial']:.2f}× |")
    lines += ['', '## What is timed', '',
              '- Whole-run wall time: executable/interpreter startup, imports, empirical input loading, '
              'tree generation and transformations, AliSim calls, duplicate checking/retries, output file I/O, '
              'worker startup/shutdown, and final summary. File validation, hashing and deletion are outside timing.',
              '- Both use the validated birth-death/inherited-rate/empirical-diameter/terminal-repair tree prior '
              'and LG amino acids without Gamma or indels. There is one MSA per tree, with up to 20 attempts '
              'on that fixed tree, skipping exhausted samples without replacement. No ancestors or window selection.',
              '- Rust is the unmodified production CLI with `--preset python --model LG --max-attempts 20`. '
              'It uses numeric duplicate checking, direct tree input and one native thread per worker.',
              '- Python imports the actual `SimulateAlterRescale.py.simulateTree` and '
              '`alisim.py.simulate_alignment` functions. The reference combines both stages in one interpreter '
              'rather than paying Snakemake and two separate interpreter startups; this is a conservative reference.',
              '- `python_threads` mirrors the script: sequential per-tree work and N native AliSim threads '
              '(the script’s `--processes N`). `python_pool` is a benchmark-only parallel scheduler over '
              'the same functions, with N persistent Python workers and one AliSim thread each. '
              'The original script has no such pool. Rust uses N persistent process workers for N > 1.',
              '- Python retains its normal shell launch, per-attempt standalone startup, FASTA I/O and parsing '
              'for duplicate checks, and default log generation/deletion. Launch instrumentation only pins the '
              'standalone executable, deterministic attempt seeds and overwrite behavior. Default CPU kernels '
              'are used. No production Python source is modified.',
              '- Python uses standalone IQ-TREE 2.4.0 from the pinned source, matching the linked library version. '
              'The older installed 2.2.2.6 is excluded to avoid mixing version/build differences with bridge speed. '
              'Standalone OpenMP and linked sequential builds remain different builds.',
              '- Per-sample Python/NumPy/DendroPy RNG seeds are reset to make Python pool scheduling deterministic. '
              'Python and Rust tree/RNG algorithms differ, so datasets and retry counts need not match across '
              'languages. These are matched prior/policy comparisons, not identical sequence trajectories. '
              'Written counts and accepted throughput are reported rather than assuming identical work.',
              '- Rejected Python tree files are removed so both outputs contain complete accepted pairs only. '
              'The original tree script would leave these trees on disk. Python attempt counts and summed '
              'tree/MSA function times are in runs.csv; Rust does not expose an attempt total.',
              '- After each run, every accepted MSA is checked for tip count/IDs, length and unique sequences, '
              'and tree/MSA sample IDs must match. Rust hashes must match across worker counts; Python '
              'one-thread/pool hashes must match across worker counts. Multi-thread AliSim streams may differ.',
              '- Runs are sequential, with rotating case/backend order, one unmeasured warmup per setting, '
              'normal buffered filesystem I/O and no fsync. This is warm local-filesystem CLI completion time, '
              'not durable or network-storage throughput. Small batches are dominated by startup costs.', '',
              '## Reproduce', '', '```sh', 'cargo build --release',
              '/path/to/legacy/.venv/bin/python benchmarks/benchmark_python.py --output benchmarks/python-new',
              '```', '', 'Raw measurements: [runs.csv](runs.csv), [summary.csv](summary.csv), '
              '[metadata.json](metadata.json). Raw distribution validation: '
              '[report](../../validation/python_raw/REPORT.md).', '']
    (output / 'REPORT.md').write_text('\n'.join(lines))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'benchmarks/python')
    parser.add_argument('--binary', type=Path, default=ROOT / 'target/release/pf-sim')
    parser.add_argument('--reference', type=Path, default=ROOT / 'target/alisim-reference/build/iqtree/iqtree2')
    parser.add_argument('--legacy-root', type=Path, default=Path('/Users/lucblassel/Work/bayesnj_simulation/scripts'))
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--cases', nargs='+', choices=[case[0] for case in CASES])
    parser.add_argument('--summarize-only', action='store_true')
    parser.add_argument('--reference-run', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--data', type=Path, help=argparse.SUPPRESS)
    parser.add_argument('--samples', type=int, help=argparse.SUPPRESS)
    parser.add_argument('--tips', type=int, help=argparse.SUPPRESS)
    parser.add_argument('--length', type=int, help=argparse.SUPPRESS)
    parser.add_argument('--seed', type=int, help=argparse.SUPPRESS)
    parser.add_argument('--backend', choices=['python_threads', 'python_pool'], help=argparse.SUPPRESS)
    parser.add_argument('--jobs', type=int, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.reference_run:
        reference_run(args)
        return
    if args.summarize_only:
        summarize(args.output)
        return
    if args.repeats < 2:
        parser.error('--repeats must be >=2')
    cases = [case for case in CASES if not args.cases or case[0] in args.cases]
    args.output.mkdir(parents=True, exist_ok=False)
    metadata = dict(started_utc=datetime.now(timezone.utc).isoformat(), platform=platform.platform(),
                    logical_cpus=os.cpu_count(), python=platform.python_version(), repeats=args.repeats,
                    cases=cases, settings=SETTINGS, jobs=[1, 4, 8], seed_base=20261102,
                    rng_streams='Python random: sample_seed; NumPy: sample_seed mod 2^32; DendroPy: sample_seed + 1',
                    binaries={str(path): sha(path) for path in [args.binary, args.reference]},
                    reference_version=subprocess.check_output([str(args.reference), '--version'], text=True).strip(),
                    sources={str(path): sha(path) for path in [Path(__file__),
                             args.legacy_root / 'SimulateAlterRescale.py', args.legacy_root / 'alisim.py',
                             args.legacy_root / 'hogenom_diams.txt', args.legacy_root / 'raxml_diams.txt']})
    (args.output / 'metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    for backend, jobs in SETTINGS:
        measure(args, ('warmup', 2, 50, 500), backend, jobs, 19937)
    reference_hashes = {}
    fields = ['case', 'round', 'backend', 'jobs', 'seed', 'requested', 'written', 'failed',
              'wall_seconds', 'cpu_seconds', 'output_bytes', 'dataset_sha256', 'trees_sha256',
              'attempts', 'tree_seconds_sum', 'msa_seconds_sum']
    with (args.output / 'runs.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for repeat in range(args.repeats):
            case_order = cases[repeat % len(cases):] + cases[:repeat % len(cases)]
            for case in case_order:
                offset = (repeat + CASES.index(case)) % len(SETTINGS)
                order = SETTINGS[offset:] + SETTINGS[:offset]
                for backend, jobs in order:
                    row = measure(args, case, backend, jobs, 20261102 + repeat)
                    row['round'] = repeat
                    if backend == 'rust' or backend == 'python_pool' or (backend == 'python_threads' and jobs == 1):
                        key = (case[0], repeat, 'rust' if backend == 'rust' else 'python_serial_pool')
                        sig = (row['dataset_sha256'], row['written'], row['failed'], row['attempts'])
                        if key in reference_hashes and reference_hashes[key] != sig:
                            raise RuntimeError(f'Worker-dependent output: {key}')
                        reference_hashes[key] = sig
                    writer.writerow(row)
                    stream.flush()
                    print(f"round {repeat+1}/{args.repeats} {case[0]} {backend}/{jobs}: "
                          f"{row['wall_seconds']:.3f}s, {row['written']} accepted/{row['failed']} failed", flush=True)
    metadata['completed_utc'] = datetime.now(timezone.utc).isoformat()
    (args.output / 'metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    summarize(args.output)


if __name__ == '__main__':
    main()
