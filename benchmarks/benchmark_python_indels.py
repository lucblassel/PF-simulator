"""Time raw indel simulations: actual Python functions versus the Rust CLI.

No duplicate checks, retries, or find_unduplicated. The primary Python reference
includes ancestral export and the script's first-L truncation plus raw files;
a tip-only control separates the extra ancestral-output cost. All use LG (no
Gamma), Python BD tree priors, and POW indels at rates 0.002/0.002.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import statistics
import subprocess
import sys
import tempfile
import time

import benchmark_python as base

ROOT = base.ROOT
CASES = [('indels_50', 128, 50, 500), ('indels_200', 128, 200, 5000)]
SETTINGS = [(backend, jobs) for jobs in [1, 4, 8]
            for backend in ['rust', 'python_ancestral', 'python_tips']]
POW = 'POW{1.7/100},POW{1.7/100}'


def initialize(legacy_root, ancestral):
    base.initialize(legacy_root)
    base.MSA_MODULE.args.allow_duplicate_sequences = True
    simulate = base.MSA_MODULE.simulate_alignment

    def with_indels(*args):
        args = list(args)
        args[3] = ancestral
        args[6] = True
        return simulate(*args)

    # Change only the invocation options, never the legacy implementation.
    base.MSA_MODULE.simulate_alignment = with_indels


def python_sample(task):
    result = base.python_sample(task)
    assert result['attempts'] == 1 and result['written'] == 1
    name = f"sample_{task['sample']:06}.fa"
    out = Path(task['output'])
    aligned = out / 'sequences' / name
    aligned.rename(out / 'trimmed_sequences' / name)
    Path(str(aligned) + '.untrimmed').rename(aligned)
    (out / 'sequences' / name.replace('.fa', '.unaligned.fa')).rename(out / 'unaligned_sequences' / name)
    return result


def reference_run(args):
    args.data.mkdir()
    for directory in ['trees', 'sequences', 'unaligned_sequences', 'trimmed_sequences']:
        (args.data / directory).mkdir()
    ancestral = args.backend == 'python_ancestral'
    tasks = [dict(sample=i, seed=args.seed, tips=args.tips, length=args.length,
                  output=str(args.data), reference=str(args.reference), native_threads=1)
             for i in range(args.samples)]
    if args.jobs == 1:
        initialize(str(args.legacy_root), ancestral)
        rows = [python_sample(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=min(args.jobs, args.samples), initializer=initialize,
                                 initargs=(str(args.legacy_root), ancestral)) as pool:
            rows = list(pool.map(python_sample, tasks, chunksize=1))
    summary = dict(requested_samples=args.samples, written_samples=sum(r['written'] for r in rows),
                   failed_duplicate_samples=0, attempts=sum(r['attempts'] for r in rows),
                   tree_seconds_sum=sum(r['tree_seconds'] for r in rows),
                   msa_seconds_sum=sum(r['msa_seconds'] for r in rows))
    (args.data / 'simulation_summary.json').write_text(json.dumps(summary) + '\n')


def validate_dataset(output, samples, tips, length, backend):
    """Outside timing: raw reconstruction, row/ID counts, trimming and hashes."""
    from Bio import SeqIO
    trees = sorted((output / 'trees').glob('*.nwk'))
    aligned_dir = 'msa_untrimmed' if backend == 'rust' else 'sequences'
    unaligned_dir = 'unaligned' if backend == 'rust' else 'unaligned_sequences'
    aligned = sorted((output / aligned_dir).glob('*.fa'))
    unaligned = sorted((output / unaligned_dir).glob('*.fa'))
    assert len(trees) == len(aligned) == len(unaligned) == samples
    assert {p.stem for p in trees} == {p.stem for p in aligned} == {p.stem for p in unaligned}
    digest, tip_digest = hashlib.sha256(), hashlib.sha256()
    total_records = duplicate_outputs = size = 0
    expected_names = {f'T{i}' for i in range(1, tips + 1)}

    def read(path):
        records = list(SeqIO.parse(path, 'fasta'))
        result = {r.id: str(r.seq) for r in records}
        assert len(result) == len(records), f'Duplicate record ID: {path}'
        return result

    for p in aligned:
        a = read(p)
        u = read(output / unaligned_dir / p.name)
        assert a.keys() == u.keys() and expected_names <= a.keys()
        expected_records = 2 * tips - 1 if backend == 'python_ancestral' else tips
        assert len(a) == expected_records
        assert len(set(map(len, a.values()))) == 1
        assert len(a['T1']) >= length
        for name in a:
            assert a[name].replace('-', '') == u[name], f'Gap reconstruction: {p}/{name}'
        tip_rows = [a[f'T{i}'] for i in range(1, tips + 1)]
        duplicate_outputs += len(set(tip_rows)) != tips  # Observe, do not reject.
        total_records += len(a)
        if backend != 'rust':
            truncated = read(output / 'trimmed_sequences' / p.name)
            assert truncated == {name: row[:length] for name, row in a.items()}
        tip_digest.update(p.name.encode() + b'\0')
        for name in sorted(expected_names):
            tip_digest.update(name.encode() + b'\0' + a[name].encode() + b'\0' + u[name].encode() + b'\0')
    for p in sorted(p for p in output.rglob('*') if p.is_file() and p.name != 'simulation_summary.json'):
        data = p.read_bytes()
        digest.update(str(p.relative_to(output)).encode() + b'\0' + data)
        size += len(data)
    return dict(written=len(aligned), dataset_sha256=digest.hexdigest(), tips_sha256=tip_digest.hexdigest(),
                output_bytes=size, aligned_records=total_records, duplicate_outputs=duplicate_outputs,
                reconstruction_checks=total_records)


def measure(args, case, backend, jobs, seed):
    name, samples, tips, length = case
    with tempfile.TemporaryDirectory(prefix='python-indel-bench-') as temp:
        out = Path(temp) / 'data'
        if backend == 'rust':
            command = [str(args.binary), '--out', str(out), '--preset', 'python', '--model', 'LG',
                       '--samples', str(samples), '--tips', str(tips), '--length', str(length),
                       '--seed', str(seed), '--jobs', str(jobs), '--allow-duplicate-sequences',
                       '--indel', '0.002,0.002', '--indel-size', POW, '--materialize', 'unaligned']
        else:
            command = [sys.executable, str(Path(__file__).resolve()), '--reference-run',
                       '--data', str(out), '--samples', str(samples), '--tips', str(tips),
                       '--length', str(length), '--seed', str(seed), '--backend', backend,
                       '--jobs', str(jobs), '--reference', str(args.reference),
                       '--legacy-root', str(args.legacy_root)]
        before = resource.getrusage(resource.RUSAGE_CHILDREN)
        start = time.perf_counter()
        run = subprocess.run(command, capture_output=True, text=True, timeout=600)
        elapsed = time.perf_counter() - start
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        if run.returncode:
            raise RuntimeError(f'{command}\n{run.stdout}\n{run.stderr}')
        summary = json.loads((out / 'simulation_summary.json').read_text())
        assert summary['written_samples'] == samples and summary['failed_duplicate_samples'] == 0
        if backend != 'rust':
            assert summary['attempts'] == samples
        return dict(case=name, backend=backend, jobs=jobs, seed=seed, requested=samples, failed=0,
                    wall_seconds=elapsed,
                    cpu_seconds=after.ru_utime + after.ru_stime - before.ru_utime - before.ru_stime,
                    **validate_dataset(out, samples, tips, length, backend),
                    attempts=summary.get('attempts', ''),
                    tree_seconds_sum=summary.get('tree_seconds_sum', ''),
                    msa_seconds_sum=summary.get('msa_seconds_sum', ''))


def summarize(output):
    meta = json.loads((output / 'metadata.json').read_text())
    rows = list(csv.DictReader((output / 'runs.csv').open()))
    summary = []
    lines = ['# Python versus Rust: indels, no duplicate checking', '',
             'Full raw simulation before `find_unduplicated`: LG without Gamma, Python birth-death '
             'tree prior, insertion/deletion rates 0.002/0.002 and POW{1.7/100} sizes. '
             'Every tree produces one MSA; no retry, rejection, window search or cropping by Rust.', '',
             f"{meta['repeats']} rounds per case and setting, {meta['logical_cpus']} logical CPUs, "
             f"release Rust binary; Python {meta['python']}. Standalone and linked AliSim are both 2.4.0.", '',
             '| Case | Pairs | Tips × root sites | Workers | Rust (s) | Python with ancestors (s) | Python tips only (s) | Ancestor Python / Rust | Tip-only Python / Rust |',
             '| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |']
    for case, samples, tips, length in meta['cases']:
        for jobs in [1, 4, 8]:
            groups = {backend: [r for r in rows if r['case'] == case and int(r['jobs']) == jobs and r['backend'] == backend]
                      for backend in ['rust', 'python_ancestral', 'python_tips']}
            assert all(len(g) == meta['repeats'] for g in groups.values())
            timings = {backend: statistics.median(float(r['wall_seconds']) for r in group)
                       for backend, group in groups.items()}
            baseline = {r['round']: float(r['wall_seconds']) for r in groups['rust']}
            ratios = {}
            for backend in ['python_ancestral', 'python_tips']:
                values = [float(r['wall_seconds']) / baseline[r['round']] for r in groups[backend]]
                q1, q3 = base.quartiles(values)
                ratios[backend] = statistics.median(values)
                summary.append(dict(case=case, jobs=jobs, backend=backend,
                                    rust_median_seconds=timings['rust'], python_median_seconds=timings[backend],
                                    paired_speedup=ratios[backend], speedup_q1=q1, speedup_q3=q3,
                                    python_pairs_per_second=statistics.median(samples / float(r['wall_seconds']) for r in groups[backend]),
                                    rust_pairs_per_second=statistics.median(samples / float(r['wall_seconds']) for r in groups['rust'])))
            lines.append(f"| {case} | {samples} | {tips} × {length} | {jobs} | {timings['rust']:.3f} | "
                         f"{timings['python_ancestral']:.3f} | {timings['python_tips']:.3f} | "
                         f"{ratios['python_ancestral']:.2f}× | {ratios['python_tips']:.2f}× |")
    lines += ['', 'Times are medians of whole-run wall time; speedups are medians of same-round ratios. '
              'All settings retain all requested pairs. IQRs and throughput are in [summary.csv](summary.csv).', '',
              '## Reference and output scope', '',
              '- The reference calls the actual `SimulateAlterRescale.py.simulateTree` and '
              '`alisim.py.simulate_alignment` functions with `add_indels=True` and '
              '`allow_duplicate_sequences=True`. It stops before `find_unduplicated`; Gamma is '
              'disabled explicitly to match the preceding LG benchmark and the current Rust CLI.',
              '- The primary reference also sets `ancestral=True`, as the Snakefile does. The '
              'tip-only Python control changes only this flag. Python tip trajectories are checked '
              'for exact equality across these two flags; otherwise no cross-language byte identity is assumed.',
              '- Python still creates its first-L-column `.fa`, full `.fa.untrimmed`, and '
              '`.unaligned.fa` files. All that I/O and ancestral export is included. Raw and unaligned '
              'files are moved into the same directory layout as Rust, and the truncated file is '
              'retained under `trimmed_sequences/`; these renames are timed. No record filtering or '
              'sequence rewriting is added by the harness.',
              '- Rust emits full aligned and unaligned **tips only**, with duplicate checking '
              'disabled. Consequently the ancestor-enabled Python comparison includes extra output '
              'work; use the tip-only control to assess that contribution. Even this control retains '
              'Python’s truncated file and default log generation/deletion.',
              '- AliSim is forced to one native thread by `alisim.py` whenever indels are enabled. '
              'Thus 1 means sequential execution; 4 and 8 mean persistent process workers on both '
              'sides. The Python pool is a benchmark-only scheduler, not the original script’s '
              '`--processes` flag. Both versions of the Python reference share this scheduler.', '',
              '## Timing and validation', '',
              '- Includes process/interpreter startup, imports, empirical tree-input loading, tree '
              'generation/transformation, AliSim simulation, all native/Python output writes, and '
              'worker shutdown. Python stages are combined in one interpreter, excluding Snakemake '
              'and a second interpreter startup. Filesystem writes are buffered, without fsync.',
              '- One unmeasured two-sample warmup per setting; three rounds by default, rotating '
              'case/setting order. Runs are sequential, on fresh temporary output directories. '
              'Validation, hashing and directory deletion are outside the timed region.',
              '- Each Python simulation must launch exactly once; every backend must write every '
              'requested sample with zero duplicate failures. Duplicated tip sequences are counted '
              'but never rejected. Output IDs, row counts and lengths are checked. Every raw aligned '
              'row (including Python ancestors) must exactly reconstruct its corresponding '
              'unaligned sequence after gap removal; Python’s truncated rows must match the first L columns.',
              '- Sample seeds use the same independent Python/NumPy/DendroPy stream setup as the '
              'earlier benchmark. Data hashes must be invariant across worker counts per backend; '
              'Python ancestral/tip-only hashes for tip rows must match. Tree/RNG trajectories still '
              'differ between Rust and Python, while priors/settings match. Rust simulation behavior is unchanged.', '',
              f"All {len(rows)} runs completed. The raw files contain "
              f"{sum(int(r['reconstruction_checks']) for r in rows):,} exact aligned-to-unaligned row checks; "
              f"{sum(int(r['duplicate_outputs']) for r in rows):,} output MSAs contained duplicate tip sequences "
              '(counts include repeated seeds across worker settings).', '',
              '## Reproduce', '', '```sh',
              '/path/to/legacy/.venv/bin/python benchmarks/benchmark_python_indels.py --output benchmarks/python_indels-new',
              '```', '', 'Raw measurements: [runs.csv](runs.csv), [summary.csv](summary.csv), '
              '[metadata.json](metadata.json). [No-indel benchmark](../python/REPORT.md).', '']
    with (output / 'summary.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)
    (output / 'REPORT.md').write_text('\n'.join(lines))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'benchmarks/python_indels')
    parser.add_argument('--binary', type=Path, default=ROOT / 'target/release/pf-sim')
    parser.add_argument('--reference', type=Path, default=ROOT / 'target/alisim-reference/build/iqtree/iqtree2')
    parser.add_argument('--legacy-root', type=Path, default=Path('/Users/lucblassel/Work/bayesnj_simulation/scripts'))
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--cases', nargs='+', choices=[case[0] for case in CASES])
    parser.add_argument('--summarize-only', action='store_true')
    parser.add_argument('--reference-run', action='store_true', help=argparse.SUPPRESS)
    for flag, typ in [('data', Path), ('samples', int), ('tips', int), ('length', int),
                      ('seed', int), ('jobs', int), ('backend', str)]:
        parser.add_argument('--' + flag, type=typ, help=argparse.SUPPRESS)
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
    meta = dict(started_utc=datetime.now(timezone.utc).isoformat(), platform=platform.platform(),
                logical_cpus=os.cpu_count(), python=platform.python_version(), repeats=args.repeats,
                cases=cases, settings=SETTINGS, seed_base=20261102,
                model='LG', gamma=None, indel='0.002,0.002', indel_size=POW, duplicate_checks=False,
                binaries={str(path): base.sha(path) for path in [args.binary, args.reference]},
                reference_version=subprocess.check_output([str(args.reference), '--version'], text=True).strip(),
                sources={str(path): base.sha(path) for path in [Path(__file__), Path(base.__file__),
                         args.legacy_root / 'SimulateAlterRescale.py', args.legacy_root / 'alisim.py',
                         args.legacy_root / 'hogenom_diams.txt', args.legacy_root / 'raxml_diams.txt']})
    (args.output / 'metadata.json').write_text(json.dumps(meta, indent=2) + '\n')
    for backend, jobs in SETTINGS:
        measure(args, ('warmup', 2, 50, 500), backend, jobs, 19937)
    reference_hashes, tip_hashes = {}, {}
    with (args.output / 'runs.csv').open('w', newline='') as stream:
        writer = None
        for repeat in range(args.repeats):
            order = cases[repeat % len(cases):] + cases[:repeat % len(cases)]
            for case in order:
                offset = (repeat + CASES.index(case)) % len(SETTINGS)
                settings = SETTINGS[offset:] + SETTINGS[:offset]
                for backend, jobs in settings:
                    row = measure(args, case, backend, jobs, 20261102 + repeat)
                    row['round'] = repeat
                    key = case[0], repeat, backend
                    sig = row['dataset_sha256']
                    assert key not in reference_hashes or reference_hashes[key] == sig, ('Worker-dependent data', key)
                    reference_hashes[key] = sig
                    if backend.startswith('python'):
                        key = case[0], repeat
                        assert key not in tip_hashes or tip_hashes[key] == row['tips_sha256'], ('Ancestor flag changes tips', key)
                        tip_hashes[key] = row['tips_sha256']
                    if writer is None:
                        writer = csv.DictWriter(stream, fieldnames=list(row))
                        writer.writeheader()
                    writer.writerow(row)
                    stream.flush()
                    print(f"round {repeat+1}/{args.repeats} {case[0]} {backend}/{jobs}: {row['wall_seconds']:.3f}s", flush=True)
    meta['completed_utc'] = datetime.now(timezone.utc).isoformat()
    (args.output / 'metadata.json').write_text(json.dumps(meta, indent=2) + '\n')
    summarize(args.output)


if __name__ == '__main__':
    main()
