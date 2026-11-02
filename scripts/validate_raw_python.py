"""Compare actual alisim.py raw outputs with unconditional bridge output.

Matches the empirical Gamma prior and POW indels before any window selection,
trimming, or duplicate rejection. Uses independent ensembles for KS comparisons.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
import gzip
import json
import math
from pathlib import Path
import random
import shlex
import shutil
import subprocess
import tempfile
from types import SimpleNamespace

import numpy as np

import validate_msa_priors as base

CASES = [('GC', 'bd50'), ('G8', 'bd50'), ('GC', 'bd200'), ('G8', 'bd200')]


def initialize(legacy_root):
    base.initialize(legacy_root)
    base.BANKS['bd200'] = (200, ['--preset', 'python'])


def sample_alpha(seed):
    random.seed(seed)
    np.random.seed(seed & 0xffffffff)
    return float(base.LEGACY.sample_scale(base.LEGACY.alphas))


def unaligned_metrics(rows, tips):
    strings = [rows[f'T{i}'] for i in range(1, tips + 1)]
    assert all(set(row) <= set(base.AA) for row in strings), 'Unexpected gap-free alphabet'
    lengths = np.array(list(map(len, strings)))
    total = int(lengths.sum())
    assert total > 0
    result = dict(records=len(strings), mean_length=float(lengths.mean()), sd_length=float(lengths.std()),
                  min_length=int(lengths.min()), max_length=int(lengths.max()),
                  first_tip_length=len(strings[0]), total_residues=total,
                  duplicate=int(len(set(strings)) != tips), unique_fraction=len(set(strings)) / tips)
    joined = ''.join(strings)
    result.update({f'freq_{state}': joined.count(state) / total for state in base.AA})
    return result


def measure(task, aligned, unaligned):
    tips = base.BANKS[task['bank']][0]
    names = {f'T{i}' for i in range(1, tips + 1)}
    assert names <= aligned.keys() and names <= unaligned.keys()
    a = {name: aligned[name] for name in names}
    u = {name: unaligned[name] for name in names}
    for name in names:
        assert a[name].replace('-', '') == u[name], f'{name}: unaligned reconstruction mismatch'
    info = dict(sample=task['sample'], tree_sha256=base.sha(task['tree']),
                model=task['model'], seed=task['seed'])
    msa = {**base.metrics(a, tips), **info, 'alignment_sha256': base.canonical_hash(a),
           'gamma_shape': task['alpha']}
    gap_free = {**unaligned_metrics(u, tips), **info, 'alignment_sha256': base.canonical_hash(u),
                'gamma_shape': task['alpha']}
    # Measure padding separately: ancestral output may include columns in which
    # every tip is a gap. Do not silently drop these from the raw comparison.
    encoded = np.array([np.frombuffer(a[f'T{i}'].encode(), dtype=np.uint8)
                        for i in range(1, tips + 1)])
    msa['all_gap_columns'] = int((encoded == ord('-')).all(axis=0).sum())
    for result in [msa, gap_free]:
        assert all(math.isfinite(float(value)) for key, value in result.items()
                   if key not in {'model', 'tree_sha256', 'alignment_sha256'})
    return dict(aligned=msa, unaligned=gap_free,
                structure=dict(sample=task['sample'], source_aligned_records=len(aligned),
                               source_unaligned_records=len(unaligned), compared_tips=tips))


def keep_fixture(task, aligned_path, unaligned_path, command=None):
    if not task.get('fixture'):
        return
    out = Path(task['fixture'])
    out.mkdir(parents=True)
    shutil.copyfile(task['tree'], out / 'tree.nwk')
    shutil.copyfile(aligned_path, out / ('sample.fa.untrimmed' if command else 'sample.fa'))
    shutil.copyfile(unaligned_path, out / 'sample.unaligned.fa')
    data = {key: task[key] for key in ['sample', 'model', 'alpha', 'seed', 'alpha_seed', 'length']}
    data['command'] = command
    (out / 'fixture.json').write_text(json.dumps(data, indent=2) + '\n')


def python_raw(task):
    # Only process launch is instrumented. The actual function samples alpha,
    # invokes AliSim with --write-all/-mwopt and preserves the .untrimmed file.
    random.seed(task['alpha_seed'])
    np.random.seed(task['alpha_seed'] & 0xffffffff)
    legacy = base.LEGACY
    legacy.args = SimpleNamespace(gamma=task['gamma'], allow_duplicate_sequences=True, keep_logfiles=True)
    captured = {}
    with tempfile.TemporaryDirectory(prefix='python-raw-') as temp:
        def launch(command, **kwargs):
            words = shlex.split(command)
            words[0] = task['binary']
            prefix = words[words.index('--alisim') + 1]
            model = words[words.index('-m') + 1]
            assert model == task['model'], (model, task['model'])
            words += ['--seed', str(task['seed']), '-pre', prefix, '-redo']
            captured['command'] = words
            kwargs.pop('shell', None)
            captured['process'] = subprocess.Popen(words, **kwargs)
            return captured['process']
        legacy.subprocess = SimpleNamespace(Popen=launch, PIPE=subprocess.PIPE)
        error, _ = legacy.simulate_alignment(task['tree'], 'LG', task['gamma'], True,
                                             '', '', True, temp, task['length'], 20, 1, '', None)
        if error is not None or captured['process'].returncode:
            raise RuntimeError(f'Legacy simulation failed: {error}; {captured.get("command")}')
        prefix = Path(temp) / Path(task['tree']).stem
        aligned_path = Path(str(prefix) + '.fa.untrimmed')
        unaligned_path = Path(str(prefix) + '.unaligned.fa')
        result = measure(task, base.fasta(aligned_path), base.fasta(unaligned_path))
        keep_fixture(task, aligned_path, unaligned_path, captured['command'])
        return result


def bridge_raw(task):
    prefix = Path(task['output']) / task['id']
    a, u = Path(str(prefix) + '.fa'), Path(str(prefix) + '.unaligned.fa')
    result = measure(task, base.fasta(a), base.fasta(u))
    keep_fixture(task, a, u)
    return result


def make_tasks(case_index, gamma, bank, trees, backend, binary, length, output):
    rng = random.Random(3517921 + case_index * 770003 + backend * 310009)
    tasks = []
    for i, tree in enumerate(trees):
        alpha_seed = rng.randrange(1, 2**32)
        alpha = sample_alpha(alpha_seed)
        task = dict(id=f'raw_{gamma}_{bank}_{i:06}', sample=i, gamma=gamma, bank=bank,
                    tree=str(tree), alphabet='AA', alpha_seed=alpha_seed, alpha=alpha,
                    model=f'LG+{gamma}{{{alpha}}}', length=length, seed=rng.randrange(1, 2147483646),
                    binary=str(binary), indel='0.002,0.002', size=base.POW)
        if i < 3:
            task['fixture'] = str(output / 'fixtures' / f'{gamma}_{bank}' / str(backend) / f'{i:06}')
        tasks.append(task)
    return tasks


def write_report(output, metadata, comparisons, data):
    lines = ['# Raw Python / bridge alignment comparison', '',
             '**This comparison precedes duplicate rejection and window selection.** '
             'Each backend uses the actual empirical Gamma prior, LG, POW indels '
             f"(rates 0.002/0.002), and {metadata['length']} root sites. "
             'All attempts are retained. Production simulation behavior is unchanged.', '',
             '| Case | Python executable | Output | Variable checks | Detected differences | Largest KS D |',
             '| --- | --- | --- | ---: | ---: | ---: |']
    for gamma, bank in CASES:
        for label in ['python_240', 'python_2226']:
            for kind in ['aligned', 'unaligned']:
                name = f'{gamma}_{bank}/{label}/{kind}'
                rows = [r for r in comparisons if r['comparison'] == name]
                lines.append(f"| {gamma}_{bank} | {label} | {kind} | {len(rows)} | "
                             f"{sum(r['detected_difference'] for r in rows)} | {max(r['ks_d'] for r in rows):.4f} |")
    lines += ['', '## What was compared', '',
              '- Python: the actual `alisim.py.simulate_alignment` function with empirical alpha sampling, '
              '`--write-all`, `-mwopt`, POW indels and duplicates allowed. Read `.fa.untrimmed`, '
              'never the trimmed `.fa`, and the actual `.unaligned.fa`.',
              '- Bridge: the production direct-tree/numeric-output functions through `validate_alisim`, '
              'with the same Gamma sampling law and indel/model settings, with no rejection or trimming. '
              'The CLI does not currently sample empirical Gamma shapes automatically; this harness supplies '
              'one sampled model string per MSA.',
              '- There are 50 or 200 tips (`bd50`/`bd200`), using the Python tree preset. '
              'Independent tree banks, alpha draws and sequence seeds are used per backend. '
              'The raw-tree prior was validated separately; this is a comparison of the induced raw sequence prior.',
              '- Python emits ancestors as well as tips; our output contains tips. Tip IDs are matched '
              'before comparing. Raw tip-projected alignments retain all columns, including all-gap tip columns. '
              'Source row counts and all-gap columns are recorded explicitly.',
              '- Every unaligned tip sequence is checked, residue for residue, against the matching raw '
              'aligned row after removing gaps, for every simulation in every backend. Unaligned distributions '
              'are compared separately using lengths, composition and duplicate/unique-sequence fractions.',
              '- Python process launch only substitutes the executable and adds an explicit seed/log prefix. '
              'It uses one thread and the standalone default CPU kernel; no portable-kernel override is added. '
              'All other arguments and the alpha sampling/trim implementation are the actual Python code.',
              '- The bridge and standalone use different RNG initialization; matching seed values do not '
              'guarantee byte identity across executables. Distributional comparisons therefore use independent ensembles.', '',
              '## Method and scope', '',
              f"{metadata['samples']:,} independent whole MSAs per backend/case. KS comparisons use "
              'the equal-size DKWM upper bound `2 exp(-n D²)`, Bonferroni family alpha 0.01. '
              'Constant checks remain in the correction denominator. Sites and tips are not counted as independent replicates.', '',
              f"{metadata['comparisons']} variable comparisons + {metadata['constant_checks']} constant checks; "
              f"family size {metadata['family_size']}; detection distance {comparisons[0]['critical_distance']:.4f}.", '',
              'No detected difference is evidence at this resolution, not proof of equality of all joint '
              'distributions or every possible AliSim option. Custom matrices and supplied root sequences are not covered.', '',
              '## Raw defaults versus Python, without any filtering', '']
    baseline_rows = [r for r in comparisons if r['family'] == 'defaults']
    for label in ['python_240', 'python_2226']:
        for gamma in ['GC', 'G8']:
            selected = [r for r in baseline_rows if r['comparison'] == f'default_LG_vs_{gamma}/{label}']
            if selected:
                lines.append(f"- Default LG (no indels) versus raw Python {gamma}, {label}, 50 tips: "
                             f"{sum(r['detected_difference'] for r in selected)} / {len(selected)} variable "
                             'features differ **before either deduplication or cropping**.')
    lines += ['', 'This baseline reuses the earlier unconditional `lg_bd50-bridge.csv.gz` ensemble '
              '(never duplicate-filtered). It is separate from the fair, matched-prior comparisons above. '
              'Thus the earlier full-pipeline differences cannot be attributed to window deduplication alone.', '',
              '## Files and reproduction', '',
              'Raw aligned/unaligned per-MSA measurements are gzip CSVs. '
              '`*-structure.csv.gz` records original row counts. `fixtures/` preserves the first three '
              'tree/raw-alignment/unaligned files per backend/case, plus exact model strings and seeds; '
              'backend 0 is the bridge, 1 is Python + 2.4.0, 2 is Python + 2.2.2.6. '
              'Python fixture commands contain temporary paths; replace them with the retained tree and a new output prefix.', '',
              '```sh', 'cargo build --release --example validate_alisim',
              '/path/to/legacy/.venv/bin/python scripts/validate_raw_python.py --output validation/python_raw-new',
              '```', '',
              'See [comparisons.csv](comparisons.csv), [metadata.json](metadata.json), and '
              '[constants.json](constants.json). Reference executables are built as described in '
              '[the earlier report](../msa_priors/REPORT.md).', '']
    (output / 'REPORT.md').write_text('\n'.join(lines))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--samples', type=int, default=3000)
    parser.add_argument('--length', type=int, default=200)
    parser.add_argument('--jobs', type=int, default=8)
    parser.add_argument('--binary', type=Path, default=base.ROOT / 'target/release/pf-sim')
    parser.add_argument('--harness', type=Path, default=base.ROOT / 'target/release/examples/validate_alisim')
    parser.add_argument('--reference', type=Path, default=base.ROOT / 'target/alisim-reference/build/iqtree/iqtree2')
    parser.add_argument('--installed', type=Path, default=Path(shutil.which('iqtree2') or 'iqtree2'))
    parser.add_argument('--legacy-root', type=Path, default=Path('/Users/lucblassel/Work/bayesnj_simulation/scripts'))
    parser.add_argument('--baseline', type=Path, default=base.ROOT / 'validation/msa_priors/lg_bd50-bridge.csv.gz')
    parser.add_argument('--output', type=Path, default=base.ROOT / 'validation/python_raw')
    parser.add_argument('--pilot', action='store_true', help='Allow small descriptive-only smoke runs')
    args = parser.parse_args()
    if args.samples < 458 and not args.pilot:
        parser.error('Use >=458 samples for DKWM, or --pilot for a smoke run')
    if args.length != 200:
        parser.error('The reused raw-default baseline has 200 root sites; use --length 200')
    if args.jobs < 1 or args.samples < 1 or args.length < 1:
        parser.error('samples, jobs and length must be positive')
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    initialize(str(args.legacy_root))
    metadata = dict(samples=args.samples, length=args.length, jobs=args.jobs, cases=CASES,
                    alpha_count=len(base.LEGACY.alphas), pilot=args.pilot,
                    tree_bank_backend_offset=10,
                    binaries={str(p): base.sha(p) for p in [args.binary, args.harness, args.reference, args.installed]},
                    versions={str(p): subprocess.check_output([str(p), '--version'], text=True).strip()
                              for p in [args.reference, args.installed]},
                    sources={str(p): base.sha(p) for p in [Path(__file__), Path(base.__file__),
                             args.legacy_root / 'alisim.py', args.legacy_root / 'hogenom_alphas.txt',
                             base.ROOT / 'examples/validate_alisim.rs', base.ROOT / 'native/alisim_bridge.cpp']},
                    alpha_prior='actual alisim.py sample_scale: uniform empirical choice + N(0, alpha/10), floor 0.05',
                    baseline_hash=base.sha(args.baseline))
    data, comparisons, constants = {}, [], []
    with tempfile.TemporaryDirectory(prefix='python-raw-validation-') as temp, ProcessPoolExecutor(
            max_workers=args.jobs, initializer=initialize, initargs=(str(args.legacy_root),)) as pool:
        scratch = Path(temp)
        banks = {(bank, backend): base.tree_bank(args.binary, scratch, bank, backend + 10, args.samples)
                 for bank in ['bd50', 'bd200'] for backend in range(3)}
        for case_index, (gamma, bank) in enumerate(CASES):
            name = f'{gamma}_{bank}'
            for backend, label, binary in [(0, 'bridge', args.harness),
                                           (1, 'python_240', args.reference),
                                           (2, 'python_2226', args.installed)]:
                tasks = make_tasks(case_index, gamma, bank, banks[bank, backend], backend, binary,
                                   args.length, args.output)
                if backend == 0:
                    folder = scratch / name
                    folder.mkdir()
                    produced = base.generate_bridge(tasks, args.harness, folder, min(args.jobs, args.samples))
                    rows = list(pool.map(bridge_raw, produced, chunksize=8))
                    shutil.rmtree(folder)
                else:
                    rows = list(pool.map(python_raw, tasks, chunksize=8))
                for kind in ['aligned', 'unaligned', 'structure']:
                    values = [r[kind] for r in rows]
                    base.save_rows(args.output / f'{name}-{label}-{kind}.csv.gz', values)
                    data[name, label, kind] = values
                print(f'{name}: {label}, {len(rows)} raw aligned + unaligned outputs', flush=True)
            for label in ['python_240', 'python_2226']:
                for kind in ['aligned', 'unaligned']:
                    base.compare_groups(data[name, 'bridge', kind], data[name, label, kind],
                                        f'{name}/{label}/{kind}', 'matched_raw', comparisons, constants)
        # The earlier no-indel LG ensemble is unconditional too: no filtering
        # phase can explain differences already visible in these raw outputs.
        with gzip.open(args.baseline, 'rt') as stream:
            old = [{k: (v if k in {'model', 'tree_sha256', 'alignment_sha256', 'seed'} else float(v))
                    for k, v in row.items()} for row in csv.DictReader(stream)]
        for gamma in ['GC', 'G8']:
            for label in ['python_240', 'python_2226']:
                base.compare_groups(old, data[f'{gamma}_bd50', label, 'aligned'],
                                    f'default_LG_vs_{gamma}/{label}', 'defaults', comparisons, constants)
    family_size = len(comparisons) + len(constants)
    for row in comparisons:
        row['p_upper_adjusted'] = min(1.0, family_size * row['p_upper'])
        row['detected_difference'] = row['p_upper_adjusted'] < 0.01
        row['critical_distance'] = (math.sqrt(math.log(2 * family_size / 0.01) / row['n'])
                                    if row['n'] >= 458 else
                                    math.sqrt(2 * math.log(4 * family_size / 0.01) / row['n']))
    metadata.update(comparisons=len(comparisons), constant_checks=len(constants), family_size=family_size)
    with (args.output / 'comparisons.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(comparisons[0]))
        writer.writeheader()
        writer.writerows(comparisons)
    (args.output / 'constants.json').write_text(json.dumps(constants, indent=2) + '\n')
    (args.output / 'metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    write_report(args.output, metadata, comparisons, data)


if __name__ == '__main__':
    main()
