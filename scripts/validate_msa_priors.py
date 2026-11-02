"""Compare linked numeric/direct AliSim with independent standalone executables.

No simulator behavior is changed. Uses NumPy/BioPython from the legacy env.
Retains per-MSA metrics, hashes and provenance, not pooled-site pseudo-samples.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import csv
import gzip
import hashlib
import importlib.util
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
from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
AA = 'ACDEFGHIKLMNPQRSTVWY'
POW = 'POW{1.7/100},POW{1.7/100}'
CASES = [
    ('lg_u10', 'u10', 'LG', 'AA', 200, '', ''),
    ('lg_bd50', 'bd50', 'LG', 'AA', 200, '', ''),
    ('lg_u200', 'u200', 'LG', 'AA', 200, '', ''),
    ('g8_bd50', 'bd50', 'LG+G8{0.5}', 'AA', 200, '', ''),
    ('gc_bd50', 'bd50', 'LG+GC{0.5}', 'AA', 200, '', ''),
    ('invariant_bd50', 'bd50', 'LG+I{0.2}', 'AA', 200, '', ''),
    ('pow_bd50', 'bd50', 'LG', 'AA', 200, '0.002,0.002', POW),
    ('gc_pow_bd50', 'bd50', 'LG+GC{0.5}', 'AA', 200, '0.002,0.002', POW),
    ('jc_u10', 'u10', 'JC', 'DNA', 200, '', ''),
    ('duplicates_u4', 'u4', 'JC', 'DNA', 2, '', ''),
]
BANKS = {'u10': (10, []), 'bd50': (50, ['--preset', 'python']),
         'u200': (200, []), 'u4': (4, ['--branch-mean', '0.5'])}
LEGACY = None
FILTER = None


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def initialize(legacy_root):
    global LEGACY, FILTER
    LEGACY = load_module('legacy_alisim_validation', Path(legacy_root) / 'alisim.py')
    FILTER = load_module('legacy_filter_validation', Path(legacy_root) / 'find_unduplicated.py')
    LEGACY.alphas = LEGACY.load_list(LEGACY.ALPHAS_PATH)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fasta(path):
    result = {}
    for record in SeqIO.parse(str(path), 'fasta'):
        if record.id in result:
            raise RuntimeError(f'Duplicate FASTA ID {record.id}')
        result[record.id] = str(record.seq)
    return result


def canonical_hash(rows):
    return hashlib.sha256(''.join(f'>{key}\n{rows[key]}\n' for key in sorted(rows)).encode()).hexdigest()


def metrics(rows, tips):
    names = [f'T{i}' for i in range(1, tips + 1)]
    assert set(names) <= rows.keys()
    strings = [rows[name] for name in names]
    assert len(set(map(len, strings))) == 1
    length = len(strings[0])
    a = np.array([np.frombuffer(s.encode(), dtype=np.uint8) for s in strings])
    gaps = a == ord('-')
    valid = (~gaps).sum(axis=0)
    counts = np.array([(a == ord(state)).sum(axis=0) for state in AA])
    assert (counts.sum(axis=0) == valid).all(), 'Unexpected alphabet'
    probs = counts / np.maximum(valid, 1)
    entropy = -(probs * np.log(np.maximum(probs, 1e-300))).sum(axis=0)
    invariant = (counts.max(axis=0) == valid) & (valid > 0)
    complete = valid == tips
    identity = (counts * (counts - 1)).sum(axis=0) / np.maximum(valid * (valid - 1), 1)
    pairs = sorted({(i, j) for i, j in [(0, 1), (0, tips // 2), (0, tips - 1),
                                      (tips // 2, tips - 1), (tips // 3, 2 * tips // 3)] if i < j})
    distances = []
    for i, j in pairs:
        keep = ~(gaps[i] | gaps[j])
        distances.append(float((a[i, keep] != a[j, keep]).mean()) if keep.any() else 0.0)
    ungapped = [s.replace('-', '') for s in strings]
    gap_runs = []
    for row in gaps:
        padded = np.r_[False, row, False].astype(np.int8)
        edges = np.diff(padded)
        gap_runs.extend((np.where(edges == -1)[0] - np.where(edges == 1)[0]).tolist())
    result = dict(records=len(rows), tips=tips, aligned_length=length,
                  ungapped_mean=float(np.mean(list(map(len, ungapped)))),
                  ungapped_sd=float(np.std(list(map(len, ungapped)))),
                  gap_fraction=float(gaps.mean()), gapped_columns=float(gaps.any(axis=0).mean()),
                  mean_gap_run=float(np.mean(gap_runs)) if gap_runs else 0.0,
                  max_gap_run=max(gap_runs, default=0),
                  invariant_fraction=float(invariant.mean()),
                  complete_invariant_fraction=float(invariant[complete].mean()) if complete.any() else 0.0,
                  entropy_mean=float(entropy[valid > 0].mean()),
                  entropy_sd=float(entropy[valid > 0].std()),
                  identity_variance=float(identity[valid > 1].var()),
                  pairwise_mean=float(np.mean(distances)), pairwise_sd=float(np.std(distances)),
                  first_pair_distance=distances[0],
                  duplicate=float(len(set(strings)) != tips),
                  unique_fraction=len(set(strings)) / tips,
                  ungapped_duplicate=float(len(set(ungapped)) != tips))
    total = counts.sum()
    result.update({f'freq_{state}': float(counts[i].sum() / total) for i, state in enumerate(AA)})
    return result


def describe(task, rows, ungapped=None):
    tips = BANKS[task['bank']][0]
    result = metrics(rows, tips)
    if ungapped is not None:
        for i in range(1, tips + 1):
            name = f'T{i}'
            assert ungapped[name] == rows[name].replace('-', ''), 'Gap-free reconstruction mismatch'
    result.update(sample=task['sample'], tree_sha256=sha(task['tree']),
                  alignment_sha256=canonical_hash(rows), seed=task['seed'], model=task['model'])
    return result


def standalone(task):
    with tempfile.TemporaryDirectory(prefix='alisim-reference-') as temp:
        prefix = Path(temp) / 'msa'
        command = [task['binary'], '--alisim', str(prefix), '-t', task['tree'],
                   '-m', task['model'], '--seqtype', task['alphabet'], '--length', str(task['length']),
                   '--seed', str(task['seed']), '--threads', '1', '-af', 'fasta', '-lk', 'x86',
                   '-pre', str(prefix), '-redo']
        if task['indel']:
            command += ['--indel', task['indel'], '--indel-size', task['size']]
        run = subprocess.run(command, capture_output=True, text=True, timeout=60)
        if run.returncode:
            raise RuntimeError(f'{command}\n{run.stdout}\n{run.stderr}')
        rows = fasta(prefix.with_suffix('.fa'))
        ungapped_path = prefix.with_suffix('.unaligned.fa')
        ungapped = fasta(ungapped_path) if ungapped_path.exists() else None
        assert (ungapped is not None) == bool(task['indel'])
        return describe(task, rows, ungapped)


def bridge_metric(task):
    prefix = Path(task['output']) / task['id']
    rows = fasta(prefix.with_suffix('.fa'))
    ungapped_path = prefix.with_suffix('.unaligned.fa')
    result = describe(task, rows, fasta(ungapped_path) if ungapped_path.exists() else None)
    return result


def legacy_policy(task):
    # Run the actual Python function, including empirical alpha sampling,
    # trimming, and the actual downstream unique-window selection. Only process
    # launch is instrumented to pin executable/seed/kernel/log location.
    random.seed(task['alpha_seed'])
    np.random.seed(task['alpha_seed'] & 0xffffffff)
    LEGACY.args = SimpleNamespace(gamma=task['gamma'], allow_duplicate_sequences=True, keep_logfiles=True)
    captured = {}
    with tempfile.TemporaryDirectory(prefix='alisim-python-') as temp:
        def launch(command, **kwargs):
            words = shlex.split(command)
            words[0] = task['binary']
            prefix = words[words.index('--alisim') + 1]
            captured['model'] = words[words.index('-m') + 1]
            words += ['--seed', str(task['seed']), '-lk', 'x86', '-pre', prefix, '-redo']
            kwargs.pop('shell', None)
            return subprocess.Popen(words, **kwargs)
        LEGACY.subprocess = SimpleNamespace(Popen=launch, PIPE=subprocess.PIPE)
        error, _ = LEGACY.simulate_alignment(task['tree'], 'LG', task['gamma'], True,
                                             '', '', True, temp, task['length'], 20, 1, '', None)
        if error is not None:
            raise RuntimeError(f'Legacy script failed: {error}')
        path = Path(temp) / (Path(task['tree']).stem + '.fa')
        records = list(SeqIO.parse(str(path) + '.untrimmed', 'fasta'))
        tips = [record for record in records if FILTER.seq_is_tip(record)]
        offset = None if FILTER.has_duplicates(tips) else FILTER.find_unduplicated(tips, task['length'])
        accepted = offset is not None
        if accepted:
            selected = {record.id: str(record.seq[int(offset):int(offset) + task['length']]) for record in tips}
            result = describe({**task, 'model': captured['model']}, selected)
        else:
            result = dict(sample=task['sample'], tree_sha256=sha(task['tree']), seed=task['seed'], model=captured['model'])
        result['accepted'] = int(accepted)
        result['full_aligned_length'] = len(records[0])
        result['internal_records'] = len(records) - len(tips)
        result['window_offset'] = int(offset) if accepted else -1
        result['alpha'] = float(captured['model'].split('{')[-1].rstrip('}'))
        return result


def ks(a, b):
    a, b = np.sort(a), np.sort(b)
    grid = np.unique(np.r_[a, b])
    return float(np.max(np.abs(np.searchsorted(a, grid, side='right') / len(a) -
                               np.searchsorted(b, grid, side='right') / len(b))))


def save_rows(path, rows):
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with gzip.open(path, 'wt', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def manifest_rows(tasks, path):
    text = ''.join('\t'.join(map(str, [task['id'], task['tree'], task['model'], task['alphabet'],
                                      task['length'], task['seed'], task['indel'], task['size']])) + '\n'
                   for task in tasks)
    path.write_text(text)


def generate_bridge(tasks, harness, scratch, jobs):
    batches = [tasks[i::jobs] for i in range(jobs)]
    def run(item):
        i, batch = item
        manifest = scratch / f'manifest-{i}.tsv'
        output = scratch / f'bridge-{i}'
        manifest_rows(batch, manifest)
        subprocess.run([str(harness), '--manifest', str(manifest), '--output', str(output)], check=True,
                       stdout=subprocess.DEVNULL, timeout=600)
        return [{**task, 'output': str(output)} for task in batch]
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        produced = [task for batch in pool.map(run, enumerate(batches)) for task in batch]
    return sorted(produced, key=lambda task: task['sample'])


def tree_bank(binary, scratch, name, backend, samples):
    tips, flags = BANKS[name]
    index = list(BANKS).index(name)
    seed = 2601101 + backend * 1000003 + index * 170003
    output = scratch / f'trees-{name}-{backend}'
    subprocess.run([str(binary), '--out', str(output), '--samples', str(samples), '--tips', str(tips),
                    '--seed', str(seed), '--jobs', '4', '--trees-only', *flags], check=True,
                   stdout=subprocess.DEVNULL, timeout=600)
    return sorted((output / 'trees').glob('*.nwk'))


def tasks_for(case, trees, backend, binary):
    name, bank, model, alphabet, length, indel, size = case
    rng = random.Random(210011 + backend * 700001 + CASES.index(case) * 110003)
    return [dict(id=f'{name}_{i:06}', sample=i, bank=bank, tree=str(tree), model=model,
                 alphabet=alphabet, length=length, seed=rng.randrange(1, 2147483646),
                 indel=indel, size=size, binary=str(binary)) for i, tree in enumerate(trees)]


def compare_groups(group_a, group_b, name, family, rows, constants):
    skip = {'sample', 'seed', 'tree_sha256', 'alignment_sha256', 'model', 'accepted', 'alpha',
            'full_aligned_length', 'internal_records', 'window_offset'}
    a = [row for row in group_a if row.get('accepted', 1)]
    b = [row for row in group_b if row.get('accepted', 1)]
    n = min(len(a), len(b))
    if n == 0:
        raise RuntimeError(f'{name}: no retained samples to compare')
    for metric in sorted((a[0].keys() & b[0].keys()) - skip):
        x = np.array([row[metric] for row in a[:n]], dtype=float)
        y = np.array([row[metric] for row in b[:n]], dtype=float)
        if np.ptp(np.r_[x, y]) == 0:
            constants.append(dict(comparison=name, metric=metric, value=float(x[0]), n=n))
            continue
        d = ks(x, y)
        rows.append(dict(family=family, comparison=name, metric=metric, n=n, ks_d=d,
                         mean_a=float(x.mean()), mean_b=float(y.mean()),
                         p_upper=min(1.0, (2 * math.exp(-n * d * d)) if n >= 458 else
                                     (4 * math.exp(-n * d * d / 2))),
                         bound='two_sample_DKWM' if n >= 458 else 'union_of_one_sample_DKW'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--samples', type=int, default=3000)
    parser.add_argument('--jobs', type=int, default=8)
    parser.add_argument('--binary', type=Path, default=ROOT / 'target/release/pf-sim')
    parser.add_argument('--harness', type=Path, default=ROOT / 'target/release/examples/validate_alisim')
    parser.add_argument('--reference', type=Path, default=ROOT / 'target/alisim-reference/build/iqtree/iqtree2')
    parser.add_argument('--installed', type=Path, default=Path(shutil.which('iqtree2') or 'iqtree2'))
    parser.add_argument('--legacy-root', type=Path, default=Path('/Users/lucblassel/Work/bayesnj_simulation/scripts'))
    parser.add_argument('--output', type=Path, default=ROOT / 'validation/msa_priors')
    args = parser.parse_args()
    if args.samples < 458:
        parser.error('--samples must be at least 458 for the equal-size DKWM bound')
    args.output.mkdir(parents=True, exist_ok=False)
    metadata = dict(samples=args.samples, jobs=args.jobs, binaries={str(p): sha(p) for p in
                    [args.binary, args.harness, args.reference, args.installed]},
                    versions={str(p): subprocess.check_output([str(p), '--version'], text=True).strip()
                              for p in [args.reference, args.installed]},
                    source_revision=subprocess.check_output(['git', '-C', str(ROOT / 'dependencies/iqtree2'), 'rev-parse', 'HEAD'], text=True).strip(),
                    source_hashes={str(p): sha(p) for p in [Path(__file__), ROOT / 'native/alisim_bridge.cpp',
                     ROOT / 'native/memory_io.cmake', ROOT / 'examples/validate_alisim.rs',
                     ROOT / 'scripts/build_alisim_reference.py', ROOT / 'src/tree.rs',
                     ROOT / 'src/alignment.rs', args.legacy_root / 'alisim.py',
                     args.legacy_root / 'find_unduplicated.py', args.legacy_root / 'hogenom_alphas.txt',
                     args.legacy_root.parent / 'Snakefile']},
                    cases=CASES, family_alpha=0.01,
                    reference_build='Unmodified pinned source with OpenMP, 1 runtime thread, -lk x86; no memory adapters')
    (args.output / 'metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    comparisons, constants, data = [], [], {}
    with tempfile.TemporaryDirectory(prefix='pf-msa-validation-') as temp, ProcessPoolExecutor(
            max_workers=args.jobs, initializer=initialize, initargs=(str(args.legacy_root),)) as pool:
        scratch = Path(temp)
        banks = {(name, backend): tree_bank(args.binary, scratch, name, backend, args.samples)
                 for name in BANKS for backend in range(3)}
        print('Generated independent tree banks', flush=True)
        for case in CASES:
            name = case[0]
            for backend, label, executable in [(0, 'bridge', args.harness),
                                                (1, 'standalone_240', args.reference),
                                                (2, 'standalone_2226', args.installed)]:
                tasks = tasks_for(case, banks[case[1], backend], backend, executable)
                if backend == 0:
                    folder = scratch / name
                    folder.mkdir()
                    produced = generate_bridge(tasks, args.harness, folder, args.jobs)
                    rows = list(pool.map(bridge_metric, produced, chunksize=20))
                    shutil.rmtree(folder)
                else:
                    rows = list(pool.map(standalone, tasks, chunksize=8))
                data[name, label] = rows
                save_rows(args.output / f'{name}-{label}.csv.gz', rows)
                print(f'{name}: {label}, {len(rows)} independent MSAs', flush=True)
            for label in ['standalone_240', 'standalone_2226']:
                compare_groups(data[name, 'bridge'], data[name, label], f'{name}/{label}', 'engine', comparisons, constants)
        # A changed gamma prior is a substantive negative control.
        control = ('gc_bd50', 'bd50', 'LG+GC{5}', 'AA', 200, '', '')
        tasks = tasks_for(CASES[4], banks['bd50', 1], 1, args.reference)
        tasks = [{**task, 'model': control[2], 'seed': (task['seed'] + 49157) % 2147483645 + 1} for task in tasks]
        rows = list(pool.map(standalone, tasks, chunksize=8))
        save_rows(args.output / 'control_gamma5.csv.gz', rows)
        compare_groups(data['gc_bd50', 'bridge'], rows, 'changed_gamma_0.5_vs_5', 'control', comparisons, constants)
        print('Completed gamma negative control', flush=True)

        # Audit the actually emitted default Python-preset dataset independently.
        pf_output = scratch / 'pf_emitted'
        subprocess.run([str(args.binary), '--out', str(pf_output), '--preset', 'python',
                        '--samples', str(args.samples), '--tips', '50', '--length', '200',
                        '--jobs', str(args.jobs), '--seed', '871731'], check=True, stdout=subprocess.DEVNULL)
        emitted = []
        for path in sorted((pf_output / 'sequences').glob('*.fa')):
            sample = int(path.stem.split('_')[-1])
            task = dict(sample=sample, bank='bd50', tree=str(pf_output / 'trees' / f'{path.stem}.nwk'),
                        seed='master:871731', model='LG')
            emitted.append(describe(task, fasta(path)))
        summary = json.loads((pf_output / 'simulation_summary.json').read_text())
        metadata['pf_emitted_counts'] = {key: summary[key] for key in
                                        ['requested_samples', 'written_samples', 'failed_duplicate_samples']}
        save_rows(args.output / 'policy-pf_emitted.csv.gz', emitted)
        for index, gamma in enumerate(['GC', 'G8']):
            tasks = [dict(sample=i, bank='bd50', tree=str(tree), binary=str(args.reference),
                          length=200, seed=1300001 + i * 101 + index * 100003,
                          gamma=gamma, alpha_seed=661923 + i * 131 + index * 170003)
                     for i, tree in enumerate(banks['bd50', 1])]
            rows = list(pool.map(legacy_policy, tasks, chunksize=8))
            save_rows(args.output / f'policy-python_{gamma}.csv.gz', rows)
            compare_groups(emitted, rows, f'current_python_preset_vs_legacy_{gamma}', 'policy', comparisons, constants)
            metadata[f'python_{gamma}_accepted'] = sum(row['accepted'] for row in rows)
            print(f'Actual Python {gamma} pipeline: {metadata[f"python_{gamma}_accepted"]}/{args.samples} retained', flush=True)

    # Bonferroni controls false detections across engine/policy/control metrics.
    m = len(comparisons) + len(constants)
    for row in comparisons:
        row['p_upper_adjusted'] = min(1.0, m * row['p_upper'])
        row['detected_difference'] = row['p_upper_adjusted'] < metadata['family_alpha']
        row['critical_distance'] = (math.sqrt(math.log(2 * m / metadata['family_alpha']) / row['n'])
                                    if row['n'] >= 458 else
                                    math.sqrt(2 * math.log(4 * m / metadata['family_alpha']) / row['n']))
    with (args.output / 'comparisons.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(comparisons[0]))
        writer.writeheader()
        writer.writerows(comparisons)
    (args.output / 'structural_constants.json').write_text(json.dumps(constants, indent=2) + '\n')
    metadata['comparisons'] = len(comparisons)
    metadata['bonferroni_family_size'] = m
    metadata['structural_constant_checks'] = len(constants)
    (args.output / 'metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    write_report(args.output, metadata, comparisons)


def write_report(output, metadata, comparisons):
    lines = ['# Standalone AliSim and Python MSA-prior validation', '',
             f"{metadata['samples']:,} independent MSAs per engine/case, using independent tree banks and sequence seeds. "
             'No production simulation behavior was changed.', '',
             '## Engine comparisons', '',
             '| Case | Standalone | Variable metric checks | Detected differences | Largest KS D |',
             '| --- | --- | ---: | ---: | ---: |']
    for name, *_ in CASES:
        for label in ['standalone_240', 'standalone_2226']:
            rows = [r for r in comparisons if r['comparison'] == f'{name}/{label}']
            lines.append(f"| {name} | {label} | {len(rows)} | {sum(r['detected_difference'] for r in rows)} | {max(r['ks_d'] for r in rows):.4f} |")
    lines += ['', '## Actual Python-pipeline output audit', '',
              f"PF `--preset python`, default LG, 200 root sites: {metadata['pf_emitted_counts']}", '',
              '| Python policy | Retained / requested | Detected metric differences | Worst KS D |',
              '| --- | ---: | ---: | ---: |']
    for gamma in ['GC', 'G8']:
        rows = [r for r in comparisons if r['comparison'].endswith(f'legacy_{gamma}')]
        lines.append(f"| {gamma} + POW indels + actual window selection | {metadata[f'python_{gamma}_accepted']} / {metadata['samples']} | {sum(r['detected_difference'] for r in rows)} / {len(rows)} | {max(r['ks_d'] for r in rows):.4f} |")
    controls = [r for r in comparisons if r['family'] == 'control']
    lines += ['', f"Negative control: changing Gamma shape 0.5 → 5 detected {sum(r['detected_difference'] for r in controls)} / {len(controls)} metric differences.", '',
              '## What is and is not established', '',
              '- Engine tests hold the model/length/indel settings and tree prior fixed across implementations, '
              'while using independent draws. No duplicate rejection or alignment cropping is applied to these unconditional MSAs.',
              '- Cases use LG, LG+G8{0.5}, LG+GC{0.5}, LG+I{0.2}, JC DNA, and the Python POW{1.7/100} indel law. '
              'Trees cover 4, 10, 50, and 200 tips, rooted Python BD and unrooted uniform trees.',
              '- The pinned standalone 2.4.0 build has the original source writers/tree loading. '
              'It uses OpenMP with one runtime thread and the portable kernel. The bridge uses a sequential build '
              'with explicit per-call RNG initialization; identical numeric seeds need not produce identical sequences.',
              '- The separately installed 2.2.2.6 binary is also tested, so version/build differences are not hidden.',
              '- Each data point is a whole MSA. Pair distances are summarized within each MSA over fixed tip pairs; '
              'sites, taxa, and pairs are not treated as independent replicate samples.',
              '- Measured features include lengths, gap fraction/columns/run lengths, invariant sites, entropy, '
              'between-site identity variation, tip distances, duplicate incidence/unique-tip fraction, and 20 state frequencies. '
              'Constant structural values are checked separately instead of inflating the statistical test count.',
              f"- There are {metadata['comparisons']} variable-metric comparisons and "
              f"{metadata['structural_constant_checks']} constant-metric checks; all "
              f"{metadata['bonferroni_family_size']} count in the Bonferroni denominator.",
              '- Statistical comparisons use two-sample KS with the equal-size DKWM bound `2 exp(-n D²)` '
              '(n >= 458); for smaller retained ensembles, the conservative union of one-sample DKW bounds '
              '`4 exp(-n D²/2)` is used. Bonferroni family alpha is 0.01 over all listed comparisons. '
              'Retained policy ensembles are truncated to equal sample counts. No detected difference is evidence '
              'at this finite resolution, not proof of equality of every aspect of the joint distribution.',
              '- The Python-policy audit runs the actual imported `alisim.py` and `find_unduplicated.py` functions. '
              'Only process launching is instrumented to pin executable, explicit independent seeds, one thread, '
              'portable kernel and log prefix. Alpha sampling, POW indels, ancestral simulation, trimming and '
              'the downstream unique-window selection are the actual legacy code.',
              '- The current preset only selects the Python **tree** pipeline. Its default MSA is LG without '
              'rate heterogeneity or indels. It retries duplicates on a fixed tree and retains full alignments. '
              'The Snakefile uses empirical gamma shapes (GC or configured G8), POW indels, and downstream '
              'unique-window selection/cropping to a fixed length, then filters ancestral records to tips.',
              '- Therefore agreeing raw tree priors and conditional engines does not establish full emitted-dataset '
              'prior equality. Different retry/filter/window policies can also change which trees are retained.',
              '- Optional user-supplied root sequences and custom model-definition files are not covered by these cases; '
              'they are not currently exposed by the bridge CLI. No universal claim over every AliSim option is made.', '',
              '## Reproduce', '', '```sh', 'python3 scripts/build_alisim_reference.py',
              'cargo build --release --example validate_alisim',
              'cargo build --release',
              '/path/to/legacy/.venv/bin/python scripts/validate_msa_priors.py --output validation/msa_priors-new', '```', '',
              'Raw per-MSA metrics/hashes are stored as gzip CSVs. See [comparisons.csv](comparisons.csv), '
              '[metadata.json](metadata.json), and [structural_constants.json](structural_constants.json).', '',
              'References: [Wei & Dudley, two-sample DKWM](https://arxiv.org/abs/1107.5356), '
              '[AliSim documentation](https://iqtree.github.io/doc/AliSim).', '']
    (output / 'REPORT.md').write_text('\n'.join(lines))


if __name__ == '__main__':
    main()
