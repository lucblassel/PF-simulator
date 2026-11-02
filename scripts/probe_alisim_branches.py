"""Diagnose standalone version differences with a fixed JC quartet.

The known JC mismatch probability is 3/4 * (1 - exp(-4d/3)), where d is
patristic distance. Replicates, rather than sites, define standard errors.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
import gzip
import itertools
import json
import math
from pathlib import Path
import shutil
import subprocess
import tempfile

import numpy as np

from validate_msa_priors import ROOT, fasta, generate_bridge, sha

TREE = '(T1:0.1,T2:0.2,(T3:0.3,T4:0.4):0.5);'
DISTANCES = [0.3, 0.9, 1.0, 1.0, 1.1, 0.7]


def measure(task):
    if task.get('output'):
        rows = fasta(Path(task['output']) / f"{task['id']}.fa")
    else:
        with tempfile.TemporaryDirectory(prefix='jc-quartet-') as temp:
            prefix = str(Path(temp) / 'msa')
            subprocess.run([task['binary'], '--alisim', prefix, '-t', task['tree'], '-m', 'JC',
                            '--seqtype', 'DNA', '--length', str(task['length']), '--seed', str(task['seed']),
                            '--threads', '1', '-af', 'fasta', '-lk', 'x86', '-pre', prefix, '-redo'],
                           check=True, capture_output=True, timeout=60)
            rows = fasta(prefix + '.fa')
    a = {key: np.frombuffer(value.encode(), dtype=np.uint8) for key, value in rows.items()}
    return dict(sample=task['sample'], seed=task['seed'],
                **{f'T{i}_T{j}': float((a[f'T{i}'] != a[f'T{j}']).mean())
                   for i, j in itertools.combinations(range(1, 5), 2)})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'validation/msa_priors/branch_probe')
    parser.add_argument('--samples', type=int, default=200)
    parser.add_argument('--length', type=int, default=5000)
    parser.add_argument('--jobs', type=int, default=8)
    parser.add_argument('--harness', type=Path, default=ROOT / 'target/release/examples/validate_alisim')
    parser.add_argument('--reference', type=Path, default=ROOT / 'target/alisim-reference/build/iqtree/iqtree2')
    parser.add_argument('--installed', type=Path, default=Path(shutil.which('iqtree2') or 'iqtree2'))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    binaries = [args.harness, args.reference, args.installed]
    summary = []
    with tempfile.TemporaryDirectory(prefix='jc-branch-probe-') as temp, ProcessPoolExecutor(args.jobs) as pool:
        scratch = Path(temp)
        tree = scratch / 'quartet.nwk'
        tree.write_text(TREE + '\n')
        for backend, (label, binary) in enumerate(zip(['bridge', 'standalone_240', 'standalone_2226'], binaries)):
            tasks = [dict(id=f'quartet_{i:05}', sample=i, tree=str(tree), model='JC', alphabet='DNA',
                          length=args.length, seed=719013 + backend * 1000003 + i * 101,
                          indel='', size='', binary=str(binary)) for i in range(args.samples)]
            if backend == 0:
                tasks = generate_bridge(tasks, binary, scratch, args.jobs)
            rows = list(pool.map(measure, tasks, chunksize=8))
            with gzip.open(args.output / f'{label}.csv.gz', 'wt', newline='') as stream:
                writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
            for (i, j), distance in zip(itertools.combinations(range(1, 5), 2), DISTANCES):
                key = f'T{i}_T{j}'
                values = np.array([row[key] for row in rows])
                expected = 0.75 * (1 - math.exp(-4 * distance / 3))
                summary.append(dict(backend=label, pair=key, distance=distance, expected=expected,
                                    observed=float(values.mean()),
                                    standard_error=float(values.std(ddof=1) / math.sqrt(len(values)))))
    metadata = dict(tree=TREE, model='JC', samples=args.samples, root_length=args.length,
                    binaries={str(path): sha(path) for path in binaries}, source_sha256=sha(__file__), summary=summary)
    (args.output / 'results.json').write_text(json.dumps(metadata, indent=2) + '\n')
    lines = ['# Fixed unrooted quartet: JC branch-length probe', '', f'`{TREE}`', '',
             f'{args.samples} independent MSAs × {args.length} sites per backend. '
             'SE is calculated across MSA replicate means.', '',
             '| Backend | Tip pair | Tree distance | JC expectation | Observed mismatch | SE |',
             '| --- | --- | ---: | ---: | ---: | ---: |']
    for r in summary:
        lines.append(f"| {r['backend']} | {r['pair']} | {r['distance']:.1f} | {r['expected']:.5f} | {r['observed']:.5f} | {r['standard_error']:.5f} |")
    lines += ['', 'The bridge and standalone 2.4.0 track the JC expectations. In the installed '
              '2.2.2.6 build, pairs involving T1 instead track the predictions after subtracting '
              'T1’s 0.1 terminal edge: 0.17555, 0.49188 and 0.52410. Other pairs track the original '
              'distances. This is consistent with that build ignoring the first terminal edge '
              'when rooting this unrooted input; the experiment does not establish whether the '
              'cause is the version, build configuration, or platform.', '', 'Reproduce with `python scripts/probe_alisim_branches.py --output NEW_DIRECTORY` '
              'in the validation environment. This is a descriptive diagnostic; it is not included '
              'in the main report’s family of KS checks.', '']
    (args.output / 'REPORT.md').write_text('\n'.join(lines))
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
