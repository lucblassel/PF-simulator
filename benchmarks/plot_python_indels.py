"""Plot worker scaling for the no-deduplication indel benchmark."""
import argparse
import csv
import json
import os
from pathlib import Path
import statistics
import tempfile

def main():
    cache = Path(tempfile.gettempdir()) / 'pf-simulator-plot-cache'
    os.environ.setdefault('MPLCONFIGDIR', str(cache / 'matplotlib'))
    os.environ.setdefault('XDG_CACHE_HOME', str(cache / 'xdg'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=Path('benchmarks/python_indels'))
    args = parser.parse_args()
    meta = json.loads((args.input / 'metadata.json').read_text())
    rows = list(csv.DictReader((args.input / 'runs.csv').open()))
    fig, axes = plt.subplots(1, len(meta['cases']), figsize=(10, 4.5), squeeze=False)
    for ax, (name, samples, tips, length) in zip(axes[0], meta['cases']):
        for backend, label, color in [('rust', 'Rust: tips', '#167647'),
                                       ('python_ancestral', 'Python: tips + ancestors', '#be5d12'),
                                       ('python_tips', 'Python: tips only', '#3a66a6')]:
            ys, lows, highs = [], [], []
            for jobs in [1, 4, 8]:
                times = [float(r['wall_seconds']) for r in rows
                         if r['case'] == name and r['backend'] == backend and int(r['jobs']) == jobs]
                med = statistics.median(times)
                q1, _, q3 = statistics.quantiles(times, n=4, method='inclusive')
                ys.append(med)
                lows.append(med-q1)
                highs.append(q3-med)
            ax.errorbar([1, 4, 8], ys, yerr=[lows, highs], marker='o', capsize=4,
                        label=label, color=color, linewidth=1.8)
        ax.set_xscale('log', base=2)
        ax.set_yscale('log')
        ax.set_xticks([1, 4, 8], ['1', '4', '8'])
        ax.set_xlabel('Process workers')
        ax.set_ylabel('Whole-run wall time (seconds)')
        ax.set_title(f'{samples} pairs · {tips} tips × {length:,} root sites', fontsize=11)
        ax.grid(which='major', alpha=0.2)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=3, frameon=False)
    fig.suptitle('Indels · LG · duplicate checks disabled · raw outputs', fontsize=13)
    fig.text(0.5, 0.91, f"Median and IQR of {meta['repeats']} runs; includes trees, startup and all output I/O.",
             ha='center', fontsize=9, color='#444444')
    fig.tight_layout(rect=(0, 0.09, 1, 0.88))
    fig.savefig(args.input / 'scaling.png', dpi=180)
    plt.close(fig)
    path = args.input / 'REPORT.md'
    report = path.read_text()
    marker = '![Indel worker scaling](scaling.png)'
    if marker not in report:
        report = report.replace('| Case | Pairs |', marker + '\n\n| Case | Pairs |', 1)
    command = '/path/to/legacy/.venv/bin/python benchmarks/plot_python_indels.py --input benchmarks/python_indels-new'
    if command not in report:
        report = report.replace('--output benchmarks/python_indels-new\n', '--output benchmarks/python_indels-new\n' + command + '\n')
    path.write_text(report)


if __name__ == '__main__':
    main()
