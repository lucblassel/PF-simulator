"""Plot the paired Python/Rust benchmark medians and interquartile ranges."""
import argparse
import csv
import json
from pathlib import Path
import statistics

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=Path('benchmarks/python'))
    args = parser.parse_args()
    meta = json.loads((args.input / 'metadata.json').read_text())
    rows = list(csv.DictReader((args.input / 'runs.csv').open()))
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharex='row', sharey='row')
    colors = {'rust': '#167647', 'python_threads': '#be5d12', 'python_pool': '#3a66a6'}
    labels = {'rust': 'PF simulator', 'python_threads': 'Python: native AliSim threads',
              'python_pool': 'Python: process pool'}
    for row_index, tips in enumerate([50, 200]):
        cases = sorted((case for case in meta['cases'] if case[2] == tips), key=lambda case: case[1])
        for column, jobs in enumerate(meta['jobs']):
            ax = axes[row_index, column]
            for backend in ['rust', 'python_threads', 'python_pool']:
                if backend == 'python_pool' and jobs == 1:
                    continue
                xs, ys, lows, highs = [], [], [], []
                for name, samples, _, length in cases:
                    times = [float(r['wall_seconds']) for r in rows
                             if r['case'] == name and r['backend'] == backend and int(r['jobs']) == jobs]
                    med = statistics.median(times)
                    q1, _, q3 = statistics.quantiles(times, n=4, method='inclusive')
                    xs.append(samples)
                    ys.append(med)
                    lows.append(med - q1)
                    highs.append(q3 - med)
                ax.errorbar(xs, ys, yerr=[lows, highs], marker='o', markersize=4,
                            capsize=4, linewidth=1.8, label=labels[backend], color=colors[backend])
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xticks([case[1] for case in cases], [str(case[1]) for case in cases])
            ax.grid(which='major', alpha=0.2)
            ax.set_title(f'{jobs} worker/thread budget · {tips} tips × {cases[0][3]:,} sites', fontsize=10)
            if column == 0:
                ax.set_ylabel('Whole-run wall time (seconds)')
            if row_index == 1:
                ax.set_xlabel('Tree/MSA pairs requested')
    handles, legend_labels = axes[0, 2].get_legend_handles_labels()
    fig.legend(handles, legend_labels, loc='lower center', ncol=3, frameon=False)
    fig.suptitle('No indels · LG · duplicate checks · up to 20 attempts per tree', fontsize=14)
    fig.text(0.5, 0.935, f"Median and IQR of {meta['repeats']} runs; includes startup, trees, retries and file I/O. Both axes logarithmic.",
             ha='center', fontsize=10, color='#444444')
    fig.tight_layout(rect=(0, 0.055, 1, 0.91))
    fig.savefig(args.input / 'scaling.png', dpi=180)
    plt.close(fig)
    report = args.input / 'REPORT.md'
    text = report.read_text()
    marker = '![Whole-run scaling](scaling.png)'
    if marker not in text:
        text = text.replace('## Matched single-threaded and parallel wall times',
                            marker + '\n\n## Matched single-threaded and parallel wall times')
        report.write_text(text)


if __name__ == '__main__':
    main()
