#!/usr/bin/env bash
#SBATCH --job-name=pf-buffer-bench
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --array=0-7%1
#SBATCH --output=pf-buffer-%A_%a.log

# Eight jobs total: four buffers × WORK/SCRATCH. One active job avoids self-contention.
# Submit with your partition/account, and optionally --exclusive:
# sbatch benchmark_buffers_slurm.sh /path/pf-sim "$WORK/pf-buffer-bench" "$SCRATCH/pf-buffer-bench"
set -euo pipefail
binary=${1:?Usage: benchmark_buffers_slurm.sh /path/pf-sim WORK_RESULTS SCRATCH_RESULTS}
work_root=${2:?Supply a results directory on WORK}
scratch_root=${3:?Supply a results directory on SCRATCH}
task=${SLURM_ARRAY_TASK_ID:?Submit as an array, or set SLURM_ARRAY_TASK_ID for a smoke test}
[[ "$task" =~ ^[0-9]+$ ]] && ((task < 8)) || { echo 'Array index must be 0..7' >&2; exit 2; }
[[ -x "$binary" ]] || { echo "Binary is not executable: $binary" >&2; exit 2; }
help=$("$binary" --help)
[[ "$help" == *--shard-buffer-bytes* ]] || {
    echo 'Rebuild pf-sim: this benchmark requires --shard-buffer-bytes.' >&2; exit 2;
}
binary=$(python3 -c 'import os,sys; print(os.path.realpath(sys.argv[1]))' "$binary")
buffers=(262144 1048576 4194304 16777216)
buffer=${buffers[task / 2]}
filesystem=work
root=$work_root
if ((task % 2)); then filesystem=scratch; root=$scratch_root; fi
workers=${PF_BENCH_WORKERS:-${SLURM_CPUS_PER_TASK:-16}}
samples=${PF_BENCH_SAMPLES:-10000}
repeats=${PF_BENCH_REPEATS:-1}
for number in "$workers" "$samples" "$repeats"; do
    [[ "$number" =~ ^[1-9][0-9]*$ ]] || { echo 'Workers, samples and repeats must be positive integers' >&2; exit 2; }
done
((workers <= ${SLURM_CPUS_PER_TASK:-16})) || { echo 'Workers exceed allocated CPUs' >&2; exit 2; }
label="buffer-${SLURM_ARRAY_JOB_ID:-manual}-${task}"
run_dir="$root/$label"
mkdir -p -- "$root"
mkdir -- "$run_dir" # Never overwrite a previous measurement.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
launch=()
if [[ -n ${SLURM_JOB_ID:-} ]]; then launch=(srun --cpu-bind=cores); fi
"${launch[@]}" python3 - "$binary" "$run_dir" "$task" "$filesystem" "$buffer" "$workers" "$samples" "$repeats" "$label" <<'PY'
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import shutil
import socket
import subprocess
import sys
import time

binary, destination, task, filesystem, buffer, workers, samples, repeats, label = sys.argv[1:]
root = Path(destination)
buffer, workers, samples, repeats = map(int, (buffer, workers, samples, repeats))
digest = hashlib.sha256(Path(binary).read_bytes()).hexdigest()
affinity = sorted(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else []
cpu_model = platform.processor()
if Path('/proc/cpuinfo').exists():
    cpu_model = next((line.split(':', 1)[1].strip() for line in Path('/proc/cpuinfo').read_text().splitlines()
                      if line.startswith('model name')), cpu_model)
(root / 'environment.json').write_text(json.dumps({
    'node': socket.gethostname(), 'platform': platform.platform(), 'cpu_model': cpu_model,
    'cpu_affinity': affinity, 'binary': binary, 'binary_sha256': digest,
    'slurm_job_id': os.environ.get('SLURM_JOB_ID'), 'partition': os.environ.get('SLURM_JOB_PARTITION'),
    'filesystem_label': filesystem, 'results_path': str(root.resolve()),
    'buffer_bytes_per_worker': buffer, 'workers': workers, 'samples': samples, 'repeats': repeats,
    'timing_scope': 'Executable wall time including simulation, buffered output and index finalization; excludes Slurm launch, enumeration and cleanup.'
}, indent=2))
for repeat in range(repeats):
    output = root / f'data-{repeat:02}'
    identity = f'{label}-r{repeat}'
    # All eight observation points use the same seed and simulation parameters.
    seed = 20261103 + repeat
    command = [binary, '--out', str(output), '--samples', str(samples), '--ensure-samples',
               '--tips', '200', '--length', '5000', '--jobs', str(workers), '--seed', str(seed),
               '--preset', 'python', '--gamma', 'G8', '--ancestral', '--indel', '0.002,0.002',
               '--indel-size', 'POW{1.7/100},POW{1.7/100}', '--require-distinct-gapped',
               '--require-distinct-ungapped=false', '--require-unique-window', '--window-policy', 'python',
               '--max-attempts', '1', '--output-format', 'packed', '--shard-size-bytes', str(1 << 30),
               '--shard-buffer-bytes', str(buffer), '--batch-id', identity, '--file-prefix', identity]
    (root / f'command-{repeat:02}.json').write_text(json.dumps(command, indent=2))
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    unix_start, started = time.time(), time.perf_counter()
    with (root / f'simulator-{repeat:02}.log').open('w') as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
    wall, unix_end = time.perf_counter() - started, time.time()
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    cpu = usage.ru_utime + usage.ru_stime - before.ru_utime - before.ru_stime
    summary_file = output / 'simulation_summary.json'
    summary = json.loads(summary_file.read_text()) if summary_file.exists() else {}
    paths = [path for path in output.rglob('*') if path.is_file()]
    size = sum(path.stat().st_size for path in paths)
    written, completed = summary.get('written_samples') or 0, summary.get('completed_samples') or 0
    success = result.returncode == 0 and summary.get('status') == 'completed' and written == samples
    row = dict(array_job=os.environ.get('SLURM_ARRAY_JOB_ID', 'manual'), task=int(task), repeat=repeat,
               filesystem=filesystem, node=socket.gethostname(), binary_sha256=digest,
               buffer_bytes=buffer, buffer_mib=buffer / 2**20, workers=workers,
               effective_workers=summary.get('effective_jobs'), cpus_per_task=os.environ.get('SLURM_CPUS_PER_TASK', '16'),
               affinity_cpus=len(affinity), seed=seed, tips=200, length=5000, gamma='G8',
               requested=samples, completed=completed, written=written, failed=summary.get('failed_samples'),
               failed_duplicates=summary.get('failed_duplicate_samples'), failed_windows=summary.get('failed_window_samples'),
               started_unix_seconds=unix_start, ended_unix_seconds=unix_end,
               wall_s=wall, simulator_wall_s=summary.get('wall_time_seconds'), cpu_s=cpu, avg_cpu_cores=cpu / wall,
               candidates_per_s=completed / wall, accepted_per_s=written / wall,
               output_bytes=size, output_files=len(paths), shards=sum(path.suffix == '.tar' for path in paths),
               output_mib_per_s=size / 2**20 / wall, exit_code=result.returncode,
               status='completed' if success else 'failed')
    with (root / 'metrics.csv').open('a', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=list(row))
        if repeat == 0:
            writer.writeheader()
        writer.writerow(row)
    if summary_file.exists():
        shutil.copy2(summary_file, root / f'summary-{repeat:02}.json')
    print(f'{filesystem}, {buffer / 2**20:g} MiB/worker, repeat {repeat}: '
          f'{written}/{samples} outputs; {wall:.3f}s; {written / wall:.2f} examples/s; '
          f'{size / 2**20 / wall:.2f} logical MiB/s', flush=True)
    if not success:
        sys.exit(result.returncode or 1) # Preserve failed output and stop this job.
    if os.environ.get('PF_BENCH_KEEP_DATA', '0') != '1':
        shutil.rmtree(output)
PY
