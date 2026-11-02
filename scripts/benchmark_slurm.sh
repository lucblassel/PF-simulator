#!/usr/bin/env bash
#SBATCH --job-name=pf-bench
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=00:30:00
#SBATCH --array=0-95%4
#SBATCH --output=pf-bench-%A_%a.log

# Submit: sbatch [--partition=...] [--exclusive] benchmark_slurm.sh /path/pf-sim /scratch/results
# Grid: 2 workloads × 4 worker counts × 2 rate models × 2 formats × 3 seeds = 96 tasks.
# %1 avoids concurrent benchmark I/O; increase it to measure concurrent array load.
# Requires Bash + Python 3. Set PF_BENCH_SAMPLES to override the batch size,
# or PF_BENCH_KEEP_DATA=1 to retain datasets (logs, summaries and metrics always remain).
set -euo pipefail
binary=${1:?Usage: benchmark_slurm.sh /path/to/pf-sim /path/to/results}
root=${2:?Usage: benchmark_slurm.sh /path/to/pf-sim /path/to/results}
task=${SLURM_ARRAY_TASK_ID:?Submit as a Slurm array, or set SLURM_ARRAY_TASK_ID for a local smoke test}
[[ -x "$binary" ]] || {
  echo "Binary is not executable: $binary" >&2
  exit 2
}
[[ "$task" =~ ^[0-9]+$ ]] && ((task < 96)) || {
  echo "Array index must be 0..95" >&2
  exit 2
}

workers=(1 4 8 16)
gammas=(G8 GC)
formats=(packed files)
tip_counts=(50 200)
lengths=(500 5000)
sample_counts=(2048 64)
index=$task
repeat=$((index % 3))
index=$((index / 3))
format=${formats[index % 2]}
index=$((index / 2))
gamma=${gammas[index % 2]}
index=$((index / 2))
jobs=${workers[index % 4]}
case_index=$((index / 4))
tips=${tip_counts[case_index]}
length=${lengths[case_index]}
samples=${PF_BENCH_SAMPLES:-${sample_counts[case_index]}}
seed=$((20261102 + repeat)) # Same candidate stream across worker/format comparisons.
((jobs <= ${SLURM_CPUS_PER_TASK:-16})) || {
  echo "Requested workers exceed allocated CPUs" >&2
  exit 2
}

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
label="bench-${SLURM_ARRAY_JOB_ID:-manual}-${task}"
run_dir="$root/$label"
mkdir -p -- "$root"
mkdir -- "$run_dir" # Refuse to overwrite an earlier run.
command=("$binary" --out "$run_dir/data" --preset python --gamma "$gamma"
  --samples "$samples" --tips "$tips" --length "$length" --jobs "$jobs" --seed "$seed"
  --ancestral --indel 0.002,0.002 --indel-size 'POW{1.7/100},POW{1.7/100}'
  --require-distinct-gapped --require-distinct-ungapped=false --require-unique-window
  --window-policy python --max-attempts 1 --output-format "$format"
  --shard-size-bytes 1073741824 --batch-id "$label" --file-prefix "$label")
printf 'Task %s: tips=%s length=%s samples=%s workers=%s gamma=%s format=%s seed=%s\n' \
  "$task" "$tips" "$length" "$samples" "$jobs" "$gamma" "$format" "$seed"
launch=()
if [[ -n ${SLURM_JOB_ID:-} ]]; then launch=(srun --cpu-bind=cores); fi
"${launch[@]}" python3 - "$run_dir" "$task" "$repeat" "${command[@]}" <<'PY'
import csv, hashlib, json, os, resource, shutil, socket, subprocess, sys, time
from pathlib import Path

run = Path(sys.argv[1])
args = sys.argv[4:]
before = resource.getrusage(resource.RUSAGE_CHILDREN)
start = time.perf_counter()
with (run / "simulator.log").open("w") as log:
    result = subprocess.run(args, stdout=log, stderr=subprocess.STDOUT)
wall = time.perf_counter() - start
usage = resource.getrusage(resource.RUSAGE_CHILDREN)
cpu = usage.ru_utime + usage.ru_stime - before.ru_utime - before.ru_stime
summary_path = run / "data" / "simulation_summary.json"
summary = json.loads(summary_path.read_text()) if summary_path.exists() else {}
file_count = size = 0
for path in (run / "data").rglob("*"):  # Logical bytes, outside the timed interval.
    if path.is_file():
        file_count += 1
        size += path.stat().st_size

def option(name):
    return args[args.index(name) + 1]

completed, written = summary.get("completed_samples") or 0, summary.get("written_samples") or 0
row = dict(array_job=os.environ.get("SLURM_ARRAY_JOB_ID", "manual"),
           task=int(sys.argv[2]), repeat=int(sys.argv[3]), node=socket.gethostname(),
           binary_sha256=hashlib.sha256(Path(args[0]).read_bytes()).hexdigest(),
           cpus_per_task=os.environ.get("SLURM_CPUS_PER_TASK", "16"),
           effective_workers=summary.get("effective_jobs"),
           seed=option("--seed"), workers=option("--jobs"), tips=option("--tips"),
           length=option("--length"), gamma=option("--gamma"), format=option("--output-format"),
           requested=option("--samples"), completed=completed, written=written,
           failed=summary.get("failed_samples"), failed_duplicates=summary.get("failed_duplicate_samples"),
           failed_windows=summary.get("failed_window_samples"), wall_s=wall,
           simulator_wall_s=summary.get("wall_time_seconds"), cpu_s=cpu, avg_cpu_cores=cpu / wall,
           candidates_per_s=completed / wall, accepted_per_s=written / wall,
           output_bytes=size, output_files=file_count, output_mib_per_s=size / 2**20 / wall,
           exit_code=result.returncode, status=summary.get("status", "failed"))
with (run / "metrics.tsv").open("w", newline="") as file:
    writer = csv.DictWriter(file, fieldnames=list(row), delimiter="\t")
    writer.writeheader()
    writer.writerow(row)
if summary_path.exists():
    shutil.copy2(summary_path, run / "simulation_summary.json")
if os.environ.get("PF_BENCH_KEEP_DATA", "0") != "1" and result.returncode == 0:
    shutil.rmtree(run / "data")
print(f"{written}/{completed} accepted; {wall:.3f}s; {written / wall:.3f} accepted/s; metrics: {run / 'metrics.tsv'}")
sys.exit(result.returncode)
PY
