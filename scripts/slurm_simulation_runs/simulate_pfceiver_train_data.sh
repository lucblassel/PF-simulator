#!/bin/bash

#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1
#SBATCH --array=0-20%21
#SBATCH --job-name=pf-sim
#SBATCH --output=/lustre/fsn1/projects/rech/rdm/commun/Uniform_Exp0.2_50tips_LGG8_indels_2M/slurm_logs/pf-sim-%A-%a.log

set -euo pipefail

module purge
module load gcc/15.2.0

# Important because the regular tmpdir on job nodes is tiny
export TMPDIR="$JOBSCRATCH"

DATASET_ROOT="/lustre/fsn1/projects/rech/rdm/commun/Uniform_Exp0.2_50tips_LGG8_indels_2M/"
PF_SIM="/lustre/fswork/projects/rech/qwh/utj35ee/repos/PF-simulator/target-jeanzay/release/pf-sim"

if [ "$SLURM_ARRAY_TASK_ID" == 20 ]; then
  dataset="${DATASET_ROOT}/val"
  n_samples=5000
else
  dataset="${DATASET_ROOT}/train"
  n_samples=100000
fi


batch="batch-${SLURM_ARRAY_TASK_ID}"
seed=$((20261103 + SLURM_ARRAY_TASK_ID))

"$PF_SIM" --out "$dataset/batches/$batch" \
  --output-format "packed" \
  --samples "$n_samples" --ensure-samples --max-attempts 1\
  --seed "$seed" \
  --batch-id "$batch" --file-prefix "$batch" \
  --jobs "$SLURM_CPUS_PER_TASK" \
  --tree-prior "uniform" --branch-prior "exponential" --branch-mean 0.2 --tips 50 \
  --model "LG" --gamma "G8" --length 500 --indel "0.002,0.002" --indel-size "POW{1.7/100},POW{1.7/100}"

