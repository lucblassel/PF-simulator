#!/bin/bash

#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --job-name=pf-merge
#SBATCH --output=/lustre/fsn1/projects/rech/rdm/commun/Uniform_Exp0.2_50tips_LGG8_indels_2M/slurm_logs/pf-merge-%A-%a.log

set -euo pipefail

# Important because the regular tmpdir on job nodes is tiny
export TMPDIR="$JOBSCRATCH"

DATASET="/lustre/fsn1/projects/rech/rdm/commun/Uniform_Exp0.2_50tips_LGG8_indels_2M/train"
SCRIPT="/lustre/fswork/projects/rech/qwh/utj35ee/repos/PF-simulator/scripts/merge_packed.py"

python3 "$SCRIPT" \
  "$DATASET"/batches/batch-* \
  --out "$DATASET"

