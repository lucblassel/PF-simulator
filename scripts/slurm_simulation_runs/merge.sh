
#!/bin/bash

#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --job-name=pf-merge
#SBATCH --output=path/to/dataset/slurm_logs/pf-merge-%A-%a.log

set -euo pipefail

# Important because the regular tmpdir on job nodes is tiny
export TMPDIR="$JOBSCRATCH"

python3 path/to/merge_packed.py  /batches/batch-*   --out 
