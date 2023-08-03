#!/bin/bash

# slurm arguments
#SBATCH -D /work/remelgad/
#SBATCH -J mine_xy
#SBATCH -t 0-24:00:00
#SBATCH --mem-per-cpu=50G

# load required modules
module load foss/2018b Python
source /home/remelgad/test/bin/activate

# recover directory set during job submissing (base project directory)
wdir="${SLURM_SUBMIT_DIR}"

# execute task
python "$wdir"/02_code/sample.py "$wdir"/config.yml "$id"
