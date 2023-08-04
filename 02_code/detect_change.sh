#!/bin/bash

# ----------------------------------------------------------------------
# slurm arguments
# ----------------------------------------------------------------------

#SBATCH -J mine_year
#SBATCH -t 0-24:00:00
#SBATCH --mem-per-cpu=50G

# --------------------------------------------------------------------------#
# load required modules
# --------------------------------------------------------------------------#

module load foss/2018b Python
source /home/remelgad/test/bin/activate

# --------------------------------------------------------------------------#
# execute task
# --------------------------------------------------------------------------#

id="($SLURM_ARRAY_TASK_ID-1)"
wdir="${SLURM_SUBMIT_DIR}" # data paths
python "$wdir"/02_code/detect_change.py "$wdir"/config.yml
