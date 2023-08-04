#!/bin/bash

# ----------------------------------------------------------------------
# slurm arguments
# ----------------------------------------------------------------------

#SBATCH -J gee_lt_ts
#SBATCH -t 0-00:30:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-44929%50

# --------------------------------------------------------------------------#
# load required modules
# --------------------------------------------------------------------------#

module load foss/2018b Python
source /home/remelgad/test/bin/activate

# --------------------------------------------------------------------------#
# execute task
# --------------------------------------------------------------------------#

id="($SLURM_ARRAY_TASK_ID-1)"
wdir="${SLURM_SUBMIT_DIR}
python "$wdir"/02_code/extract_refletances.py "$wdir"/config.yml "$id"

