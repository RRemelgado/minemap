#!/bin/bash
#SBATCH -D /work/remelgad/

# ----------------------------------------------------------------------
# slurm arguments
# ----------------------------------------------------------------------

#SBATCH -J gee_lt_ts
#SBATCH -t 0-00:30:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-44929%50

# ----------------------------------------------------------------------
# setup job report
# ----------------------------------------------------------------------

#SBATCH -o /work/remelgad/%x-%j-%a_log.txt

# --------------------------------------------------------------------------#
# load required modules
# --------------------------------------------------------------------------#

module load foss/2018b Python
source /home/remelgad/test/bin/activate

# --------------------------------------------------------------------------#
# execute task
# --------------------------------------------------------------------------#

id="$(printf "%05d" $(($SLURM_ARRAY_TASK_ID-1)))"

wdir=/data/idiv_meyer/01_projects/Ruben/GlobES/tmp/mining_sites

python "$wdir"/02_code/extract_refletances.py "$wdir"/config.yml "$id"

