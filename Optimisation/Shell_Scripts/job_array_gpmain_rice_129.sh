#!/bin/bash

#SBATCH -p serial
#SBATCH --mem=8g
#SBATCH -J job_array_gpmain_rice_129
#SBATCH -a 1-15
#SBATCH --output=Logs/job_array_gpmain_rice_129.%A_%a.out
#SBATCH --error=Logs/job_array_gpmain_rice_129.%A_%a.err

mkdir -p Logs

source /etc/profile
module add matlab/2022a

matlab -nodisplay -nodesktop -r "run(fullfile('$(pwd)', '../gpmain/gpmain_rice_129_new.m')); quit"

echo This is job task ${SLURM_ARRAY_TASK_ID} running on compute node `uname -n`