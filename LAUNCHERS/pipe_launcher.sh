#!/bin/bash

#SBATCH --job-name="pipe"
#SBATCH --output=pipe.out
#SBATCH --mem=500G
#SBATCH --cpus-per-task=16
#SBATCH --partition=sihp # cluster partition to use
#SBATCH --mail-user=<user-email> # user email
#SBATCH --mail-type=ALL # receive emails for all updates
#SBATCH --chdir=<output-dir>

input_seurat=.rds
output_seurat=.rds

#--------------------------------------------------------------------------------
# Do not modify below
#--------------------------------------------------------------------------------

# check if script is started via SLURM or bash
# if with SLURM: there variable '$SLURM_JOB_ID' will exist
# `if [ -n $SLURM_JOB_ID ]` checks if $SLURM_JOB_ID is not an empty string
if [ -n $SLURM_JOB_ID ];  then
    # check the original location through scontrol and $SLURM_JOB_ID
    CURRPATH=$(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2}')
else
    # otherwise: started with bash. Get the real location.
    CURRPATH=$(realpath $0)
fi

BASEPATH="${CURRPATH%/*/*}"
SCRIPTPATH="$BASEPATH/SCRIPTS/pipe.R"

. /local/env/envconda.sh # load conda
conda activate /home/genouest/irset/privaud/.conda/envs/seurat4 # activate R environment
Rscript $SCRIPTPATH $input_seurat $output_seurat # launcher pre-pipe script
