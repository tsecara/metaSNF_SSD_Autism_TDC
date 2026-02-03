#!/bin/bash
#SBATCH --job-name=snf_clusters_%A
#SBATCH --time=24:00:00
#SBATCH --array=2-8
#SBATCH --mem-per-cpu=9000
#SBATCH --chdir=/projects/tsecara/metaSNF/code
#SBATCH --output=.../SNF_output1/logs/SNF_%j.out
#SBATCH --error=.../SNF_output1/logs/SNF_%j.err



# script to run SNF across parameters in parallel by cluster number 

# add R module
module load R

Rscript SNF_standard.R ${SLURM_ARRAY_TASK_ID}
