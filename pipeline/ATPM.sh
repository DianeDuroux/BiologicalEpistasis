#!/bin/bash
#SBATCH --ntasks=1 #each job has one task
#SBATCH --cpus-per-task=1 # each task uses 1 cpu
#SBATCH --partition=urtgen_24hrs 
#SBATCH --mem-per-cpu=128000 #128GB

module load R/3.5.1

cd $1/pvalues

########
# ATPM #
########

Rscript $1/ATPM_step1.R $1 $2 $3
Rscript $1/ATPM_step2.R $1 $2 $3
