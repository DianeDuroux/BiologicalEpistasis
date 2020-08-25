#!/bin/bash
#SBATCH --ntasks=1 #each job has one task
#SBATCH --cpus-per-task=1 # each task uses 1 cpu
#SBATCH --partition=urtgen_24hrs 
#SBATCH --mem-per-cpu=60000 #60GB

module load R/3.5.1

cd $1/pathway

########
# ATPM #
########

Rscript $1/pathway.R $1 $2 $3 $4 $5 $6 $7 $8
