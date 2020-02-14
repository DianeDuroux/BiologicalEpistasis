#!/bin/bash
#SBATCH --partition=urtgen_24hrs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=4500

##################
# Input for ATPM #
##################

Rscript $7/output_format.R $1 $2 $3 $4 $5 $6