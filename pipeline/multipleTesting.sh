#!/bin/bash
#SBATCH --partition=urtgen_5days
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4500

mkdir $1/noFilter
cd $1/noFilter

##################
# Input for ATPM #
##################

Rscript $1/multipleTesting.R $1 $2 $3

