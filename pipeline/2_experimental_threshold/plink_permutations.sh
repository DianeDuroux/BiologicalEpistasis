#!/bin/bash
#SBATCH --partition=urtgen_24hrs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=4500
#SBATCH --array=1

options="--continuous"

Rscript $7/subset_snps.R $1 $2 $3 $5 $7


