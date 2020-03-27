#!/bin/bash
#SBATCH --partition=urtgen_24hrs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=4500

##################
# Input for ATPM #
##################
#only keep SNP pairs existing in "true phenopype epistasis detection analysis" (file signSNPpairs.txt)
Rscript $1/output_format.R $1 $2 $3 $4 $5