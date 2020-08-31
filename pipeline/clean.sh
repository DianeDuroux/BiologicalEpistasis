#!/bin/bash
#SBATCH --ntasks=1 #each job has one task
#SBATCH --cpus-per-task=1 # each task uses 1 cpu
#SBATCH --partition=urtgen_24hrs
#SBATCH --mem-per-cpu=30000 #30GB

dir=/home/mass/ifilesets/URT/UG_STG/PRIV/Team/Diane/consensus_network/September/MBMDR/noFilter/

#rm $1/hla/hla*
rm $1/perm/perm_1/data*
rm $1/perm/perm_1/listSNPs_pval*
rm $1/perm/perm_1/plink*