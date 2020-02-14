#!/bin/bash
#SBATCH --ntasks=1 #each job has one task
#SBATCH --cpus-per-task=1 # each task uses 1 cpu
#SBATCH --partition=kosmos

dir=/home/mass/ifilesets/URT/UG_STG/PRIV/Team/Diane/consensus_network/September/PLINK_response/Epistasis_detection/permutations/pheno2
data=/home/mass/ifilesets/URT/UG_STG/PRIV/Data/IBD/dataset_response_residuals/CD_UC_CON_QCed_rel1_without_relatives_maf0.05_hwe0.001_Liu2015_232SNPs_LD0.75_chromatine_response
module load R/3.5.1

###############################
# Create 99 random phenotypes #
###############################

Rscript $dir/random_pheno.R "$dir" "$data"
