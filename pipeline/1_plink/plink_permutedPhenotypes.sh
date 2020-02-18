#!/bin/bash
#SBATCH --partition=urtgen_24hrs
#SBATCH --array=1-200
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=4500

softwares=/home/mass/ifilesets/URT/UG_STG/PRIV/Team/Diane/softwares
dir=/home/mass/ifilesets/URT/UG_STG/PRIV/Team/Diane/SNPtoGene/1_plink/permutations/chromatine
phenotypes=/home/mass/ifilesets/URT/UG_STG/PRIV/Team/Diane/SNPtoGene/1_plink/pheno
data=/home/mass/ifilesets/URT/UG_STG/PRIV/Data/IBD/dataset_response_residuals/CD_UC_CON_QCed_rel1_without_relatives_maf0.05_hwe0.001_Liu2015_232SNPs_LD0.75_chromatine_response

j=$((1200+$SLURM_ARRAY_TASK_ID)) #To be adapted to perform 1400 permutations

mkdir $dir/plink_$j
cd $dir/plink_$j

$softwares/plink_17_oct --bfile $data --pheno $phenotypes/pheno_$j.txt --epistasis --allow-no-sex --threads 6