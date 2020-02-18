#!/bin/bash
#SBATCH --partition=urtgen_24hrs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=4500

softwares=/home/mass/ifilesets/URT/UG_STG/PRIV/Team/Diane/softwares
dir=/home/mass/ifilesets/URT/UG_STG/PRIV/Team/Diane/SNPtoGene/1_plink/true_phenotype/withoutBiofilter/eqtl
data=/home/mass/ifilesets/URT/UG_STG/PRIV/Data/IBD/dataset_response_residuals/withoutBiofilter/CD_UC_CON_QCed_rel1_without_relatives_maf0.05_hwe0.001_Liu2015_232SNPs_LD0.75_eqtl_continuous_withoutBiof


$softwares/plink_17_oct --bfile $data --epistasis --allow-no-sex --threads 6

