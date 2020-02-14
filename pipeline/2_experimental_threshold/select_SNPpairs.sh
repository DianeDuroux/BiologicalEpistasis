#!/bin/bash
#SBATCH --ntasks=1 #each job has one task
#SBATCH --cpus-per-task=1 # each task uses 1 cpu
#SBATCH --partition=urtgen_24hrs
#SBATCH --mem-per-cpu=30000 #30GB

mkdir $1/perm/perm_1
cd $1/perm/perm_1

############################################################################################################
# Select pairs having pvalue below the corrected threshold and remove pairs where 2 SNPs are in hla region #
############################################################################################################

##Select pairs
Rscript $6/hla_filtering.R $1 $2 $3 $5

#Extract sign SNPs
$5/plink_17_oct --bfile $3 --noweb --allow-no-sex --make-bed --extract listSNPs.txt --out data

###############################
# Find LD of significant snps #
###############################

#Pairwise LD measures for multiple SNPs
$5/plink_17_oct --bfile data --allow-no-sex --r2

Rscript $6/LD_signSNPpairs.R $1 $2 $3 $4 $5

