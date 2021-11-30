#!/bin/bash
#SBATCH --ntasks=1 #each job has one task
#SBATCH --cpus-per-task=1 # each task uses 1 cpu
#SBATCH --partition=urtgen_24hrs
#SBATCH --mem-per-cpu=80000 #80GB

mkdir $1/perm/perm_1
cd $1/perm/perm_1

############################################################################################################
# Select pairs having pvalue below the corrected threshold and remove pairs where 2 SNPs are in hla region #
############################################################################################################

##Select pairs
Rscript $1/hla_filtering.R $1 $2 $3 $6

#Extract sign SNPs
$6/plink_17_oct --bfile $3 --noweb --allow-no-sex --make-bed --extract listSNPs.txt --out data

###############################
# Find LD of significant snps #
###############################

#Pairwise LD measures for multiple SNPs
$6/plink_17_oct --bfile data --allow-no-sex --r2

Rscript $1/LD_signSNPpairs.R $1 $2 $3 $4 $5 $6

Rscript $1/subset_snps.R $1 $3 $7 $6 $5 $4
