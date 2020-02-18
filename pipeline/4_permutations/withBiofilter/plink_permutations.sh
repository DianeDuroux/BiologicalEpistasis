#!/bin/bash
#SBATCH --partition=urtgen_24hrs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=4500
#SBATCH --array=1-10

options="--continuous"

for i in {1..100}
 do

j=$(($i+($SLURM_ARRAY_TASK_ID-1)*100))

  mkdir $1/perm/perm_$j
  cd $1/perm/perm_$j

  $4/plink_17_oct --bfile $2 --pheno $7/pheno_$j.txt --extract $1/perm/subset_SNPpairs.txt --epistasis --allow-no-sex

 done


