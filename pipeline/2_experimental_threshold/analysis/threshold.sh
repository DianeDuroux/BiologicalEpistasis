#!/bin/bash
#SBATCH --ntasks=1 #each job has one task
#SBATCH --cpus-per-task=1 # each task uses 1 cpu
#SBATCH --partition=urtgen_24hrs
#SBATCH --mem-per-cpu=40000 #40GB


dir_scripts=/home/mass/ifilesets/URT/UG_STG/PRIV/Team/Diane/SNPtoGene/2_experimentalThreshold/eqtl/analysis
for a in {1..400}
 do


dir_output=/home/mass/ifilesets/URT/UG_STG/PRIV/Team/Diane/SNPtoGene/2_experimentalThreshold/eqtl
dir=/home/mass/ifilesets/URT/UG_STG/PRIV/Team/Diane/SNPtoGene/2_experimentalThreshold/eqtl/analysis

cd $dir

cp $dir_output/perm_$a/perm/mappable_SNPpairs.txt $dir/mappable_SNPpairs_$a.txt

done



module load R/3.5.1

Rscript $dir_scripts/threshold.R "$dir_scripts"
