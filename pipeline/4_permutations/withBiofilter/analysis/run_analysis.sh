#!/bin/bash
#SBATCH --ntasks=1 #each job has one task
#SBATCH --cpus-per-task=1 # each task uses 1 cpu
#SBATCH --partition=urtgen_24hrs
#SBATCH --mem-per-cpu=30000 #30GB


#for a in {401..1400}
# do

#dir_output=/home/mass/ifilesets/URT/UG_STG/PRIV/Team/Diane/SNPtoGene/4_permutations/chromatine
#dir=/home/mass/ifilesets/URT/UG_STG/PRIV/Team/Diane/SNPtoGene/4_permutations/chromatine/analysis
#cd $dir
#cp $dir_output/perm_$a/pvalues/sign_GenePairs.txt $dir/sign_GenePairs_$a.txt

#done

cd /home/mass/ifilesets/URT/UG_STG/PRIV/Team/Diane/SNPtoGene/4_permutations/withBiofilter/chromatine/analysis
find . -maxdepth 1 -type f -name 'sign_*.txt' -print0 | sort -zV | xargs -0 cat >concatGenePairs_tmp.txt
sed '/genePairs_names,MinP,SNP1,SNP2,k/d' ./concatGenePairs_tmp.txt >concatGenePairs.txt
rm concatGenePairs_tmp.txt
