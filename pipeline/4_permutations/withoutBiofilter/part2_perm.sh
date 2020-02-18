#!/bin/bash
#SBATCH --ntasks=1 #each job has one task
#SBATCH --cpus-per-task=1 # each task uses 1 cpu
#SBATCH --partition=urtgen_24hrs
#SBATCH --mem-per-cpu=30000

dir_scripts=/home/mass/ifilesets/URT/UG_STG/PRIV/Team/Diane/SNPtoGene/4_permutations/withoutBiofilter/chromatine
softwares=/home/mass/ifilesets/URT/UG_STG/PRIV/Team/Diane/softwares
data=/home/mass/ifilesets/URT/UG_STG/PRIV/Data/IBD/dataset_response_residuals/withoutBiofilter/CD_UC_CON_QCed_rel1_without_relatives_maf0.05_hwe0.001_Liu2015_232SNPs_LD0.75_chromatine_continuous_withoutBiof
dir_epistasis_output=/home/mass/ifilesets/URT/UG_STG/PRIV/Team/Diane/SNPtoGene/1_plink/permutations/withoutBiofilter/chromatine/plink_$a
pathway=/home/mass/ifilesets/URT/UG_STG/PRIV/Team/Diane/SNPtoGene/0_data/pathway
B=999 #Number of permutations
mapFile=/home/mass/ifilesets/URT/UG_STG/PRIV/Data/IBD/FUMA/noFilter_chromatine_FUMA.tsv
pheno_ATPM=/home/mass/ifilesets/URT/UG_STG/PRIV/Team/Diane/SNPtoGene/1_plink/pheno
threshold="0.00000000009011"


cd $dir_scripts/perm_$1
dir=$dir_scripts/perm_$1


if [ -s $dir/perm/subset_SNPpairs.txt ] 

then
	echo "file has some data."
	jid2=$(sbatch  $dir_scripts/plink_permutations.sh "$dir" "$data" "$B" "$softwares" "$threshold" "$mapFile" "$pheno_ATPM" "$dir_scripts")
	jid2_cut=$(echo "$jid2" | cut -d ' ' -f4)

	jid3=$(sbatch  --dependency=afterany:$jid2_cut  $dir_scripts/output_format.sh "$dir" "$data" "$B" "$softwares" "$threshold" "$mapFile" "$dir_scripts")
	jid3_cut=$(echo "$jid3" | cut -d ' ' -f4)

	###################################################################################
	# Calculate adjusted P-values, select significant pairs and visualize in networks #
	###################################################################################
	jid4=$(sbatch  --dependency=afterany:$jid3_cut  $dir_scripts/ATPM.sh "$dir" "$mapFile" "$B" "$dir_scripts")
	jid4_cut=$(echo "$jid4" | cut -d ' ' -f4)

	###############################
	# Pathway enrichment analysis #
	###############################
	jid5=$(sbatch  --dependency=afterany:$jid4_cut  $dir_scripts/pathway.sh "$dir" "$mapFile" "$B" "$pathway" "$dir_scripts")
	jid5_cut=$(echo "$jid5" | cut -d ' ' -f4)

	#################
	# Visualization #
	#################
	jid6=$(sbatch  --dependency=afterany:$jid5_cut  $dir_scripts/network_visualization.sh "$dir" "$dir_scripts")
	jid6_cut=$(echo "$jid6" | cut -d ' ' -f4)

	#########
	# Clean #
	#########
	jid7=$(sbatch  --dependency=afterany:$jid6_cut  $dir_scripts/clean.sh "$dir" "$dir_scripts")
	
else
	echo "file is empty."
        jid2=$(sbatch  $dir_scripts/clean.sh "$dir" "$dir_scripts")

fi
