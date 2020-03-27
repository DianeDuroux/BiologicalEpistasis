#!/bin/bash
#SBATCH --ntasks=1 #each job has one task
#SBATCH --cpus-per-task=1 # each task uses 1 cpu
#SBATCH --partition=urtgen_24hrs
#SBATCH --mem-per-cpu=40000 #40GB

dir=path_currentDirectory
softwares=Path_softwaresDirectory
data=path/data
dir_epistasis_output=path_EpistasisoutputDirectory
pathway=path_pathwayDirectory
B=999 #Number of permutations
mapFile=path_SNPpairsToGenepairs_mapping
pheno_ATPM=path_permutedPhenotypes
threshold="experimental_threshold"
snpToGene=path_SNPtoGene_mappingFile
genes=path_converterGeneName_HUGO_Ensembl
biofModels=path_geneModels


module load R/3.5.1

mkdir $dir/hla
mkdir $dir/perm
mkdir $dir/pvalues
mkdir $dir/pathway

###########################################
# Extract SNPs in hla region from dataset #
###########################################
$softwares/plink_17_oct --bfile $data  --noweb --make-bed --allow-no-sex --chr 6 --from-mb 25 --to-mb 34 --out $dir/hla/hla
cut -f2 $dir/hla/hla.bim > $dir/hla/hla.txt
sed -i '1s/^/snp\n/' $dir/hla/hla.txt

######################################################################
# Remove SNPs in hla region and in too high LD from epistasis output #
######################################################################
jid1=$(sbatch $dir/select_SNPpairs.sh "$dir" "$dir_epistasis_output" "$data" "$mapFile" "$threshold" "$softwares" "$B")
jid1_cut=$(echo "$jid1" | cut -d ' ' -f4)

################
# Permutations #
################
jid2=$(sbatch  --dependency=afterany:$jid1_cut  plink_permutations.sh "$dir" "$data" "$B" "$softwares" "$threshold" "$mapFile" "$pheno_ATPM")
jid2_cut=$(echo "$jid2" | cut -d ' ' -f4)

jid3=$(sbatch  --dependency=afterany:$jid2_cut  output_format.sh "$dir" "$data" "$B" "$softwares" "$threshold" "$mapFile")
jid3_cut=$(echo "$jid3" | cut -d ' ' -f4)


###################################################################################
# Calculate adjusted P-values, select significant pairs and visualize in networks #
###################################################################################
jid4=$(sbatch  --dependency=afterany:$jid3_cut  ATPM.sh "$dir" "$mapFile" "$B")
jid4_cut=$(echo "$jid4" | cut -d ' ' -f4)

###############################
# Pathway enrichment analysis #
###############################
jid5=$(sbatch  --dependency=afterany:$jid4_cut  pathway.sh "$dir" "$mapFile" "$B" "$pathway")
jid5_cut=$(echo "$jid5" | cut -d ' ' -f4)

#################
# Visualization #
#################
jid6=$(sbatch  --dependency=afterany:$jid5_cut  network_visualization.sh "$dir")
jid6_cut=$(echo "$jid6" | cut -d ' ' -f4)

#########
# Clean #
#########
jid7=$(sbatch  --dependency=afterany:$jid6_cut  clean.sh "$dir")


