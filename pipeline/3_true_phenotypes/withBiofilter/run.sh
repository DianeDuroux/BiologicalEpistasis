#!/bin/bash
#SBATCH --ntasks=1 #each job has one task
#SBATCH --cpus-per-task=1 # each task uses 1 cpu
#SBATCH --partition=urtgen_24hrs
#SBATCH --mem-per-cpu=40000 #40GB

dir=/home/mass/ifilesets/URT/UG_STG/PRIV/Team/Diane/SNPtoGene/3_truePhenotypes/withBiofilter/chromatine
softwares=/home/mass/ifilesets/URT/UG_STG/PRIV/Team/Diane/softwares
data=/home/mass/ifilesets/URT/UG_STG/PRIV/Data/IBD/dataset_response_residuals/CD_UC_CON_QCed_rel1_without_relatives_maf0.05_hwe0.001_Liu2015_232SNPs_LD0.75_chromatine_response
#dir_epistasis_output=/home/mass/ifilesets/URT/UG_STG/PRIV/Team/Diane/SNPtoGene/1_plink/true_phenotype/chromatine
dir_epistasis_output=/home/mass/ifilesets/URT/UG_STG/PRIV/Team/Diane/consensus_network/September/PLINK_response/Epistasis_detection/true_phenotypes/chromatine
B=999 #Number of permutations
mapFile=/home/mass/ifilesets/URT/UG_STG/PRIV/Data/IBD/mapping/mapping_3d.txt
pheno_ATPM=/home/mass/ifilesets/URT/UG_STG/PRIV/Team/Diane/SNPtoGene/1_plink/pheno
threshold="0.0000000102"
genes=/home/mass/ifilesets/URT/UG_STG/PRIV/Data/IBD/genes.tsv
biofModels=/home/mass/ifilesets/URT/UG_STG/PRIV/Data/IBD/biofilter_models_ensembl.tsv
pathway=/home/mass/ifilesets/URT/UG_STG/PRIV/Team/Diane/SNPtoGene/0_data/pathway/kegg_go_biocarta_canonical_header.txt
snpToGene=/home/mass/ifilesets/URT/UG_STG/PRIV/Data/IBD/FUMA/chromatin_interaction_mapping.tsv


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
jid5=$(sbatch  --dependency=afterany:$jid4_cut  pathway.sh "$dir" "$mapFile" "$B" "$snpToGene" "$genes" "$biofModels" "$pathway")
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


