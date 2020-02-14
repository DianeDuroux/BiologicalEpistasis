#!/bin/bash
#SBATCH --ntasks=1 #each job has one task
#SBATCH --cpus-per-task=1 # each task uses 1 cpu
#SBATCH --partition=urtgen_24hrs
#SBATCH --mem-per-cpu=50000 #40GB

for a in {1..400}
 do


dir_scripts=/home/mass/ifilesets/URT/UG_STG/PRIV/Team/Diane/SNPtoGene/2_experimentalThreshold/eqtl
softwares=/home/mass/ifilesets/URT/UG_STG/PRIV/Team/Diane/softwares
data=/home/mass/ifilesets/URT/UG_STG/PRIV/Data/IBD/dataset_response_residuals/CD_UC_CON_QCed_rel1_without_relatives_maf0.05_hwe0.001_Liu2015_232SNPs_LD0.75_eqtl_response
dir_epistasis_output=/home/mass/ifilesets/URT/UG_STG/PRIV/Team/Diane/SNPtoGene/1_plink/permutations/eqtl/plink_$a
pathway=/home/mass/ifilesets/URT/UG_STG/PRIV/Team/Diane/SNPtoGene/0_data/pathway
B=999 #Number of permutations
mapFile=/home/mass/ifilesets/URT/UG_STG/PRIV/Data/IBD/mapping/mapping_eqtl.txt
pheno_ATPM=/home/mass/ifilesets/URT/UG_STG/PRIV/Team/Diane/SNPtoGene/1_plink/pheno

module load R/3.5.1

mkdir $dir_scripts/perm_$a
cd $dir_scripts/perm_$a
dir=$dir_scripts/perm_$a


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
jid1=$(sbatch $dir_scripts//select_SNPpairs.sh "$dir" "$dir_epistasis_output" "$data" "$mapFile" "$softwares" "$dir_scripts")
jid1_cut=$(echo "$jid1" | cut -d ' ' -f4)

################
# Permutations #
################
jid2=$(sbatch  --dependency=afterany:$jid1_cut  $dir_scripts/plink_permutations.sh "$dir" "$data" "$B" "$softwares" "$mapFile" "$pheno_ATPM" "$dir_scripts")
jid2_cut=$(echo "$jid2" | cut -d ' ' -f4)

#########
# Clean #
#########
jid3=$(sbatch  --dependency=afterany:$jid2_cut  $dir_scripts/clean.sh "$dir" "$dir_scripts")

done
