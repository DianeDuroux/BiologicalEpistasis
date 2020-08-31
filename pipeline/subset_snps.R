library(data.table)
library(plyr)
library(sensitivitymv)

###############
# Set options #
###############
args=(commandArgs(TRUE))
print(args)

################################
# Map SNP pairs to gene pairs  #
################################

mapping=fread(args[6])
colnames(mapping)=c("gene1", "gene2", "SNP_1", "SNP_2")

#Import summary statistics at the SNP-level
output <-list()
  output[[1]] = data.frame(fread(paste(args[1], "/perm/perm_1/final_signSNPpairs.txt",sep="")))
  colnames(output[[1]])=c("SNP_1", "SNP_2", "pvalue")
  output_SNPpairs=output[[1]][,1:2]
  output_SNPpairs <-  cbind(do.call(pmax, output_SNPpairs), do.call(pmin, output_SNPpairs)) #order SNPs for each SNP-pair
  output[[1]]=data.frame(output_SNPpairs, output[[1]]$pvalue)
  colnames(output[[1]])=c("SNP_1", "SNP_2", "pvalue")
  output[[1]]=merge(output[[1]], mapping, by=c("SNP_1", "SNP_2")) #Map SNP pairs to gene-pairs

  gene_pairs=output[[1]][,4:5]
  gene_pairs=unique(gene_pairs)
  colnames(gene_pairs)=c("gene1", "gene2")
  subset=merge(gene_pairs, mapping, by=c("gene1", "gene2"))
  head(subset)
  fwrite(subset, paste(args[1],"/perm/subset_mapping.txt", sep=""))

  snp_pairs=subset[,3:4]
  snp_pairs=unique(snp_pairs)
  fwrite(snp_pairs, paste(args[1],"/perm/subset_SNPpairs.txt", sep=""), sep=" ")

