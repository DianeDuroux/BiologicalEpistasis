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

mapping=fread(args[4])
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
  fwrite(output[[1]], paste(args[1],"/perm/mappable_SNPpairs.txt", sep=""))

