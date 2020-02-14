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

#Import summary statistics at the SNP-level
output <-list()
  output[[1]] = data.frame(fread(paste(args[1], "/perm/perm_1/final_signSNPpairs.txt",sep="")))
  colnames(output[[1]])=c("SNP_1", "SNP_2", "pvalue")
  colnames(mapping)=c("SNP_1", "gene1")
  mapped_1=merge(output[[1]], mapping, by="SNP_1")
  colnames(mapping)=c("SNP_2", "gene2")
  mapped_2=merge(output[[1]], mapping, by="SNP_2")

  mapped_1=mapped_1[,c(1,4)]
  mapped_2=mapped_2[,c(1,4)]
  colnames(mapped_1)=c("SNP", "gene")
  colnames(mapped_2)=c("SNP", "gene")
  mapped=rbind(mapped_1, mapped_2)
  mapped=mapped[,2]
  mapped=data.frame(unique(mapped))
  colnames(mapped)=c("gene")
  colnames(mapping)=c("SNP", "gene")

head(mapped)
head(mapping)

  subset=merge(mapped, mapping, by="gene")
  fwrite(subset, paste(args[1],"/perm/subset_mapping.txt", sep=""))

  snp=data.frame(subset[,2])
  fwrite(snp, paste(args[1],"/perm/subset_SNPpairs.txt", sep=""), sep=" ")

