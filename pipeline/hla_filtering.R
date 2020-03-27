library(data.table)

args=(commandArgs(TRUE))
print(args)
a = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

output=fread(paste(args[2],"/plink.epi.qt", sep=""))
colnames(output)=c("CHR1","SNP1","CHR2","SNP2", "OR_INT", "STAT", "pvalue")
output=output[,c(2,4,7)]
head(output)
dim(output)

##############################################
# Remove pairs when 2 SNPs are in hla region #
##############################################

hla=fread(paste(args[1], "/hla/hla.txt", sep=""))

#Create all possible combination between SNPs in hla region
expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))
all_possible_pairs <- expand.grid.df(x = hla$snp, y = hla$snp)
all_possible_pairs$indicator=rep(0,nrow(all_possible_pairs))
colnames(all_possible_pairs)=c("SNP1", "SNP2", "indicator")

#Merge with PLINK output
merge_hla=merge(all_possible_pairs, output, by=c("SNP1","SNP2"), all.y=T)
head(merge_hla)
pairs_to_keep_pval=merge_hla[is.na(merge_hla$indicator),c(1,2,4)]
dim(pairs_to_keep_pval)
pairs_to_keep=merge_hla[is.na(merge_hla$indicator),1:2]
dim(pairs_to_keep)

fwrite(pairs_to_keep, paste(args[1], "/perm/perm_1/listSNPs.txt",sep=""), sep=" ", col.names = T)
fwrite(pairs_to_keep_pval, paste(args[1], "/perm/perm_1/listSNPs_pval.txt",sep=""), sep=" ", col.names = T)


