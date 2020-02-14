library(data.table)

args=(commandArgs(TRUE))
print(args)
nbSNPs=as.integer(args[2])

#Load data
SNPs_notLD_epistasis=fread(paste(args[3],"/perm_1/SNPs_notLD_epistasis.txt", sep=""))

#Multiple testing correction using FDR
nb_possiblePairs=choose(nbSNPs,2) #nb of possible pairs in the filtered dataset
SNPs_notLD_epistasis=data.frame(SNPs_notLD_epistasis$SNP_A, SNPs_notLD_epistasis$SNP_B, SNPs_notLD_epistasis$pval*nb_possiblePairs)
#SNPs_notLD_epistasis=data.frame(SNPs_notLD_epistasis$SNP_A, SNPs_notLD_epistasis$SNP_B, p.adjust(SNPs_notLD_epistasis$pval, method =c("bonferroni"), n = nb_possiblePairs)) #we could also use this p.adjust function
colnames(SNPs_notLD_epistasis)=c("SNP_A","SNP_B", "pval_adj")
SNPs_notLD_epistasis=SNPs_notLD_epistasis[SNPs_notLD_epistasis$pval_adj<=0.1,]
dim(SNPs_notLD_epistasis)
head(SNPs_notLD_epistasis)
fwrite(SNPs_notLD_epistasis, paste(args[1], "/pvalues/final_signSNPpairs.txt",sep=""), col.names=T, sep=" ")


