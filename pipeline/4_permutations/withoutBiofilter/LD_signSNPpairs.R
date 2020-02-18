library(data.table)

args=(commandArgs(TRUE))
print(args)
a = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
threshold=as.numeric(args[5])

#Load PLINK LD output and the selected pairs
LD=fread(paste(args[1], "/perm/perm_1/plink.ld",sep=""))
signSNPpairs=fread(paste(args[1], "/perm/perm_1/listSNPs_pval.txt",sep=""))
colnames(signSNPpairs)=c("SNP_A","SNP_B", "pval")

#Re-order SNPs in each pair
LD=LD[,c(3,6,7)]
LD_SNPpair=LD[,1:2]
LD_SNPpair <-  cbind(do.call(pmax, LD_SNPpair), do.call(pmin, LD_SNPpair)) #order SNPs for each SNP-pair
LD=data.frame(LD_SNPpair, LD[,3])
colnames(LD)=c("SNP_A","SNP_B", "R2")

signSNPpairs_SNPpair=signSNPpairs[,1:2]
signSNPpairs_SNPpair <-  cbind(do.call(pmax, signSNPpairs_SNPpair), do.call(pmin, signSNPpairs_SNPpair)) #order SNPs for each SNP-pair
signSNPpairs=data.frame(signSNPpairs_SNPpair, signSNPpairs[,3])
colnames(signSNPpairs)=c("SNP_A","SNP_B", "Pvalue")

#Combine the two datasets
merged=merge(LD,signSNPpairs, by=c("SNP_A", "SNP_B"), all.y=T)
merged=merged[order(-merged$R2),]
merged[is.na(merged$R2),3] <- 0

#Pairs in LD
SNPsLD=merged[merged$R2>0.75,c(1,2,4)]
print("Dim SNPs in LD")
dim(SNPsLD)

##########################
# Select pairs not in LD #
##########################

SNPs_notLD=merged[merged$R2<0.75,c(1,2,4)]
SNPs_notLD_epistasis=merge(SNPs_notLD, signSNPpairs, by=c("SNP_A", "SNP_B"))
SNPs_notLD_epistasis=SNPs_notLD_epistasis[,c(1,2,3)]
class(SNPs_notLD_epistasis$pval)
colnames(SNPs_notLD_epistasis)=c("SNP_A","SNP_B", "pval")
print("head SNP pairs not in LD")
head(SNPs_notLD_epistasis)

###########################################
# Reduce searching space for permutations #
###########################################

#Multiple testing correction 
SNPs_notLD_epistasis=data.frame(SNPs_notLD_epistasis$SNP_A, SNPs_notLD_epistasis$SNP_B, SNPs_notLD_epistasis$pval*(0.05/threshold))
colnames(SNPs_notLD_epistasis)=c("SNP_A","SNP_B", "pval_adj")
SNPs_notLD_epistasis=SNPs_notLD_epistasis[SNPs_notLD_epistasis$pval_adj<=0.05,]
dim(SNPs_notLD_epistasis)
head(SNPs_notLD_epistasis)
fwrite(SNPs_notLD_epistasis, paste(args[1], "/perm/perm_1/final_signSNPpairs.txt",sep=""), col.names=T, sep=" ")
fwrite(SNPs_notLD_epistasis, paste(args[1], "/perm/perm_1/isEmpty.txt",sep=""), col.names=F, sep=" ")

data=data.frame()
fwrite(data, paste(args[1], "/perm/subset_SNPpairs.txt",sep=""), col.names=F, sep=" ")

