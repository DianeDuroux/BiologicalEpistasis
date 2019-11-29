library(data.table)

args=(commandArgs(TRUE))
print(args)
threshold=as.integer(args[5])

for (a in 2:1000){
	output=fread(paste(args[1], "/perm/perm_",a,"/plink.epi.qt", sep=""))
	output=output[,c(2,4,7)]
	colnames(output)=c("SNP1", "SNP2", "pval")

	#Multiple testing correction
	output$pval=output$pval*(0.05/threshold)
	output=output[output$pval<=0.05,]

	fwrite(output, paste(args[1], "/perm/perm_",a,"/final_signSNPpairs.txt", sep=""), col.names=T)
}