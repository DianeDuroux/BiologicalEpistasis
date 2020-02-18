library(data.table)

set.seed(2019)

args=(commandArgs(TRUE))
print(args)

nb_perm=as.integer(args[3])+1
nb_perm

true_pheno=fread(paste(args[2],".fam", sep=""))
pheno=true_pheno[,6]

for(i in 1:1400){
  new_pheno <- pheno[sample(nrow(pheno)),]
  data_newpheno=data.frame(true_pheno[,1:2], new_pheno)
  fwrite(data_newpheno, paste(args[1], "/pheno_",i,".txt",sep=""), sep="\t", row.names=FALSE, quote=FALSE, col.names = FALSE)
}
