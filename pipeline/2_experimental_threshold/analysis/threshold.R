library(data.table)
args=(commandArgs(TRUE))
print(args)


minp=data.frame()
for(i in 1:400){
  data=fread(paste(args[1], "/mappable_SNPpairs_", i, ".txt", sep=""))
  data=data[order(data$pvalue),]
  minp=rbind(minp, data[1,])
}

minp=minp[order(minp$pvalue),]
minp=minp[as.character(minp$gene1)!=as.character(minp$gene2),]

minp[1:30,]
nrow(minp)

threshold=minp[(5*nrow(minp))/100,3]


print("threshold T*")
threshold

print("N*")
0.05/threshold
