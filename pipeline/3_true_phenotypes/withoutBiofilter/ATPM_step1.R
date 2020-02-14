library(data.table)
library(plyr)
library(sensitivitymv)

###############
# Set options #
###############
args=(commandArgs(TRUE))
print(args)
K=c(0.001, 0.01, 0.05)
B=as.integer(args[3])

################################
# Map SNP pairs to gene pairs  #
################################

mapping=fread(paste(args[1],"/perm/subset_mapping.txt", sep=""))
mapping=data.frame(mapping[,2], mapping[,1])

#Import summary statistics at the SNP-level
output <-list()

for (i in 1:(B+1)) {
  output[[i]] = data.frame(fread(paste(args[1],"/perm/perm_", i, "/final_signSNPpairs.txt", sep="")))
  if(nrow(output[[i]])==0){
  output[[i]]=data.frame(0,0,0,0,0,0)
  colnames(output[[i]])=c("SNP_1", "SNP_2", "pvalue", "gene1", "gene2", "combined ")
  fwrite(output[[i]], paste(args[1],"/pvalues/output",i,".txt", sep="")) } else { 
  colnames(output[[i]])=c("SNP_1", "SNP_2", "pvalue")
  colnames(mapping)=c("SNP_1", "gene1")
  mapped_1=merge(output[[i]], mapping, by="SNP_1")
  colnames(mapping)=c("SNP_2", "gene2")
  mapped_2=merge(mapped_1, mapping, by="SNP_2")
  mapped_2=mapped_2[as.character(mapped_2$gene1)!=as.character(mapped_2$gene2),]
  gene_pair=mapped_2[,4:5]
  gene_pair=cbind(do.call(pmax, gene_pair), do.call(pmin, gene_pair)) #order SNPs for each SNP-pair
  mapped_2=data.frame(mapped_2$SNP_1, mapped_2$SNP_2, mapped_2$pvalue, gene_pair)
  colnames(mapped_2)=c("SNP_1", "SNP_2", "pvalue", "gene1", "gene2")
  output[[i]]=mapped_2
output[[i]]=unique(output[[i]])
  output[[i]]$combined <- factor(paste(output[[i]]$gene1,output[[i]]$gene2))
  fwrite(output[[i]], paste(args[1],"/pvalues/output",i,".txt", sep=""))
 }
}


length(output)


##############################################
# Calculate the truncated product statistics #
##############################################
output[[1]]$combined=as.factor(output[[1]]$combined)
TP=matrix(nrow=length(levels(output[[1]]$combined))*length(output), ncol=length(K)*3) #create output matrix for the truncated product
genePairs_names=levels(output[[1]]$combined)
fwrite(as.data.frame(genePairs_names), paste(args[1],"/pvalues/genePairs_names.txt", sep=""))
rownames(TP)=rep(genePairs_names, length(output)) #First column = gene-pairs

#Truncated product
row=1
for (i in 1:(B+1)){ #for each output (observed + permuted)
  for(j in levels(as.factor(genePairs_names))){ #for each observed gene-pair
    col=1
    sub_data=output[[i]][output[[i]]$combined==j,]
    for(k in K){
      if(nrow(sub_data[sub_data$pvalue<=k,])==0){
        TP[row,col]=1 #1 if no pvalue <Kj in permuted dataset for SNP-pairs associated to gene-pair i
        TP[row,col+1]=NA
      }
      else { TP[row,col]=prod(sub_data[sub_data$pvalue<=k,3]) #calculate truncated product
      TP[row,col+1]=paste(unique(c(as.character(sub_data$SNP_1))),collapse=" ") #save original SNPs
      TP[row,col+2]=paste(unique(c(as.character(sub_data$SNP_2))),collapse=" ") #save original SNPs
      }
      col=col+3
    }
    row=row+1
  }
}
TP=data.frame(rep(genePairs_names, length(output)), TP)
fwrite(TP, paste(args[1],"/pvalues/truncatedProduct.txt", sep=""))

####################
# Estimated Pvalue #
####################

S_kb=TP
for(i in 1:length(K)){
  S_kb[,((i-1)*3+2)]=0
}
nb_pairs=length(levels(S_kb[,1]))

for(i in 1:length(levels(S_kb[,1]))){ #for each gene
  genePair=levels(S_kb[,1])[i]
  sub_data=TP[TP[,1]==genePair,]
  for (j in 1:(B+1)){ #for each output
    for(k in 1:length(K)){ #for each cutoff
      S_kb[nb_pairs*(j-1)+i,3*(k-1)+2]=sum(as.numeric(as.character(sub_data[j,3*(k-1)+2]))>=as.numeric(as.character(sub_data[,3*(k-1)+2])))/(B+1)
    }
  }
}
fwrite(S_kb, paste(args[1],"/pvalues/estimated_Pvalue.txt", sep=""))

