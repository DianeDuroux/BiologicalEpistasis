library(data.table)
library(plyr)
library(sensitivitymv)

###############
# Set options #
###############
args=(commandArgs(TRUE))
print(args)
B=as.integer(args[3])
K=c(0.001, 0.01, 0.05)

S_kb=fread(paste(args[1],"/pvalues/estimated_Pvalue.txt", sep=""))

#########
# Min P #
#########
S_kb=as.data.frame(S_kb)
S_kb[,1]=as.factor(S_kb[,1])
nb_pairs=length(levels(S_kb[,1]))

cutoff_index=c(seq(2, (length(K)*3)-1, 3))
subset_S_kb=S_kb[,cutoff_index]
min_S_kb=list()

for(i in 1:nrow(S_kb)){
  min_K=ifelse( !all(is.na(subset_S_kb[i,])), which(S_kb[i,]==min(subset_S_kb[i,], na.rm=T))[1], 2) #allow NA, even for all cutoffs + when several pvalues are equal for a specific gene-pair and a specific output, we set k as max k corresponding to these pvalues
  min_S_kb[[i]]=data.frame(S_kb[i,c(min_K, min_K+1, min_K+2)],K[(min_K+1)/3]) #when several pvalues are equal for a specific gene-pair and a specific output, we set k as max k corresponding to these pvalues
  colnames(min_S_kb[[i]])=c("MinP", "SNP1", "SNP2", "k")
}
min_S_kb=data.frame(ldply(t(min_S_kb), rbind))
min_S_kb=data.frame(S_kb[,1], min_S_kb)
colnames(min_S_kb)=c("genePairs_names", "MinP", "SNP1", "SNP2", "k")
fwrite(min_S_kb[,1:5], paste(args[1],"/pvalues/min_S_kb.txt", sep=""))

############################################################################
# Overall adjusted P-value for the adaptative truncated product statistic  #
############################################################################
Pvalues=min_S_kb[,1:5]
Pvalues$MinP=0

for(i in 1:length(levels(min_S_kb$genePairs_names))){ #for each gene
  genePair=levels(min_S_kb$genePairs_names)[i]
  sub_data=min_S_kb[min_S_kb[,1]==genePair,]
  for (j in 1:(B+1)){ #for each output
    if(sub_data[j,2]==1){
      Pvalues[nb_pairs*(j-1)+i,2]=1
      Pvalues[nb_pairs*(j-1)+i,5]=NA
    } else { if(all(sub_data[j,2]<=sub_data[,2])){
      Pvalues[nb_pairs*(j-1)+i,2]=1/(B+1)
    } else { Pvalues[nb_pairs*(j-1)+i,2]=sum(sub_data[,2]<=sub_data[j,2])/(B+1)
    }
    }
  }
}

##################
# Export Pvalues #
##################

#All pvalues
fwrite(Pvalues, paste(args[1],"/pvalues/Pvalues.txt", sep=""), col.names=T)

#Pvalues per output
n=length(levels(min_S_kb$genePairs_names))
for(i in 1:(B+1)){
  sub_output=Pvalues[as.numeric(((i-1)*n)+1):as.numeric(((i-1)*n)+length(levels(min_S_kb$genePairs_names))),]
  fwrite(sub_output, paste(args[1],"/pvalues/Pvalues",i,".txt", sep=""), col.names=T)
}

