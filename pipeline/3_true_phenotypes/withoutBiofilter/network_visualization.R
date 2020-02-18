library(data.table)
library(stringr)
library(data.table)
library(igraph)
library(RColorBrewer)


###############
# Set options #
###############
args=(commandArgs(TRUE))
print(args)

#########################
# Significant SNP-pairs #
#########################
Pvalue_SNPs=fread(paste(args[1], "/pvalues/sign_SNPpairs.txt", sep=""))
data_SNP=Pvalue_SNPs[,1:3]
data_SNP=unique(data_SNP)
head(Pvalue_SNPs)

##########################
# Significant gene-pairs #
##########################
Pvalue_genes=fread(paste(args[1], "/pvalues/sign_GenePairs_withoutThreshold.txt", sep=""))
genePair=Pvalue_genes[,1]
genePair=str_split_fixed(genePair$genePairs_names, " ", 2)
data=data.frame(genePair, Pvalue_genes[,2])
head(data)

########################################################

########
# SNP  #
########
node1=data.frame(data_SNP[,1])
node2=data.frame(data_SNP[,2])
dim(node1)
colnames(node1)=colnames(node2)
nodes_snp=rbind(node1,node2)
nodes_snp=data.frame(nodes_snp[!duplicated(nodes_snp),])
dim(nodes_snp)
data_SNP$width=-log10(as.numeric(data_SNP$pvalue))
links_snp=data_SNP[,c(1,2,4)]
colnames(nodes_snp)=c("id")
colnames(links_snp)=c("id1","id2","width")
net_SNP <- graph_from_data_frame(d=links_snp, vertices=nodes_snp)

#Color cluster in Gene-based exhaustive and visualize them in SNP based exhaustive
comp=components(net_SNP, mode = c("weak", "strong"))
nodes_snp=data.frame(comp$membership)
nodes_snp$gene=rownames(nodes_snp)
colnames(nodes_snp)=c("group","id")
nodes_snp=data.frame(nodes_snp$id, nodes_snp$group)
colnames(nodes_snp)=c("id","group")


########
# Gene #
########
node1=data.frame(data[,1])
node2=data.frame(data[,2])
colnames(node1)=colnames(node2)
nodes_gene=rbind(node1,node2)
nodes_gene=data.frame(nodes_gene[!duplicated(nodes_gene),])
data$width=data$MinP
links_gene=data[,c(1,2,4)]
colnames(nodes_gene)=c("id")
colnames(links_gene)=c("id1","id2","width")
net_gene <- graph_from_data_frame(d=links_gene, vertices=nodes_gene)

#Color cluster in Gene-based exhaustive and visualize them in SNP based exhaustive
#Create community for the gene-based network
comp=components(net_gene, mode = c("weak", "strong"))
nodes_gene=data.frame(comp$membership)
nodes_gene$gene=rownames(nodes_gene)
colnames(nodes_gene)=c("group","id")
nodes_gene=data.frame(nodes_gene$id, nodes_gene$group)
colnames(nodes_gene)=c("id","group")


#####################
# Visualization SNP #
#####################

V(net_SNP)$size <- 2
V(net_SNP)$frame.color <- "white"
V(net_SNP)$color <- nodes_snp$group
#V(net_SNP)$label <- ""
E(net_SNP)$arrow.mode <- 0
E(net_SNP)$width <- links_snp$width*0.7
E(net_SNP)$label <- ""
plot(net_SNP, main=paste("Chromatin FUMA mapping without biofilter models", sep=""), cex.main=50, margin=c(0,0,0,0))


######################
# Visualization Gene #
######################

#Parameters
V(net_gene)$size <- 1
V(net_gene)$frame.color <- "white"
V(net_gene)$color <- nodes_gene$group
V(net_gene)$label.cex = 0.5
#V(net_gene)$label <- ""
E(net_gene)$arrow.mode <- 0
E(net_gene)$width <- links_gene$width*0.5
E(net_gene)$label <- ""
plot(net_gene, main=paste("Chromatin FUMA mapping without biofilter models", sep=""), cex.main=50)



############
# Analysis #
############

#Largest component
comp_snp=components(net_SNP, mode = c("weak", "strong"))
max(comp_snp$csize)

comp_gene=components(net_gene, mode = c("weak", "strong"))
max(comp_gene$csize)


#degree
degSNP=degree(net_SNP,mode="all")
averageDegreeSNP=mean(degSNP)
medDegreeSNP=median(degSNP)
averageDegreeSNP
medDegreeSNP

degGene=degree(net_gene,mode="all")
averageDegreeGene=mean(degGene)
medDegreeGene=median(degGene)
averageDegreeGene
medDegreeGene


