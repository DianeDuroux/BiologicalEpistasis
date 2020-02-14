library(data.table)
library(stringr)
library(plyr)

###############
# Set options #
###############
args=(commandArgs(TRUE))
print(args)
K=c(0.001, 0.01, 0.05, 0.1, 0.2)
B=as.integer(args[3])

################################
# Export significant SNP pairs  #
################################
#pval<0.05 with Bonferroni correction and mappable to gene-pairs with biofilter

output_true=fread(paste(args[1],"/pvalues/output1.txt", sep=""))
output_true=output_true[output_true$pvalue<0.05,]
dim(output_true)
fwrite(output_true, paste(args[1],"/pvalues/sign_SNPpairs.txt", sep=""))

############################
# Import gene pairs output #
############################

truePhenoOutput=data.frame(fread(paste(args[1],"/pvalues/Pvalues1.txt", sep=""))) #Load gene-pairs pvalues associated to the true phenotype
nb_pairs=as.integer(nrow(truePhenoOutput))
Pvalues=fread(paste(args[1],"/pvalues/Pvalues.txt", sep="")) #Load all gene-pairs pvalues
Pvalues=Pvalues[,1:4]
perm=c()
for (i in 1:(B+1)){
  temp=rep(i,nb_pairs)
  perm=c(perm,temp)
}
Pvalues=data.frame(Pvalues,perm)
colnames(Pvalues)=c("genePair", "P", "SNP1", "SNP2", "perm")

# Corect for multiple testing (fdr nb gene pairs tested)
truePhenoOutput$MinP=truePhenoOutput$MinP*nb_pairs
Pvalues$P=Pvalues$P*nb_pairs

fwrite(truePhenoOutput, paste(args[1],"/pvalues/sign_GenePairs_withoutThreshold.txt", sep=""))

#Select only gene pairs when the pairs are significant for the true phenotype
signPairs=truePhenoOutput[truePhenoOutput$MinP<0.05,] #significant gene pairs
fwrite(signPairs, paste(args[1],"/pvalues/sign_GenePairs.txt", sep=""))

truePhenoOutput=truePhenoOutput[truePhenoOutput$MinP<0.1,] #allow less significant pairs for ATPM
genePair=data.frame(truePhenoOutput$genePair)
colnames(genePair)="genePair"
Pvalues=merge(Pvalues, genePair, by="genePair", all.y=T)
Pvalues=Pvalues[order(Pvalues$perm),]

####################
# Pathway analysis #
####################

###########################
# STEP 1: Load IBD results #
###########################
genes_IBD=fread(paste(args[1],"/pvalues/sign_GenePairs.txt", sep="")) #Load gene-pairs pvalues associated to the true phenotype
genePair=genes_IBD[,1]
genes_IBD=str_split_fixed(genePair$genePairs_names, " ", 2)


##################################################
# STEP 2: Create the combined functional network #
##################################################
library(igraph)
ppi=fread(paste(args[4],"/ppi", sep=""))
coexpr=fread(paste(args[4],"/coexpr", sep=""))

func_network=rbind(coexpr, ppi)
func_network=unique(func_network)
func_network$int=rep(1,nrow(func_network))

nodes_func_network=c(func_network$gene1, func_network$gene2)
nodes_func_network=unique(nodes_func_network)

net_func_network <- graph_from_data_frame(d=func_network, vertices=nodes_func_network)


#####################################
# STEP 3: create clusters of genes  #
#####################################
##For each epistatic gene pair, create a cluster of genes that are densely connected with each other in terms of epistatic interactions
#and functional (PPI and co-expression) interactions

clusters_allPairs=list()
for(i in 1:nrow(genes_IBD)){ #for each gene pair
  gene1=genes_IBD[i,1]
  gene2=genes_IBD[i,2]
  if((gene1 %in% nodes_func_network)==T & (gene2 %in% nodes_func_network)==T){
    sp=get.all.shortest.paths(graph=net_func_network, #find ALL the shortest paths between the two gene in the functional network (length>1)
                              from=gene1, 
                              to = gene2)
    df <- as.data.frame(t(sapply(sp$res, as_ids)))
    
    genes=data.frame(c(levels(unlist(df))))
    if(length(genes)!=0){
      #retain only genes with a loosely epistatic interaction with at least one other gene on these paths (??)
      expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))#create all possible pairs
      all_possible_pairs <- expand.grid.df(x = genes, y = genes)
      colnames(all_possible_pairs)=c("V1", "V2")
      cluster=list(merge(all_possible_pairs, genes_IBD, by=c("V1", "V2")))
      clusters_allPairs=c(clusters_allPairs, cluster)
    }
  }
}

clusters_allPairs

###############################################
# STEP 4: enrichment analysis of each cluster #
###############################################
#Load GSEA pathways
pathway=fread(paste(args[4],"/kegg_go_biocarta_canonical_header.txt", sep=""))
pathway=pathway[,-2]
pathway=unique(pathway)
pathway=data.frame(t(pathway))

nbGenes_PPI_coexpr=length(nodes_func_network)
nbGenes_pathway=nrow(pathway)*(ncol(pathway)-1)-sum(is.na(pathway))
nbGenes_total=nbGenes_PPI_coexpr+nbGenes_pathway

pathway_analysis=data.frame()
for (i in 1:length(clusters_allPairs)){ #For each cluster
  for(j in 1:ncol(pathway)){ # and each pathway
    genesInCluster=data.frame(t(clusters_allPairs[[i]]))
    colnames(genesInCluster)=c("gene")
    genesInPath=data.frame(pathway[complete.cases(pathway[,j]), j])
    colnames(genesInPath)=c("gene")
    genesInBoth=merge(genesInCluster, genesInPath, by="gene")
    if(nrow(genesInBoth)!=0){
      contingency=matrix(nrow=2, ncol=2) #2 × 2 contingency table for the genes in the cluster or not, and in the pathway or not
      contingency[1,1]=nrow(genesInBoth)
      contingency[1,2]=nrow(genesInPath)-nrow(genesInBoth)
      contingency[2,1]=nrow(genesInCluster)-nrow(genesInBoth)
      contingency[2,2]=nbGenes_total-contingency[1,1]-contingency[1,2]-contingency[2,1] #the background set contains all genes contained in the PPI network, co-expression network or GSEA pathways
      chi2=chisq.test(contingency,correct = TRUE) #Based on this contingency table, we computed a corrected chi-square statistic
      output=data.frame(t(levels(genesInBoth$gene)), genesInPath[1,1], chi2$p.value) #and the corresponding P-value
      pathway_analysis=rbind(pathway_analysis, output)
    }
  }
}


#These raw P-values were then corrected by Bonferroni correction, based on the total number of GSEA pathway terms and the number
#of clusters identified based on the respective SNP–gene association method.
pathway_analysis$correcter_pval=pathway_analysis$chi2.p.value*(ncol(pathway)+length(clusters_allPairs))
pathway_analysis=pathway_analysis[pathway_analysis$correcter_pval<=0.05,]
pathway_analysis=pathway_analysis[order(pathway_analysis$correcter_pval),]

unique_pathway=data.frame(unique(pathway_analysis$genesInPath.1..1.))


fwrite(pathway_analysis, paste(args[1],"/pathway/pathwayEnrichment.txt", sep=""), col.names=T)
fwrite(unique_pathway, paste(args[1],"/pathway/unique_pathway.txt", sep=""), col.names=T)
