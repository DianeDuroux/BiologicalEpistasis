library(data.table)
library(stringr)
library(plyr)
library(igraph)

args=(commandArgs(TRUE))
print(args)

###############
# Set options #
###############
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

##############################
# Create SNP to gene mapping #
##############################

snpToGene=fread(args[4]) #import snp to gene mapping
snpToGene=snpToGene[,1:2]
genes=fread(args[5]) #import genes ensg -> HUGO
genes=genes[,c(1, 4)]
colnames(genes)=c("ensg", "gene")
snpToGene=merge(snpToGene, genes, by="ensg")
snpToGene=snpToGene[,2:3]
snpToGene=unique(snpToGene)
nodes_func_network=data.frame(unique(snpToGene$gene))
colnames(nodes_func_network)="gene"

######################
# Create gene weight #
######################

nbGenes=length(levels(as.factor(snpToGene$gene)))
nbSNPs=nrow(snpToGene)
geneWeight=data.frame(unique(snpToGene$gene))
geneWeight$weight=rep(NA, nrow(geneWeight))

a=1
for(i in levels(geneWeight$unique.snpToGene.gene.)){
  geneWeight[a,2]=nrow(snpToGene[snpToGene$gene==i,])
  a=a+1
}

geneWeight=na.omit(geneWeight)
geneWeight$weight2=(1/geneWeight$weight)+(1-mean(1/geneWeight$weight))
summary(geneWeight$weight2)
geneWeight=geneWeight[,c(1,3)]
colnames(geneWeight)=c("gene", "weight")

############################
# STEP 1: Load IBD results #
############################
gene_epi=fread(paste(args[1],"/pvalues/sign_GenePairs_withoutThreshold.txt", sep="")) #Load gene-pairs pvalues associated to the true phenotype
genePair=gene_epi[,1]
gene_epi=str_split_fixed(genePair$genePairs_names, " ", 2)

genes_network=data.frame(c(levels(unlist(data.frame(gene_epi)))))
colnames(genes_network)="genes"

########################################################################
# STEP 2: Create the combined functional network, ie biofilter network #
########################################################################
biof=fread(args[6])
biof=biof[,1:2]
func_network <-  data.frame(cbind(do.call(pmax, biof), do.call(pmin, biof))) #order gene for each pair
colnames(func_network)=c("gene1", "gene2")
nodes_func_network=c(as.character(func_network$gene1), as.character(func_network$gene2))
nodes_func_network=as.factor(unique(nodes_func_network))


#####################################
# STEP 3: create clusters of genes  #
#####################################
##For each epistatic gene pair, create a cluster of genes that are densely connected with each other in terms of epistatic interactions
clusters_allPairs=list()

for(i in 1:nrow(gene_epi)){ #for each gene pair
  genea=gene_epi[i,1]
  geneb=gene_epi[i,2]
  if((genea %in% nodes_func_network)==T & (geneb %in% nodes_func_network)==T){
    sub_network=func_network[!(func_network$gene1==genea & func_network$gene2==geneb) | (func_network$gene1==geneb & func_network$gene2==genea), ] #remove direct link
    sub_net_func_network <- graph_from_data_frame(d=sub_network, vertices=nodes_func_network)
    sp=get.all.shortest.paths(graph=sub_net_func_network, #find ALL the shortest paths between the two gene in the functional network (length>1)
                              from=genea, 
                              to = geneb)
    df <- as.data.frame(t(sapply(sp$res, as_ids)))
    genes=data.frame(c(levels(unlist(df))))
    if(ncol(genes)!=0 & nrow(genes)!=0){
      colnames(genes)="genes"
      cluster=merge(genes_network, genes, by="genes")
      if(nrow(cluster)>2){
        clusters_allPairs=c(clusters_allPairs, cluster)
      }
    }
  }
}
clusters_allPairs

###############################################
# STEP 4: enrichment analysis of each cluster #
###############################################
#Load GSEA pathways
pathway=fread(args[7])

pathway=pathway[,-2]
pathway=unique(pathway)
pathway=data.frame(t(pathway))

pathway_analysis=data.frame()
for (i in 1:length(clusters_allPairs)){ #For each cluster
  for(j in 1:ncol(pathway)){ # and each pathway
    genesInCluster=data.frame(clusters_allPairs[[i]])
    colnames(genesInCluster)=c("gene")
    genesInPath=data.frame(pathway[complete.cases(pathway[,j]), j])
    colnames(genesInPath)="gene"
    path=genesInPath[1,1]
    genesInPath=merge(genesInPath, geneWeight, by="gene")
    genesInPath=data.frame(genesInPath[,1])
    colnames(genesInPath)="gene"
    genesInBoth=merge(genesInCluster, genesInPath, by="gene")
    if(nrow(genesInBoth)>2){
      #not weighted
      contingency=matrix(nrow=2, ncol=2) #2 × 2 contingency table for the genes in the cluster or not, and in the pathway or not
      contingency[1,1]=nrow(genesInBoth)
      contingency[1,2]=nrow(genesInPath)-nrow(genesInBoth)
      contingency[2,1]=nrow(genesInCluster)-nrow(genesInBoth)
      contingency[2,2]=nbGenes-contingency[1,1]-contingency[1,2]-contingency[2,1] #the background set contains all genes contained in the PPI network, co-expression network or GSEA pathways
      chi2=chisq.test(contingency,correct = TRUE) #Based on this contingency table, we computed a corrected chi-square statistic
      fish=fisher.test(contingency) 
      hyperGeom=phyper(nrow(genesInBoth)-1, nrow(genesInCluster), nbGenes-nrow(genesInCluster), 
                           nrow(genesInPath),lower.tail = FALSE)
      #weighted
      genesInCluster=data.frame(clusters_allPairs[[i]])
      colnames(genesInCluster)="gene"
      genesInCluster=merge(genesInCluster, geneWeight, by="gene")
      genesInPath=data.frame(pathway[complete.cases(pathway[,j]), j])
      colnames(genesInPath)="gene"
      path=genesInPath[1,1]
      genesInPath=merge(genesInPath, geneWeight, by="gene")
      genesInBoth=merge(genesInCluster, genesInPath, by="gene")
      genesInBoth=genesInBoth[,1:2]
      colnames(genesInBoth)=c("gene", "weight")
      
      contingency=matrix(nrow=2, ncol=2) #2 × 2 contingency table for the genes in the cluster or not, and in the pathway or not
      contingency[1,1]=sum(genesInBoth$weight)
      contingency[1,2]=sum(genesInPath$weight)-sum(genesInBoth$weight)
      contingency[2,1]=sum(genesInCluster$weight)-sum(genesInBoth$weight)
      contingency[2,2]=nbGenes-contingency[1,1]-contingency[1,2]-contingency[2,1] #the background set contains all genes contained in the PPI network, co-expression network or GSEA pathways
      chi2_weighted=chisq.test(contingency,correct = TRUE) #Based on this contingency table, we computed a corrected chi-square statistic
      fish_weighted=fisher.test(contingency)
      hyperGeom_Weighted=phyper(sum(genesInBoth$weight)-1, sum(genesInCluster$weight), nbGenes-sum(genesInCluster$weight), 
                       sum(genesInPath$weight),lower.tail = FALSE)
      output=data.frame(paste(unique(c(as.character(genesInBoth$gene))),collapse=" "), path, chi2$statistic,chi2$p.value, fish$p.value, hyperGeom, chi2_weighted$statistic, chi2_weighted$p.value, fish_weighted$p.value, hyperGeom_Weighted) #and the corresponding P-value
      pathway_analysis=rbind(pathway_analysis, output)
    }
  }
}



#These raw P-values were then corrected by Bonferroni correction, based on the total number of GSEA pathway terms and the number
#of clusters identified based on the respective SNP-gene association method.
pathway_analysis[,3:9]=pathway_analysis[,3:9]*(ncol(pathway)+length(clusters_allPairs))
unique_pathway=data.frame(unique(pathway_analysis$path))

#fwrite(pathway_analysis, paste(args[1],"/pathway/pathwayEnrichment.txt", sep=""), col.names=T)
#fwrite(unique_pathway, paste(args[1],"/pathway/unique_pathway.txt", sep=""), col.names=T)
