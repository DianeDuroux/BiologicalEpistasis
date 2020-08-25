library(data.table)
library(stringr)
library(plyr)
library(igraph)
library(dplyr)

args=(commandArgs(TRUE))
print(args)

setwd(paste(args[1], "/pathway/", sep=""))

###############
# Set options #
###############
K=c(0.001, 0.01, 0.05, 0.1, 0.2)
B=as.integer(args[3])
exp_threshold=as.numeric(args[8])
exp_threshold

################################
# Export significant SNP pairs  #
################################
#pval<0.05 with Bonferroni correction and mappable to gene-pairs with biofilter
output_true=fread(paste(args[1],"/pvalues/output1.txt", sep=""))
output_true=output_true[output_true$pvalue<0.05,]
output_true$pvalue_uncorrected_MT=(output_true$pvalue*exp_threshold)/0.05
output_true=data.frame(output_true$SNP_1, output_true$SNP_2, output_true$pvalue_uncorrected_MT, output_true$pvalue, output_true$gene1, output_true$gene2)
colnames(output_true)=c("SNP_1", "SNP_2", "pvalue_uncorrected_MT", "pvalue_corrected_MT", "gene1", "gene2")
head(output_true)
fwrite(output_true, paste(args[1],"/pvalues/sign_SNPpairs.txt", sep=""))

############################
# Import gene pairs output #
############################

truePhenoOutput=data.frame(fread(paste(args[1],"/pvalues/Pvalues1.txt", sep=""))) #Load gene-pairs pvalues associated to the true phenotype
Pvalues=fread(paste(args[1],"/pvalues/Pvalues.txt", sep="")) #Load all gene-pairs pvalues

#Select only gene pairs when the pairs are significant for the true phenotype
signPairs=truePhenoOutput[truePhenoOutput$MinP<0.05,] #significant gene pairs
fwrite(signPairs, paste(args[1],"/pvalues/sign_GenePairs.txt", sep=""))

####################
# Pathway analysis #
####################
path=args[1]

################################################################################################
# STEP 1: Create gene sets/clusters from significant gene pairs detected and biofilter network #
################################################################################################

#Load significant gene-pairs
gene_epi=fread(paste(path, "/pvalues/sign_GenePairs.txt", sep=""))
genePair=gene_epi[,1]
gene_epi=data.frame(str_split_fixed(genePair$genePairs_names, " ", 2))
gene_epi=t(apply(gene_epi, 1, sort)) #order genes for each pair

genes_network=data.frame(c(levels(unlist(data.frame(gene_epi))))) #create list of significant genes
colnames(genes_network)="genes"

#Load functional/combined network = biofilter network
biofilter=fread(args[6])
func_network=biofilter[,1:2]
func_network=data.frame(t(apply(func_network, 1, sort))) #order genes for each pair
colnames(func_network)=c("gene1", "gene2")

nodes_func_network=c(as.character(func_network$gene1), as.character(func_network$gene2)) #create list of biofilter genes
nodes_func_network=as.factor(unique(nodes_func_network))

#For each epistatic gene pair, create a cluster of genes that are densely connected with each other in terms of epistatic interactions
expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...)) #function to create all possible pairs

clusters_allPairs=list()
for(i in 1:nrow(gene_epi)){ #for each gene pair
  genea=gene_epi[i,1]
  geneb=gene_epi[i,2]
  if((genea %in% nodes_func_network)==T & (geneb %in% nodes_func_network)==T){ #if the two genes are present in biofilter network
    sub_network=func_network[!(func_network$gene1==genea & func_network$gene2==geneb) | (func_network$gene1==geneb & func_network$gene2==genea), ] #remove direct link
    sub_net_func_network <- graph_from_data_frame(d=sub_network, vertices=nodes_func_network)
    sp=get.all.shortest.paths(graph=sub_net_func_network, #find ALL the shortest paths between the two gene in the functional network (length>1)
                              from=genea, 
                              to = geneb)
    df <- as.data.frame(t(sapply(sp$res, as_ids)))
    genes=data.frame(c(levels(unlist(df))))
    if(ncol(genes)!=0 & nrow(genes)!=0){
      all_possible_pairs <- expand.grid.df(x = genes, y = genes) #create all possible pairs with 2 genes in the shortest path
      colnames(all_possible_pairs)=c("X1", "X2")
      epi_net=data.frame(gene_epi)
      genes=data.frame(unique(unlist(merge(gene_epi, all_possible_pairs, by=c("X1", "X2"))))) #keep pairs present in the epistasis network
      colnames(genes)="gene"
      cluster=data.frame(unique(c(as.character(genes$gene),genea, geneb))) #add the originaly investigated gene pair
      #cluster=merge(genes_network, genes, by="genes") #select genes that are also involved in a significant gene-level epistasis signal
      if(nrow(cluster)>2){
        clusters_allPairs=c(clusters_allPairs, cluster)
      }
    }
  }
}
clusters_allPairs


#########################################################################################################################################################
# STEP 2: Create gene universe: mappable gene with genes in biofilter. Reduce pathway to pathways containing these genes and derive number of pathways. #
#########################################################################################################################################################

#Import SNP to gene mapping using HUGO names
snpToGene=fread(args[4]) #import snp to gene mapping (dataset specific)
snpToGene=snpToGene[,1:2]
genes=fread(args[5]) #import genes ensg -> HUGO for name conversion
genes=genes[,c(1, 4)]
colnames(genes)=c("ensg", "gene")
snpToGene=merge(snpToGene, genes, by="ensg") #Create SNP -> HUGO mapping
snpToGene=snpToGene[,2:3]
snpToGene=unique(snpToGene)
colnames(snpToGene)=c("rsID", "gene")
snpToGene=na.omit(snpToGene)
mappableGenes=data.frame(unique(snpToGene$gene))
colnames(mappableGenes)="gene"

#Reduce biofilter genes to mappable genes
biofGenes=data.frame(nodes_func_network)
colnames(biofGenes)="gene"
biofMap=merge(biofGenes, mappableGenes, by="gene")
#nbGenes=nrow(data.frame(unique(biofMap$gene))) #6455 instead of 14550 with genes in pathways

#Reduce genes in pathways to genes in biofilter and mapplable
pathway=fread(args[7], na.strings=c("","NA")) #Load pathways
pathway=pathway[-1,-2] #Remove unnecessary columns
pathway=unique(pathway)
pathway$V1=make.unique(pathway$V1) #add extension to duplicated pathway names (avoid to remove dupplicates)
pathNames=pathway$V1
pathway=pathway[,-1] #remove pathway names

genes_tot=unlist(pathway)
genes_tot=data.frame(unique(genes_tot))
colnames(genes_tot)="gene"
genes_tot=merge(genes_tot,biofMap, by="gene" ) 
nbGenes=nrow(data.frame(unique(genes_tot$gene))) #6399 genes

#Reduce genes in pathways to genes in Biofilter and mappable
BiofPath=list()
for(i in 1:nrow(pathway)){
  tmp=na.omit(data.frame(t(pathway[i,])))
  colnames(tmp)="gene"
  BiofPath[[i]]=t(merge(tmp, biofMap, by="gene"))
}
df <- ldply(BiofPath, data.frame)
df=data.frame(pathNames, df)
write.table(df, "biofPathways.csv", col.names = F)

########################################################
# STEP 3: enrichment analysis of each gene set/cluster #
########################################################
#Load GSEA pathways
#pathway=fread("biofPathways.csv")
pathway=df
#pathway=pathway[,-1]
pathway=unique(pathway)
pathway=data.frame(t(pathway))

pathway_analysis=data.frame()
for (i in 1:length(clusters_allPairs)){ #For each cluster
  for(j in 1:ncol(pathway)){ # and for each pathway
    genesInCluster=data.frame(clusters_allPairs[[i]]) #Select the genes of the cluster
    colnames(genesInCluster)=c("gene")
    genesInPath=data.frame(pathway[complete.cases(pathway[,j]), j]) #Select the genes of the pathway
    path=genesInPath[1,1] #save pathway name wich corresponds to the first line of the table
    genesInPath=data.frame(genesInPath[-1,]) #remove the name of the pathway to only keep the gene names
    colnames(genesInPath)="gene"
    genesInBoth=merge(genesInCluster, genesInPath, by="gene") #Select the genes present in the cluster AND in the pathway
    if(nrow(genesInBoth)>=1){ #if the overlap contains at least 1 gene
      contingency=matrix(nrow=2, ncol=2) #2 × 2 contingency table for the genes in the cluster or not, and in the pathway or not
      contingency[1,1]=nrow(genesInBoth)
      contingency[1,2]=nrow(genesInPath)-nrow(genesInBoth)
      contingency[2,1]=nrow(genesInCluster)-nrow(genesInBoth)
      contingency[2,2]=nbGenes-contingency[1,1]-contingency[1,2]-contingency[2,1] #the background set contains all genes contained in the PPI network, co-expression network or GSEA pathways
      chi2=chisq.test(contingency,correct = TRUE) #Based on this contingency table, we computed a corrected chi-square statistic
      fish=fisher.test(contingency) 
      hyperGeom=phyper(nrow(genesInBoth)-1, nrow(genesInCluster), nbGenes-nrow(genesInCluster), 
                       nrow(genesInPath),lower.tail = FALSE)
      output=data.frame(paste(unique(c(as.character(genesInBoth$gene))),collapse=" "), path, chi2$statistic,chi2$p.value, fish$p.value, hyperGeom) #and the corresponding P-value
      pathway_analysis=rbind(pathway_analysis, output) #save results
    }
   }
}

#Correct the P-values with Bonferroni correction, based on the total number of pathways and the number of gene clusters identified
pathway_analysis[,4:6]=pathway_analysis[,4:6]*(ncol(pathway)+length(clusters_allPairs))
pathway_analysis_thresholdChi2=pathway_analysis[pathway_analysis$chi2.p.value<0.05,]
pathway_analysis_thresholdHyper=pathway_analysis[pathway_analysis$hyperGeom<0.05,]

fwrite(pathway_analysis_thresholdHyper, paste(args[1],"/pathway/pathwayEnrichment.txt", sep=""), col.names=T)




