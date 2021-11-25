To run this protocol on your dataset you need to first adapt the `run.sh` file. 

# Pre-requisites

Before launching `run.sh`, you need to prepare the dataset and the permutations:

1. Perform the quality control and other biological filters to the dataset (e.g. remove SNPs not mapped to any Biofilter interaction).
2. Run the epistasis detection using PLINK's `--epistasis` for the observed phenotype and for 400 permuted phenotypes.
3. Compute the experimental threshold of significance from these permutations.
4. Compute 999 permutations of the phenotypes for the ATPM procedure.

# Launching `run.sh`

First adapt the header by inserting paths to your specific files, then execute the script using the `sbatch` command. The main results will populate the folders `pvalues` (`sign_SNPpairs.txt` and `sign_GenePairs.txt`) and the `pathway` (`pathway_analysis_thresholdHyper.txt` and `Rplots.pdf`).

Specifically, the following variables nee to be set up in `run.sh`'s header:

- data: bed/bim/fam files
- mapFile: tsv containing the Biofilter models. It consists of four columns: Gene1, Gene2, SNP1 and SNP2.
- snpToGene: tsv containing the SNP-gene mapping. It consists of two columns: SNP and gene.
- genes: tsv containing the gene information. It consists of five columns: ensg, symbol, entrezID, HUGO and type.
- biofModels: tsv containing the Biofilter models. It consists of six columns: gene_2, gene_1, num_sources, num_instances, ensembl_1 and ensembl_2.
- pathway: file in [GMT format](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29) containing the pathways.
