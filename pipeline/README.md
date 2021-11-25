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

- data: [bed](https://www.cog-genomics.org/plink2/formats#bed)/[bim](https://www.cog-genomics.org/plink2/formats#bim)/[fam](https://www.cog-genomics.org/plink2/formats#fam) files
- mapFile: [CSV](https://en.wikipedia.org/wiki/Comma-separated_values) file containing the Biofilter models. It consists of four columns: Gene1, Gene2, SNP1 and SNP2.
- snpToGene: [TSV](https://en.wikipedia.org/wiki/Tab-separated_values) file containing the SNP-gene mapping. It consists of two columns: SNP and gene.
- genes: [TSV](https://en.wikipedia.org/wiki/Tab-separated_values) containing the gene information. It consists of five columns: ensg, symbol, entrezID, HUGO and type.
- biofModels: [CSV](https://en.wikipedia.org/wiki/Comma-separated_values) containing the Biofilter models. It consists of six columns: gene_2, gene_1, num_sources, num_instances, ensembl_1 and ensembl_2.
- pathway: [GMT](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29) file containing the pathways.
