# BiologicalEpistasis

We propose a procedure to detect epistatic interactions at the gene level along the edges of a gene-gene co-function network:

> Duroux, D., Climente-González, H., Azencott, C.-A., & Van Steen, K. (2020). **Interpretable network-guided epistasis detection.** bioRxiv. https://doi.org/10.1101/2020.09.24.310136

This repository includes the code used to generate the results presented in this article. It consists of three subdirectories:

- `results/` contains the results presented in the article. The main script of the pipeline is [run.sh](pipeline/run.sh), and includes the following steps:

  1. Quality controls: hla and LD
  2. Gene-pairs p-values computation via the ATPM procedure
  3. Pathway enrichment analysis
  4. Visualization of the results

- `pipeline/` contains the code used to produce the aforementioned results.
- `doc/` contains the code used to produce the figures and analyses in the article.

## Apply it to your dataset

If you are interested in applying this protocol to your own dataset, we provide two options. The first one is the `network_epistasis.nf` script in [hclimente/gwas-tools](https://github.com/hclimente/gwas-tools#network_epistasis).

The second one is cloning this repository, and modifying and running the `run.sh` file. First adapt the header by inserting paths to your specific files, then execute the script using the `sbatch` command. The main results will populate the folders `pvalues` (`sign_SNPpairs.txt` and `sign_GenePairs.txt`) and the `pathway` (`pathway_analysis_thresholdHyper.txt` and `Rplots.pdf`).

Specifically, the following variables nee to be set up in `run.sh`'s header:

- data: bed/bim/fam files
- mapFile: tsv containing the Biofilter models. It consists of four columns: Gene1, Gene2, SNP1 and SNP2.
- snpToGene: tsv containing the SNP-gene mapping. It consists of two columns: SNP and gene.
- genes: tsv containing the gene information. It consists of five columns: ensg, symbol, entrezID, HUGO and type.
- biofModels: tsv containing the Biofilter models. It consists of six columns: gene_2, gene_1, num_sources, num_instances, ensembl_1 and ensembl_2.
- pathway: file in [GMT format](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29) containing the pathways.

## Citation

If you use this code in your projects please cite our work:

```
@article{duroux2020interpretable,
  title={Interpretable network-guided epistasis detection},
  author={Duroux, Diane and Climente-Gonz{\'a}lez, H{\'e}ctor and Azencott, Chlo{\'e}-Agathe and Van Steen, Kristel},
  journal={bioRxiv},
  year={2020},
  publisher={Cold Spring Harbor Laboratory}
}
```
