**BiologicalEpistasis**

Code and results from "Duroux, D., Climente-Gonzáles, H., Azencott, C. A., & Van Steen, K. (2020). Interpretable network-guided epistasis detection. bioRxiv.", to detect epistatic interactions at the gene level along the edges of a gene-gene co-function network.

**Usage**

The code necessary to replicate the analysis is in the folder "pipeline".

Steps that are included in the pipeline:

\* Quality controls: hla and LD.

\* Gene-pairs p-values computation via the ATPM procedure.

\* Pathway enrichment analysis.

\* Visualization of the results.

The main script of the pipeline is “run.sh”. To apply the analysis on a new dataset, this is the only file to modify and to run.  The header of the file run.sh needs to be adapted by inserting paths to your specific files. To start the analysis, execute the run.sh script (via the command sbatch). 

The main results will be uploaded in the folder "pvalues" (sign\_SNPpairs.txt and sign\_GenePairs.txt) and the folder "pathway" (pathway\_analysis\_thresholdHyper.txt and Rplots.pdf)

**Input files format (in the header of the run.sh file)**

data: bed/bim/fam (https://www.cog-genomics.org/plink/1.9/formats)

mapFile: Gene1,Gene2,SNP1,SNP2

snpToGene: SNP	gene (tab-separator)

genes: ensg	symbol	entrezID	HUGO	type (tab-separator)

biofModels: gene\_2,gene\_1,num\_sources,num\_instances,ensembl\_1,ensembl\_2 (from biofilter)

pathway: "pathway" "gene1" "gene2"... (from the R package msigdbr). 

**Citation**

If you use this code in your projects please cite our work:

@article{duroux2020interpretable,

`  `title={Interpretable network-guided epistasis detection},

`  `author={Duroux, Diane and Climente-Gonz{\'a}les, H{\'e}ctor and Azencott, Chlo{\'e}-Agathe and Van Steen, Kristel},

`  `journal={bioRxiv},

`  `year={2020},

`  `publisher={Cold Spring Harbor Laboratory}

}
