To do before the pipeline:
* Perform quality controls and apply biological filters on the dataset.
* Run the epistasis detection (for example PLINK --epistasis) for the observed phenotype and for 400 permuted phenotypes
* From these results, compute the experimental threshold
* Create 999 new permuted phenotypes for the ATPM procedure.

The presented scripts will perform additional quality controls (hla, LD), compute gene-pairs pvalues (ATPM), pathway enrichment analysis and visualization of the results.
To run the pipeline, adapt the header of the file run.sh first, by inserting paths to your specific files. Then, only execute the run.sh script (via the command sbatch). It will start the complete analysis and automaticaly call the other scripts.
Main results will be created in the folder "pvalues" (sign_SNPpairs.txt and sign_GenePairs.txt) and the folder "pathway" (pathway_analysis_thresholdHyper.txt and Rplots.pdf)