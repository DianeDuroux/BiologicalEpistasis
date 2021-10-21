The results of the proposed epistasis detection pipeline are available in the [withBiofilter_filter](withBiofilter_filter) folder. Our pipeline uses gene-gene interactions obtained from the Biofilter database as candidates for epistasis detection. There are multiple ways of mapping genotypes to these gene-gene interactions. Each subfolder contains the results for one such way; for instance, [withBiofilter_filter/eqtl_chrom_phys](withBiofilter_filter/eqtl_chrom_phys) uses all mappings (eQTL, 3D structure of the chromatin and physical mapping). Within each of these subfolders, the following files can be found:

- pathway_analysis_thresholdHyper.txt: pathway enrichment analyses of the significant genes.
- sign_GenePairs.txt: significant gene pairs.
- sign_SNPpairs.txt: significant SNP pairs.
- sign_SNPpairs_add.txt: P-values of the SNP pairs (from sign_SNPpairs.txt), including the polygenic risk score in the regression to adjust for main effects. SNPs are encoded using an additive encoding.
- sign_SNPpairs_codom.txt: P-values of the SNP pairs (from sign_SNPpairs.txt), including the polygenic risk score in the regression to adjust for main effects. SNPs are encoded using a one-hot encoding.

We also include the results of experiments to characterize the pipeline, in folders which have a similar organization. In [withoutBiofilter_filter](withoutBiofilter_filter) we run the pipeline without considering Biofilter interactions (all SNP-SNP interactions are possible), and different SNP-gene mappings to produce gene-gene results.
