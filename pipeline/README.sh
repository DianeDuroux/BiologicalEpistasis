Folders 1_plink (epistasis detection), 2_experimental_threshold, 3_true_phenotypes (analysis) and 4_permutations (type I error) are the ordered steps of the pipeline.

  * 1_plink:

random_pheno.sh: to compute 1400 random phenotypes in R
plink_truePhenotypes.sh: to detect epistasis with the true dataset (to be performed for each process: eqtl with biofilter, eqtl without biofilter, eqtl+chrom with biofilter, eqtl+chrom without biofilter, physical with biofilter...)
plink_permutedPhenotypes.sh: to detect epistasis with the 1400 permuted datasets (to be performed for each process)

  * 2_experimental_threshold

run.sh: pipeline to select snp-pairs (remove hla, LD...) from 400 epistasis detection with permuted phenotypes. Adapt the header of this run file according to the process tested. (to be performed for each process)
threshold.sh (in folder "analysis"): to find the experimental threshold. (to be performed for each process)

  * 3_true_phenotypes

One pipeline when using Biofilter models, one pipeline without using it.
run.sh: to run entire pipeline (SNP pairs selection, SNP to gene conversion with ATPM, pathway analysis, visualisation). Adapt the header of the run file according to the analysis. (to be performed for each process)

  * 4_permutations

run.sh: to run the pipeline from 1000 epistasis detections with permuted phenotypes. Adapt the header of the run AND the part2_perm.sh files according to the analysis (same header in the two files). (to be performed for each process)
Analysis folder and scripts to summarize the outputs


In the header of a run file:
- dir= current directory
- softwares= location of mbmdr and plink softwares
- data=bed, bim, fam dataset
- dir_epistasis_output=location of the epistasis detection runs plink.epi.qt
- pathway= location of file kegg_go_biocarta_canonical_header.txt
- B= number of permutations in the ATPM procedure
- mapFile= mapping of pairs for example mapping_eqtl.txt
	gene1,gene2,SNP_1,SNP_2
	ABO,A3GALT2,rs505922,rs10429857
	...

pheno_ATPM=location of the permuted phenotypes for ATPM
threshold= threshold obtained in the log file of folder 2_experimental_threshold/analysis
snpToGene=mapping of snps, for example eqtl_mapping.tsv
	rsID    ensg    tissue
	rs1240708       ENSG00000235098 xQTLServer_eQTLs
	rs1240708       ENSG00000215915 xQTLServer_eQTLs


genes=location of genes.tsv
	ensg    symbol  entrezID        HUGO    type
	ENSG00000237094 RP4-669L17.10   NA      NA      lincRNA
	ENSG00000230021 RP5-857K21.4    101928626       LOC101928626    lincRNA

biofModels=location of biofilter_models_ensembl.tsv
	gene_2,gene_1,num_sources,num_instances,ensembl_1,ensembl_2
	A1CF,APOBEC2,2,3,ENSG00000124701,ENSG00000148584
	A1CF,APOBEC3B,2,3,ENSG00000179750,ENSG00000148584
