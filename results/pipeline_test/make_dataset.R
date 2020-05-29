#!/usr/bin/env Rscript
library(tidyverse)

eligible_snps <- read_tsv('CD_UC_CON_QCed_rel1_without_relatives_maf0.05_hwe0.001_Liu2015_232SNPs_LD0.75_noFilter_continuous.bim', col_names = F) %>% filter(X1 != 6) %>% .$X2
pos <- read_csv('~/projects/BiologicalEpistasis/results/withBiofilter_filter/eqtl/sign_SNPpairs.txt')
pos <- unique(c(pos$SNP_1, pos$SNP_2))

eqtl <- read_tsv('eqtl_mapping.tsv')

other_snps <- filter(eqtl, (rsID %in% eligible_snps) & (!rsID %in% pos)) %>%
    .$rsID %>%
    unique %>%
    sample(5000)
snps = c(pos, other_snps) 

tibble(snp = snps) %>%
    write_tsv('synthetic_snps.txt')

eqtl <- read_tsv('eqtl_mapping.tsv') %>% filter(rsID %in% snps)
info <- read_tsv('genes.tsv') %>% select(ensg, HUGO)

inner_join(eqtl, info) %>% select(rsID, HUGO, tissue) %>% rename(snp = rsID, gene = HUGO) %>% write_tsv('eqtl_benchmark.tsv')

genes <- read_tsv('eqtl_benchmark.tsv')$gene
read_tsv('biofilter.gene.models') %>% filter(`#gene1` %in% genes & gene2 %in% genes)  %>% rename(`Official Symbol Interactor A` = `#gene1`, `Official Symbol Interactor B` = gene2) %>% write_tsv('models.tab2')
