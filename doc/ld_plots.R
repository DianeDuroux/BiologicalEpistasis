#!/usr/bin/env Rscript

library(tidyverse)
library(cowplot)

rslt <- '../results/'

bg <- read_tsv(paste0(rslt, 'plink.ld'), col_types = 'iicdiicddd') %>%
        mutate(uniq_snp_id = cbind(SNP_A, SNP_B) %>% apply(1, sort) %>% apply(2, paste, collapse = '_'),
               distance = abs(BP_A - BP_B),
               maf_diff = abs(MAF_A - MAF_B)) %>%
        select(uniq_snp_id, R2, DP, distance, maf_diff) %>%
        mutate(where = 'All SNPs')

plot_significant <- function(folder, tag) {
    
    pairs <- read_csv(paste0(rslt, 'withBiofilter_filter/', folder, '/sign_SNPpairs.txt'), col_types = 'ccdccc') %>%
        mutate(uniq_snp_id = cbind(SNP_1, SNP_2) %>% apply(1, sort) %>% apply(2, paste, collapse = '_')) %>%
        select(uniq_snp_id, pvalue)
    read_tsv(paste0(rslt, 'withBiofilter_filter/', folder, '/plink.ld'), col_types = 'iicdiicddd') %>%
        mutate(uniq_snp_id = cbind(SNP_A, SNP_B) %>% apply(1, sort) %>% apply(2, paste, collapse = '_'),
               distance = abs(BP_A - BP_B),
               maf_diff = abs(MAF_A - MAF_B)) %>%
        select(uniq_snp_id, R2, DP, distance, maf_diff) %>%
        inner_join(pairs, by = "uniq_snp_id") %>%
        mutate(where = tag)
    
}

snp_info <- bind_rows(bg,
                      plot_significant('../noFilter', 'Standard'),
                      plot_significant('eqtl', 'eQTL'),
                      plot_significant('chromatin', 'Chromatin'),
                      plot_significant('eqtl_chrom', 'eQTL +\nChromatin'),
                      plot_significant('eqtl_chrom_phys', 'eQTL +\nChromatin +\nPhysical'),) %>%
    mutate(where = factor(where, levels = c('All SNPs', 'Standard', 'eQTL', 'Chromatin', 
                                            'eQTL +\nChromatin', 'eQTL +\nChromatin +\nPhysical')))

snp_dist <- ggplot(snp_info, aes(x = where, y = distance / 1000)) +
        geom_boxplot(outlier.shape = NA) +
        labs(x = 'Protocol', y = 'Distance (kb)') +
        theme_bw() +
        theme(text = element_text(size = 20))

snp_maf <- ggplot(snp_info, aes(x = where, y = maf_diff)) +
        geom_boxplot(outlier.shape = NA) +
        labs(x = 'Protocol', y = 'MAF difference') +
        theme_bw() +
        theme(text = element_text(size = 20))

snp_r2 <- ggplot(snp_info, aes(x = where, y = R2)) +
        geom_boxplot(outlier.shape = NA) +
        labs(x = 'Protocol', y = expression(R^2)) +
        theme_bw() +
        theme(text = element_text(size = 20))

snp_dp <- ggplot(snp_info, aes(x = where, y = DP)) +
        geom_boxplot(outlier.shape = NA) +
        labs(x = 'Protocol', y = expression(DP)) +
        theme_bw() +
        theme(text = element_text(size = 20))

plt <- plot_grid(snp_r2, snp_dp, snp_dist, snp_maf)
save(plt, file = 'ld_plots.RData')
ggsave('ld_plots.png', plt, width=19, height=12, bg = "transparent")
