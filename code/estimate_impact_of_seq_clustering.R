# If we compute V gene frequencies and clone sizes with and without sequence clustering, how strongly are they correlated?
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

source('gene_frequency_functions.R')

# Basic info for each clone (germline genes, CDR lenght, naive CDR seq)
clone_info <- read_csv('../processed_data/clone_info.csv') %>%
  dplyr::rename(clone_id = clone_id_partis, v_gene = v_segment_partis, j_gene = j_segment_partis,
                d_gene = d_segment_partis)

# Number of unique productive sequences in each clone, by tissue and cell type ** WITH sequence clustering **
seq_counts_clustered <- read_csv('../processed_data/unique_seq_counts.csv')

# ** WITHOUT sequence clustering **
seq_counts_unclustered <- read_csv('../processed_data/seq_counts_unclustered.csv') %>%
  # (code for computing gene freqs. requires count variable to be named uniq_prod_seqs)
  dplyr::rename(uniq_prod_seqs = prod_seqs)

seq_counts_clustered <- left_join(seq_counts_clustered, clone_info)
seq_counts_clustered <- get_info_from_mouse_id(seq_counts_clustered)

seq_counts_unclustered <- left_join(seq_counts_unclustered, clone_info)
seq_counts_unclustered <- get_info_from_mouse_id(seq_counts_unclustered)

# ======== Are gene frequencies with and without sequence clustering correlated? =========

# Get naive frequencies excluding the LN, tissue-specific frequencies for other cell types
naive_from_tissue <- c('spleen','BM')

naive_freqs_clustered <- (calc_gene_freqs(seq_counts_clustered, long_format = F, by_tissue = F, tissue_subset = naive_from_tissue))$naive_freqs
exp_freqs_clustered <- (calc_gene_freqs(seq_counts_clustered, long_format = F, by_tissue = T))$exp_freqs

naive_freqs_unclustered <- (calc_gene_freqs(seq_counts_unclustered, long_format = F, by_tissue = F, tissue_subset = naive_from_tissue))$naive_freqs
exp_freqs_unclustered <- (calc_gene_freqs(seq_counts_unclustered, long_format = F, by_tissue = T))$exp_freqs


naive_freqs <- full_join(naive_freqs_clustered, naive_freqs_unclustered, 
                         by = c('mouse_id','day','infection_status', 'group',
                                'group_controls_pooled','v_gene'),
                         suffix = c('_clustered','_unclustered')) 
  
exp_freqs <- full_join(exp_freqs_clustered, exp_freqs_unclustered, 
                  by = c('mouse_id','day','infection_status', 'group',
                         'group_controls_pooled','v_gene','tissue','cell_type'),
                  suffix = c('_clustered','_unclustered')) %>%
  mutate(tissue = factor(tissue, levels = c('LN','spleen','BM')),
         cell_type = factor(cell_type, levels = c('GC','PC','mem','experienced')))


v_gene_freq_correlation_naive <- naive_freqs %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled) %>%
  summarise(cor_coef = cor.test(naive_vgene_seq_freq_clustered, naive_vgene_seq_freq_unclustered, method = 'spearman')$estimate) %>%
  ungroup() %>%
  mutate(group_controls_pooled = factor(group_controls_pooled, levels = group_controls_pooled_factor_levels))
  
v_gene_freq_correlation_exp <- exp_freqs %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, cell_type, tissue) %>%
  mutate(n_valid_pairs = sum(!is.na(vgene_seq_freq_clustered) & !is.na(vgene_seq_freq_unclustered))) %>%
  filter(n_valid_pairs > 0) %>%
  summarise(cor_coef = cor.test(vgene_seq_freq_clustered, vgene_seq_freq_unclustered, method = 'spearman')$estimate) %>%
  mutate(group_controls_pooled = factor(group_controls_pooled, levels = group_controls_pooled_factor_levels))


v_gene_freq_correlation_naive %>%
  ggplot(aes(x = group_controls_pooled, y = cor_coef, color = infection_status)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point() +
  geom_hline(yintercept = 0.95, linetype = 2) +
  xlab('Group') +
  ylab('Correlation between naive V gene\nfrequencies with and without clustering') +
  background_grid() +
  theme(legend.position = 'top') +
  scale_color_discrete(name = 'Infection status')

v_gene_freq_correlation_exp %>%
  ggplot(aes(x = group_controls_pooled, y = cor_coef, color = infection_status)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point() +
  geom_hline(yintercept = 0.95, linetype = 2) +
  facet_grid(tissue~cell_type) +
  xlab('Group') +
  ylab('Correlation between V gene frequencies with and without clustering') +
  background_grid() +
  theme(legend.position = 'top') +
  scale_color_discrete(name = 'Infection status') +
  scale_x_discrete(labels = function(x){str_replace(x,'-','\n')})

# ======== Are clone frequencies with and without sequence clustering correlated? =========
clone_size_dist_clustered <- get_clone_size_distribution(seq_counts_clustered)
clone_size_dist_unclustered <- get_clone_size_distribution(seq_counts_unclustered)

clone_size_dist <- full_join(clone_size_dist_clustered, clone_size_dist_unclustered, 
          by = c('mouse_id','day','infection_status', 'group',
                 'group_controls_pooled','clone_id','tissue','cell_type'),
          suffix = c('_clustered','_unclustered'))

clone_size_correlation <- clone_size_dist %>%
  group_by(mouse_id, day, infection_status, group, group_controls_pooled, cell_type, tissue) %>%
  mutate(n_valid_pairs = sum(!is.na(clone_freq_clustered) & !is.na(clone_freq_unclustered))) %>%
  filter(n_valid_pairs >= 3) %>%
  summarise(cor_coef = cor.test(clone_freq_clustered, clone_freq_unclustered, method = 'spearman')$estimate) %>%
  mutate(tissue = factor(tissue, levels = c('LN','spleen','BM')),
         cell_type = factor(cell_type, levels = c('naive','GC','PC','mem','experienced')),
         group_controls_pooled = factor(group_controls_pooled, levels = group_controls_pooled_factor_levels))
  
clone_size_correlation %>%
  ggplot(aes(x = group_controls_pooled, y = cor_coef, color = infection_status)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point() +
  geom_hline(yintercept = 0.95, linetype = 2) +
  facet_grid(tissue~cell_type) +
  xlab('Group') +
  ylab('Correlation between clone frequencies with and without clustering') +
  background_grid() +
  theme(legend.position = 'top') +
  scale_color_discrete(name = 'Infection status') +
  scale_x_discrete(labels = function(x){str_replace(x,'-','\n')}) +
  ylim(c(0,1))

