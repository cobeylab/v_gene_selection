library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggdendro)
library(cowplot)
library(scales)
library(viridis)
#library(ComplexHeatmap)
#library(circlize)

theme_set(theme_cowplot())

source('gene_frequency_functions.R')
source('plot_options.R')

# Load precomputed gene frequencies, neutral realizations, pairwise correlations 

load('../results/precomputed_gene_freqs.RData')
#load('~/Desktop/v_gene_selection_files/precomputed_gene_freqs.RData')

seq_counts <- bind_rows(exp_seq_counts, naive_seq_counts)


# Basic info for each clone (germline genes, CDR lenght, naive CDR seq)
clone_info <- read_csv('../processed_data/clone_info.csv')
# clone_info <- read_csv('~/Desktop/v_gene_selection_files/clone_info.csv')


figures_output_dir <- '../figures/exported_ggplot_objects/'
#figures_output_dir <- '~/Desktop/v_gene_selection_files/figures/exported_ggplot_objects/'

# ======= Example of gene rank frequency plot ======================================
LN_PC_freqs_primary8 <- gene_freqs %>% 
  filter(group_controls_pooled == 'primary-8', tissue == 'LN', cell_type == 'PC',
         total_compartment_seqs >= 100)

vgene_order <- LN_PC_freqs_primary8 %>%
  group_by(v_gene) %>%
  summarise(median_freq_across_mice = median(vgene_seq_freq),
            median_naive_freq_across_mice = median(naive_vgene_seq_freq)) %>%
  arrange(desc(median_freq_across_mice)) %>%
  pull(v_gene)

LN_PC_freqs_primary8  %>%
  mutate(v_gene = factor(v_gene, levels = vgene_order)) %>%
  ggplot(aes(x = v_gene)) +
  geom_point(aes(y = vgene_seq_freq)) +
  facet_wrap('mouse_id', nrow = 2) +
  scale_y_log10(breaks = c(1e-4,1e-3,1e-2,1e-1)) +
  background_grid(minor = 'none', major = 'y') +
  theme(axis.text.x = element_blank()) +
  xlab('V genes\n(ordered by median frequency across mice)') +
  ylab('Frequency in lymph node plasma cells')

# ===== CORRELATION BETWEEN NAIVE AND EXPERIENCED FREQUENCIES ==========

# -- Select scatterplots showing exp. vs. naive freq correlations ------

plot_naive_freq_corr <- function(group_controls_pooled, tissue, cell_type){
  stopifnot(tissue == 'LN') # If looking at other tissues have to change y axis
  
  cell_type_label <- switch(cell_type,
                            'GC' = 'germinal center cells',
                            'PC' = 'plasma cells',
                            'mem' = 'memory cells',
                            'nonnaive_IgD+B220+' = 'non-naive IgD+B220+ cells')
  
  gene_freqs %>% 
    filter(group_controls_pooled == !!group_controls_pooled, tissue == !!tissue, cell_type == !!cell_type) %>%
    #filter(total_compartment_seqs >= 100, total_mouse_naive_seqs >= 100) %>%
    ggplot(aes(x = naive_vgene_seq_freq, vgene_seq_freq)) +
    geom_point() +
    facet_wrap('mouse_id', nrow = 2) +
    scale_y_log10() +
    scale_x_log10() + 
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    background_grid() +
    xlab('Frequency in naive repertoire') +
    ylab(paste0('Frequency in lymph node ', cell_type_label))
}
plot_naive_freq_corr('primary-8', 'LN','PC')
plot_naive_freq_corr('primary-16', 'LN','GC')
plot_naive_freq_corr('secondary-40', 'LN','GC')

plot_naive_freq_corr('primary-8', 'LN','nonnaive_IgD+B220+')
plot_naive_freq_corr('primary-16', 'LN','nonnaive_IgD+B220+')
plot_naive_freq_corr('secondary-40', 'LN','nonnaive_IgD+B220+')


# Spearman correlation coefficients
naive_exp_correlations_obs <- get_naive_exp_correlations(gene_freqs) %>%
  filter(total_compartment_seqs > 100, total_mouse_naive_seqs > 100) %>%
  mutate(group_controls_pooled = factor(group_controls_pooled, levels = group_controls_pooled_factor_levels)) %>%
  mutate(day = as.integer(as.character(day)))

naive_exp_correlations_neutral <- lapply(neutral_realizations %>% group_by(replicate) %>% group_split(),
                                         FUN = get_naive_exp_correlations)

naive_exp_correlations_neutral <- bind_rows(naive_exp_correlations_neutral, .id = 'replicate') %>%
  mutate(replicate = as.numeric(replicate)) %>%
  filter(total_compartment_seqs > 100, total_mouse_naive_seqs > 100) %>%
  mutate(group_controls_pooled = factor(group_controls_pooled, levels = group_controls_pooled_factor_levels)) %>%
  mutate(day = as.integer(as.character(day)))

naive_exp_corr_weighted_means <- naive_exp_correlations_obs %>%
  group_by(group_controls_pooled, infection_status, cell_type, tissue) %>%
  mutate(weights = total_compartment_seqs/sum(total_compartment_seqs)) %>%
  summarise(naive_exp_corr_weighted_mean = sum(naive_exp_corr*weights)) %>%
  ungroup()

naive_exp_correlations_plot <- naive_exp_correlations_obs %>%
  filter(tissue == 'LN', group_controls_pooled != 'control') %>%
  filter(cell_type %in% c('GC','PC','mem')) %>%
  mutate(cell_type = case_when(
    cell_type == 'GC' ~ 'Lymph node GC cells',
    cell_type == 'PC' ~ 'Lymph node plasma cells',
    cell_type == 'mem' ~ 'Lymph node memory cells'
  )) %>%
  mutate(cell_type = factor(cell_type, levels = c('Lymph node GC cells',
                                                  'Lymph node plasma cells',
                                                  'Lymph node memory cells'))) %>%
  ggplot(aes(x = day, y = naive_exp_corr, color = infection_status, group = day)) +
  geom_boxplot(data = naive_exp_correlations_neutral %>%
                 filter(tissue == 'LN', group_controls_pooled != 'control',
                        cell_type %in% c('GC','PC','mem')) %>%
                 mutate(cell_type = case_when(
                   cell_type == 'GC' ~ 'Lymph node GC cells',
                   cell_type == 'PC' ~ 'Lymph node plasma cells',
                   cell_type == 'mem' ~ 'Lymph node memory cells'
                 )) %>%
                 mutate(cell_type = factor(cell_type, levels = c('Lymph node GC cells',
                                                                 'Lymph node plasma cells',
                                                                 'Lymph node memory cells'))),
               outlier.alpha = 0, linetype = 2) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(aes(size = total_compartment_seqs),
             position = position_jitter(width = 0.1),
             alpha = point_alpha) +
  facet_grid(.~cell_type, scales = 'free') +
  theme(legend.position = 'top') +
  xlab('Days after primary infection') +
  ylab('Spearman correlation in V gene frequencies between\nnaive repertoire and influenza induced populations') +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_color_manual(values = c('green3','dodgerblue2')) +
  scale_size_continuous(name = 'Number of sequences',
                        breaks = c(1000,10000,20000)) +
  guides(color = 'none') +
  background_grid() +
  scale_x_continuous(breaks = unique(naive_exp_correlations_obs$day))

plot(naive_exp_correlations_plot)

naive_exp_correlations_obs %>%
  filter(tissue == 'spleen') %>%
  ggplot(aes(x = group_controls_pooled, y = naive_exp_corr, color = infection_status)) +
  geom_boxplot(data = naive_exp_correlations_neutral %>% filter(tissue == 'spleen'),
               outlier.alpha = 0) +
  geom_point(aes(size = total_compartment_seqs), shape = 1,
             position = position_jitter(width = 0.1)) +
  geom_point(data = naive_exp_corr_weighted_means %>%
               filter(tissue == 'spleen'),
             aes(y = naive_exp_corr_weighted_mean),
             shape = 4, size = 4, stroke = 2, show.legend = F) +
  facet_grid(.~cell_type, scales = 'free') +
  theme(axis.text.x = element_text(angle = 20, vjust = 0.5), legend.position = 'top') +
  xlab('Group') +
  ylab('Correlation in gene frequencies\nin naive repertoire vs. spleen cells') +
  geom_hline(yintercept = 0, linetype = 2) +
  #scale_color_manual(values = c('green3','dodgerblue2')) +
  scale_size_continuous(name = 'Number of unique sequences',
                        breaks = c(1000,10000,20000)) +
  guides(color = 'none') +
  background_grid()


# Rho: ratio of experienced to naive frequency
obs_rhos <- gene_freqs_adj_naive_zeros %>%
  mutate(log_rho = log(vgene_seq_freq) - log(naive_vgene_seq_freq)) %>%
  mutate(obs_rho = exp(log_rho)) %>%
  rename(obs_n_vgene_seqs = n_vgene_seqs) %>%
  select(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type, v_gene,
         obs_n_vgene_seqs,obs_rho)

neutral_rhos <- neutral_realizations %>%
  mutate(log_rho = log(vgene_seq_freq) - log(naive_vgene_seq_freq)) %>%
  mutate(rho = exp(log_rho)) %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type, v_gene) %>%
  summarise(mean_sim_rho = mean(rho),
            lbound_sim_rho = quantile(rho, 0.025),
            ubound_sim_rho = quantile(rho, 0.975)) %>%
  ungroup()

deviation_from_naive <- left_join(obs_rhos, neutral_rhos) %>%
  mutate(deviation_from_naive = case_when(
    obs_rho > ubound_sim_rho ~ 'positive',
    obs_rho < lbound_sim_rho ~ 'negative',
    (obs_rho >= lbound_sim_rho) & (obs_rho <= ubound_sim_rho) ~ 'neutral'
  )) %>%
  mutate(group_controls_pooled = factor(group_controls_pooled, levels = group_controls_pooled_factor_levels)) %>%
  # deviation from naive is NA for mice lacking naive sequences in the tissues in naive_from_tissue
  filter(!is.na(deviation_from_naive))

n_genes_by_deviation <- deviation_from_naive %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type) %>%
  mutate(total_genes = length(unique(v_gene))) %>%
  ungroup() %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type, total_genes, deviation_from_naive) %>%
  summarise(n_genes = length(unique(v_gene))) %>%
  ungroup() %>%
  mutate(fraction_genes = n_genes / total_genes)

n_genes_by_deviation %>% filter(group_controls_pooled != 'control', tissue == 'LN') %>%
  #filter(!is.na(deviation_from_naive)) %>%
  ggplot(aes(x = group_controls_pooled, y = fraction_genes, color = infection_status)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point() +
  facet_grid(deviation_from_naive~cell_type) +
  background_grid() + 
  scale_color_manual(values = c('green3','dodgerblue2')) +
  theme(axis.text.x = element_text(angle = 20, vjust = 0.5), legend.position = 'top') +
  xlab('Group') + ylab('Fraction of genes')


# Genes deviating somewhat consistently across mice 
genes_by_n_mice_called <- 
  left_join(deviation_from_naive %>% filter(!is.na(deviation_from_naive)), 
            gene_freqs %>% select(mouse_id, cell_type, tissue, total_compartment_seqs, total_mouse_naive_seqs) %>%
              unique()) %>%
  filter(total_compartment_seqs >= 100, total_mouse_naive_seqs >= 100) %>%
  group_by(group_controls_pooled, v_gene) %>% 
  mutate(n_mice_gene_occurs_in_this_group = length(unique(mouse_id))) %>%
  group_by(group_controls_pooled, v_gene, n_mice_gene_occurs_in_this_group, cell_type, tissue, deviation_from_naive) %>%
  summarise(n_mice_with_call = length(unique(mouse_id))) %>% ungroup()

get_consistent_genes_table <- function(candidate_genes_tissue, candidate_genes_cell_type, candidate_genes_group){
  candidate_genes <- genes_by_n_mice_called %>% 
    filter(tissue == candidate_genes_tissue, group_controls_pooled == candidate_genes_group, cell_type == candidate_genes_cell_type,
           deviation_from_naive == 'positive') %>%
    arrange(desc(n_mice_with_call)) %>%
    filter(n_mice_with_call >= 2) %>% pull(v_gene)
  
  output_table <- genes_by_n_mice_called %>% 
    filter(tissue == candidate_genes_tissue , group_controls_pooled == candidate_genes_group, cell_type == candidate_genes_cell_type) %>%
    filter(v_gene %in% candidate_genes) %>% pivot_wider(names_from = deviation_from_naive, values_from = n_mice_with_call) %>%
    arrange(desc(positive)) %>% select(v_gene, positive, matches('neutral'), negative) 
  
  if(('neutral' %in% colnames(output_table)) == F){
    output_table <- output_table %>% mutate(neutral = 0)
  }
  
  output_table <- output_table %>%
    mutate(positive = ifelse(is.na(positive), 0, positive),
           negative = ifelse(is.na(negative), 0, negative),
           neutral = ifelse(is.na(neutral), 0, neutral)
    ) %>%
    rename(nonsignificant = neutral)
  return(output_table)
}

get_consistent_genes_table('LN','PC', 'primary-8') %>%
  filter(positive >= 3)
get_consistent_genes_table('LN','PC', 'primary-16') %>%
  filter(positive >= 3)
get_consistent_genes_table('LN','PC', 'primary-24') %>%
  filter(positive >= 2)

get_consistent_genes_table('LN','PC', 'secondary-56') %>%
  filter(positive >= 3)


get_consistent_genes_table('LN','GC', 'primary-8') %>%
  filter(positive >= 2)
get_consistent_genes_table('LN','GC', 'primary-16') %>%
  filter(positive >= 3)
get_consistent_genes_table('LN','GC', 'primary-24') %>%
  filter(positive >= 2)

get_consistent_genes_table('LN','nonnaive_IgD+B220+', 'primary-8') 
get_consistent_genes_table('LN','nonnaive_IgD+B220+', 'control') 

# Plots showing deviation from naive freq for major genes
gene_freqs <- left_join(gene_freqs, 
                        deviation_from_naive %>%
                          select(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type, v_gene, matches('rho'), deviation_from_naive))

gene_freqs <- gene_freqs %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type) %>%
  mutate(v_gene_rank = rank(-vgene_seq_freq, ties.method = 'first')) %>%
  ungroup()

plot_most_common_genes <- function(plot_cell_type, plot_tissue, plot_group){
  gene_freqs %>% 
    filter(cell_type == plot_cell_type, tissue == plot_tissue, group_controls_pooled == plot_group, v_gene_rank <= 20) %>%
    filter(total_compartment_seqs >= 100, total_mouse_naive_seqs >= 100) %>%
    rowwise() %>%
    mutate(label_position = ifelse(vgene_seq_freq > naive_vgene_seq_freq, 1.05*vgene_seq_freq, 1.05*naive_vgene_seq_freq)) %>%
    ggplot() +
    geom_point(aes(x = v_gene_rank, y = naive_vgene_seq_freq, color = deviation_from_naive)) +
    geom_segment(aes(x = v_gene_rank, xend = v_gene_rank, y = naive_vgene_seq_freq, yend = vgene_seq_freq,
                     color = deviation_from_naive),
                 arrow = arrow(ends = 'last', length = unit(10, "pt"), type = 'closed')) +
    geom_text(aes(label = str_remove(v_gene, 'IGHV'),
                  x = v_gene_rank, y = label_position), angle = 20, size = 3, alpha = 0.8) +
    facet_wrap('mouse_id', scales = 'free') +
    xlab(paste0('V gene rank in ', plot_group, ' ', plot_tissue, ' ',
                switch(plot_cell_type, 'PC' = 'plasma cells', 'GC' = 'germinal center cells', 'mem' = 'memory cells',
                       'nonnaive_IgD+B220+' = 'Non-naive IgD+B220+ cells'), ' (top 20 genes only)')) +
    ylab('V gene frequency') +
    theme(legend.position = c(0.4,0.25)) + 
    scale_x_continuous(expand = c(0.15,0)) +
    scale_color_discrete(name = 'Deviation from naive\nfrequency (bootstrap)', labels = c('negative','non-significant','positive'))
  
}

top_genes_LN_PC_day8_plot <- plot_most_common_genes('PC','LN','primary-8') +
  theme(legend.position = c(0.87,0.35)) + background_grid()
plot(top_genes_LN_PC_day8_plot)

top_genes_LN_PC_day16_plot <- plot_most_common_genes('PC','LN','primary-16') +
  theme(legend.position = c(0.67,0.25))
plot(top_genes_LN_PC_day16_plot)

plot_most_common_genes('PC','LN','primary-24') +
  theme(legend.position = c(0.85,0.30))

plot_most_common_genes('PC','LN','secondary-56') +
  theme(legend.position = c(0.87,0.35))

top_genes_LN_GC_day8_plot <- plot_most_common_genes('GC','LN','primary-8') +
  theme(legend.position = c(0.8,0.3))

plot(top_genes_LN_GC_day8_plot)

top_genes_LN_GC_day16_plot <- plot_most_common_genes('GC','LN','primary-16') +
  theme(legend.position = c(0.33,0.2))
plot(top_genes_LN_GC_day16_plot)

plot_most_common_genes('GC','LN','primary-24') +
  theme(legend.position = c(0.85,0.35))

plot_most_common_genes('nonnaive_IgD+B220+','LN','control') +
  theme(legend.position = c(0.85,0.3))

plot_most_common_genes('nonnaive_IgD+B220+','LN','primary-8') +
  theme(legend.position = c(0.8,0.2))
plot_most_common_genes('nonnaive_IgD+B220+','LN','primary-16') +
  theme(legend.position = c(0.85,0.35))
plot_most_common_genes('nonnaive_IgD+B220+', 'LN','primary-24') +
  theme(legend.position = c(0.85,0.35))

plot_most_common_genes('nonnaive_IgD+B220+','spleen','primary-8') +
  theme(legend.position = c(0.8,0.2))

plot_most_common_genes('nonnaive_IgD+B220+','spleen','primary-8') +
  theme(legend.position = c(0.8,0.2))

plot_most_common_genes('mem','LN','primary-24') +
  theme(legend.position = c(0.85,0.30))

# ============== CLONE SIZE DISTRIBUTIONS


# Add general clone information to object with clone size distributions
clone_freqs_by_tissue_and_cell_type <- left_join(clone_freqs_by_tissue_and_cell_type, clone_info)


# Add information on whether the V gene used by each clone has evidence of gene-level selection
clone_freqs_by_tissue_and_cell_type <- left_join(clone_freqs_by_tissue_and_cell_type,
                                                 deviation_from_naive %>%
                                                   select(mouse_id, day, infection_status, group_controls_pooled,
                                                          tissue, cell_type, v_gene, matches('rho'), deviation_from_naive) %>%
                                                   dplyr::rename(compartment_tissue = tissue, compartment_cell_type = cell_type))

clone_freqs_by_tissue_and_cell_type <- left_join(clone_freqs_by_tissue_and_cell_type,
                                                 gene_freqs %>% select(mouse_id, total_mouse_naive_seqs) %>% unique()) %>%
  mutate(compartment_tissue = factor(compartment_tissue, levels = c('LN','spleen','BM')),
         compartment_cell_type = factor(compartment_cell_type, levels = c('naive', 'nonnaive_IgD+B220+','GC','PC','mem')))

# For each clone, include a list of mutations above a certain frequency threshold
# (For now, frequency is calculated relative to the number of productive sequences from a clone in a particular cell type and tissue combination)

mutations_above_threshold <- list_clone_mutations_above_threshold(mutation_freqs_within_clones_by_tissue_and_cell_type,
                                                                  threshold = 0.5)

clone_freqs_by_tissue_and_cell_type <- left_join(clone_freqs_by_tissue_and_cell_type,
                                                 mutations_above_threshold %>% select(mouse_id, clone_id, tissue, cell_type,
                                                                                      mutations_above_threshold) %>%
                                                   dplyr::rename(compartment_tissue = tissue, compartment_cell_type = cell_type)) %>%
  mutate(mutations_above_threshold = ifelse(is.na(mutations_above_threshold), '', mutations_above_threshold))


# Fraction of sequences in the largest 10 clones
fraction_in_top_10_clones_plot <- clone_freqs_by_tissue_and_cell_type %>% 
  filter(compartment_tissue == 'LN', group_controls_pooled != 'control') %>%
  filter(clone_rank_in_compartment <=10) %>% 
  filter(compartment_cell_type %in% c('GC','PC','mem')) %>%
  mutate(compartment_cell_type = factor(compartment_cell_type, levels = c('GC','PC','mem'))) %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, compartment_cell_type, compartment_tissue, total_seqs_in_compartment) %>%
  summarise(seqs_in_top_clones = sum(n_clone_seqs_in_compartment)) %>%
  mutate(fraction_seqs_in_top_clones = seqs_in_top_clones / total_seqs_in_compartment) %>%
  ungroup() %>%
  mutate(day = as.integer(as.character(day))) %>%
  mutate(compartment_cell_type = case_when(
    compartment_cell_type == 'GC' ~ 'Lymph node GC cells',
    compartment_cell_type == 'PC' ~ 'Lymph node plasma cells',
    compartment_cell_type == 'mem' ~ 'Lymph node memory cells'
  )) %>%
  mutate(compartment_cell_type = factor(compartment_cell_type,
                                        levels = c('Lymph node GC cells',
                                                  'Lymph node plasma cells',
                                                  'Lymph node memory cells'))) %>%
  ggplot(aes(x = day, y = fraction_seqs_in_top_clones, color = infection_status, group = day)) +
  geom_boxplot(outlier.alpha =  F, show.legend = F) +
  geom_point(aes(size = total_seqs_in_compartment), alpha = 0.8) +
  facet_grid(.~compartment_cell_type) +
  background_grid() +
  scale_color_manual(values = c('green3','dodgerblue2'), name = 'Infection') +
  #ggtitle() +
  theme(legend.position = 'top',
        legend.key.size = unit(50,'pt')) +
  xlab("Days after primary infection") +
  ylab("Fraction of sequences in the 10 largest clones") +
  scale_x_continuous(breaks = sort(as.integer(unique(clone_freqs_by_tissue_and_cell_type$day)))) +
  scale_size_continuous(breaks = c(100,1000,10000,50000), name = ' Number of sequences')

plot(fraction_in_top_10_clones_plot)

plot_clone_size_dist <- function(plot_cell_type, plot_tissue, plot_group, plot_abs_size = F,
                                 annotation = 'v_gene'){
  
  
  if(plot_abs_size){
    y_axis_var <- 'n_clone_seqs_in_compartment'
    y_axis_label <- 'Number of sequences'
  }else{
    y_axis_var <- 'clone_freq_in_compartment'
    y_axis_label <- 'Clone frequency'
  }
  
  plotting_data <- clone_freqs_by_tissue_and_cell_type %>%
    mutate(across(c('v_gene','d_gene','j_gene'),
                  function(x){str_remove(str_remove(x,'IGH'), '\\*[0-9]+')})) %>%
    unite('annotation', annotation, sep = ' ; ') 
  
  pl <- plotting_data %>%
    filter(total_seqs_in_compartment >= 100) %>%
    filter(compartment_cell_type == plot_cell_type, compartment_tissue == plot_tissue, group_controls_pooled == plot_group,
           clone_rank_in_compartment <= 20) %>%
    ungroup() %>%
    ggplot(aes_string(x = 'clone_rank_in_compartment', y = y_axis_var, group = 'mouse_id')) +
    geom_line() +
    facet_wrap('mouse_id', scales = 'free') +
    xlab('Clone rank (top 20 clones only)') +
    ylab(y_axis_label) +
    theme(legend.position = c(0.8,0.2)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.3 * length(annotation)))) +
    scale_x_continuous(expand = expansion(mult = c(0.05,0.2 * length(annotation)))) +
    scale_color_discrete(name = 'V gene deviation from naive frequency',
                         labels = c('negative','non-significant','positive'))
  
  if(length(annotation) == 1 & all(annotation == 'v_gene')){
    pl <- pl + 
      geom_point(aes(color = deviation_from_naive)) +
      geom_text(aes(x = clone_rank_in_compartment + 0.3, color = deviation_from_naive, angle = 30, hjust = 0,
                    label = annotation), show.legend = F, size = 3.5) 
  }else{
    pl <- pl + geom_point() +
      geom_text(aes(x = clone_rank_in_compartment + 0.3, angle = 35, hjust = 0,
                    label = annotation), show.legend = F, size = 3.5)
  }
  
  return(pl)
  
  
}
plot_clone_size_dist('PC','LN', 'primary-8', plot_abs_size = F) + theme(legend.position = c(0.7,0.2))
plot_clone_size_dist('PC','LN', 'primary-8', plot_abs_size = F, annotation = 'clone_consensus_cdr3_partis')
plot_clone_size_dist('PC','LN', 'primary-8', plot_abs_size = F, annotation = 'd_gene')
plot_clone_size_dist('PC','LN', 'primary-8', plot_abs_size = F, annotation = 'j_gene')
plot_clone_size_dist('PC','LN', 'primary-8', plot_abs_size = F, annotation = c('v_gene','j_gene',
                                                                               'clone_consensus_cdr3_partis'))
plot_clone_size_dist('PC','LN', 'primary-8', plot_abs_size = F, annotation = c('d_gene','j_gene',
                                                                               'clone_consensus_cdr3_partis'))
plot_clone_size_dist('PC','LN', 'primary-8', plot_abs_size = F, annotation = c('v_gene', 'mutations_above_threshold'))


plot_clone_size_dist('PC','LN', 'primary-16', plot_abs_size = F) + theme(legend.position = c(0.7,0.2))
plot_clone_size_dist('PC','LN', 'primary-16', plot_abs_size = F, annotation = 'clone_consensus_cdr3_partis')
plot_clone_size_dist('PC','LN', 'primary-16', plot_abs_size = F, annotation = 'd_gene')
plot_clone_size_dist('PC','LN', 'primary-16', plot_abs_size = F, annotation = 'j_gene')


plot_clone_size_dist('PC','LN', 'primary-24', plot_abs_size = F) + theme(legend.position = c(0.81,0.35))
plot_clone_size_dist('PC','LN', 'primary-24', plot_abs_size = F, annotation = 'clone_consensus_cdr3_partis')
plot_clone_size_dist('PC','LN', 'primary-24', plot_abs_size = F, annotation = c('v_gene', 'mutations_above_threshold'))


plot_clone_size_dist('GC','LN', 'primary-8', plot_abs_size = F) + theme(legend.position = c(0.81,0.35))
plot_clone_size_dist('GC','LN', 'primary-8', plot_abs_size = T) + theme(legend.position = c(0.81,0.35))

plot_clone_size_dist('GC','LN', 'primary-16', plot_abs_size = F) + theme(legend.position = c(0.7,0.2))
plot_clone_size_dist('GC','LN', 'primary-16', plot_abs_size = T) + theme(legend.position = c(0.7,0.2))
plot_clone_size_dist('GC','LN', 'primary-16', plot_abs_size = T, annotation = c('v_gene', 'mutations_above_threshold')) 

plot_clone_size_dist('GC','LN', 'primary-24', plot_abs_size = F) + theme(legend.position = c(0.81,0.35))
plot_clone_size_dist('GC','LN', 'primary-24', plot_abs_size = T) + theme(legend.position = c(0.81,0.35))
plot_clone_size_dist('GC','LN', 'primary-24', plot_abs_size = T, annotation = c('v_gene', 'mutations_above_threshold'))

plot_clone_size_dist('PC','BM', 'primary-16', plot_abs_size = F, c('v_gene', 'mutations_above_threshold')) 


plot_mutations_top_clones <- function(plot_cell_type, plot_tissue, plot_group, cdr3_only){
  
  if(cdr3_only){
    y_axis_var <- 'mean_cdr3_mutations_partis_aa'
    ymin <- 'min_cdr3_mutations_partis_aa'
    ymax <- 'max_cdr3_mutations_partis_aa'
    y_axis_label <- 'Mean number of amino acid mutations in CDR3'
  }else{
    y_axis_var <- 'mean_n_mutations_partis_aa'
    ymin <- 'min_n_mutations_partis_aa'
    ymax <- 'max_n_mutations_partis_aa'
    y_axis_label <- 'Mean number of amino acid mutations'
  }
  
  clone_freqs_by_tissue_and_cell_type %>%
    filter(compartment_tissue == plot_tissue, compartment_cell_type == plot_cell_type,
           group_controls_pooled == plot_group) %>%
    filter(total_seqs_in_compartment >= 100) %>%
    filter(clone_rank_in_compartment <= 100) %>%
    ggplot(aes_string(x = 'clone_rank_in_compartment', y = y_axis_var)) +
    geom_hline(yintercept = c(1), linetype = 2) +
    geom_linerange(aes_string(ymin = ymin, ymax = ymax)) +
    geom_point(aes(size = n_clone_seqs_in_compartment), color = 'dodgerblue') +
    #geom_text(aes(label = clone_id)) +
    facet_wrap('mouse_id', scales = 'free') +
    scale_x_log10() +
    xlab('Clone rank (top 100 clones only)') +
    ylab(paste0(y_axis_label,'\n(whisker = min, max)')) +
    scale_size_continuous(name = 'Clone size\n(n sequences)') +
    scale_y_continuous(limits = c(0,NA))
}
plot_mutations_top_clones('PC','LN','primary-8', cdr3_only = F) + theme(legend.position = c(0.75,0.3))
plot_mutations_top_clones('PC','LN','primary-16', cdr3_only = F) + theme(legend.position = c(0.75,0.3))
plot_mutations_top_clones('PC','LN','primary-24', cdr3_only = F) + theme(legend.position = c(0.01,0.94))
plot_mutations_top_clones('PC','LN','secondary-56', cdr3_only = F) + theme(legend.position = c(0.01,0.91))

plot_mutations_top_clones('PC','LN','primary-16', cdr3_only = T) + theme(legend.position = c(0.75,0.3))

plot_mutations_top_clones('GC','LN','primary-8', cdr3_only = F) + theme(legend.position = c(0.05,0.94))
plot_mutations_top_clones('GC','LN','primary-16', cdr3_only = F) + theme(legend.position = c(0.75,0.15))
plot_mutations_top_clones('GC','LN','primary-24', cdr3_only = F) + theme(legend.position = c(0.01,0.93))

# ------------ PAIRWISE CORRELATIONS BETWEEN MICE ----------------

# mean_neutral_pw_cor_freqs <- neutral_pairwise_correlations$freqs %>%
#   group_by(mouse_pair, pair_type, mouse_id_i, mouse_id_j, day_i, day_j, cell_type, tissue, total_compartment_seqs_i, total_compartment_seqs_j) %>%
#   summarise(mean_neutral_cor_coef = mean(cor_coef_freqs),
#             llim_neutral_cor_coef = quantile(cor_coef_freqs,0.025, na.rm = T),
#             ulim_neutral_cor_coef = quantile(cor_coef_freqs, 0.975, na.rm = T))
# 
# mean_neutral_pw_cor_freq_ratios <- neutral_pairwise_correlations$freq_ratios %>%
#   group_by(mouse_pair, pair_type, mouse_id_i, mouse_id_j, day_i, day_j, cell_type, tissue, total_compartment_seqs_i, total_compartment_seqs_j) %>%
#   summarise(mean_neutral_cor_coef = mean(cor_coef_freq_ratios),
#             median_neutral_cor_coef = median(cor_coef_freq_ratios),
#             llim_neutral_cor_coef = quantile(cor_coef_freq_ratios,0.025, na.rm = T),
#             ulim_neutral_cor_coef = quantile(cor_coef_freq_ratios, 0.975, na.rm = T))

# Get null percentile for observed correlation values
# null_percentiles_of_obs_values_freqs <- left_join(neutral_pairwise_correlations$freqs, pairwise_correlations$freqs %>%
#                                                     select(mouse_pair, cell_type, tissue, cor_coef_freqs) %>% dplyr::rename(obs_cor_coef_freqs = cor_coef_freqs)) %>%
#   group_by(mouse_pair, cell_type, tissue) %>%
#   summarise(null_percentile_of_obs_value = sum(obs_cor_coef_freqs >= cor_coef_freqs) / length(cor_coef_freqs)) %>%
#   ungroup()
# 
# null_percentiles_of_obs_values_freq_ratios <- left_join(neutral_pairwise_correlations$freq_ratios, pairwise_correlations$freq_ratios %>%
#                                                           select(mouse_pair, cell_type, tissue, cor_coef_freq_ratios) %>% dplyr::rename(obs_cor_coef_freq_ratios = cor_coef_freq_ratios)) %>%
#   group_by(mouse_pair, cell_type, tissue) %>%
#   summarise(null_percentile_of_obs_value = sum(obs_cor_coef_freq_ratios >= cor_coef_freq_ratios) / length(cor_coef_freq_ratios)) %>%
#   ungroup()


# Add neutral summary statistics to observed pairwise correlations
# pairwise_correlations$freqs <- left_join(pairwise_correlations$freqs, mean_neutral_pw_cor_freqs) %>%
#   mutate(neutrality_test = case_when(
#     cor_coef_freqs > ulim_neutral_cor_coef ~ 'Higher than neutral',
#     cor_coef_freqs < llim_neutral_cor_coef ~ 'Lower than neutral',
#     (cor_coef_freqs >= llim_neutral_cor_coef) & (cor_coef_freqs <= ulim_neutral_cor_coef) ~ 'Neutral'
#   ))
# 
# pairwise_correlations$freq_ratios <- left_join(pairwise_correlations$freq_ratios, mean_neutral_pw_cor_freq_ratios) %>%
#   mutate(neutrality_test = case_when(
#     cor_coef_freq_ratios > ulim_neutral_cor_coef ~ 'Higher than neutral',
#     cor_coef_freq_ratios < llim_neutral_cor_coef ~ 'Lower than neutral',
#     (cor_coef_freq_ratios >= llim_neutral_cor_coef) & (cor_coef_freq_ratios <= ulim_neutral_cor_coef) ~ 'Neutral'
#   ))

#pairwise_correlations$freqs$neutrality_test[pairwise_correlations$freqs$total_mouse_naive_seqs_i < 100 | pairwise_correlations$freqs$total_mouse_naive_seqs_j < 100] <- 'Insufficient naive seqs in one or both mice'
#pairwise_correlations$freq_ratios$neutrality_test[pairwise_correlations$freq_ratios$total_mouse_naive_seqs_i < 100 | pairwise_correlations$freq_ratios$total_mouse_naive_seqs_j < 100] <- 'Insufficient naive seqs in one or both mice'

#pairwise_correlations$freqs <- left_join(pairwise_correlations$freqs, null_percentiles_of_obs_values_freqs)
#pairwise_correlations$freq_ratios <- left_join(pairwise_correlations$freq_ratios, null_percentiles_of_obs_values_freq_ratios)


#pairwise_correlations$freqs <- pairwise_correlations$freqs %>%
#  mutate(neutrality_test = factor(neutrality_test, levels = c('Lower than neutral','Neutral','Higher than neutral','Insufficient naive seqs in one or both mice')))

#pairwise_correlations$freq_ratios <- pairwise_correlations$freq_ratios %>%
#  mutate(neutrality_test = factor(neutrality_test, levels = c('Lower than neutral','Neutral','Higher than neutral','Insufficient naive seqs in one or both mice')))


pairwise_naive_correlations_plot <- pairwise_correlations$freqs %>%
  filter(cell_type == 'naive') %>%
  filter(total_compartment_seqs_i >= 100, total_compartment_seqs_j >= 100) %>%
  ggplot(aes(x = pair_type, y = cor_coef_freqs, color = pair_type)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position = position_jitter(width = 0.1),
             alpha = point_alpha) +
  scale_y_continuous(limits = c(0,1)) +
  xlab('Type of pair') +
  ylab('Spearman correlation in naive V gene frequencies\nbetween mouse pairs (excluding mice with < 100 sequences)') +
  theme(legend.position = 'none') +
  background_grid()

plot(pairwise_naive_correlations_plot)

pairwise_correlations$freqs %>%
  filter(cell_type != 'naive') %>%
  mutate(tissue = factor(tissue, levels = c('LN','spleen','BM')),
         cell_type = factor(cell_type, levels = c('experienced','nonnaive_IgD+B220+','GC','PC','mem'))) %>%
  filter(total_compartment_seqs_i >= 100, total_compartment_seqs_j >= 100) %>%
  ggplot(aes(x = pair_type, y = cor_coef_freqs, color = pair_type)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position = position_jitter(width = 0.1),
             alpha = 0.5) +
  scale_y_continuous(limits = c(0,1)) +
  theme(legend.position = 'top') +
  xlab('Type of pair') +
  ylab('Correlation in gene frequencies between mouse pairs') +
  theme(legend.position = 'none',
        axis.text.x = element_text(size = 9)) +
  facet_grid(tissue~cell_type) +
  scale_x_discrete(labels = function(x){str_replace(x, '/','\n&\n')}) +
  background_grid()


pairwise_freq_correlations_plot <- pairwise_correlations$freqs %>%
  filter(cell_type %in% c('GC','PC','mem')) %>%
  mutate(cell_type = factor(cell_type, levels = c('GC','PC','mem'))) %>%
  filter(day_i == day_j, total_compartment_seqs_i >= 100, total_compartment_seqs_j >= 100,
         tissue == 'LN', pair_type %in% c('primary','secondary')) %>%
  mutate(day_i = as.integer(as.character(day_i))) %>%
  mutate(cell_type = case_when(
    cell_type == 'GC' ~ 'Lymph node GC cells',
    cell_type == 'PC' ~ 'Lymph node plasma cells',
    cell_type == 'mem' ~ 'Lymph node memory cells'
  )) %>%
  mutate(cell_type = factor(cell_type, levels = c('Lymph node GC cells',
                                                  'Lymph node plasma cells',
                                                  'Lymph node memory cells'))) %>%
  ggplot(aes(x = day_i, y = cor_coef_freqs, group = day_i)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position = position_jitter(width =0.1), alpha = 0.8, size = 4,
             aes(color = pair_type),
  ) +
  facet_wrap('cell_type') +
  background_grid() +
  scale_color_manual(values = c('green3','dodgerblue2'), guide = 'none') +
  xlab('Days after primary infection') +
  ylab('Correlation in V gene frequencies\nbetween mouse pairs (excluding mice with < 100 seqs.)') +
  theme(legend.position = 'top') +
  scale_x_continuous(breaks = c(8,16,24,40,56)) +
  geom_hline(yintercept = 0, linetype = 2)
plot(pairwise_freq_correlations_plot)

pairwise_freq_deviations_plot <- pairwise_correlations$freq_ratios %>%
  filter(cell_type %in% c('GC','PC','mem')) %>%
  mutate(cell_type = factor(cell_type, levels = c('GC','PC','mem'))) %>%
  filter(day_i == day_j, total_compartment_seqs_i >= 100, total_compartment_seqs_j >= 100,
         tissue == 'LN', pair_type %in% c('primary','secondary')) %>%
  mutate(day_i = as.integer(as.character(day_i))) %>%
  mutate(cell_type = case_when(
    cell_type == 'GC' ~ 'Lymph node GC cells',
    cell_type == 'PC' ~ 'Lymph node plasma cells',
    cell_type == 'mem' ~ 'Lymph node memory cells'
  )) %>%
  mutate(cell_type = factor(cell_type, levels = c('Lymph node GC cells',
                                                  'Lymph node plasma cells',
                                                  'Lymph node memory cells'))) %>%
  ggplot(aes(x = day_i, y = cor_coef_freq_ratios, group = day_i)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position = position_jitter(width =0.1), alpha = 0.8, size = 4,
             aes(color = pair_type),
  ) +
  facet_wrap('cell_type') +
  background_grid() +
  scale_color_manual(values = c('green3','dodgerblue2'), guide = 'none') +
  xlab('Days after primary infection') +
  ylab('Correlation in V gene frequency deviations\nbetween mouse pairs (excluding mice with < 100 seqs.)') +
  theme(legend.position = 'top') +
  scale_x_continuous(breaks = c(8,16,24,40,56)) +
  geom_hline(yintercept = 0, linetype = 2) 

plot(pairwise_freq_deviations_plot)

# ===== Linear models for correlation vs. time in LN plasma cells ======

# First, we subset values for the LN only, and for pairs of mice from the same time point and with a min. number of seqs.

# For observed data
freq_correlations_LN_OBS <- pairwise_correlations$freqs %>% filter(cell_type %in% c('GC','PC','mem')) %>%
  filter(day_i == day_j, total_compartment_seqs_i >= 100, total_compartment_seqs_j >= 100,
         tissue == 'LN', pair_type %in% c('primary','secondary')) %>%
  mutate(time = as.integer(as.character(day_i)))

deviation_correlations_LN_OBS <- pairwise_correlations$freq_ratios %>%
  filter(cell_type %in% c('GC','PC','mem')) %>%
  filter(day_i == day_j, total_compartment_seqs_i >= 100, total_compartment_seqs_j >= 100,
         tissue == 'LN', pair_type %in% c('primary','secondary')) %>%
  mutate(time = as.integer(as.character(day_i)))

# And for data with randomized noncontrol mice
freq_correlations_LN_randomized <- pairwise_correlations_randomized_noncontrol_groups$freqs %>%
  filter(cell_type %in% c('GC','PC','mem')) %>%
  filter(day_i == day_j, total_compartment_seqs_i >= 100, total_compartment_seqs_j >= 100,
         tissue == 'LN', pair_type %in% c('primary','secondary')) %>%
  mutate(time = as.integer(as.character(day_i)))

deviation_correlations_LN_randomized <- pairwise_correlations_randomized_noncontrol_groups$freq_ratios %>%
  filter(cell_type %in% c('GC','PC','mem')) %>%
  filter(day_i == day_j, total_compartment_seqs_i >= 100, total_compartment_seqs_j >= 100,
         tissue == 'LN', pair_type %in% c('primary','secondary')) %>%
  mutate(time = as.integer(as.character(day_i)))
  
# Now we compute the regression coefficient separately for each cell type, first for observed data
beta_freqs_OBS <- freq_correlations_LN_OBS  %>%
  group_by(cell_type) %>%
  summarise(beta_obs = lm(cor_coef_freqs~time)$coefficients[2])
 
beta_deviations_OBS <- deviation_correlations_LN_OBS  %>%
  group_by(cell_type) %>%
  summarise(beta_obs = lm(cor_coef_freq_ratios~time)$coefficients[2])

beta_freqs_primary_only_OBS <- freq_correlations_LN_OBS  %>%
  filter(pair_type == 'primary') %>%
  group_by(cell_type) %>%
  summarise(beta_obs = lm(cor_coef_freqs~time)$coefficients[2])

beta_deviations_primary_only_OBS <- deviation_correlations_LN_OBS  %>%
  filter(pair_type == 'primary') %>%
  group_by(cell_type) %>%
  summarise(beta_obs = lm(cor_coef_freq_ratios~time)$coefficients[2])


# Then for the randomizations

rdm_beta_freqs <- freq_correlations_LN_randomized %>%
  group_by(replicate, cell_type) %>%
  summarise(beta_rdm = lm(cor_coef_freqs~time)$coefficients[2]) %>%
  group_by(cell_type) 


rdm_beta_freqs_primary_only <- freq_correlations_LN_randomized %>%
  group_by(replicate, cell_type) %>%
  filter(pair_type == 'primary') %>%
  summarise(beta_rdm = lm(cor_coef_freqs~time)$coefficients[2]) %>%
  group_by(cell_type) 

rdm_beta_deviations <- deviation_correlations_LN_randomized %>%
  group_by(replicate, cell_type) %>%
  summarise(beta_rdm = lm(cor_coef_freq_ratios~time)$coefficients[2]) %>%
  group_by(cell_type) 

rdm_beta_deviations_primary_only<- deviation_correlations_LN_randomized %>%
  group_by(replicate, cell_type) %>%
  filter(pair_type == 'primary') %>%
  summarise(beta_rdm = lm(cor_coef_freq_ratios~time)$coefficients[2]) %>%
  group_by(cell_type) 

# Now we ask what % of randomizations equal to or more negative than observed values
  
left_join(rdm_beta_freq, beta_freqs_OBS) %>%
  group_by(cell_type) %>%
  summarise(n_replicates = n(),
            n_replicates_smaller_than_or_equal_to_obs = sum(beta_rdm <= beta_obs))

left_join(rdm_beta_freqs_primary_only, beta_freqs_primary_only_OBS) %>%
  group_by(cell_type) %>%
  summarise(n_replicates = n(),
            n_replicates_smaller_than_or_equal_to_obs = sum(beta_rdm <= beta_obs))

left_join(rdm_beta_deviations, beta_deviations_OBS) %>%
  group_by(cell_type) %>%
  summarise(n_replicates = n(),
            n_replicates_smaller_than_or_equal_to_obs = sum(beta_rdm <= beta_obs))

left_join(rdm_beta_deviations_primary_only, beta_deviations_primary_only_OBS) %>%
  group_by(cell_type) %>%
  summarise(n_replicates = n(),
            n_replicates_smaller_than_or_equal_to_obs = sum(beta_rdm <= beta_obs))


###### HIERARCHICAL CLUSTERING ANALYSIS

# Exclude mice with fewer than
min_seqs <- 100


# ----- Dendrograms and heatmaps based on correlation in V gene frequencies -------
pdf(file = paste0(figures_output_dir, 'freqs_heatmap_LN_PCs.pdf'), height = 7,
    width = 7.5)
get_vgene_freq_correlation_clustering(pairwise_correlations = pairwise_correlations,
                                      metric = 'freqs',
                                      cell_type = 'PC',
                                      tissue = 'LN',
                                      min_seqs = min_seqs)
dev.off()

pdf(file = paste0(figures_output_dir, 'freqs_heatmap_LN_GCs.pdf'), height = 7,
    width = 7.5)
get_vgene_freq_correlation_clustering(pairwise_correlations = pairwise_correlations,
                                      metric = 'freqs',
                                      cell_type = 'GC',
                                      tissue = 'LN',
                                      min_seqs = min_seqs)
dev.off()


get_vgene_freq_correlation_clustering(pairwise_correlations = pairwise_correlations,
                                      metric = 'freqs',
                                      cell_type = 'GC',
                                      tissue = 'spleen',
                                      min_seqs = min_seqs)

get_vgene_freq_correlation_clustering(pairwise_correlations = pairwise_correlations,
                                      metric = 'freqs',
                                      cell_type = 'GC',
                                      tissue = 'BM')

get_vgene_freq_correlation_clustering(pairwise_correlations = pairwise_correlations,
                                      metric = 'freqs',
                                      cell_type = 'PC',
                                      tissue = 'BM')

get_vgene_freq_correlation_clustering(pairwise_correlations = pairwise_correlations,
                                      metric = 'freqs',
                                      cell_type = 'PC',
                                      tissue = 'spleen')

pdf(file = paste0(figures_output_dir, 'freqs_heatmap_LN_mem.pdf'), height = 7,
    width = 7.5)
get_vgene_freq_correlation_clustering(pairwise_correlations = pairwise_correlations,
                                      metric = 'freqs',
                                      cell_type = 'mem',
                                      tissue = 'LN',
                                      min_seqs = min_seqs)
dev.off()

get_vgene_freq_correlation_clustering(pairwise_correlations = pairwise_correlations,
                                      metric = 'freqs',
                                      cell_type = 'mem',
                                      tissue = 'spleen')


get_vgene_freq_correlation_clustering(pairwise_correlations = pairwise_correlations,
                                      metric = 'freqs',
                                      cell_type = 'nonnaive_IgD+B220+',
                                      tissue = 'LN')

get_vgene_freq_correlation_clustering(pairwise_correlations = pairwise_correlations,
                                      metric = 'freqs',
                                      cell_type = 'nonnaive_IgD+B220+',
                                      tissue = 'BM')

get_vgene_freq_correlation_clustering(pairwise_correlations = pairwise_correlations,
                                      metric = 'freqs',
                                      cell_type = 'nonnaive_IgD+B220+',
                                      tissue = 'spleen')

# ----- Dendrograms and heatmaps based on correlation in V gene frequency deviations from naive repertoire -------

pdf(file = paste0(figures_output_dir, 'deviations_heatmap_LN_PCs.pdf'), height = 7,
    width = 7.5)
get_vgene_freq_correlation_clustering(pairwise_correlations = pairwise_correlations,
                                      metric = 'freq_ratios',
                                      cell_type = 'PC',
                                      tissue = 'LN',
                                      min_seqs = min_seqs)
dev.off()

pdf(file = paste0(figures_output_dir, 'deviations_heatmap_LN_GCs.pdf'), height = 7,
    width = 7.5)
get_vgene_freq_correlation_clustering(pairwise_correlations = pairwise_correlations,
                                      metric = 'freq_ratios',
                                      cell_type = 'GC',
                                      tissue = 'LN',
                                      min_seqs = min_seqs)
dev.off()

pdf(file = paste0(figures_output_dir, 'deviations_heatmap_LN_mem.pdf'), height = 7,
    width = 7.5)
get_vgene_freq_correlation_clustering(pairwise_correlations = pairwise_correlations,
                                      metric = 'freq_ratios',
                                      cell_type = 'mem',
                                      tissue = 'LN',
                                      min_seqs = min_seqs)
dev.off()

get_vgene_freq_correlation_clustering(pairwise_correlations = pairwise_correlations,
                                      metric = 'freq_ratios',
                                      cell_type = 'nonnaive_IgD+B220+',
                                      tissue = 'LN')



get_vgene_freq_correlation_clustering(pairwise_correlations = pairwise_correlations,
                                      metric = 'freq_ratios',
                                      cell_type = 'PC',
                                      tissue = 'spleen')

# PLOTS TO EXPORT
save(naive_exp_correlations_plot,
     pairwise_naive_correlations_plot,
     pairwise_freq_correlations_plot,
     pairwise_freq_deviations_plot,
     file = paste0(figures_output_dir, 'correlation_plots.RData')
     )

save(top_genes_LN_PC_day8_plot,
     top_genes_LN_PC_day16_plot,
     top_genes_LN_GC_day8_plot,
     top_genes_LN_GC_day16_plot,
     file = paste0(figures_output_dir, 'top_genes_plots.RData'))

save(fraction_in_top_10_clones_plot,
     file = paste0(figures_output_dir, 'fraction_in_top_10_clones_plot.RData'))




# Clustering based on Euclidean distance of V gene frequencies
# get_vgene_freq_clustering <- function(gene_freqs, cell_type, tissue, n_top_genes = NULL){
#   
#   if(!is.null(n_top_genes)){
#     wide_format_freqs <- gene_freqs %>% 
#       filter(group_controls_pooled != 'control', v_gene_rank <= n_top_genes) 
#   }else{
#     wide_format_freqs <- gene_freqs %>% 
#       filter(group_controls_pooled != 'control')
#   }
#   
#   wide_format_freqs <- wide_format_freqs %>%
#     filter(total_compartment_seqs >= 100, total_mouse_naive_seqs >= 100) %>%
#     filter(cell_type == !!cell_type, tissue == !!tissue) %>%
#     select(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type, total_compartment_seqs, v_gene,
#            vgene_seq_freq) %>%
#     pivot_wider(names_from = v_gene, values_from = vgene_seq_freq, values_fill = 0)
#   
#   gene_freq_matrix <- as.matrix(wide_format_freqs[grepl('IGHV', colnames(wide_format_freqs))])
#   rownames(gene_freq_matrix) <- wide_format_freqs$mouse_id
#   
#   # Remove all-zero columns
#   gene_freq_matrix <- gene_freq_matrix[,colSums(gene_freq_matrix) != 0]
#   
# 
#   z_score_matrix <- apply(gene_freq_matrix, 2,
#                             function(x){(x - mean(x))/sd(x)})
#   
#   
#   z_score_dist_matrix <- dist(z_score_matrix, method = 'euclidean')
#   cluster <- hclust(z_score_dist_matrix, method = 'complete')
#   
#   cluster_for_ggplot <- dendro_data(cluster, type="rectangle")
#   
#   cluster_for_ggplot$labels <-
#     left_join(cluster_for_ggplot$labels, 
#               wide_format_freqs %>% select(mouse_id, group_controls_pooled) %>% dplyr::rename(label = mouse_id))
#    
#   ggplot() +
#     geom_segment(data = segment(cluster_for_ggplot),
#                  aes(x = x, y = y, xend = xend, yend = yend)) +
#     geom_point(data = label(cluster_for_ggplot), 
#               aes(x = x, y = y, color = group_controls_pooled), 
#               size = 3
#     ) +
#     geom_text(data = label(cluster_for_ggplot), 
#               aes(x = x, y = y-1.5, label = label, color = group_controls_pooled), 
#               size = 3
#     ) +
#     coord_flip() +
#     scale_y_reverse(expand = c(0.2, 0)) +
#     theme(axis.line = element_blank(),
#           axis.ticks = element_blank(),
#           axis.text = element_blank(),
#           axis.title = element_blank(),
#           legend.position = c(0.05,0.7),
#           plot.title = element_text(size = 12)) +
#     scale_color_discrete(name = 'Group', type = 'qual') +
#     ggtitle('Clustering based on V gene frequencies in lymph node plasma cells')
#   
#     
# }






