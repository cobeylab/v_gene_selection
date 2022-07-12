library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(scales)
library(viridis)

theme_set(theme_cowplot())

source('gene_frequency_functions.R') 
source('plot_options.R')

args <- commandArgs(trailingOnly = T)
frequency_type <- args[1] # frequency_type <- 'all_seqs'
use_Greiff2017_naive_freqs <- as.logical(args[2]) # use_Greiff2017_naive_freqs <- F
collapsed_novel_alleles <- as.logical(args[3]) # collapsed_novel_alleles <- F

# For some analyses, exclude mice with fewer than min_compartment_size reads
min_compartment_size = 100

# ===== LOAD PRE-COMPUTED OBJECTS AND DEFINE OUTPUT DIRECTORIES =====

results_directory <- '../results/'

processed_data_directory <- '../processed_data/'

precomputed_freqs_file <- paste0('precomputed_gene_freqs_', frequency_type, '.RData')

germline_mutability_by_region <- read_csv('../results/germline_mutability_by_region.csv')
germline_mutability_by_region_type <- read_csv('../results/germline_mutability_by_region_type.csv')


figure_directory <- paste0('../figures/', frequency_type, '_freqs/')
deviations_by_allele_results_path <- paste0(results_directory, 'deviations_by_allele_', frequency_type, '_freqs.csv')


if(use_Greiff2017_naive_freqs){
  stopifnot(frequency_type == 'all_seqs')
  stopifnot(collapsed_novel_alleles == F)
  figure_directory <- '../figures/all_seqs_freqs_Greiff2017_naive_freqs/'
  precomputed_freqs_file <- 'precomputed_gene_freqs_all_seqs_Greiff2017_naive_freqs.RData'
  deviations_by_allele_results_path  <-paste0(results_directory, 'deviations_by_allele_', frequency_type, 'freqs_Greiff2017_naive_freqs.csv')
}
if(collapsed_novel_alleles){
  stopifnot(frequency_type == 'all_seqs')
  stopifnot(use_Greiff2017_naive_freqs == F)
  figure_directory <- '../figures/all_seqs_freqs_collapsed_novel_alleles/'
  precomputed_freqs_file <- 'precomputed_gene_freqs_all_seqs_collapsed_novel_alleles.RData'
  deviations_by_allele_results_path  <-paste0(results_directory, 'deviations_by_allele_', frequency_type, '_collapsed_novel_alleles.csv')
  
}

# Load precomputed gene frequencies, neutral realizations, pairwise correlations 
load(paste0(results_directory, precomputed_freqs_file))

# Create necessary directories
#dir.create(figure_directory, showWarnings = F, recursive = T)
exported_figure_objects_dir <- paste0(figure_directory,'exported_ggplot_objects/')
dir.create(exported_figure_objects_dir, showWarnings = F, recursive = T)


# ===== DEFINE SOME ADDITIONAL FUNCTIONS =====
plot_naive_exp_correlations <- function(naive_exp_correlations, method){
  naive_exp_correlations %>%
    filter(tissue == 'LN', method == !!method) %>%
    filter(cell_type %in% c('GC','PC','mem')) %>%
    cell_type_facet_labeller() %>%
    mutate(day = ifelse(group_controls_pooled == 'control', 0, day)) %>%
    ggplot(aes(x = day, y = naive_exp_corr, color = infection_status, group = day)) +
    geom_boxplot(outlier.alpha = 0) +
    geom_point(aes(size = total_compartment_seqs),
               position = position_jitter(width = 0.1, height = 0),
               alpha = point_alpha) +
    facet_grid(.~cell_type, scales = 'free') +
    theme(legend.position = 'top') +
    xlab('Days after primary infection') +
    ylab('Correlation in V gene frequencies between\nnaive repertoire and influenza-induced populations') +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_size_continuous(name = 'Number of sequences',
                          breaks = c(1000,10000,20000)) +
    guides(color = 'none') +
    background_grid() +
    groups_color_scale(name = 'Infection') +
    label_controls_as_day_0
    #scale_color_discrete(name = 'Infection')
}

plot_most_common_genes <- function(plot_group, gene_freqs, plot_cell_type, plot_tissue, max_rank,
                                   min_compartment_size){
  gene_freqs %>% 
    filter(cell_type == plot_cell_type, tissue == plot_tissue, group_controls_pooled == plot_group, v_gene_rank <= max_rank) %>%
    filter(total_compartment_seqs >= min_compartment_size, total_mouse_naive_seqs >= min_compartment_size) %>%
    rowwise() %>%
    mutate(label_position = ifelse(vgene_seq_freq > naive_vgene_seq_freq, 1.05*vgene_seq_freq, 1.05*naive_vgene_seq_freq)) %>%
    ggplot() +
    geom_point(aes(x = v_gene_rank, y = naive_vgene_seq_freq, color = deviation_from_naive)) +
    geom_segment(aes(x = v_gene_rank, xend = v_gene_rank, y = naive_vgene_seq_freq, yend = vgene_seq_freq,
                     color = deviation_from_naive),
                 arrow = arrow(ends = 'last', length = unit(10, "pt"), type = 'closed')) +
    geom_text(aes(label = str_remove(v_gene, 'IGHV'),
                  x = v_gene_rank, y = label_position), angle = 20, size = 3, alpha = 0.8) +
    facet_wrap('mouse_id', scales = 'free', ncol = 1) +
    xlab(paste0('V gene rank in ', plot_group, ' ', plot_tissue, ' ',
                switch(plot_cell_type, 'PC' = 'plasma cells', 'GC' = 'germinal center cells', 'mem' = 'memory cells',
                       'nonnaive_IgD+B220+' = 'Non-naive IgD+B220+ cells'))) +
    ylab('V gene frequency') +
    theme(legend.position = 'bottom') + 
    scale_x_continuous(expand = c(0.15,0), breaks = 1:max_rank) +
    #scale_color_discrete(name = 'Deviation from naive\nfrequency (bootstrap)', labels = c('negative','non-significant','positive'))
    deviations_color_scale(name = 'Deviation from naive\nfrequency (bootstrap)')
  
}

plot_focal_genes <- function(gene_freqs, clone_freqs_by_tissue_and_cell_type,
                             focal_genes, min_compartment_size){
  
  base_function <- function(gene_freqs, clone_freqs_by_tissue_and_cell_type,
                            tissue, min_compartment_size){
    base_gene_freq_data <- gene_freqs %>% filter(tissue == !!tissue,
                                                 cell_type %in% c('GC','mem','PC')) %>%
      filter(total_compartment_seqs >= min_compartment_size, total_mouse_naive_seqs >= min_compartment_size,
             v_gene %in% focal_genes) %>%
      mutate(v_gene = factor(v_gene, levels = focal_genes)) 
    
    base_gene_freq_data <- cell_type_facet_labeller(base_gene_freq_data )
    base_gene_freq_data <- set_controls_as_day_0(base_gene_freq_data)
    
    clone_ranks_data <- clone_freqs_by_tissue_and_cell_type %>%
      filter(compartment_tissue == !!tissue,
             compartment_cell_type %in% c('GC','mem','PC')) %>%
      filter(total_seqs_in_compartment >= min_compartment_size,
             total_mouse_naive_seqs >= min_compartment_size,
             v_gene %in% focal_genes) %>%
      mutate(v_gene = factor(v_gene, levels = focal_genes)) %>%
      group_by(mouse_id, day, group_controls_pooled, infection_status, compartment_cell_type, compartment_tissue,
               v_gene) %>%
      summarise(rank_biggest_clone = min(clone_rank_in_compartment)) %>%
      ungroup()
    
    clone_ranks_data <- cell_type_facet_labeller(clone_ranks_data)
    clone_ranks_data <- set_controls_as_day_0(clone_ranks_data)
    
    x_axis_breaks <- sort(as.integer(unique(base_gene_freq_data$day)))
    x_axis_labels <- c('control', sort(as.integer(unique(gene_freqs$day))))
    
    freq_ratios_pl <- base_gene_freq_data %>%
      ggplot(aes(x = day, y = obs_rho + 1e-4)) +
      geom_boxplot(aes(group = day), outlier.alpha = 0) +
      geom_point(aes(color = deviation_from_naive), position = position_jitter(width = 0.5, height = 0), size = 4, alpha = 0.5) +
      facet_grid(cell_type~v_gene, scales = 'free') +
      geom_hline(yintercept = 1, linetype =2) +
      theme(legend.position = 'top',
            panel.border = element_rect(color = 'black')) +
      scale_y_log10() +
      deviations_color_scale(name =  'Deviation from naive repertoire (bootstrap)') +
      xlab("Days after primary infection") +
      ylab("Ratio of experienced-to-naive frequencies \n (log10 + 1e-4)") +
      label_controls_as_day_0 + 
      background_grid()
    
    freqs_pl <- base_gene_freq_data %>%
      ggplot(aes(x = day, y = vgene_seq_freq + 1e-4)) +
      geom_boxplot(aes(group = day), outlier.alpha = 0) +
      geom_point(aes(color = deviation_from_naive), position = position_jitter(width = 0.5, height = 0), size = 4, alpha = 0.5) +
      facet_grid(cell_type~v_gene, scales = 'free') +
      theme(legend.position = 'top',
            panel.border = element_rect(color = 'black')) +
      scale_y_log10() +
      deviations_color_scale(name =  'Deviation from naive repertoire (bootstrap)') +
      xlab("Days after primary infection") +
      label_controls_as_day_0 + 
      ylab("Germline allele frequency \n (log10 + 1e-4)") +
      background_grid()
    
    clone_ranks_pl <- clone_ranks_data %>%
      ggplot(aes(x = day, y = rank_biggest_clone)) +
      geom_boxplot(aes(group = day), outlier.alpha = 0) +
      geom_point(aes(color = infection_status), position = position_jitter(width = 0.5, height = 0), size = 4, alpha = 0.5) +
      facet_grid(compartment_cell_type~v_gene, scales = 'free') +
      theme(legend.position = 'top',
            panel.border = element_rect(color = 'black')) +
      scale_y_log10() +
      deviations_color_scale(name =  'Deviation from naive repertoire (bootstrap)') +
      xlab("Days after primary infection") +
      ylab("Rank of biggest clone") +
      label_controls_as_day_0 +
      background_grid()
    
    return(list(freq_ratios = freq_ratios_pl, freqs = freqs_pl, clone_ranks = clone_ranks_pl))
    
  }
  
  plots_list <- lapply(as.list(unique(gene_freqs$tissue)),
                       FUN = base_function,
                       gene_freqs = gene_freqs,
                       clone_freqs_by_tissue_and_cell_type = clone_freqs_by_tissue_and_cell_type,
                       min_compartment_size = min_compartment_size)
  names(plots_list) <- unique(gene_freqs$tissue)
  return(plots_list)
  
}

# ===== CUMULATIVE NAIVE FREQUENCIES =====

cumulative_naive_freqs <- naive_freqs %>%
  group_by(mouse_id) %>%
  filter(total_mouse_naive_seqs >= min_compartment_size) %>%
  arrange(mouse_id, naive_vgene_seq_freq) %>%
  mutate(cumulative_naive_freq = cumsum(naive_vgene_seq_freq)) %>%
  ungroup() %>%
  mutate(group_controls_pooled = factor(group_controls_pooled, levels = group_controls_pooled_factor_levels))

cumulative_naive_freqs_pl <- cumulative_naive_freqs %>% 
  ggplot(aes(x = naive_vgene_seq_freq, y = cumulative_naive_freq, group = mouse_id)) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = c(0.02,0.03), linetype = 2)  +
  geom_hline(yintercept = 0.5, linetype = 2) +
  facet_wrap('group_controls_pooled') +
  xlab('Allele frequency in the naive repertoire') +
  ylab('Cumulative frequency')


# ===== CORRELATION BETWEEN NAIVE AND EXPERIENCED FREQUENCIES ==========

naive_exp_correlations_obs <- get_naive_exp_correlations(gene_freqs) %>%
  filter(total_compartment_seqs > min_compartment_size, total_mouse_naive_seqs > min_compartment_size) %>%
  mutate(group_controls_pooled = factor(group_controls_pooled, levels = group_controls_pooled_factor_levels)) %>%
  mutate(day = as.integer(as.character(day)))

naive_exp_pearson_corr_plot <- plot_naive_exp_correlations(naive_exp_correlations_obs, method = 'pearson')
naive_exp_spearman_corr_plot <- plot_naive_exp_correlations(naive_exp_correlations_obs, method = 'spearman')

# ====== CREATE TABLE WITH NUMBER OF DEVIATION CALLS FOR EACH ALLELE =======
deviations_by_allele <- gene_freqs %>%
  filter(total_compartment_seqs >= min_compartment_size, total_mouse_naive_seqs >= min_compartment_size) %>%
  group_by(group_controls_pooled, tissue, cell_type, v_gene, deviation_from_naive) %>%
  count() %>%
  pivot_wider(names_from = deviation_from_naive, values_from = n,values_fill = 0) %>%
  arrange(group_controls_pooled, tissue, cell_type, desc(positive))

write_csv(deviations_by_allele, deviations_by_allele_results_path)

# ====== PLOT THE MOST ABUNDANT ALLELES =======
# Plots showing deviation from naive freq for major genes
infected_groups <- c('primary-8', 'primary-16', 'primary-24', 'secondary-40', 'secondary-56')

top_genes_LN_GC <- sapply(infected_groups,
                          FUN = plot_most_common_genes, gene_freqs = gene_freqs,
                          plot_cell_type = 'GC', plot_tissue = 'LN', max_rank = 10,
                          min_compartment_size = min_compartment_size,
                          USE.NAMES = T, simplify = F)

top_genes_LN_PC <- sapply(infected_groups,
                          FUN = plot_most_common_genes, gene_freqs = gene_freqs,
                          plot_cell_type = 'PC', plot_tissue = 'LN', max_rank = 10,
                          min_compartment_size = min_compartment_size,
                          USE.NAMES = T, simplify = F)

top_genes_LN_mem <- sapply(infected_groups,
                          FUN = plot_most_common_genes, gene_freqs = gene_freqs,
                          plot_cell_type = 'mem', plot_tissue = 'LN', max_rank = 10,
                          min_compartment_size = min_compartment_size,
                          USE.NAMES = T, simplify = F)

# More detailed plot for day 8 LN PCs (for main text fig.)
top_20_genes_day8_LN_PC <- plot_most_common_genes(plot_group = 'primary-8', gene_freqs = gene_freqs,
                                                  plot_cell_type = 'PC', plot_tissue = 'LN', max_rank = 20,
                                                  min_compartment_size = min_compartment_size)

# ======= MAKE DETAILED PLOTS TRACKING THE FATE OF FOCAL GENES =======
# How genes consistently overrepresented on day 8 LN PCs (14-4, 1-69, 1-82) fare at later time points and other cell types
focal_alleles_plot <- plot_focal_genes(gene_freqs = gene_freqs,
                                             clone_freqs_by_tissue_and_cell_type = clone_freqs_by_tissue_and_cell_type,
                                             focal_genes = c('IGHV14-4*01', 'IGHV1-69*01', 'IGHV1-82*01'),
                                             min_compartment_size = min_compartment_size)

# ======= PAIRWISE CORRELATIONS BETWEEN MICE =======

pairwise_naive_correlations_plot <- pairwise_correlations$freqs %>%
  filter(cell_type == 'naive', method == 'pearson') %>% 
  filter(total_compartment_seqs_i >= min_compartment_size, total_compartment_seqs_j >= min_compartment_size) %>%
  ggplot(aes(x = pair_type, y = cor_coef_freqs, color = pair_type)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position = position_jitter(width = 0.1, height = 0),
             alpha = point_alpha) +
  scale_y_continuous(limits = c(0,1)) +
  xlab('Type of pair') +
  ylab('Correlation in naive V gene frequencies\nbetween mouse pairs (excluding mice with < 100 sequences)') +
  theme(legend.position = 'none') +
  background_grid() +
  pair_types_color_scale(name = '')

pairwise_freq_correlations_plot <- pairwise_correlations$freqs %>%
  filter(cell_type %in% c('GC','PC','mem'), method == 'pearson') %>%
  filter((day_i == day_j | pair_type == 'control'),
         total_compartment_seqs_i >= min_compartment_size,
         total_compartment_seqs_j >= min_compartment_size,
         tissue == 'LN', pair_type %in% c('control','primary','secondary')) %>%
  mutate(day_i = ifelse(pair_type == 'control', 0, as.integer(as.character(day_i)))) %>%
  cell_type_facet_labeller() %>%
    ggplot(aes(x = day_i, y = cor_coef_freqs, group = day_i)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position = position_jitter(width =0.1, height = 0), alpha = 0.8, size = 4,
             aes(color = pair_type),
  ) +
  facet_grid(method ~ cell_type) +
  background_grid() +
  xlab('Days after primary infection') +
  ylab('Correlation in V gene frequencies\nbetween mouse pairs (excluding mice with < 100 seqs.)') +
  theme(legend.position = 'top') +
  scale_x_continuous(breaks = c(0,8,16,24,40,56),
                     labels = c('control', '8','16','24','40','56')) +
  geom_hline(yintercept = 0, linetype = 2) +
  groups_color_scale(name = 'Infection')

pairwise_freq_deviations_plot <- pairwise_correlations$freq_ratios %>%
  filter(cell_type %in% c('GC','PC','mem'), method == 'pearson') %>%
  filter((day_i == day_j | pair_type == 'control'),
         total_compartment_seqs_i >= min_compartment_size,
         total_compartment_seqs_j >= min_compartment_size,
         total_mouse_naive_seqs_i >= min_compartment_size, # Since estimating deviations requires naive freqs,
         total_mouse_naive_seqs_j >= min_compartment_size, # ...remove mice with too few naive seqs.
         tissue == 'LN', pair_type %in% c('control','primary','secondary')) %>%
  mutate(day_i = ifelse(pair_type == 'control', 0, as.integer(as.character(day_i)))) %>%
  cell_type_facet_labeller() %>%
  ggplot(aes(x = day_i, y = cor_coef_freq_ratios, group = day_i)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position = position_jitter(width =0.1, height = 0), alpha = 0.8, size = 4,
             aes(color = pair_type),
  ) +
  facet_grid(method~cell_type) +
  background_grid() +
  #scale_color_manual(values = c('green3','dodgerblue2'), guide = 'none') +
  xlab('Days after primary infection') +
  ylab('Correlation in V gene frequency deviations\nbetween mouse pairs (excluding mice with < 100 seqs.)') +
  theme(legend.position = 'top') +
  scale_x_continuous(breaks = c(0,8,16,24,40,56),
                     labels = c('control', '8','16','24','40','56')) +
  geom_hline(yintercept = 0, linetype = 2) +
  groups_color_scale(name = 'Infection') 

# ======= FREQS AND FREQ. DEVIATIONS VS. MUTABILITY =======

freq_ratio_mutability_correlations <- get_freq_ratio_mutability_correlations(
  gene_freqs,
  germline_mutability_by_region_type,
  min_compartment_size = min_compartment_size,
  method = 'pearson') 

freq_ratio_mutability_correlations_pl <- freq_ratio_mutability_correlations %>%
  filter(cell_type %in% c('GC','PC', 'mem'), tissue == 'LN') %>%
  filter(mutability_metric %in% c('average_RS5NF_mutability_cdr', 'average_RS5NF_mutability_fwr')) %>%
  mutate(mutability_metric = case_when(
    mutability_metric == 'average_RS5NF_mutability_cdr' ~ 'CDRs',
    mutability_metric == 'average_RS5NF_mutability_fwr' ~ 'FRs'
  )) %>%
  set_controls_as_day_0() %>%
  cell_type_facet_labeller() %>%
  ggplot(aes(x = day, y = correlation, color = infection_status)) +
  geom_boxplot(outlier.alpha = 0,aes(group = day)) +
  geom_point(aes(size = total_compartment_seqs), alpha = 0.5) +
  facet_grid(mutability_metric~cell_type) +
  geom_hline(yintercept = 0, linetype = 2) +
  theme(legend.position = 'top') + 
  scale_size_continuous(name = 'Number of sequences', breaks = c(1000,10000,50000)) +
  label_controls_as_day_0 +
  xlab('Days after primary infection') + 
  ylab('Correlation between average RS5NF mutability\nand experienced-to-naive frequency ratio') +
  background_grid() +
  groups_color_scale(name = 'Infection')

germline_mutability_by_region <- germline_mutability_by_region %>%
  mutate(region = case_when(
    region == 'cdr1' ~ 'CDR1',
    region == 'cdr2' ~ 'CDR2',
    region == 'fwr1' ~ 'FR1',
    region == 'fwr2' ~ 'FR2',
    region == 'fwr3' ~ 'FR3',
    region == 'whole_sequence' ~ 'whole sequence'
  ))

germline_mutability_by_region_pl <- germline_mutability_by_region  %>%
  ggplot(aes(x = average_RS5NF_mutability)) +
  geom_histogram() +
  facet_wrap('region') +
  geom_vline(data = germline_mutability_by_region %>% filter(v_gene %in% c('IGHV14-4*01',
                                                                           'IGHV1-82*01',
                                                                           'IGHV1-69*01')),
             aes(xintercept = average_RS5NF_mutability, color = v_gene), size = 1.5) +
  theme(legend.position = 'top') +
  xlab('Average RS5NF mutability') +
  ylab('Number of alleles') +
  scale_color_brewer(name = 'V allele', type = 'qual') +
  background_grid()


# Null model with fixed lineage sizes but randomized V alleles

# Compute summary statistics for randomizations
summary_pwcorr_randomized_lineage_V_alleles <-
  list(freqs = pairwise_correlations_randomized_lineage_V_alleles$freqs %>%
         filter(total_compartment_seqs_i >= min_compartment_size, total_compartment_seqs_i >= min_compartment_size,
                total_mouse_naive_seqs_i >= min_compartment_size, total_mouse_naive_seqs_j >= min_compartment_size) %>%
         group_by(pair_type, day_i, tissue, cell_type, method) %>%
         summarise(cor_coef_lower = quantile(cor_coef_freqs, 0.25, na.rm = T),
                   cor_coef_upper = quantile(cor_coef_freqs, 0.75, na.rm = T),
                   cor_coef_freqs = median(cor_coef_freqs, na.rm = T)),
       freq_ratios = pairwise_correlations_randomized_lineage_V_alleles$freq_ratios %>%
         filter(total_compartment_seqs_i >= min_compartment_size, total_compartment_seqs_i >= min_compartment_size,
                total_mouse_naive_seqs_i >= min_compartment_size, total_mouse_naive_seqs_j >= min_compartment_size) %>%
         group_by(pair_type, day_i, tissue, cell_type, method) %>%
         summarise(cor_coef_lower = quantile(cor_coef_freq_ratios, 0.25, na.rm = T),
                   cor_coef_upper = quantile(cor_coef_freq_ratios, 0.75, na.rm = T),
                   cor_coef_freq_ratios = median(cor_coef_freq_ratios, na.rm = T))
  )



pw_freq_cors_randomized_lineage_V_alleles <- pairwise_freq_correlations_plot +
  geom_pointrange(data = summary_pwcorr_randomized_lineage_V_alleles$freqs %>% 
                    filter(cell_type %in% c('GC','PC','mem'),
                           method == 'pearson') %>%
                    cell_type_facet_labeller() %>%
                    mutate(day_i = ifelse(pair_type == 'control', 0, as.integer(day_i))),
                  color = '#FF3333',
                  shape = 15,
                  aes(ymin = cor_coef_lower,
                      ymax = cor_coef_upper),
                  size = 1.2
  )


pw_freq_deviations_randomized_lineage_V_alleles <- pairwise_freq_deviations_plot +
  geom_pointrange(data = summary_pwcorr_randomized_lineage_V_alleles$freq_ratios %>% 
                    filter(cell_type %in% c('GC','PC','mem'),
                           method == 'pearson') %>%
                    cell_type_facet_labeller() %>%
                    mutate(day_i = ifelse(pair_type == 'control', 0, as.integer(day_i))),
                  color = '#FF3333',
                  shape = 15,
                  aes(ymin = cor_coef_lower,
                      ymax = cor_coef_upper),
                  size = 1.2
  ) 

# PLOTS TO EXPORT
save(cumulative_naive_freqs_pl,
     file = paste0(exported_figure_objects_dir, 'cumulative_naive_freqs.RData'))

save(naive_exp_pearson_corr_plot,
     naive_exp_spearman_corr_plot,
     pairwise_naive_correlations_plot,
     pairwise_freq_correlations_plot,
     pairwise_freq_deviations_plot,
     freq_ratio_mutability_correlations_pl,
     germline_mutability_by_region_pl,
     pw_freq_cors_randomized_lineage_V_alleles,
     pw_freq_deviations_randomized_lineage_V_alleles,
     file = paste0(exported_figure_objects_dir, 'correlation_plots.RData')
     )

save(top_genes_LN_GC,
     top_genes_LN_PC,
     top_genes_LN_mem,
     top_20_genes_day8_LN_PC,
     focal_alleles_plot,
     file = paste0(exported_figure_objects_dir, 'top_genes_plots.RData'))
