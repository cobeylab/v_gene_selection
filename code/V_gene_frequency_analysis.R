library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggdendro)
library(cowplot)
library(scales)
library(viridis)

theme_set(theme_cowplot())

source('gene_frequency_functions.R')
source('plot_options.R')

args <- commandArgs(trailingOnly = T)


frequency_type <- args[1] #frequency_type <- 'all_seqs'
use_Greiff2017_naive_freqs <- as.logical(args[2]) # use_Greiff2017_naive_freqs <- F


results_directory <- '../results/'
#results_directory <- '~/Desktop/v_gene_selection/results/'

processed_data_directory <- '../processed_data/'
#processed_data_directory <- '~/Desktop/v_gene_selection/processed_data/'

figure_directory <- paste0('../figures/', frequency_type, '_freqs/')
#figure_directory <- paste0('~/Desktop/v_gene_selection/figures/',frequency_type, '_freqs/')

precomputed_freqs_file <- paste0('precomputed_gene_freqs_', frequency_type, '.RData')

correlation_rdm_test_results <- paste0(results_directory, 'correlation_rdm_tests_', frequency_type, '_freqs.csv')
deviations_by_allele_results_path <- paste0(results_directory, 'deviations_by_allele_', frequency_type, '_freqs.csv')

if(use_Greiff2017_naive_freqs){
  stopifnot(frequency_type == 'all_seqs')
  figure_directory <- '../figures/all_seqs_freqs_Greiff2017_naive_freqs/'
  precomputed_freqs_file <- 'precomputed_gene_freqs_all_seqs_Greiff2017_naive_freqs.RData'
  correlation_rdm_test_results <- paste0(results_directory, 'correlation_rdm_tests_', frequency_type, 'freqs_Greiff2017_naive_freqs.csv')
  deviations_by_allele_results_path  <-paste0(results_directory, 'deviations_by_allele_', frequency_type, 'freqs_Greiff2017_naive_freqs.csv')
}


# Load precomputed gene frequencies, neutral realizations, pairwise correlations 
load(paste0(results_directory, precomputed_freqs_file))

# These count objects are from the precomputed file, so they use either all sequences or unique sequences
seq_counts <- bind_rows(exp_seq_counts, naive_seq_counts)

# Basic info for each clone (germline genes, CDR lenght, naive CDR seq)
clone_info <- read_csv(paste0(processed_data_directory,'clone_info.csv'))

exported_figure_objects_dir <- paste0(figure_directory,'exported_ggplot_objects/')
dir.create(exported_figure_objects_dir, showWarnings = F)

# For some analyses, exclude mice with fewer than min_compartment_size reads
min_compartment_size = 100

# ===== CUMULATIVE NAIVE FREQUENCIES

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

plot_naive_exp_correlations <- function(naive_exp_correlations_obs, naive_exp_correlations_neutral, method){
  naive_exp_correlations_obs %>%
    filter(tissue == 'LN', group_controls_pooled != 'control', method == !!method) %>%
    filter(cell_type %in% c('GC','PC','mem')) %>%
    cell_type_facet_labeller() %>%
    ggplot(aes(x = day, y = naive_exp_corr, color = infection_status, group = day)) +
    geom_boxplot(data = naive_exp_correlations_neutral %>%
                   filter(tissue == 'LN', group_controls_pooled != 'control',
                          cell_type %in% c('GC','PC','mem'), method == !!method) %>%
                   cell_type_facet_labeller(),
                 outlier.alpha = 0, linetype = 2) +
    geom_boxplot(outlier.alpha = 0) +
    geom_point(aes(size = total_compartment_seqs),
               position = position_jitter(width = 0.1, height = 0),
               alpha = point_alpha) +
    facet_grid(.~cell_type, scales = 'free') +
    theme(legend.position = 'top') +
    xlab('Days after primary infection') +
    ylab('Correlation in V gene frequencies between\nnaive repertoire and influenza-induced populations') +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_color_manual(values = c('green3','dodgerblue2')) +
    scale_size_continuous(name = 'Number of sequences',
                          breaks = c(1000,10000,20000)) +
    guides(color = 'none') +
    background_grid() +
    scale_x_continuous(breaks = unique(naive_exp_correlations_obs$day))
}


naive_exp_correlations_obs <- get_naive_exp_correlations(gene_freqs) %>%
  filter(total_compartment_seqs > min_compartment_size, total_mouse_naive_seqs > min_compartment_size) %>%
  mutate(group_controls_pooled = factor(group_controls_pooled, levels = group_controls_pooled_factor_levels)) %>%
  mutate(day = as.integer(as.character(day)))

naive_exp_correlations_neutral <- lapply(neutral_realizations %>% group_by(replicate) %>% group_split(),
                                         FUN = get_naive_exp_correlations)

naive_exp_correlations_neutral <- bind_rows(naive_exp_correlations_neutral, .id = 'replicate') %>%
  mutate(replicate = as.numeric(replicate)) %>%
  filter(total_compartment_seqs > min_compartment_size, total_mouse_naive_seqs > min_compartment_size) %>%
  mutate(group_controls_pooled = factor(group_controls_pooled, levels = group_controls_pooled_factor_levels)) %>%
  mutate(day = as.integer(as.character(day)))

naive_exp_pearson_corr_plot <- plot_naive_exp_correlations(naive_exp_correlations_obs, naive_exp_correlations_neutral,
                                                           method = 'pearson')

naive_exp_spearman_corr_plot <- plot_naive_exp_correlations(naive_exp_correlations_obs, naive_exp_correlations_neutral,
                                                            method = 'spearman')

# ====== FIND GENES DEVIATING CONSISTENTLY BETWEEN MICE =======
genes_by_n_mice_called <- 
  left_join(deviation_from_naive %>% filter(!is.na(deviation_from_naive)), 
            gene_freqs %>% select(mouse_id, cell_type, tissue, total_compartment_seqs, total_mouse_naive_seqs) %>%
              unique()) %>%
  filter(total_compartment_seqs >= min_compartment_size, total_mouse_naive_seqs >= min_compartment_size) %>%
  group_by(group_controls_pooled, v_gene) %>% 
  mutate(n_mice_gene_occurs_in_this_group = length(unique(mouse_id))) %>%
  group_by(group_controls_pooled, v_gene, n_mice_gene_occurs_in_this_group, cell_type, tissue, deviation_from_naive) %>%
  summarise(n_mice_with_call = length(unique(mouse_id))) %>% ungroup()

count_deviations_by_gene <- function(candidate_genes_tissue, candidate_genes_cell_type, candidate_genes_group){
  candidate_genes <- genes_by_n_mice_called %>% 
    filter(tissue == candidate_genes_tissue, group_controls_pooled == candidate_genes_group, cell_type == candidate_genes_cell_type,
           deviation_from_naive == 'positive') %>%
    arrange(desc(n_mice_with_call)) %>%
    #filter(n_mice_with_call >= 2) %>%
    pull(v_gene)
  
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
    rename(nonsignificant = neutral) %>%
    mutate(tissue = candidate_genes_tissue, cell_type = candidate_genes_cell_type, group_controls_pooled = candidate_genes_group) %>%
    select(tissue, cell_type, group_controls_pooled, everything())
  return(output_table)
}


flu_induced_pops <- expand_grid(tissue = 'LN', cell_type = c('GC','PC','mem'),
                                group = c('primary-8', 'primary-16', 'primary-24', 'secondary-40', 'secondary-56'))

deviations_by_allele <- bind_rows(mapply(count_deviations_by_gene, candidate_genes_tissue = flu_induced_pops$tissue,
       candidate_genes_cell_type = flu_induced_pops$cell_type, candidate_genes_group = flu_induced_pops$group, SIMPLIFY = F))
write_csv(deviations_by_allele, deviations_by_allele_results_path)

# Plots showing deviation from naive freq for major genes

plot_most_common_genes <- function(plot_cell_type, plot_tissue, plot_group){
  gene_freqs %>% 
    filter(cell_type == plot_cell_type, tissue == plot_tissue, group_controls_pooled == plot_group, v_gene_rank <= 20) %>%
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
top_genes_LN_PC_day16_plot <- plot_most_common_genes('PC','LN','primary-16') +
  theme(legend.position = c(0.67,0.25)) + background_grid()


top_genes_LN_GC_day8_plot <- plot_most_common_genes('GC','LN','primary-8') +
  theme(legend.position = c(0.8,0.3))

top_genes_LN_GC_day16_plot <- plot_most_common_genes('GC','LN','primary-16') +
  theme(legend.position = c(0.33,0.2))

top_genes_LN_GC_day24_plot <- plot_most_common_genes('GC','LN','primary-24') +
  theme(legend.position = c(0.33,0.2))

top_genes_LN_mem_day8_plot <- plot_most_common_genes('mem','LN','primary-8') +
  theme(legend.position = c(0.85,0.30))

top_genes_LN_mem_day16_plot <- plot_most_common_genes('mem','LN','primary-16') +
  theme(legend.position = c(0.85,0.30))

top_genes_LN_mem_day24_plot <- plot_most_common_genes('mem','LN','primary-24') +
  theme(legend.position = c(0.85,0.30))


# ======= MAKE DETAILED PLOTS TRACKING THE FATE OF FOCAL GENES =======
# (e.g., what do consistently overrepresented genes on day 8 LN PCs do at other time points / cell types)


# How genes consistently overrepresented on day 8 LN PCs (14-4, 1-69, 1-82) fare at later time points, etc.

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
      #scale_color_manual(values = c('green3','dodgerblue2'), name = 'Infection') +
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
      #scale_color_manual(values = c('green3','dodgerblue2'), name = 'Infection') +
      xlab("Days after primary infection") +
      label_controls_as_day_0 + 
      ylab("Experienced frequency \n (log10 + 1e-4)") +
      
      scale_x_continuous(breaks = sort(as.integer(unique(gene_freqs$day)))) + background_grid()
    
    clone_ranks_pl <- clone_ranks_data %>%
      ggplot(aes(x = day, y = rank_biggest_clone)) +
      geom_boxplot(aes(group = day), outlier.alpha = 0) +
      geom_point(aes(color = infection_status), position = position_jitter(width = 0.5, height = 0), size = 4, alpha = 0.5) +
      facet_grid(compartment_cell_type~v_gene, scales = 'free') +
      theme(legend.position = 'top',
            panel.border = element_rect(color = 'black')) +
      scale_y_log10() +
      #scale_color_manual(values = c('green3','dodgerblue2'), name = 'Infection') +
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

day8_LN_PC_overrep_genes <- plot_focal_genes(gene_freqs = gene_freqs,
                                             clone_freqs_by_tissue_and_cell_type = clone_freqs_by_tissue_and_cell_type,
                                             focal_genes = c('IGHV14-4*01', 'IGHV1-69*01', 'IGHV1-82*01'),
                                             min_compartment_size = min_compartment_size)



# ============== CLONE SIZE PLOTS ==========

# Fraction of sequences in the ten most common genes 
fraction_in_top_10_genes_plot <- gene_freqs %>% 
  filter(tissue == 'LN', group_controls_pooled != 'control',
         cell_type %in% c('GC','PC','mem')) %>%
  filter(v_gene_rank <= 10) %>% 
  group_by(mouse_id, day, infection_status, group_controls_pooled, cell_type, tissue, total_compartment_seqs) %>%
  summarise(fraction_seqs_in_top_genes = sum(vgene_seq_freq)) %>%
  ungroup() %>%
  mutate(day = as.integer(as.character(day))) %>%
  mutate(cell_type = case_when(
    cell_type == 'GC' ~ 'Lymph node GC cells',
    cell_type == 'PC' ~ 'Lymph node plasma cells',
    cell_type == 'mem' ~ 'Lymph node memory cells'
  )) %>%
  mutate(cell_type = factor(cell_type,
                            levels = c('Lymph node GC cells',
                                       'Lymph node plasma cells',
                                       'Lymph node memory cells'))) %>%
  ggplot(aes(x = day, y = fraction_seqs_in_top_genes, color = infection_status, group = day)) +
  geom_boxplot(outlier.alpha =  F, show.legend = F) +
  geom_point(aes(size = total_compartment_seqs), alpha = 0.8) +
  facet_grid(.~cell_type) +
  background_grid() +
  scale_color_manual(values = c('green3','dodgerblue2'), name = 'Infection') +
  #ggtitle() +
  theme(legend.position = 'top',
        legend.key.size = unit(50,'pt')) +
  xlab("Days after primary infection") +
  ylab("Fraction of sequences in the top 10 genes") +
  scale_x_continuous(breaks = sort(as.integer(unique(clone_freqs_by_tissue_and_cell_type$day)))) +
  scale_size_continuous(breaks = c(100,1000,10000,50000), name = ' Number of sequences')

plot(fraction_in_top_10_genes_plot)


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


# Ranked clones plot only implemented if using all-sequence frequencies
# (because clones are annotated with high frequency mutations)
# if(frequency_type == 'all_seqs'){
#   plot_top_clones <- function(plot_cell_type, plot_tissue, plot_group, plot_abs_size = F,
#                                    annotation = 'v_gene'){
#     
#     
#     if(plot_abs_size){
#       y_axis_var <- 'n_clone_seqs_in_compartment'
#       y_axis_label <- 'Number of sequences'
#     }else{
#       y_axis_var <- 'clone_freq_in_compartment'
#       y_axis_label <- 'Clone frequency'
#     }
#     
#     if(length(annotation) == 2){
#       plotting_data <- clone_freqs_by_tissue_and_cell_type %>%
#         mutate(across(c('v_gene','d_gene','j_gene'),
#                       function(x){str_remove(str_remove(x,'IGH'), '\\*[0-9]+')})) %>%
#         unite('annotation', annotation, sep = ' (') %>%
#         mutate(annotation = paste0(annotation,')')) %>%
#         mutate(annotation = str_remove(annotation, ' \\(\\)'))
#     }else{
#       plotting_data <- clone_freqs_by_tissue_and_cell_type %>%
#         mutate(across(c('v_gene','d_gene','j_gene'),
#                       function(x){str_remove(str_remove(x,'IGH'), '\\*[0-9]+')})) %>%
#         unite('annotation', annotation, sep = ' ; ') 
#     }
#     
#     pl <- plotting_data %>%
#       filter(total_seqs_in_compartment >= min_compartment_size, total_mouse_naive_seqs >= min_compartment_size) %>%
#       filter(compartment_cell_type == plot_cell_type, compartment_tissue == plot_tissue, group_controls_pooled == plot_group,
#              clone_rank_in_compartment <= 20) %>%
#       ungroup() %>%
#       ggplot(aes_string(x = 'clone_rank_in_compartment', y = y_axis_var, group = 'mouse_id')) +
#       geom_line() +
#       facet_wrap('mouse_id', scales = 'free') +
#       xlab('Clone rank (top 20 clones only)') +
#       ylab(y_axis_label) +
#       theme(legend.position = c(0.8,0.2)) +
#       scale_y_continuous(expand = expansion(mult = c(0.05, 0.3 * length(annotation)))) +
#       scale_x_continuous(expand = expansion(mult = c(0.05,0.2 * length(annotation)))) +
#       scale_color_discrete(name = 'V gene deviation from naive frequency',
#                            labels = c('negative','non-significant','positive'))
#     
#     if(length(annotation) == 1 & all(annotation == 'v_gene')){
#       pl <- pl + 
#         geom_point(aes(color = deviation_from_naive)) +
#         geom_text(aes(x = clone_rank_in_compartment + 0.3, color = deviation_from_naive, angle = 30, hjust = 0,
#                       label = annotation), show.legend = F, size = 3.5) 
#     }else{
#       pl <- pl + geom_point() +
#         geom_text(aes(x = clone_rank_in_compartment + 0.3, angle = 35, hjust = 0,
#                       label = annotation), show.legend = F, size = 3.5)
#     }
#     
#     return(pl)
#     
#     
#   }
# 
#   top_clones_LN_PC_day8 <- plot_top_clones('PC','LN', 'primary-8', plot_abs_size = T, annotation = c('v_gene', 'mutations_above_threshold'))
#   top_clones_LN_PC_day16 <- plot_top_clones('PC','LN', 'primary-16', plot_abs_size = T, annotation = c('v_gene', 'mutations_above_threshold'))
#   top_clones_LN_PC_day24 <- plot_top_clones('PC','LN', 'primary-24', plot_abs_size = T, annotation = c('v_gene', 'mutations_above_threshold'))
#   top_clones_LN_PC_day40 <- plot_top_clones('PC','LN', 'secondary-40', plot_abs_size = T, annotation = c('v_gene', 'mutations_above_threshold'))
#   top_clones_LN_PC_day56 <- plot_top_clones('PC','LN', 'secondary-56', plot_abs_size = T, annotation = c('v_gene', 'mutations_above_threshold'))
#   
#   top_clones_LN_PC_day56_CDR3_annotated <- plot_top_clones('PC','LN', 'secondary-56', plot_abs_size = F, annotation = c('clone_consensus_cdr3_partis')) + ylim(0,1.1)
#   
#   top_clones_LN_GC_day8 <- plot_top_clones('GC','LN', 'primary-8', plot_abs_size = T, annotation = c('v_gene', 'mutations_above_threshold'))
#   top_clones_LN_GC_day16 <- plot_top_clones('GC','LN', 'primary-16', plot_abs_size = T, annotation = c('v_gene', 'mutations_above_threshold'))
#   top_clones_LN_GC_day24 <- plot_top_clones('GC','LN', 'primary-24', plot_abs_size = T, annotation = c('v_gene', 'mutations_above_threshold'))
#   top_clones_LN_GC_day40 <- plot_top_clones('GC','LN', 'secondary-40', plot_abs_size = T, annotation = c('v_gene', 'mutations_above_threshold'))
#   top_clones_LN_GC_day56 <- plot_top_clones('GC','LN', 'secondary-56', plot_abs_size = T, annotation = c('v_gene', 'mutations_above_threshold'))
#   
#     
#   top_clones_LN_mem_day8 <- plot_top_clones('mem','LN', 'primary-8', plot_abs_size = T, annotation = c('v_gene', 'mutations_above_threshold'))
#   top_clones_LN_mem_day16 <- plot_top_clones('mem','LN', 'primary-16', plot_abs_size = T, annotation = c('v_gene', 'mutations_above_threshold'))
#   top_clones_LN_mem_day24 <- plot_top_clones('mem','LN', 'primary-24', plot_abs_size = T, annotation = c('v_gene', 'mutations_above_threshold'))
#   top_clones_LN_mem_day40 <- plot_top_clones('mem','LN', 'secondary-40', plot_abs_size = T, annotation = c('v_gene', 'mutations_above_threshold'))
#   top_clones_LN_mem_day56 <- plot_top_clones('mem','LN', 'secondary-56', plot_abs_size = T, annotation = c('v_gene', 'mutations_above_threshold'))
#   
# }


# ------------ PAIRWISE CORRELATIONS BETWEEN MICE ----------------

pairwise_naive_correlations_plot <- pairwise_correlations$freqs %>%
  filter(cell_type == 'naive', method == 'pearson') %>% # Results look very similar between pearson / spearman
  filter(total_compartment_seqs_i >= min_compartment_size, total_compartment_seqs_j >= min_compartment_size) %>%
  ggplot(aes(x = pair_type, y = cor_coef_freqs, color = pair_type)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position = position_jitter(width = 0.1, height = 0),
             alpha = point_alpha) +
  scale_y_continuous(limits = c(0,1)) +
  xlab('Type of pair') +
  ylab('Correlation in naive V gene frequencies\nbetween mouse pairs (excluding mice with < 100 sequences)') +
  theme(legend.position = 'none') +
  background_grid() 

pairwise_freq_correlations_plot <- pairwise_correlations$freqs %>%
  filter(cell_type %in% c('GC','PC','mem')) %>%
  mutate(cell_type = factor(cell_type, levels = c('GC','PC','mem'))) %>%
  filter(day_i == day_j, total_compartment_seqs_i >= min_compartment_size, total_compartment_seqs_j >= min_compartment_size,
         tissue == 'LN', pair_type %in% c('primary','secondary')) %>%
  mutate(day_i = as.integer(as.character(day_i))) %>%
  cell_type_facet_labeller() %>%
    ggplot(aes(x = day_i, y = cor_coef_freqs, group = day_i)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position = position_jitter(width =0.1, height = 0), alpha = 0.8, size = 4,
             aes(color = pair_type),
  ) +
  facet_grid(method ~ cell_type) +
  background_grid() +
  scale_color_manual(values = c('green3','dodgerblue2'), guide = 'none') +
  xlab('Days after primary infection') +
  ylab('Correlation in V gene frequencies\nbetween mouse pairs (excluding mice with < 100 seqs.)') +
  theme(legend.position = 'top') +
  scale_x_continuous(breaks = c(8,16,24,40,56)) +
  geom_hline(yintercept = 0, linetype = 2)

pairwise_freq_deviations_plot <- pairwise_correlations$freq_ratios %>%
  filter(cell_type %in% c('GC','PC','mem')) %>%
  mutate(cell_type = factor(cell_type, levels = c('GC','PC','mem'))) %>%
  filter(day_i == day_j, total_compartment_seqs_i >= min_compartment_size, total_compartment_seqs_j >= min_compartment_size,
         total_mouse_naive_seqs_i >= min_compartment_size, total_mouse_naive_seqs_j >= min_compartment_size,
         tissue == 'LN', pair_type %in% c('primary','secondary')) %>%
  mutate(day_i = as.integer(as.character(day_i))) %>%
  cell_type_facet_labeller() %>%
  ggplot(aes(x = day_i, y = cor_coef_freq_ratios, group = day_i)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position = position_jitter(width =0.1, height = 0), alpha = 0.8, size = 4,
             aes(color = pair_type),
  ) +
  facet_grid(method~cell_type) +
  background_grid() +
  scale_color_manual(values = c('green3','dodgerblue2'), guide = 'none') +
  xlab('Days after primary infection') +
  ylab('Correlation in V gene frequency deviations\nbetween mouse pairs (excluding mice with < 100 seqs.)') +
  theme(legend.position = 'top') +
  scale_x_continuous(breaks = c(8,16,24,40,56)) +
  geom_hline(yintercept = 0, linetype = 2) 


# ===== Linear models for correlation vs. time in LN plasma cells ======

# First, we subset values for the LN only, and for pairs of mice from the same time point and with a min. number of seqs.
subset_LN_values <- function(data){
  data %>% filter(cell_type %in% c('GC','PC','mem')) %>%
    filter(day_i == day_j, total_compartment_seqs_i >= min_compartment_size, total_compartment_seqs_j >= min_compartment_size,
           tissue == 'LN', pair_type %in% c('primary','secondary')) %>%
    mutate(time = as.integer(as.character(day_i)))
}

# Compare obs and randomized LMs
compare_obs_vs_random_LMs <- function(data_obs, data_rdm){
  if('cor_coef_freqs' %in% names(data_obs)){
    beta_obs <- data_obs %>%
      group_by(cell_type, method) %>%
      summarise(beta_obs = lm(cor_coef_freqs~time)$coefficients[2]) %>% ungroup()
    
    beta_rdm <- data_rdm %>%
      group_by(replicate, cell_type, method) %>%
      summarise(beta_rdm = lm(cor_coef_freqs~time)$coefficients[2]) %>%
      ungroup() 
      

    
  }else{
    stopifnot('cor_coef_freq_ratios' %in% names(data_obs))
    beta_obs <- data_obs %>%
      group_by(cell_type, method) %>%
      summarise(beta_obs = lm(cor_coef_freq_ratios~time)$coefficients[2]) %>% ungroup()
    
    beta_rdm <- data_rdm %>%
      group_by(replicate, cell_type, method) %>%
      summarise(beta_rdm = lm(cor_coef_freq_ratios~time)$coefficients[2]) %>%
      ungroup() 
    
  }
  
  output <- left_join(beta_rdm, beta_obs) %>%
    group_by(cell_type, method) %>%
    summarise(beta_obs = unique(beta_obs),
              n_replicates = n(),
              n_replicates_smaller_than_or_equal_to_obs = sum(beta_rdm <= beta_obs),
              fraction = n_replicates_smaller_than_or_equal_to_obs / n_replicates)
  return(output)
}


# For observed data
freq_correlations_LN_OBS <- pairwise_correlations$freqs %>% subset_LN_values()
deviation_correlations_LN_OBS <- pairwise_correlations$freq_ratios %>% subset_LN_values()

# And for data with randomized noncontrol mice
freq_correlations_LN_randomized <- pairwise_correlations_randomized_noncontrol_groups$freqs %>% subset_LN_values()
deviation_correlations_LN_randomized <- pairwise_correlations_randomized_noncontrol_groups$freq_ratios %>% subset_LN_values()


# Now we ask what % of randomizations are equal to or more negative than observed values

correlation_randomization_tests <- bind_rows(compare_obs_vs_random_LMs(data_obs = freq_correlations_LN_OBS, data_rdm = freq_correlations_LN_randomized) %>%
            mutate(correlation_type = 'frequencies'),
          compare_obs_vs_random_LMs(data_obs = freq_correlations_LN_OBS %>% filter(pair_type == 'primary'),
                          data_rdm = freq_correlations_LN_randomized %>% filter(pair_type == 'primary')) %>%
            mutate(correlation_type = 'frequencies (primary infection only)'),
          compare_obs_vs_random_LMs(data_obs = deviation_correlations_LN_OBS, data_rdm = deviation_correlations_LN_randomized) %>%
            mutate(correlation_type = 'frequency deviations'),
          compare_obs_vs_random_LMs(data_obs = deviation_correlations_LN_OBS %>% filter(pair_type == 'primary'),
                                    data_rdm = deviation_correlations_LN_randomized %>% filter(pair_type == 'primary')) %>%
            mutate(correlation_type = 'frequency deviations (primary infection_only')
          )

write_csv(correlation_randomization_tests, file = correlation_rdm_test_results)



#=======  Number of top N genes shared over time =======
top_genes_threshold <- 10

shared_top_genes_by_freq <- pairwise_gene_freqs %>%
  group_by(mouse_pair, pair_type, day_i, day_j, tissue, cell_type, total_compartment_seqs_i,
           total_compartment_seqs_j) %>%
  mutate(rank_i = rank(-vgene_seq_freq_i, ties.method = 'first'),
         rank_j = rank(-vgene_seq_freq_j, ties.method = 'first')) %>%
  filter(rank_i <= top_genes_threshold, rank_j <= top_genes_threshold,
         total_compartment_seqs_i >= min_compartment_size,
         total_compartment_seqs_j >= min_compartment_size,
         tissue == 'LN',
         pair_type %in% c('control', 'primary', 'secondary'),
         cell_type %in% c('GC','PC','mem'),
         day_i == day_j) %>%
  summarise(genes_in_top_set = n()) %>%
  ungroup() %>%
  cell_type_facet_labeller()

shared_top_genes_by_freq_pl <- shared_top_genes_by_freq %>%
  filter(pair_type != 'control') %>%
  mutate(day_i = as.integer(as.character(day_i))) %>%
  ggplot(aes(x = day_i, y = genes_in_top_set, group = day_i)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position = position_jitter(width =0.1, height = 0), alpha = 0.8, size = 4,
             aes(color = pair_type),
  ) +
  facet_wrap('cell_type') +
  background_grid() +
  scale_color_manual(values = c('green3','dodgerblue2'), guide = 'none') +
  xlab('Days after primary infection') +
  ylab(paste0('Number of shared top-',top_genes_threshold,' genes\n(excluding mice with < 100 seqs.)')) +
  theme(legend.position = 'top') +
  scale_x_continuous(breaks = c(8,16,24,40,56)) +
  scale_y_continuous(limits = c(0, NA))
  

###### HIERARCHICAL CLUSTERING ANALYSIS

# ----- Dendrograms and heatmaps based on correlation in V gene frequencies -------
# Have to be exported invididually because not easily integrated with ggplot...

pdf(file = paste0(exported_figure_objects_dir, 'freqs_heatmap_LN_PCs_spearman.pdf'), height = 7,
    width = 7.5)
get_vgene_freq_correlation_clustering(pairwise_correlations = pairwise_correlations,
                                      metric = 'freqs',
                                      cell_type = 'PC',
                                      tissue = 'LN',
                                      min_seqs = min_compartment_size,
                                      cor_method = 'spearman')
dev.off()

pdf(file = paste0(exported_figure_objects_dir, 'freqs_heatmap_LN_PCs_pearson.pdf'), height = 7,
    width = 7.5)
get_vgene_freq_correlation_clustering(pairwise_correlations = pairwise_correlations,
                                      metric = 'freqs',
                                      cell_type = 'PC',
                                      tissue = 'LN',
                                      min_seqs = min_compartment_size,
                                      cor_method = 'pearson')
dev.off()



# ----- Dendrograms and heatmaps based on correlation in V gene frequency deviations from naive repertoire -------

pdf(file = paste0(exported_figure_objects_dir, 'deviations_heatmap_LN_PCs_spearman.pdf'), height = 7,
    width = 7.5)
get_vgene_freq_correlation_clustering(pairwise_correlations = pairwise_correlations,
                                      metric = 'freq_ratios',
                                      cell_type = 'PC',
                                      tissue = 'LN',
                                      min_seqs = min_compartment_size,
                                      cor_method = 'spearman')
dev.off()

pdf(file = paste0(exported_figure_objects_dir, 'deviations_heatmap_LN_PCs_pearson.pdf'), height = 7,
    width = 7.5)
get_vgene_freq_correlation_clustering(pairwise_correlations = pairwise_correlations,
                                      metric = 'freq_ratios',
                                      cell_type = 'PC',
                                      tissue = 'LN',
                                      min_seqs = min_compartment_size,
                                      cor_method = 'pearson')
dev.off()



# PLOTS TO EXPORT
save(cumulative_naive_freqs_pl,
     file = paste0(exported_figure_objects_dir, 'cumulative_naive_freqs.RData'))

save(naive_exp_pearson_corr_plot,
     naive_exp_spearman_corr_plot,
     pairwise_naive_correlations_plot,
     pairwise_freq_correlations_plot,
     pairwise_freq_deviations_plot,
     shared_top_genes_by_freq_pl,
     file = paste0(exported_figure_objects_dir, 'correlation_plots.RData')
     )

save(top_genes_LN_PC_day8_plot,
     top_genes_LN_PC_day16_plot,
     top_genes_LN_GC_day8_plot,
     top_genes_LN_GC_day16_plot,
     top_genes_LN_mem_day24_plot,
     day8_LN_PC_overrep_genes,
     file = paste0(exported_figure_objects_dir, 'top_genes_plots.RData'))


# save(top_clones_LN_PC_day8,
#      top_clones_LN_PC_day16,
#      top_clones_LN_PC_day24,
#      top_clones_LN_PC_day40,
#      top_clones_LN_PC_day56,
#      top_clones_LN_GC_day8,
#      top_clones_LN_GC_day16,
#      top_clones_LN_GC_day24,
#      top_clones_LN_GC_day40,
#      top_clones_LN_GC_day56,
#      top_clones_LN_mem_day8,
#      top_clones_LN_mem_day16,
#      top_clones_LN_mem_day24,
#      top_clones_LN_mem_day40,
#      top_clones_LN_mem_day56,
#      file = paste0(exported_figure_objects_dir, 'top_clones_plots.RData'))

save(fraction_in_top_10_clones_plot,
     fraction_in_top_10_genes_plot,
     file = paste0(exported_figure_objects_dir, 'fraction_in_top_10_clones_plot.RData'))

