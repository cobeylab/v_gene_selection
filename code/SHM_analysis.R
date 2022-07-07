library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

source('gene_frequency_functions.R')
source('plot_options.R')

# When looking at mutation frequencies within clones, restrict analysis to clones with at least this many seqs in the relevant compartment.
min_clone_size = 10
# Remove compartments with fewer than min_compartment_size seqs
min_compartment_size = 100

clone_info <- read_csv('../processed_data/clone_info.csv')
#clone_info <- read_csv('~/Desktop/v_gene_selection/processed_data/clone_info.csv')
clone_info <- get_info_from_mouse_id(clone_info) %>%
  mutate(group_controls_pooled = factor(group_controls_pooled, levels = group_controls_pooled_factor_levels))

load('../results/precomputed_gene_freqs_all_seqs.RData')
#load('~/Desktop/v_gene_selection/results/precomputed_gene_freqs_all_seqs.RData')


clone_freqs_by_tissue_and_cell_type <- clone_freqs_by_tissue_and_cell_type %>%
  filter(total_seqs_in_compartment >= min_compartment_size)

# clone_freqs_by_tissue_and_cell_type is annotated with mutations above freq. threshold.
# Use character vector of mutations above threshold to count them
clone_freqs_by_tissue_and_cell_type <- count_mutations_above_threshold(clone_freqs_by_tissue_and_cell_type)


# For each compartment, find fraction of clones with at least 1 mutation above threshold
fraction_clones_with_1plus_high_freq_muts <- get_fraction_of_clones_with_mutations_above_threshold(clone_freqs_by_tissue_and_cell_type,
                                                                                                   target_n_mutations = 1, 
                                                                                                   min_clone_size = min_clone_size)

# Plot this fraction for the lymph node
fraction_clones_with_high_freq_muts_LN_plot <- fraction_clones_with_1plus_high_freq_muts %>% 
  filter(compartment_cell_type %in% c('GC','mem','PC'), compartment_tissue == 'LN') %>%
  cell_type_facet_labeller() %>%
  set_controls_as_day_0() %>%
  ggplot(aes(x = day,
             y = fraction_clones_with_at_least_n_mutations,
             color = infection_status,
             group = day)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position = position_jitter(width = 0.1), aes(size = n_clones_denominator), alpha = 0.8) +
  background_grid() +
  theme(legend.position = 'top') +
  xlab('Days after primary infection') +
  ylab('Fraction of clones with at least 10 sequences\nthat have mutations at or above 50% frequency') +
  #scale_color_manual(values = c('green3','dodgerblue2')) +
  scale_size_continuous(name = 'Number of clones') +
  scale_y_continuous(limits = c(0,NA)) +
  facet_grid(.~compartment_cell_type)  +
  guides(color = 'none') +
  label_controls_as_day_0

# Then for all tissues:
fraction_clones_with_high_freq_muts_all_tissues_plot <- fraction_clones_with_1plus_high_freq_muts %>%
  filter(compartment_cell_type %in% c('GC','PC','mem')) %>%
  mutate(compartment_cell_type = case_when(
    compartment_cell_type == 'GC' ~ 'Germinal center cells',
    compartment_cell_type == 'PC' ~ 'Plasma cells',
    compartment_cell_type == 'mem' ~ 'Memory cells'
  ),
  compartment_tissue = case_when(
    compartment_tissue == 'LN' ~ 'Mediastinal LN',
    compartment_tissue == 'spleen' ~ 'Spleen',
    compartment_tissue == 'BM' ~ 'Bone marrow'
  )) %>%
  mutate(compartment_cell_type = factor(compartment_cell_type, levels = c('Germinal center cells',
                                                                          'Plasma cells',
                                                                          'Memory cells')),
         compartment_tissue = factor(compartment_tissue, levels = c('Mediastinal LN',
                                                                    'Spleen','Bone marrow'))) %>%
  ggplot(aes(x = group_controls_pooled,
             y = fraction_clones_with_at_least_n_mutations,
             color = infection_status)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position = position_jitter(width = 0.1), aes(size = n_clones_denominator), alpha = 0.8) +
  background_grid() +
  theme(legend.position = 'top', 
        axis.text.x = element_text(size = 8, angle = 40,
                                   vjust = 0.5)) +
  xlab('Group') +
  ylab('Fraction of clones with at least 10 sequences\nthat have mutations at or above 50% frequency') +
  scale_color_discrete(name = 'Infection') +
  scale_size_continuous(name = 'Number of clones') +
  scale_y_continuous(limits = c(0,NA)) +
  facet_grid(compartment_tissue~compartment_cell_type)

# Mean number of high-frequency mutation per clone (clones with 10+ seqs.)...
mean_n_mutations_above_threshold_by_compartment <- clone_freqs_by_tissue_and_cell_type %>%
  filter(n_clone_seqs_in_compartment >= min_clone_size) %>%
  group_by(mouse_id, day, infection_status, group, group_controls_pooled, compartment_tissue, compartment_cell_type) %>%
  dplyr::summarise(mean_n_mutations_above_threshold = mean(n_mutations_above_threshold),
                   n_clones_in_compartment = dplyr::n()) %>%
  ungroup()

#...in the lymph node
mean_n_mutations_LN_plot <- mean_n_mutations_above_threshold_by_compartment %>%
  filter(compartment_cell_type %in% c('GC','mem','PC'), compartment_tissue == 'LN') %>%
  cell_type_facet_labeller() %>%
  set_controls_as_day_0() %>%
  ggplot(aes(x = day,
             y = mean_n_mutations_above_threshold,
             color = infection_status,
             group = day)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position = position_jitter(width = 0.2), aes(size = n_clones_in_compartment), alpha = 0.9) +
  background_grid() +
  theme(legend.position = 'top') +
  xlab('Days after primary infection') +
  ylab('Average number of mutations at or\nabove 50% frequency in clones with at least 10 seqs.') +
  #scale_color_manual(values = c('green3','dodgerblue2')) +
  scale_y_continuous(limits = c(-0.05,NA)) +
  facet_grid(.~compartment_cell_type) +
  scale_size_continuous(name = 'Number of clones') +
  guides(color = 'none') +
  label_controls_as_day_0

# ...and in all tissues
mean_n_mutations_all_tissues_plot <- mean_n_mutations_above_threshold_by_compartment %>%
  filter(compartment_cell_type %in% c('GC','PC','mem')) %>%
  mutate(compartment_cell_type = case_when(
    compartment_cell_type == 'GC' ~ 'Germinal center cells',
    compartment_cell_type == 'PC' ~ 'Plasma cells',
    compartment_cell_type == 'mem' ~ 'Memory cells'
  ),
  compartment_tissue = case_when(
    compartment_tissue == 'LN' ~ 'Mediastinal LN',
    compartment_tissue == 'spleen' ~ 'Spleen',
    compartment_tissue == 'BM' ~ 'Bone marrow'
  )) %>%
  mutate(compartment_cell_type = factor(compartment_cell_type, levels = c('Germinal center cells',
                                                                          'Plasma cells',
                                                                          'Memory cells')),
         compartment_tissue = factor(compartment_tissue, levels = c('Mediastinal LN',
                                                                    'Spleen','Bone marrow'))) %>%
  ggplot(aes(x = group_controls_pooled,
             y = mean_n_mutations_above_threshold,
             color = infection_status)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position = position_jitter(width = 0.2), aes(size = n_clones_in_compartment), alpha = 0.9) +
  background_grid() +
  theme(legend.position = 'top', 
        axis.text.x = element_text(size = 8, angle = 40,
                                   vjust = 0.5)) +
  xlab('Group') +
  ylab('Average number of mutations at or\nabove 50% frequency in clones with at least 10 seqs.') +
  scale_color_discrete(name = 'Infection status') +
  scale_y_continuous(limits = c(-0.05,NA)) +
  facet_grid(compartment_tissue~compartment_cell_type) +
  scale_size_continuous(name = 'Number of clones')


# How often do LN clones of infected mice using the same V gene (from the same or different mice) share mutations?
shared_mutations_in_LN_clones <- count_mutations_shared_by_clone_pairs(clone_freqs_by_tissue_and_cell_type %>%
                                                                         filter(compartment_tissue == 'LN',
                                                                                compartment_cell_type %in% c('GC','PC','mem'),
                                                                                group_controls_pooled != 'control',
                                                                                n_clone_seqs_in_compartment >= min_clone_size))

fraction_LN_clones_sharing_mutations <- get_fraction_clones_with_shared_mutations(shared_mutations_in_LN_clones)


shared_mutations_in_LN_clones_pl <- fraction_LN_clones_sharing_mutations %>%
  set_controls_as_day_0() %>%
  cell_type_facet_labeller() %>%
  ggplot(aes(x = day, y = probability, fill = infection_status)) +
  geom_col() +
  geom_text(aes(y = probability + 0.02, label = total_n_clone_pairs)) +
  facet_wrap('compartment_cell_type', nrow = 1) +
  xlab('Days after primary infection') +
  ylab(paste0('Probability that two clones (', min_clone_size, '+ seqs.) sharing the same V allele\nhave high-frequency mutations in common')) +
  ylim(0,1) +
  theme(legend.position = 'top') +
  label_controls_as_day_0 +
  scale_fill_manual(values = c('green3','dodgerblue2'), name = 'Infection')


# Do the most abundant clones in lymph node populations tend to have more mutations at or above 50% frequency?
# Plasma cells...
clone_rank_vs_high_freq_muts_LN_PCs_plot <- clone_freqs_by_tissue_and_cell_type %>%
  filter(compartment_cell_type == 'PC', compartment_tissue == 'LN',
         group_controls_pooled != 'control',
         n_clone_seqs_in_compartment >= min_clone_size) %>%
  ggplot(aes(x = clone_rank_in_compartment, y = n_mutations_above_threshold)) +
  geom_point(aes(group = mouse_id, color = infection_status), alpha = 0.2) +
  facet_wrap('group_controls_pooled', nrow = 2, scales = 'free') +
  theme(legend.position = 'none') +
  geom_smooth(se = T, color = 'black', method = 'loess') +
  ylim(-0.05,20) +
  scale_color_manual(values = c('green3','dodgerblue2')) +
  xlab('Clone rank in lymph node plasma cells') +
  ylab('Number of amino acid mutations\nat or above 50% frequency') +
  ggtitle('(Y axis truncated at 20 mutations to improve visualization)')

# Germinal center cells
clone_rank_vs_high_freq_muts_LN_GCs_plot <- clone_freqs_by_tissue_and_cell_type %>%
  filter(compartment_cell_type == 'GC', compartment_tissue == 'LN',
         group_controls_pooled != 'control',
         n_clone_seqs_in_compartment >= min_clone_size) %>%
  ggplot(aes(x = clone_rank_in_compartment, y = n_mutations_above_threshold)) +
  geom_point(aes(group = mouse_id, color = infection_status), alpha = 0.2) +
  facet_wrap('group_controls_pooled', nrow = 2, scales = 'free') +
  theme(legend.position = 'none') +
  geom_smooth(se = T, color = 'black', method = 'loess') +
  ylim(-0.05,20) +
  scale_color_manual(values = c('green3','dodgerblue2')) +
  xlab('Clone rank in lymph node germinal center cells') +
  ylab('Number of amino acid mutations\nat or above 50% frequency') +
  ggtitle('(Y axis truncated at 20 mutations to improve visualization)')

# Export plots
save(fraction_clones_with_high_freq_muts_LN_plot,
     fraction_clones_with_high_freq_muts_all_tissues_plot,
     mean_n_mutations_LN_plot,
     mean_n_mutations_all_tissues_plot,
     clone_rank_vs_high_freq_muts_LN_PCs_plot,
     clone_rank_vs_high_freq_muts_LN_GCs_plot,
     shared_mutations_in_LN_clones_pl,
     file = '../figures/all_seqs_freqs/exported_ggplot_objects/high_frequency_mutations.RData')