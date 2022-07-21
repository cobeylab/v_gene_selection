library(readr)
source('gene_frequency_functions.R')

theme_set(theme_cowplot())

# This analysis is based on counts of all productive sequences (as opposed to unique sequences only)
frequency_type <- 'all_seqs'

# When looking at clone's tissue or cell type composition, consider only clones with at least this many seqs.
min_clone_size = 10
min_compartment_size = 100 # When looking at fraction seqs in top 10 clones, disregard compartments with fewer seqs than this.

results_directory <- '../results/'

processed_data_directory <- '../processed_data/'

figure_directory <- paste0('../figures/', frequency_type, '_freqs/')

precomputed_freqs_file <- paste0('precomputed_gene_freqs_', frequency_type, '.RData')

# Load precomputed gene frequencies, neutral realizations, pairwise correlations 
load(paste0(results_directory, precomputed_freqs_file))

if(frequency_type == 'all_seqs'){
  seq_counts <- read_csv('../processed_data/seq_counts.csv')
}else{
  stopifnot(frequency_type == 'unique_seqs')
  read_csv('../processed_data/unique_seq_counts.csv')
}

# Basic info for each clone (incl. tissue and cell type composition)
clone_info <- read_csv(paste0(processed_data_directory,'clone_info.csv'))

exported_figure_objects_dir <- paste0(figure_directory,'exported_ggplot_objects/')

# Set order of tissues and cell types for plotting
clone_freqs_by_tissue <- clone_freqs_by_tissue %>%
  mutate(compartment_tissue = factor(compartment_tissue, levels = c('LN','spleen','BM')))
clone_freqs_by_tissue_and_cell_type <- clone_freqs_by_tissue_and_cell_type %>%
  mutate(compartment_tissue = factor(compartment_tissue, levels = c('LN','spleen','BM')),
         compartment_cell_type = factor(compartment_cell_type, 
                                        levels = c('naive','nonnaive_IgD+B220+',
                                                   'GC','PC','mem')))

# Exclude clones with no productive sequences 
clone_info <- clone_info %>%
  filter(!is.na(total_clone_prod_seqs))

# Add mouse info (day, group, etc.)
clone_info <- get_info_from_mouse_id(clone_info) %>%
  mutate(group_controls_pooled = factor(group_controls_pooled, levels = group_controls_pooled_factor_levels)) %>%
  mutate(mouse_id = factor(mouse_id, levels = mouse_id_factor_levels))

# What fraction of clones has 90% or more reads in a single tissue?
fraction_clones_dominated_by_single_tissue_plot <- clone_info %>%
  filter(total_clone_prod_seqs >= min_clone_size) %>%
  mutate(size_filter = (biggest_tissue_fraction_prod_seqs >= 0.9)) %>%
  group_by(mouse_id, group_controls_pooled, infection_status, size_filter) %>%
  dplyr::count() %>%
  group_by(mouse_id, group_controls_pooled, infection_status) %>%
  mutate(clones_above_min_size = sum(n),
         fraction_dominated_by_single_tissue = n/clones_above_min_size) %>%
  filter(size_filter == T) %>%
  ggplot(aes(x = group_controls_pooled, y = fraction_dominated_by_single_tissue,
             color = infection_status)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(aes(size = clones_above_min_size )) +
  xlab('Group') +
  ylim(0,1.05) + 
  ylab('Fraction of clones with 90%\nor more sequences in a single tissue') + 
  background_grid() +
  groups_color_scale(name = 'Infection') +
  scale_size_continuous(name = paste0('Number of clones\n(', min_clone_size, '+ seqs.)')) +
  theme(legend.position = 'top')

# What fraction of clones has 90% or more reads in a single cell type?
fraction_clones_dominated_by_single_cell_type_plot <- clone_info %>%
  filter(total_clone_prod_seqs >= min_clone_size) %>%
  mutate(size_filter = (biggest_cell_type_fraction_prod_seqs >= 0.9)) %>%
  group_by(mouse_id, group_controls_pooled, infection_status, size_filter) %>%
  dplyr::count() %>%
  group_by(mouse_id, group_controls_pooled, infection_status) %>%
  mutate(clones_above_min_size = sum(n),
         fraction_dominated_by_single_type = n/clones_above_min_size) %>%
  filter(size_filter == T) %>%
  ggplot(aes(x = group_controls_pooled, y = fraction_dominated_by_single_type,
             color = infection_status)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(aes(size = clones_above_min_size )) +
  xlab('Group') +
  ylim(0,1.05) + 
  ylab('Fraction of clones with 90%\nor more sequences from a single cell type') + 
  background_grid() +
  groups_color_scale(name = 'Infection') +
  scale_size_continuous(name = paste0('Number of clones\n(', min_clone_size, '+ seqs.)')) +
  theme(legend.position = 'top')

# For each cell type in the lymph node, fraction of sequences in the largest 10 clones
fraction_in_top_10_clones_plot <- clone_freqs_by_tissue_and_cell_type %>% 
  filter(compartment_tissue == 'LN', total_seqs_in_compartment >= min_compartment_size) %>%
  # Select top 10 clones
  filter(clone_rank_in_compartment <= 10) %>% 
  filter(compartment_cell_type %in% c('GC','PC','mem')) %>%
  mutate(compartment_cell_type = factor(compartment_cell_type, levels = c('GC','PC','mem'))) %>%
  # Total number of sequences in top-10 clones of each compartment
  group_by(mouse_id, day, infection_status, group_controls_pooled, compartment_cell_type, compartment_tissue,
           total_seqs_in_compartment) %>%
  summarise(seqs_in_top_clones = sum(n_clone_seqs_in_compartment)) %>%
  # Divide by total number of sequences in each compartment
  mutate(fraction_seqs_in_top_clones = seqs_in_top_clones / total_seqs_in_compartment) %>%
  ungroup() %>%
  set_controls_as_day_0() %>%
  cell_type_facet_labeller() %>%
  ggplot(aes(x = day, y = fraction_seqs_in_top_clones, color = infection_status, group = day)) +
  geom_boxplot(outlier.alpha =  F, show.legend = F) +
  geom_point(aes(size = total_seqs_in_compartment), alpha = 0.8) +
  facet_grid(.~compartment_cell_type) +
  background_grid() +
  guides(color = 'none') +
  groups_color_scale(name = 'Infection') +
  theme(legend.position = 'top') +
  xlab("Days after primary infection") +
  ylab("Fraction of sequences in the 10 largest clones") +
  label_controls_as_day_0 +
  scale_size_continuous(breaks = c(100,1000,10000,50000), name = ' Number of sequences')

save(fraction_in_top_10_clones_plot,
     file = paste0(exported_figure_objects_dir, 'fraction_in_top_10_clones_plot.RData'))

save(fraction_clones_dominated_by_single_cell_type_plot,
     fraction_clones_dominated_by_single_tissue_plot,
     file = paste0(exported_figure_objects_dir,'clonal_composition.RData'))


