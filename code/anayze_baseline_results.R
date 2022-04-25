library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(cowplot)
library(stringr)
theme_set(theme_cowplot())
source('gene_frequency_functions.R')

min_compartment_size <- 100

RData_files <- list.files('../results/baseline_analysis/', pattern = 'RData', full.names = T)
# RData_files <- list.files('~/Desktop/v_gene_selection/results/baseline_analysis/', pattern = 'RData', full.names = T)
baseline_results_list <- list()

for(f in RData_files){
  mouse_id = str_extract(rev(str_split(f, '/')[[1]])[1], '[0-9]*-[0-9]*')
  load(f)
  baseline_results_list[[mouse_id]] <- convolved_baseline_results
}
rm(convolved_baseline_results)

# Read gene frequencies objects to filter results for compartments with 100+ seqs.
load('../results/precomputed_gene_freqs_all_seqs.RData')
#load('~/Desktop/v_gene_selection/results/precomputed_gene_freqs_all_seqs.RData')

# REVIEW THIS NOW THAT RESULTS ARE ALREADY CONVOLVED

# Combine results across mice
combined_stats <- baseline_results_list[[1]]@stats %>%
  mutate(mouse_id = names(baseline_results_list)[1])

for(i in 2:length(baseline_results_list)){
  next_object <- baseline_results_list[[i]]
  combined_stats <- bind_rows(combined_stats, next_object@stats %>%
                                mutate(mouse_id = names(baseline_results_list)[i])) 
}

combined_stats <- as_tibble(combined_stats) %>%
  mutate(tissue =  "LN")

# Add information of number of sequences per compartment
combined_stats <- left_join(combined_stats, gene_freqs %>% select(mouse_id, cell_type, tissue, total_compartment_seqs) %>% unique())



combined_stats <- get_info_from_mouse_id(combined_stats) %>%
  filter(cell_type %in% c('GC','PC','mem')) %>%
  cell_type_facet_labeller() %>%
  mutate(group_controls_pooled = factor(group_controls_pooled, levels = group_controls_pooled_factor_levels))


baseline_points_with_CIs <- combined_stats %>%
  filter(total_compartment_seqs >= 100) %>%
  ggplot(aes(x = group_controls_pooled, y = baseline_sigma, color = infection_status)) +
  geom_pointrange(aes(ymin = baseline_ci_lower, ymax = baseline_ci_upper), position = position_jitter(height = 0, width = 0.2),
                  alpha = 0.5) +
  facet_grid(cell_type~region) +
  geom_hline(yintercept = 0, linetype = 2) +
  theme(legend.position = 'top') +
  scale_color_discrete(name = 'Infection') +
  xlab('Group') +
  ylab('Strength of selection (sigma)\n(mice with 100+ seqs, lineages with 10+ seqs)') +
  scale_x_discrete(labels = function(x){str_replace(x,'-','\n')})

baseline_point_estimates_per_mouse <- combined_stats %>%
  filter(total_compartment_seqs >= 100) %>%
  ggplot(aes(x = group_controls_pooled, y = baseline_sigma, color = infection_status)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(aes(size = n_clones_in_compartment)) +
  facet_grid(cell_type~region) +
  geom_hline(yintercept = 0, linetype = 2) +
  theme(legend.position = 'top') +
  scale_color_discrete(name = 'Infection') +
  xlab('Group') +
  ylab('Strength of selection (sigma)\n(mice with 100+ seqs, lineages with 10+ seqs)') +
  scale_x_discrete(labels = function(x){str_replace(x,'-','\n')})

