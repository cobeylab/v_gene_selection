library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(cowplot)
library(stringr)
theme_set(theme_cowplot())
source('gene_frequency_functions.R')

min_compartment_size <- 100

# Read convolution analysis results

convolution_analysis_files <- list.files('../results/baseline_analysis/', pattern = 'convolution_baseline_analysis', full.names = T)
#convolution_analysis_files <- list.files('~/Desktop/v_gene_selection/results/baseline_analysis/', pattern = 'convolution_baseline_analysis', full.names = T)

convolution_analysis_list <- list()

for(f in convolution_analysis_files){
  mouse_id = str_extract(rev(str_split(f, '/')[[1]])[1], '[0-9]*-[0-9]*')
  load(f)
  convolution_analysis_list[[mouse_id]] <- convolved_baseline_results
}
rm(convolved_baseline_results) # Removes object from last imported RData file

# Read per-clone analysis results

per_clone_analysis_files <- list.files('../results/baseline_analysis/', pattern = 'baseline_analysis_per_clone', full.names = T)
#per_clone_analysis_files <- list.files('~/Desktop/v_gene_selection/results/baseline_analysis/', pattern = 'baseline_analysis_per_clone', full.names = T)

per_clone_analysis_list <- list()

for(f in per_clone_analysis_files){
  mouse_id = str_extract(rev(str_split(f, '/')[[1]])[1], '[0-9]*-[0-9]*')
  load(f)
  per_clone_analysis_list[[mouse_id]] <- baseline_results
}
rm(baseline_results) # Removes object from last imported RData file


# Read gene frequencies objects to filter results for compartments with 100+ seqs.
load('../results/precomputed_gene_freqs_all_seqs.RData')
#load('~/Desktop/v_gene_selection/results/precomputed_gene_freqs_all_seqs.RData')

# Combine results across mice
combined_convolution_analysis <- convolution_analysis_list[[1]]@stats %>%
  mutate(mouse_id = names(convolution_analysis_list)[1])

combined_per_clone_analysis <- per_clone_analysis_list[[1]]

for(i in 2:length(convolution_analysis_list)){
  print(i)
  next_convolution_object <- convolution_analysis_list[[i]]
  combined_convolution_analysis <- bind_rows(combined_convolution_analysis, next_convolution_object@stats %>%
                                mutate(mouse_id = names(convolution_analysis_list)[i])) 
  
  next_per_mouse_object <- per_clone_analysis_list[[i]]
  combined_per_clone_analysis@db <- bind_rows(combined_per_clone_analysis@db,
                                     next_per_mouse_object@db)
  combined_per_clone_analysis@numbOfSeqs <- rbind(combined_per_clone_analysis@numbOfSeqs,
                                         next_per_mouse_object@numbOfSeqs)
  combined_per_clone_analysis@binomK <- rbind(combined_per_clone_analysis@binomK,
                                     next_per_mouse_object@binomK)
  combined_per_clone_analysis@binomN <- rbind(combined_per_clone_analysis@binomN,
                                              next_per_mouse_object@binomN)
  combined_per_clone_analysis@binomP <- rbind(combined_per_clone_analysis@binomP,
                                              next_per_mouse_object@binomP)
  combined_per_clone_analysis@pdfs$fwr <- rbind(combined_per_clone_analysis@pdfs$fwr,
                                                next_per_mouse_object@pdfs$fwr)
  combined_per_clone_analysis@pdfs$cdr <- rbind(combined_per_clone_analysis@pdfs$cdr,
                                                next_per_mouse_object@pdfs$cdr)
  

}

# Get per-clone CIs for the selection strength sigma
per_clone_stats <- summarizeBaseline(combined_per_clone_analysis, returnType = 'df')
per_clone_stats <- left_join(as_tibble(per_clone_stats) %>%
                               dplyr::mutate(compartment_cell_type = cell_type,
                                             compartment_tissue = tissue),
                             clone_freqs_by_tissue_and_cell_type %>% mutate(clone_id = as.character(clone_id)))


# ======== PLots of convolution analysis ==========

combined_convolution_analysis <- as_tibble(combined_convolution_analysis) %>%
  mutate(tissue =  "LN")

# Add information of number of sequences per compartment
combined_convolution_analysis <- left_join(combined_convolution_analysis, gene_freqs %>% select(mouse_id, cell_type, tissue, total_compartment_seqs) %>% unique())



combined_convolution_analysis <- get_info_from_mouse_id(combined_convolution_analysis) %>%
  filter(cell_type %in% c('GC','PC','mem')) %>%
  cell_type_facet_labeller() %>%
  mutate(group_controls_pooled = factor(group_controls_pooled, levels = group_controls_pooled_factor_levels))


baseline_points_with_CIs <- combined_convolution_analysis %>%
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

baseline_point_estimates_per_mouse <- combined_convolution_analysis %>%
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

# ======== Plots of per-clone analysis

per_clone_stats %>%
  filter(cell_type %in% c('GC','PC','mem')) %>%
  filter(clone_rank_in_compartment <= 10, total_seqs_in_compartment >= min_compartment_size) %>%
  ggplot(aes(x = group_controls_pooled, y = baseline_sigma, color = infection_status)) +
  geom_pointrange(aes(ymin = baseline_ci_lower, ymax = baseline_ci_upper),
                  position = position_jitter(height = 0, width = 0.15)) +
  facet_grid(cell_type ~ region) +
  geom_hline(yintercept =  0, linetype = 2) +
  ylab('Strength of selection in top 10 clones\n(mice with 100+ seqs, lineages with 10+ seqs)')  +
  xlab('group') +
  scale_x_discrete(labels = function(x){str_replace(x,'-','\n')}) +
  theme(legend.position = 'top')

