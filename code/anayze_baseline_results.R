library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(cowplot)
library(stringr)
theme_set(theme_cowplot())
source('gene_frequency_functions.R')

RData_files <- list.files('../results/baseline_analysis/', pattern = 'RData', full.names = T)
# RData_files <- list.files('~/Desktop/v_gene_selection/results/baseline_analysis/', pattern = 'RData', full.names = T)
baseline_results_list <- list()

for(f in RData_files){
  mouse_id = str_extract(rev(str_split(f, '/')[[1]])[1], '[0-9]*-[0-9]*')
  load(f)
  baseline_results_list[[mouse_id]] <- convolved_baseline_results
}
rm(convolved_baseline_results)


# REVIEW THIS NOW THAT RESULTS ARE ALREADY CONVOLVED

# Combine results across mice
combined_stats <- baseline_results_list[[1]]@stats %>%
  mutate(mouse_id = names(baseline_results_list)[1])

for(i in 2:length(baseline_results_list)){
  next_object <- baseline_results_list[[i]]
  combined_stats <- bind_rows(combined_stats, next_object@stats %>%
                                mutate(mouse_id = names(baseline_results_list)[i])) 
}

combined_stats <- as_tibble(combined_stats %>% select(mouse_id, everything()))

stopifnot(("tissue" %in% names(combined_stats)) == F) # Just in case we choose to do baseline for multiple tissues later on

combined_stats <- get_info_from_mouse_id(combined_stats) %>%
  filter(cell_type %in% c('GC','PC','mem')) %>%
  mutate(tissue =  "LN") %>%
  cell_type_facet_labeller() %>%
  mutate(group_controls_pooled = factor(group_controls_pooled, levels = group_controls_pooled_factor_levels))

combined_stats %>%
  ggplot(aes(x = group_controls_pooled, y = baseline_sigma, color = infection_status)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(size = 2) +
  facet_grid(cell_type~region) +
  geom_hline(yintercept = 0, linetype = 2)

