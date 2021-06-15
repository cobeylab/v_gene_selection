library(tidyverse)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(readr)
source('gene_frequency_functions.R')

processed_data_dir <- '../processed_data/'
#processed_data_dir <- '~/Desktop/v_gene_selection_files/'


seq_counts <- read_csv(paste0(processed_data_dir, 'seq_counts.csv'))
unique_seq_counts <- read_csv(paste0(processed_data_dir, 'unique_seq_counts.csv'))

seq_counts <- bind_rows(seq_counts %>% mutate(seq_type = 'all productive sequences') %>%
                          dplyr::rename(n_seqs = prod_seqs),
                        unique_seq_counts %>% mutate(seq_type = 'unique productive sequences') %>%
                          dplyr::rename(n_seqs = unique_prod_seqs)) %>%
  select(seq_type, everything())


seq_counts <- get_info_from_mouse_id(seq_counts)

# Aggregated sequence numbers across clones, to look at n seqs per mouse / tissue / cell type
seq_counts <- seq_counts %>%
  group_by(seq_type, mouse_id, day, infection_status, group, group_controls_pooled, tissue, cell_type) %>%
  summarise(n_seqs = sum(n_seqs)) %>%
  mutate(group_controls_pooled = factor(group_controls_pooled, levels = group_controls_pooled_factor_levels)) %>%
  mutate(tissue = factor(tissue, levels = c('LN','spleen','BM'))) %>%
  mutate(cell_type = factor(cell_type, levels = c('naive','nonnaive_IgD+B220+', 'GC','PC','mem')))


base_plotting_function <- function(seq_counts_tibble, seq_type){
  if(seq_type == 'all productive sequences'){
    ylabel = 'Number of productive sequences'
  }else{
    stopifnot(seq_type == 'unique productive sequences')
    ylabel = 'Number of unique productive sequences'
  }
  seq_counts_tibble %>%
    filter(seq_type == !!seq_type) %>%
    ggplot(aes(x = group_controls_pooled, y = n_seqs, color = infection_status)) +
    geom_boxplot(outlier.alpha = 0) +
    geom_point() +
    facet_grid(tissue~cell_type, scales = 'free') +
    theme(legend.position = 'top',
          axis.text.x = element_text(size = 8, angle = 40,
                                     vjust = 0.5)) +
    #scale_x_discrete(labels = function(x){str_replace(x, '-','\n')}) +
    background_grid() +
    xlab('Group') +
    ylab(ylabel) +
    scale_color_discrete(name = 'Infection status')
    
}


n_unique_seqs_pl <- base_plotting_function(seq_counts_tibble = seq_counts,
                       seq_type = 'unique productive sequences')
save_plot('../figures/n_unique_productive_seqs.pdf', 
          n_unique_seqs_pl,
          base_width = 12, base_height = 8)

