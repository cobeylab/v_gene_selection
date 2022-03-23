library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggdendro)
library(cowplot)
theme_set(theme_cowplot())
library(scales)
library(viridis)

source('gene_frequency_functions.R')
source('plot_options.R')

min_compartment_size = 100

# These analyses are based on all productive sequences

clone_info <- read_csv('../processed_data/clone_info.csv')
# clone_info <- read_csv('~/Desktop/v_gene_selection/processed_data/clone_info.csv')

annotated_seqs <- read_csv('~/Desktop/v_gene_selection/processed_data/annotated_seqs.csv') %>%
  mutate(seq_id = as.character(seq_id),
         clone_id = as.character(clone_id))

annotated_seqs <- annotated_seqs %>%
  filter(productive_partis)

annotated_seqs <- get_info_from_mouse_id(annotated_seqs) %>%
  mutate(group_controls_pooled = factor(group_controls_pooled, levels = group_controls_pooled_factor_levels))

# Plot CDR3 length distributions for each mouse for each cell type in the LN.
# For naive B cells, using distribution across all tissues

CDR3_lengths_naive <- annotated_seqs %>%
  filter(group_controls_pooled != 'control',
         cell_type == 'naive') %>%
  select(mouse_id, day, infection_status, group, group_controls_pooled, seq_id, tissue, cell_type, cdr3_seq_partis) %>%
  group_by(mouse_id, day, infection_status, group, group_controls_pooled) %>%
  mutate(total_compartment_seqs = n(),
         tissue = 'all') %>%
  ungroup() %>%
  mutate(cdr3_length = nchar(cdr3_seq_partis))

CDR3_lengths_LN_experienced <- annotated_seqs %>%
  filter(group_controls_pooled != 'control',
         cell_type %in% c('GC','PC','mem'),
         tissue == 'LN') %>%
  select(mouse_id, day, infection_status, group, group_controls_pooled, seq_id, tissue, cell_type, cdr3_seq_partis) %>%
  group_by(mouse_id, day, infection_status, group, group_controls_pooled, tissue, cell_type) %>%
  mutate(total_compartment_seqs = n()) %>%
  ungroup() %>%
  mutate(cdr3_length = nchar(cdr3_seq_partis)) 


CDR3_lengths <- bind_rows(CDR3_lengths_naive,
                          CDR3_lengths_LN_experienced) %>%
  mutate(cell_type = case_when(
    cell_type == 'naive' ~ 'Naive cells (all tissues)',
    cell_type == 'GC' ~ 'Lymph node GC cells',
    cell_type == 'PC' ~ 'Lymph node PC cells',
    cell_type == 'mem' ~ 'Lymph node memory cells'
  )) %>%
  mutate(cell_type = factor(cell_type,
                            levels = c('Naive cells (all tissues)', 'Lymph node GC cells',
                                       'Lymph node PC cells', 'Lymph node memory cells')))


CDR3_length_distribution_pl <- CDR3_lengths %>%
  ggplot(aes(x = group_controls_pooled, y = cdr3_length, group = mouse_id, color = group_controls_pooled)) +
  geom_boxplot() +
  facet_wrap('cell_type', nrow = 1) +
  theme(legend.position = 'none',
        axis.text.x = element_text(size = 10)) +
  xlab('Group') +
  ylab('CDR3 amino acid length') +
  scale_x_discrete(labels = function(x){str_replace(x, '-','\n')}) +
  scale_y_continuous(breaks = seq(0,25,5), limits = c(0, NA))
  
save_plot('../figures/all_seqs_freqs/CDR3_length_distribution.pdf',
          CDR3_length_distribution_pl,
          base_height = 5, base_width = 16)

  




# Old stuff: review/delete
# unique(annotated_seqs$cdr3_seq_partis)
# 
# 
# 
# CDR3_freqs <- annotated_seqs %>% filter(productive_partis) %>%
#   group_by(mouse_id, tissue, cell_type, cdr3_seq_partis, cdr3_mutations_partis_aa) %>%
#   summarise(n_seqs = dplyr::n()) %>%
#   group_by(mouse_id, tissue, cell_type) %>%
#   mutate(cdr3_seq_freq = n_seqs / sum(n_seqs)) %>%
#   ungroup()
# 
# CDR3_freqs <- get_info_from_mouse_id(CDR3_freqs)
# 
# CDR3_freqs %>% 
#   group_by(mouse_id, tissue, cell_type) %>%
#   mutate(rank_cdr3_seq_freq = rank(-cdr3_seq_freq, ties.method = 'average')) %>%
#   arrange(mouse_id, rank_cdr3_seq_freq) %>%
#   filter(rank_cdr3_seq_freq <=5) %>%
#   filter(group_controls_pooled == 'primary-8', tissue == 'LN', cell_type == 'PC')
