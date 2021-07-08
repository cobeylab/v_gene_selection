library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggdendro)
library(cowplot)
library(scales)
library(viridis)

source('gene_frequency_functions.R')

clone_info <- read_csv('../processed_data/clone_info.csv')
# clone_info <- read_csv('~/Desktop/v_gene_selection_files/clone_info.csv')

annotated_seqs <- read_csv('~/Desktop/v_gene_selection_files/annotated_seqs.csv')

unique(annotated_seqs$cdr3_seq_partis)



CDR3_freqs <- annotated_seqs %>% filter(productive_partis) %>%
  group_by(mouse_id, tissue, cell_type, cdr3_seq_partis, cdr3_mutations_partis_aa) %>%
  summarise(n_seqs = dplyr::n()) %>%
  group_by(mouse_id, tissue, cell_type) %>%
  mutate(cdr3_seq_freq = n_seqs / sum(n_seqs)) %>%
  ungroup()

CDR3_freqs <- get_info_from_mouse_id(CDR3_freqs)

CDR3_freqs %>% 
  group_by(mouse_id, tissue, cell_type) %>%
  mutate(rank_cdr3_seq_freq = rank(-cdr3_seq_freq, ties.method = 'average')) %>%
  arrange(mouse_id, rank_cdr3_seq_freq) %>%
  filter(rank_cdr3_seq_freq <=5) %>%
  filter(group_controls_pooled == 'primary-8', tissue == 'LN', cell_type == 'PC')
