# Checking distribution of mutations in naive seqs. versus null model without the sequence clustering step

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

source('gene_frequency_functions.R')

estimated_seq_error_rate <- 0.0018

# Pre-computes the probability of observing a range of mutations for seq. lengths observed in the data, given an estimated sequencing error rate
generate_mutation_null_model <- function(seq_cluster_stats, estimated_seq_error_rate){
  length_set <- unique(seq_cluster_stats$seq_length_partis)
  
  # Find the range of the number of nt mutations observed in the data
  n_mutations_range <- seq(min(seq_cluster_stats$n_mutations_partis_nt), max(seq_cluster_stats$n_mutations_partis_nt))
  
  null_model_mutations <- expand_grid(length_set, n_mutations_range) %>%
    dplyr::rename(length = length_set, n_mutations = n_mutations_range) %>%
    group_by(length) %>%
    mutate(null_prob = dbinom(x = n_mutations, size = length[1], prob = estimated_seq_error_rate))
  
}

annotated_seqs <- read_csv('../processed_data/annotated_seqs.csv') 

annotated_seqs <- annotated_seqs %>%
  mutate(across(c('clone_id_partis','partis_uniq_ref_seq','seq_id'), as.character))

annotated_seqs$specimen_cell_subset[annotated_seqs$specimen_cell_subset == 'naïve'] <- 'naive'

annotated_seqs <- annotated_seqs %>% filter(!is.na(n_mutations_partis_nt))
annotated_seqs <- get_info_from_mouse_id(annotated_seqs)
annotated_seqs <- annotated_seqs %>% dplyr::rename(tissue = specimen_tissue, cell_type = specimen_cell_subset)




# Distribution of the number of nucleotide mutations by mouse and tissue
distribution_nt_mutations_by_mouse_and_tissue <- annotated_seqs %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type, n_mutations_partis_nt) %>%
  summarise(n_seqs = n()) %>%
  ungroup() %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type) %>%
  mutate(compartment_seqs = sum(n_seqs),
         obs_fraction = n_seqs / compartment_seqs) %>%
  ungroup() %>%
  mutate(cell_type = factor(cell_type, levels = c('naive','GC','PC','mem')),
         tissue = factor(tissue, levels = c('LN','spleen','BM')))

# Distribution of sequence lengths by mouse and tissue
seq_length_distribution <- annotated_seqs %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type, seq_length_partis) %>%
  summarise(n_seqs = n()) %>%
  ungroup() %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type) %>%
  mutate(compartment_seqs = sum(n_seqs),
         obs_fraction = n_seqs / compartment_seqs) %>%
  ungroup() %>%
  mutate(cell_type = factor(cell_type, levels = c('naive','GC','PC','mem')),
         tissue = factor(tissue, levels = c('LN','spleen','BM')))

mean_lengths <- seq_length_distribution %>% group_by(mouse_id, tissue, cell_type) %>%
  summarise(mean_length = sum(obs_fraction*seq_length_partis)) %>% ungroup()




null_model_mutations <- generate_mutation_null_model(annotated_seqs, estimated_seq_error_rate)

# Calculate expected null distribution of mutations given the observed distribution of sequence lengths
null_distribution_given_obs_lengths <- left_join(seq_length_distribution %>%
                                                   select(mouse_id, tissue, cell_type, seq_length_partis, obs_fraction),
                                                 null_model_mutations %>% dplyr::rename(seq_length_partis = length)) %>%
  group_by(mouse_id, tissue, cell_type, n_mutations) %>%
  # For each number of mutations, calculate null probability as 
  # a weighted average by length given the obs. freq distribution of lengths
  summarise(null_prob = sum(obs_fraction*null_prob)) %>%
  ungroup()


distribution_nt_mutations_by_mouse_and_tissue <- left_join(distribution_nt_mutations_by_mouse_and_tissue,
                                                           null_distribution_given_obs_lengths %>%
                                                             dplyr::rename(n_mutations_partis_nt = n_mutations))

lapply(list('spleen','LN','BM'),
       FUN = function(tis, distribution_nt_mutations_by_mouse_and_tissue){
         pl <- distribution_nt_mutations_by_mouse_and_tissue %>%
           filter(cell_type == 'naive') %>%
           pivot_longer(cols = c('obs_fraction', 'null_prob')) %>%
           mutate(mouse_id = factor(mouse_id, levels = mouse_id_factor_levels)) %>%
           mutate(name = factor(name, levels = c('obs_fraction', 'null_prob'))) %>%
           filter(tissue == tis) %>%
           filter(n_mutations_partis_nt <= 15) %>%
           ggplot(aes(x = n_mutations_partis_nt, y = value, color = name, group = name, shape = name)) +
           geom_point() +
           geom_line() +
           facet_wrap('mouse_id', ncol = 4, scales = 'free') +
           xlab('Number of mutations from inferred germline sequence') +
           ylab('Fraction of sequences') +
           scale_color_discrete(labels = c('Observed fraction','Null expectation from estimated error rate'),
                                name = '') +
           scale_shape_manual(name = '', values  = c(19,1)) +
           theme(legend.position = 'top') +
           guides(shape = 'none')
         save_plot(paste0('../results/n_SHM_over_time/distribution_nt_mutations_naive_vs_null_',tis,'_UNCLUSTERED.pdf'),
                   pl,
                   base_width = 20, base_height = 30)
       },
       distribution_nt_mutations_by_mouse_and_tissue = distribution_nt_mutations_by_mouse_and_tissue
)
