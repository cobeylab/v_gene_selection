library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(cowplot)
library(stringr)
theme_set(theme_cowplot())

# Directory containing simulation results 
args <- commandArgs(trailingOnly = T) 
results_directory <- args[1] # results_directory <- '../results/simulations/neutral_scenario_1/'

# Read main objects from simulation
load(paste0(results_directory, basename(results_directory),'_results.RData'))



# Read object with results from observed data
load('../results/precomputed_gene_freqs_all_seqs.RData')
#load('~/Desktop/results/precomputed_gene_freqs_all_seqs.RData')

allele_counts_non_neutral_scenario_1 %>% filter(t == max(t)) %>%
  group_by(individual, t) %>%
  filter(experienced_freq == max(experienced_freq))

gene_freqs %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type, total_compartment_seqs) %>%
  summarise(allele_diversity = 1 - sum(vgene_seq_freq^2)) %>%
  ungroup() %>%
  filter(tissue == 'LN', cell_type == 'PC') %>%
  ggplot(aes(x = group_controls_pooled, y = allele_diversity, color = infection_status)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(aes(size = total_compartment_seqs),
             position = position_jitter(width = 0.1, height = 0))

gene_freqs %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, tissue,
           cell_type, total_compartment_seqs) %>%
  filter(naive_vgene_seq_freq > 0 & vgene_seq_freq > 0) %>%
  summarise(genes_present_in_the_response = n()) %>%
  filter(tissue == 'LN', cell_type == 'PC') %>%
  ggplot(aes(x = group_controls_pooled, y = genes_present_in_the_response, color = infection_status)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(aes(size = total_compartment_seqs),
             position = position_jitter(width = 0.1, height = 0))

gene_freqs %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, tissue,
           cell_type, total_compartment_seqs) %>%
  mutate(gene_rank = rank(-vgene_seq_freq, ties.method = 'random')) %>%
  arrange(mouse_id, day, infection_status, group_controls_pooled, tissue,
          cell_type, total_compartment_seqs, gene_rank) %>%
  mutate(cumulative_freq = cumsum(vgene_seq_freq)) %>%
  ungroup() %>%
  filter(cell_type == 'PC', tissue == 'LN', group_controls_pooled != 'control') %>%
  ggplot(aes(x = gene_rank, cumulative_freq)) +
  geom_line(aes(group = mouse_id)) +
  facet_wrap('group_controls_pooled')

allele_counts_non_neutral_scenario_1 %>%
  group_by(t, individual) %>%
  mutate(gene_rank = rank(-experienced_freq, ties.method = 'random')) %>%
  arrange(t, individual, gene_rank) %>%
  mutate(cumulative_freq = cumsum(experienced_freq)) %>%
  ungroup() %>%
  filter(t %in% c(1, 100, 200, 300, 400, 500)) %>%
  ggplot(aes(x = gene_rank, cumulative_freq)) +
  geom_line(aes(group = individual)) +
  facet_wrap('t')



