library(tidyverse)
library(yaml)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

source('partis_output_functions.R')
source('gene_frequency_functions.R')

# Load and analyze partis results for the sequences from naive B cells from Vijay et al. 2021 (JEM) 
yaml_object <- read_yaml('../results/partis/seq_data_vijay2021JEM.yaml')

vijay2021JEM_annotated_seqs <- format_partis_info(yaml_object)

vijay2021JEM_annotated_seqs %>% group_by(productive_partis) %>%
  dplyr::count()

# Retain only sequences inferred to be productive.
vijay2021JEM_annotated_seqs <- vijay2021JEM_annotated_seqs %>% filter(productive_partis)

gene_freqs_seq_data_vijay2021JEM <- vijay2021JEM_annotated_seqs %>% filter(productive_partis) %>%
  group_by(v_segment_partis) %>%
  dplyr::summarise(n_seqs_vijay2021JEM = dplyr::n()) %>%
  ungroup() %>%
  mutate(freq_vijay2021JEM = n_seqs_vijay2021JEM/sum(n_seqs_vijay2021JEM)) %>%
  arrange(dplyr::desc(freq_vijay2021JEM )) %>%
  dplyr::rename(v_gene = v_segment_partis)

vgene_mutations_distribution_vijay2021JEM <- get_distribution_of_mutations(vijay2021JEM_annotated_seqs,
                                                                           'vgene_mutations_partis_nt', disable_grouping = T)

null_mutations_given_seq_lengths_vijay2021JEM <- get_null_mutation_distribution_given_length_distribution(vijay2021JEM_annotated_seqs,
                                                                                                          n_mutations_variable = 'vgene_mutations_partis_nt',
                                                                                                          seq_length_variable = 'sequenced_bases_in_vgene_region_partis',
                                                                                                          estimated_seq_error_rate = estimated_seq_error_rate,
                                                                                                          disable_grouping = T)
vgene_mutations_distribution_vijay2021JEM <- left_join(
  vgene_mutations_distribution_vijay2021JEM,
  null_mutations_given_seq_lengths_vijay2021JEM %>% dplyr::rename(vgene_mutations_partis_nt = n_mutations)
)

vgene_mutations_distribution_vijay2021JEM %>%
  pivot_longer(cols = c('obs_fraction','null_prob'),
               names_to = 'value_type') %>%
  mutate(obs_and_pred_number = value * compartment_seqs) %>%
  mutate(obs_and_pred_number = ifelse(value_type == 'obs_fraction',n_seqs, obs_and_pred_number)) %>%
  ggplot(aes(x = vgene_mutations_partis_nt, y = value)) +
  geom_point(aes(color = value_type)) +
  geom_line(aes(color = value_type)) +
  geom_text(aes(label = round(obs_and_pred_number,2)),
            position = position_dodge(width = 10)) +
  #scale_y_log10() +
  geom_hline(aes(yintercept = 1/compartment_seqs), linetype = 2) +
  xlab('Number of mutations in V gene region') +
  ylab('Frequency') +
  scale_color_discrete(name = '', 
                       labels = c('Null model','Observed')) +
  theme(legend.position = c(0.75,0.9)) +
  scale_x_continuous(limits = c(0,10),
                     breaks = 0:10) +
  ggtitle('A sequence with 179 mutations was ommitted to improve visualization')



# Now load pre-computed gene frequencies from the main (influenza infection) experiments
load('../results/precomputed_gene_freqs.RData')
# load('~/Desktop/v_gene_selection_files/precomputed_gene_freqs.RData')

naive_freqs_main_experiment <- left_join(gene_freqs %>% # Object with gene frequencies in the main experiment
                                          select(mouse_id, day, infection_status, group_controls_pooled, v_gene,
                                                 total_mouse_naive_seqs, naive_vgene_seq_freq) %>% unique(), 
                                        gene_freqs_seq_data_vijay2021JEM)


naive_freq_summary_stats <- naive_freqs_main_experiment %>% group_by(v_gene) %>%
  filter(total_mouse_naive_seqs >= 100) %>%
  dplyr::summarise(mean_naive_freq_main_experiment = mean(naive_vgene_seq_freq),
                   min_naive_freq_main_experiment = min(naive_vgene_seq_freq),
                   max_naive_freq_main_experiment = max(naive_vgene_seq_freq)) %>%
  arrange(dplyr::desc(mean_naive_freq_main_experiment))

naive_freq_summary_stats <- left_join(naive_freq_summary_stats, 
                                      gene_freqs_seq_data_vijay2021JEM) %>%
  mutate(v_gene = factor(v_gene, levels = v_gene))

naive_freq_summary_stats %>% 
  ggplot(aes(x = v_gene)) +
  geom_pointrange(aes(y = mean_naive_freq_main_experiment, ymin = min_naive_freq_main_experiment,
                      ymax = max_naive_freq_main_experiment)) +
  geom_point(aes(y = freq_vijay2021JEM), color = 'blue', size = 2)  +
  ylab('Naive frequency') +
  theme(axis.text.x = element_blank()) +
  ggtitle('Note: I excluded mice with fewer than 100 naive sequences') +
  xlab('V genes ranked by their average naive frequency in our experiment') 
  
naive_freqs_main_experiment %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled) %>%
  dplyr::summarise(cor_coef = cor.test(naive_vgene_seq_freq, freq_vijay2021JEM, method = 'spearman')$estimate) %>%
  ggplot(aes(x = cor_coef)) +
  geom_histogram(bins = 12, fill = 'white', color = 'black') +
  xlab('Spearman correlation in naive frequencies with the data from Vijay et al. 2021') +
  ylab('Number of mice in our experiments')


naive_freqs_main_experiment %>%
  ggplot(aes(x = naive_vgene_seq_freq, y = freq_vijay2021JEM)) +
  geom_point() +
  facet_wrap('mouse_id', scales = 'free')
  

