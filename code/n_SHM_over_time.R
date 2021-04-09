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


# Read sequence statistics aggregaed by sequence cluster
seq_cluster_stats <- read_csv('../processed_data/seq_cluster_stats.csv')

seq_cluster_stats <- get_info_from_mouse_id(seq_cluster_stats) %>%
  mutate(group_controls_pooled = factor(group_controls_pooled, levels = group_controls_pooled_factor_levels))

seq_cluster_stats <- seq_cluster_stats %>%
  mutate(across(matches('mutations'), round)) %>%
  mutate(seq_length_partis = round(seq_length_partis))


# Distribution of the number of nucleotide mutations by mouse andtissue
distribution_nt_mutations_by_mouse_and_tissue <- seq_cluster_stats %>%
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
seq_length_distribution <- seq_cluster_stats %>%
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
#mean(mean_lengths$mean_length)


# Plot with fraction of unmutated sequences by mouse / tissue / cell type
distribution_nt_mutations_by_mouse_and_tissue %>%
  filter(n_mutations_partis_nt == 0, compartment_seqs >= 100) %>%
  ggplot(aes(x = group_controls_pooled, y = obs_fraction, color = infection_status)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point() +
  facet_grid(tissue~cell_type) +
  xlab('Group') +
  ylab('Fraction of sequences unmutated') +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 20, vjust = 0.5)) +
  background_grid()

# Plot with fraction of sequences with more than 1 NT mutation by mouse / tissue / cell type
distribution_nt_mutations_by_mouse_and_tissue %>%
  filter(compartment_seqs >= 100) %>%
  filter(n_mutations_partis_nt <= 1) %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type) %>%
  summarise(obs_fraction = sum(obs_fraction)) %>%
  ggplot(aes(x = group_controls_pooled, y = obs_fraction, color = infection_status)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point() +
  facet_grid(tissue~cell_type) +
  xlab('Group') +
  ylab('Fraction of sequences with less than two mutations') +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 20, vjust = 0.5)) +
  background_grid()

# Plot with detailed distributions for each mouse
lapply(list('spleen','LN','BM'),
       FUN = function(tis, distribution_nt_mutations_by_mouse_and_tissue){
         pl <- distribution_nt_mutations_by_mouse_and_tissue %>%
           filter(compartment_seqs >= 100) %>%
           mutate(mouse_id = factor(mouse_id, levels = mouse_id_factor_levels)) %>%
           filter(tissue == tis) %>%
           filter(n_mutations_partis_nt <= 15) %>%
           ggplot(aes(x = n_mutations_partis_nt, y = obs_fraction, group = mouse_id)) +
           geom_point() +
           geom_line() +
           scale_x_continuous(breaks = c(0,5,10,15)) +
           facet_grid(group_controls_pooled~cell_type) +
           #facet_wrap('cell_type', ncol = 4, scales = 'free') +
           xlab('Number of mutations from inferred germline sequence') +
           ylab('Fraction of sequences') +
           background_grid()
         save_plot(paste0('../results/n_SHM_over_time/distribution_nt_mutations_by_mouse_',tis,'.pdf'),
                   pl,
                   base_width = 12, base_height = 10)
       },
       distribution_nt_mutations_by_mouse_and_tissue = distribution_nt_mutations_by_mouse_and_tissue
)



distribution_nt_mutations_by_group_and_tissue %>%
  ggplot(aes(x = group_controls_pooled, y = n_mutations_partis_nt, size = obs_fraction, alpha = obs_fraction)) +
  facet_grid(tissue~cell_type, scales = 'free') +
  geom_point()

seq_cluster_stats %>%
  mutate(tissue = factor(tissue, levels = c('LN','spleen','BM')),
         cell_type = factor(cell_type, levels = c('naive','GC','PC','mem'))) %>%
  ggplot(aes(x = group_controls_pooled, y = n_mutations_partis_nt, color = infection_status)) +
  geom_violin() +
  facet_grid(tissue~cell_type, scales = 'free') +
  scale_y_log10() +
  theme(legend.position = 'none') +
  xlab('Group') +
  ylab('Number of mutations') +
  background_grid() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 20, vjust = 0.5)) +
  background_grid()

distribution_nt_mutations_by_group_and_tissue %>%
  filter(n_mutations_partis <=10) %>%
  ggplot(aes(x = group_controls_pooled, y = obs_fraction, fill = n_mutations_partis)) +
  geom_col() +
  facet_grid(tissue~cell_type)




mean_n_mutations <- seq_level_data %>% group_by(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type) %>%
  summarise(mean_n_mutations = mean(n_mutations_partis))

mean_n_mutations %>% ggplot(aes(x = group_controls_pooled, y = mean_n_mutations, color = infection_status,)) +
  geom_point() +
  geom_boxplot() +
  facet_grid(tissue ~ cell_type, scales = 'free') +
  background_grid()



#-------  Distribution of n. of apparent mutations in naive sequences

# DEVELOP A BETTER NULL MODEL FOR THIS

# First across all mice and tissues
distribution_nt_mutations_global %>%
  filter(cell_type == 'naive') %>%
  pivot_longer(cols = c('obs_fraction', 'expected_fraction_from_error_rate')) %>%
  mutate(name = factor(name, levels = c('obs_fraction', 'expected_fraction_from_error_rate'))) %>%
  ggplot(aes(x = n_mutations_partis_nt, y = value, color = name, group = name)) +
  geom_point() +
  geom_line() +
  xlim(0, 20) +
  xlab('Number of mutations from inferred germline sequence') +
  ylab('Frequency') +
  scale_color_discrete(labels = c('Observed fraction','Binomial expectation from estimated error rate'),
                       name = '') +
  theme(legend.position = 'top')

distribution_nt_mutations_by_group_and_tissue %>%
  filter(cell_type == 'naive') %>%
  pivot_longer(cols = c('obs_fraction', 'expected_fraction_from_error_rate')) %>%
  mutate(name = factor(name, levels = c('obs_fraction', 'expected_fraction_from_error_rate'))) %>%
  ggplot(aes(x = n_mutations_partis_nt, y = value, color = name, group = name)) +
  geom_point() +
  geom_line() +
  xlim(0, 20) +
  facet_grid(group_controls_pooled~tissue) +
  xlab('Number of mutations from inferred germline sequence') +
  ylab('Frequency') +
  scale_color_discrete(labels = c('Observed fraction','Binomial expectation from estimated error rate'),
                       name = '') +
  theme(legend.position = 'top')

# =========== Null model for the number of mutations in naive sequences (which in theory should not be mutated) ===========

# Pre-computed null model distributions
null_model_mutations <- generate_mutation_null_model(seq_cluster_stats, estimated_seq_error_rate)

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
         save_plot(paste0('../results/n_SHM_over_time/distribution_nt_mutations_naive_vs_null_',tis,'.pdf'),
                   pl,
                   base_width = 20, base_height = 30)
       },
       distribution_nt_mutations_by_mouse_and_tissue = distribution_nt_mutations_by_mouse_and_tissue
)



