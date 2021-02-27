library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

source('gene_frequency_functions.R')

seq_level_data <- read_csv('../processed_data/seq_level_files/seq_level_data.csv')
# seq_level_data <- read_csv('~/Desktop/seq_level_data.csv')

estimated_seq_error_rate <- 0.0018

seq_level_data$specimen_cell_subset[seq_level_data$specimen_cell_subset == 'naïve'] <- 'naive'

mouse_info <- get_info_from_mouse_id(seq_level_data %>% select(mouse_id) %>% unique())

seq_level_data <- left_join(mouse_info, seq_level_data, by = 'mouse_id') %>%
  filter(productive_partis) %>% dplyr::rename(tissue = specimen_tissue, cell_type = specimen_cell_subset) %>%
  mutate(group_controls_pooled = factor(group_controls_pooled, levels = group_controls_pooled_factor_levels),
         cell_type = factor(cell_type, levels = c('naive','GC','PC','mem')),
         tissue = factor(tissue, levels = c('LN','spleen','BM')))

#-------  Distribution of n. of apparent mutations in naive sequences (first across all mice and tissues)
overall_naive_mutation_rate <- seq_level_data %>%
  filter(cell_type == 'naive') %>%
  summarise(total_n_mutations = sum(n_mutations_partis),
            total_n_bases = sum(seq_length_partis)) %>%
  mutate(mutation_rate = total_n_mutations/total_n_bases) %>%
  pull(mutation_rate) 

mean_length_naive_seqs <- seq_level_data %>%
  filter(cell_type == 'naive') %>% 
  summarise(mean_length_of_naive_seqs = round(mean(seq_length_partis))) %>% pull(mean_length_of_naive_seqs)
mean_length_naive_seqs

distribution_mutations_naive_global <- seq_level_data %>%
  filter(cell_type == 'naive', tissue == 'spleen') %>%
  group_by(n_mutations_partis) %>%
  dplyr::count() %>%
  ungroup() %>%
  mutate(obs_fraction = n/sum(n),
         expected_fraction_from_error_rate = dbinom(x = n_mutations_partis,
                                                    size = mean_length_naive_seqs,
                                                    prob = estimated_seq_error_rate),
         expected_fraction_from_mean_naive_mutation_rate = dbinom(x = n_mutations_partis,
                                                                  size = mean_length_naive_seqs,
                                                                  prob = overall_naive_mutation_rate))

distribution_mutations_naive_global %>%
  pivot_longer(cols = c('obs_fraction', 'expected_fraction_from_error_rate','expected_fraction_from_mean_naive_mutation_rate')) %>%
  mutate(name = factor(name, levels = c('obs_fraction', 'expected_fraction_from_error_rate','expected_fraction_from_mean_naive_mutation_rate'))) %>%
  ggplot(aes(x = n_mutations_partis, y = value, color = name, group = name)) +
  geom_point() +
  geom_line() +
  xlim(0, 20) +
  xlab('Number of mutations from inferred germline sequence') +
  ylab('Frequency') +
  scale_color_discrete(labels = c('Observed fraction','Binomial expectation from estimated error rate',
                                  'Binomial expectation from mean mutation rate in naive sequences'),
                       name = '') +
  theme(legend.position = 'top')

# Distribution pooling mice in the same group across tissues




# Distribution by mouse/tissue
distribution_mutations_naive_by_mouse_tissue <- seq_level_data %>%
  filter(cell_type == 'naive') %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, tissue, n_mutations_partis) %>%
  summarise(n_naive_seqs = n(), mean_naive_seq_length = mean(seq_length_partis)) %>%
  ungroup() %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, tissue) %>%
  mutate(total_naive_seqs = sum(n_naive_seqs),
         obs_fraction = n_naive_seqs / total_naive_seqs,
         mean_naive_seq_length = round(sum(mean_naive_seq_length*obs_fraction))) %>%
  mutate(expected_fraction_from_error_rate = dbinom(x = n_mutations_partis,
                                                    size = mean_naive_seq_length[1],
                                                    prob = estimated_seq_error_rate)) %>%
  ungroup()


lapply(list('spleen','LN','BM'),
       FUN = function(tis, distribution_mutations_naive_by_mouse_tissue){
         pl <- distribution_mutations_naive_by_mouse_tissue %>%
           pivot_longer(cols = c('obs_fraction', 'expected_fraction_from_error_rate')) %>%
           mutate(mouse_id = factor(mouse_id, levels = mouse_id_factor_levels)) %>%
           mutate(name = factor(name, levels = c('obs_fraction', 'expected_fraction_from_error_rate'))) %>%
           filter(tissue == tis) %>%
           filter(n_mutations_partis <= 10) %>%
           ggplot(aes(x = n_mutations_partis, y = value, color = name, group = name, shape = name)) +
           geom_point() +
           geom_line() +
           facet_wrap('mouse_id', ncol = 4, scales = 'free') +
           xlab('Number of mutations from inferred germline sequence') +
           ylab('Frequency') +
           scale_color_discrete(labels = c('Observed fraction','Binomial expectation from estimated error rate'),
                                name = '') +
           scale_shape_manual(name = '', values  = c(19,1)) +
           theme(legend.position = 'top') +
           guides(shape = 'none')
         save_plot(paste0('../results/n_SHM_over_time/mutations_naive_seqs_by_mouse_',tis,'.pdf'),
                   pl,
                   base_width = 20, base_height = 30)
       },
       distribution_mutations_naive_by_mouse_tissue = distribution_mutations_naive_by_mouse_tissue
       )


  






# Plot with fraction of unmutate naive sequences by mouse / tissue
distribution_mutations_naive %>% 
  filter(n_mutations_partis == 0) %>%
  ggplot(aes(x = group_controls_pooled, y = fraction_seqs, color = infection_status)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point() +
  facet_wrap('tissue') +
  background_grid()




seq_level_data %>% filter(group_controls_pooled != 'control', tissue == 'LN', cell_type == 'GC') %>%
  filter(n_mutations_partis > 0) %>%
  ggplot(aes(x = group_controls_pooled, y = n_mutations_partis)) +
  geom_violin() +
  facet_wrap('cell_type') +
  scale_y_log10()

mean_n_mutations <- seq_level_data %>% group_by(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type) %>%
  summarise(mean_n_mutations = mean(n_mutations_partis))

mean_n_mutations %>% ggplot(aes(x = group_controls_pooled, y = mean_n_mutations, color = infection_status)) +
  geom_point() +
  geom_boxplot() +
  facet_grid(tissue ~ cell_type, scales = 'free') +
  background_grid()


