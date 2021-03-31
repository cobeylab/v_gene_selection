library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

source('gene_frequency_functions.R')

seq_level_data <- read_csv('../processed_data/annotated_seqs.csv')
# seq_level_data <- read_csv('~/Desktop/seq_level_data.csv')

estimated_seq_error_rate <- 0.0018

# Seq information (isotype, n mutations in CDRs, etc.)
annotated_seqs <- read_csv('../processed_data/annotated_seqs.csv')

annotated_seqs <- annotated_seqs %>%
  mutate(across(c('clone_id_partis','partis_uniq_ref_seq','seq_id'), as.character))

annotated_seqs$specimen_cell_subset[annotated_seqs$specimen_cell_subset == 'naïve'] <- 'naive'

# Sequence clustering information to estimate unique sequences
unique_seq_clusters <- read_csv('../processed_data/unique_seq_clusters.csv')
unique_seq_clusters <- unique_seq_clusters %>%
  mutate(across(c('clone_id','cluster_ref_seq','seq_id'), as.character))


combined_tibble <- left_join(unique_seq_clusters,
                             annotated_seqs %>% dplyr::rename(clone_id = clone_id_partis,
                                                              cell_type = specimen_cell_subset,
                                                              tissue = specimen_tissue)) 

#-------  Distribution of n. of apparent mutations in naive sequences (first across all mice and tissues)
distribution_mutations_global <- seq_level_data %>%
  group_by(cell_type) %>%
  mutate(mean_seq_length = round(mean(seq_length_partis))) %>%
  group_by(n_mutations_partis, cell_type, mean_seq_length) %>%
  summarise(n_seqs = n()) %>%
  ungroup() %>%
  group_by(cell_type) %>%
  mutate(obs_fraction = n_seqs/sum(n_seqs),
         expected_fraction_from_error_rate = dbinom(x = n_mutations_partis,
                                                    size = mean_seq_length,
                                                    prob = estimated_seq_error_rate))

distribution_mutations_global %>%
  filter(cell_type == 'naive') %>%
  pivot_longer(cols = c('obs_fraction', 'expected_fraction_from_error_rate')) %>%
  mutate(name = factor(name, levels = c('obs_fraction', 'expected_fraction_from_error_rate'))) %>%
  ggplot(aes(x = n_mutations_partis, y = value, color = name, group = name)) +
  geom_point() +
  geom_line() +
  xlim(0, 20) +
  xlab('Number of mutations from inferred germline sequence') +
  ylab('Frequency') +
  scale_color_discrete(labels = c('Observed fraction','Binomial expectation from estimated error rate'),
                       name = '') +
  theme(legend.position = 'top')

# Distribution pooling mice in the same group across tissues
distribution_mutations_naive_by_group_and_tissue <- seq_level_data %>%
  group_by(group_controls_pooled, tissue, cell_type) %>%
  mutate(mean_seq_length = round(mean(seq_length_partis))) %>%
  group_by(group_controls_pooled, tissue, cell_type, mean_seq_length, n_mutations_partis) %>%
  summarise(n_seqs = n()) %>%
  ungroup() %>%
  group_by(group_controls_pooled, tissue, cell_type) %>%
  mutate(obs_fraction = n_seqs/sum(n_seqs)) %>%
  mutate(expected_fraction_from_error_rate = dbinom(x = n_mutations_partis,
                                                    size = mean_seq_length,
                                                    prob = estimated_seq_error_rate)) %>%
  ungroup()

distribution_mutations_naive_by_group_and_tissue %>%
  filter(cell_type == 'naive') %>%
  pivot_longer(cols = c('obs_fraction', 'expected_fraction_from_error_rate')) %>%
  mutate(name = factor(name, levels = c('obs_fraction', 'expected_fraction_from_error_rate'))) %>%
  ggplot(aes(x = n_mutations_partis, y = value, color = name, group = name)) +
  geom_point() +
  geom_line() +
  xlim(0, 20) +
  facet_grid(group_controls_pooled~tissue) +
  xlab('Number of mutations from inferred germline sequence') +
  ylab('Frequency') +
  scale_color_discrete(labels = c('Observed fraction','Binomial expectation from estimated error rate'),
                       name = '') +
  theme(legend.position = 'top')


# Distribution by mouse/tissue
distribution_mutations_naive_by_mouse_tissue <- seq_level_data %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type) %>%
  mutate(mean_seq_length = round(mean(seq_length_partis))) %>%
  ungroup() %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type, mean_seq_length, n_mutations_partis) %>%
  summarise(n_seqs = n()) %>%
  ungroup() %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type) %>%
  mutate(obs_fraction = n_seqs / sum(n_seqs)) %>%
  mutate(expected_fraction_from_error_rate = dbinom(x = n_mutations_partis,
                                                    size = mean_seq_length,
                                                    prob = estimated_seq_error_rate)) %>%
  ungroup()


lapply(list('spleen','LN','BM'),
       FUN = function(tis, distribution_mutations_naive_by_mouse_tissue){
         pl <- distribution_mutations_naive_by_mouse_tissue %>%
           filter(cell_type == 'naive') %>%
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


# Plot with fraction of unmutated sequences by mouse / tissue / cell type
distribution_mutations_naive_by_mouse_tissue %>%
  filter(n_mutations_partis == 0) %>%
  ggplot(aes(x = group_controls_pooled, y = obs_fraction, color = infection_status)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point() +
  facet_grid(tissue~cell_type) +
  xlab('Group') +
  ylab('Fraction of sequences unmutated') +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 20, vjust = 0.5)) +
  background_grid()


distribution_mutations_naive_by_group_and_tissue %>%
  ggplot(aes(x = group_controls_pooled, y = n_mutations_partis, size = obs_fraction, alpha = obs_fraction)) +
  facet_grid(tissue~cell_type, scales = 'free') +
  geom_point()

seq_level_data %>%
  ggplot(aes(x = group_controls_pooled, y = n_mutations_partis, color = infection_status)) +
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

distribution_mutations_naive_by_group_and_tissue %>%
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


