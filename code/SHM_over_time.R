library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

source('gene_frequency_functions.R')

estimated_seq_error_rate <- 0.0018

# Read sequence statistics aggregated by sequence cluster
seq_cluster_stats <- read_csv('../processed_data/seq_cluster_stats.csv')

seq_cluster_stats <- get_info_from_mouse_id(seq_cluster_stats) %>%
  mutate(group_controls_pooled = factor(group_controls_pooled, levels = group_controls_pooled_factor_levels))

seq_cluster_stats <- seq_cluster_stats %>%
  mutate(across(matches('mutations'), round)) %>%
  mutate(seq_length_partis = round(seq_length_partis),
         sequenced_bases_in_vgene_region_partis = round(sequenced_bases_in_vgene_region_partis))


# Distribution of the number of nucleotide mutations (whole sequence) by mouse and tissue
distribution_nt_mutations_by_mouse_and_tissue_whole_seq <- seq_cluster_stats %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type, n_mutations_partis_nt) %>%
  summarise(n_seqs = n()) %>%
  ungroup() %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type) %>%
  mutate(compartment_seqs = sum(n_seqs),
         obs_fraction = n_seqs / compartment_seqs) %>%
  ungroup() %>%
  mutate(cell_type = factor(cell_type, levels = c('naive','GC','PC','mem')),
         tissue = factor(tissue, levels = c('LN','spleen','BM')))

# Distribution of sequence lengths (whole sequence) by mouse and tissue
seq_length_distribution_whole_seq <- seq_cluster_stats %>%
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

# Distribution of the number of nucleotide mutations (V gene region only) by mouse and tissue
distribution_nt_mutations_by_mouse_and_tissue_v_gene_region <- seq_cluster_stats %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type, vgene_mutations_partis_nt) %>%
  summarise(n_seqs = n()) %>%
  ungroup() %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type) %>%
  mutate(compartment_seqs = sum(n_seqs),
         obs_fraction = n_seqs / compartment_seqs) %>%
  ungroup() %>%
  mutate(cell_type = factor(cell_type, levels = c('naive','GC','PC','mem')),
         tissue = factor(tissue, levels = c('LN','spleen','BM')))

# Distribution of sequence lengths (V gene region only) by mouse and tissue
seq_length_distribution_v_gene_region <- seq_cluster_stats %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type, sequenced_bases_in_vgene_region_partis) %>%
  summarise(n_seqs = n()) %>%
  ungroup() %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type) %>%
  mutate(compartment_seqs = sum(n_seqs),
         obs_fraction = n_seqs / compartment_seqs) %>%
  ungroup() %>%
  mutate(cell_type = factor(cell_type, levels = c('naive','GC','PC','mem')),
         tissue = factor(tissue, levels = c('LN','spleen','BM')))


# Plot with fraction of unmutated sequences by mouse / tissue / cell type
distribution_nt_mutations_by_mouse_and_tissue_whole_seq %>%
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

distribution_nt_mutations_by_mouse_and_tissue_v_gene_region %>%
  filter(vgene_mutations_partis_nt == 0, compartment_seqs >= 100) %>%
  mutate(label = paste0(n_seqs,'/',compartment_seqs)) %>%
  ggplot(aes(x = group_controls_pooled, y = obs_fraction, color = infection_status)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point() +
  #geom_text(aes(label = label), size = 3, color = 'black') +
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

# Pre-computed null model distributions (whole sequence)
null_model_mutations_whole_seq <- generate_mutation_null_model(seq_cluster_stats, estimated_seq_error_rate, n_mutations_variable = 'n_mutations_partis_nt', 
                                                     seq_length_variable = 'seq_length_partis')

null_model_mutations_v_gene_region <- generate_mutation_null_model(seq_cluster_stats, estimated_seq_error_rate, n_mutations_variable = 'vgene_mutations_partis_nt', 
                                                            seq_length_variable = 'sequenced_bases_in_vgene_region_partis')


# Calculate expected null distribution of mutations given the observed distribution of sequence lengths
null_distribution_given_obs_lengths_whole_seq <- left_join(seq_length_distribution_whole_seq %>%
                                                   select(mouse_id, tissue, cell_type, seq_length_partis, obs_fraction),
                                                 null_model_mutations_whole_seq %>% dplyr::rename(seq_length_partis = length)) %>%
  group_by(mouse_id, tissue, cell_type, n_mutations) %>%
  # For each number of mutations, calculate null probability as 
  # a weighted average by length given the obs. freq distribution of lengths
  summarise(null_prob = sum(obs_fraction*null_prob)) %>%
  ungroup()

null_distribution_given_obs_lengths_v_gene_region <- left_join(seq_length_distribution_v_gene_region %>%
                                                                 select(mouse_id, tissue, cell_type, sequenced_bases_in_vgene_region_partis, obs_fraction),
                                                               null_model_mutations_v_gene_region %>% dplyr::rename(sequenced_bases_in_vgene_region_partis = length)) %>%
  group_by(mouse_id, tissue, cell_type, n_mutations) %>%
  # For each number of mutations, calculate null probability as 
  # a weighted average by length given the obs. freq distribution of lengths
  summarise(null_prob = sum(obs_fraction*null_prob)) %>%
  ungroup()


distribution_nt_mutations_by_mouse_and_tissue_whole_seq <- left_join(distribution_nt_mutations_by_mouse_and_tissue_whole_seq,
          null_distribution_given_obs_lengths_whole_seq %>%
            dplyr::rename(n_mutations_partis_nt = n_mutations))

distribution_nt_mutations_by_mouse_and_tissue_v_gene_region <- left_join(distribution_nt_mutations_by_mouse_and_tissue_v_gene_region,
                                                                     null_distribution_given_obs_lengths_v_gene_region %>%
                                                                       dplyr::rename(vgene_mutations_partis_nt = n_mutations))


lapply(list('spleen','LN','BM'),
       FUN = function(tis, distribution_nt_mutations_by_mouse_and_tissue_whole_seq){
         pl <- distribution_nt_mutations_by_mouse_and_tissue_whole_seq %>%
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
           xlab('Number of mutations from inferred germline sequence (whole sequence)') +
           ylab('Fraction of sequences') +
           scale_color_discrete(labels = c('Observed fraction','Null expectation from estimated error rate'),
                                name = '') +
           scale_shape_manual(name = '', values  = c(19,1)) +
           theme(legend.position = 'top') +
           guides(shape = 'none')
         save_plot(paste0('../results/n_SHM_over_time/distribution_nt_mutations_naive_vs_null_',tis,'_whole_seq.pdf'),
                   pl,
                   base_width = 20, base_height = 30)
       },
       distribution_nt_mutations_by_mouse_and_tissue_whole_seq = distribution_nt_mutations_by_mouse_and_tissue_whole_seq
)

lapply(list('spleen','LN','BM'),
       FUN = function(tis, distribution_nt_mutations_by_mouse_and_tissue_v_gene_region){
         pl <- distribution_nt_mutations_by_mouse_and_tissue_v_gene_region %>%
           filter(cell_type == 'naive') %>%
           pivot_longer(cols = c('obs_fraction', 'null_prob')) %>%
           mutate(mouse_id = factor(mouse_id, levels = mouse_id_factor_levels)) %>%
           mutate(name = factor(name, levels = c('obs_fraction', 'null_prob'))) %>%
           filter(tissue == tis) %>%
           filter(vgene_mutations_partis_nt <= 15) %>%
           ggplot(aes(x = vgene_mutations_partis_nt, y = value, color = name, group = name, shape = name)) +
           geom_point() +
           geom_line() +
           facet_wrap('mouse_id', ncol = 4, scales = 'free') +
           xlab('Number of mutations from inferred germline sequence (V gene region only)') +
           ylab('Fraction of sequences') +
           scale_color_discrete(labels = c('Observed fraction','Null expectation from estimated error rate'),
                                name = '') +
           scale_shape_manual(name = '', values  = c(19,1)) +
           theme(legend.position = 'top') +
           guides(shape = 'none')
         save_plot(paste0('../results/n_SHM_over_time/distribution_nt_mutations_naive_vs_null_',tis,'_v_region_only.pdf'),
                   pl,
                   base_width = 20, base_height = 30)
       },
       distribution_nt_mutations_by_mouse_and_tissue_v_gene_region = distribution_nt_mutations_by_mouse_and_tissue_v_gene_region
)



