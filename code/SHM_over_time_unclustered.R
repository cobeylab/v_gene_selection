# Looks at the distribution of the number of mutations over time, accounting for the expected distribution due to seq/amplification error alone

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

source('gene_frequency_functions.R')

# Basic info for each clone (germline genes, CDR length, naive CDR seq)
clone_info <- read_csv('../processed_data/clone_info.csv')
# clone_info  <- read_csv('~/Desktop/clone_info.csv')

annotated_seqs <- read_csv('../processed_data/annotated_seqs.csv') 
# annotated_seqs <- read_csv('~/Desktop/annotated_seqs.csv')


annotated_seqs <- left_join(annotated_seqs, clone_info %>% select(mouse_id, clone_id, v_gene))

annotated_seqs <- get_info_from_mouse_id(annotated_seqs)

# What if we exclude naive sequences that are not IGD and IGM
#annotated_seqs <- annotated_seqs %>%
#  filter(!(cell_type == 'naive' & !(isotype %in% c('IGM','IGD'))))

# Observed distribution of the number of nucleotide mutations by mouse and tissue

# ...for the whole sequence:
distribution_nt_mutations_whole_seq <- get_distribution_of_mutations(annotated_seqs, n_mutations_variable = "n_mutations_partis_nt")

# ... and for the V gene region only:
distribution_nt_mutations_v_region <- get_distribution_of_mutations(annotated_seqs, n_mutations_variable = "vgene_mutations_partis_nt")

# Combine these tibbles with the null distribution given the observed sequence length:
null_distribution_given_obs_lengths_whole_seq <- get_null_mutation_distribution_given_length_distribution(
  annotated_seqs, n_mutations_variable = 'n_mutations_partis_nt', seq_length_variable = 'seq_length_partis',
  estimated_seq_error_rate = estimated_seq_error_rate
)

null_distribution_given_obs_lengths_v_region <- get_null_mutation_distribution_given_length_distribution(
  annotated_seqs, n_mutations_variable = 'vgene_mutations_partis_nt', seq_length_variable = 'sequenced_bases_in_vgene_region_partis',
  estimated_seq_error_rate = estimated_seq_error_rate
)

distribution_nt_mutations_whole_seq <- left_join(distribution_nt_mutations_whole_seq,
                                                 null_distribution_given_obs_lengths_whole_seq %>%
                                                   dplyr::rename(n_mutations_partis_nt = n_mutations))


distribution_nt_mutations_v_region <- left_join(distribution_nt_mutations_v_region,
                                                null_distribution_given_obs_lengths_v_region %>%
                                                                       dplyr::rename(vgene_mutations_partis_nt = n_mutations))




lapply(list('spleen','LN','BM'),
       FUN = function(tis, distribution_nt_mutations_whole_seq){
         pl <- distribution_nt_mutations_whole_seq %>%
           filter(cell_type == 'naive') %>%
           pivot_longer(cols = c('obs_fraction', 'null_prob')) %>%
           mutate(mouse_id = factor(mouse_id, levels = mouse_id_factor_levels)) %>%
           mutate(name = factor(name, levels = c('obs_fraction', 'null_prob'))) %>%
           filter(tissue == tis) %>%
           filter(n_mutations_partis_nt <= 15) %>%
           ggplot(aes(x = n_mutations_partis_nt, y = value, color = name, shape = name)) +
           geom_point() +
           geom_line() +
           scale_y_log10() +
           facet_wrap('mouse_id', ncol = 4, scales = 'free') +
           xlab('Number of mutations from inferred germline sequence (whole sequence)') +
           ylab('Fraction of sequences') +
           scale_color_discrete(labels = c('Observed fraction','Null expectation from estimated error rate'),
                                name = '') +
           scale_shape_manual(name = '', values  = c(19,1)) +
           theme(legend.position = 'top') +
           guides(shape = 'none')
         save_plot(paste0('../results/SHM_over_time/distribution_nt_mutations_naive_vs_null_',tis,'_whole_seq_UNCLUSTERED.pdf'),
                   pl,
                   base_width = 20, base_height = 30)
       },
       distribution_nt_mutations_whole_seq = distribution_nt_mutations_whole_seq
)

lapply(list('spleen','LN','BM'),
       FUN = function(tis, distribution_nt_mutations_v_region ){
         pl <- distribution_nt_mutations_v_region  %>%
           filter(cell_type == 'naive') %>%
           pivot_longer(cols = c('obs_fraction', 'null_prob')) %>%
           mutate(mouse_id = factor(mouse_id, levels = mouse_id_factor_levels)) %>%
           mutate(name = factor(name, levels = c('obs_fraction', 'null_prob'))) %>%
           filter(tissue == tis) %>%
           filter(vgene_mutations_partis_nt <= 15) %>%
           ggplot(aes(x = vgene_mutations_partis_nt, y = value, color = name, shape = name)) +
           geom_point() +
           geom_line() +
           scale_y_log10() +
           geom_hline(aes(yintercept = 1/compartment_seqs), linetype = 2) +
           facet_wrap('mouse_id', ncol = 4, scales = 'free') +
           xlab('Number of mutations from inferred germline sequence (V gene region only)') +
           ylab('Fraction of sequences') +
           scale_color_discrete(labels = c('Observed fraction','Null expectation from estimated error rate'),
                                name = '') +
           scale_shape_manual(name = '', values  = c(19,1)) +
           theme(legend.position = 'top') +
           guides(shape = 'none')
         save_plot(paste0('../results/SHM_over_time/distribution_nt_mutations_naive_vs_null_',tis,'_v_region_only_UNCLUSTERED.pdf'),
                   pl,
                   base_width = 20, base_height = 30)
       },
       distribution_nt_mutations_v_region = distribution_nt_mutations_v_region 
)
