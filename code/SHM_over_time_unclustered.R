# Looks at the distribution of the number of mutations over time, accounting for the expected distribution due to seq/amplification error alone

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

source('gene_frequency_functions.R')

estimated_seq_error_rate <- 0.0018

annotated_seqs <- read_csv('../processed_data/annotated_seqs.csv') 

annotated_seqs <- annotated_seqs %>%
  mutate(across(c('clone_id_partis','partis_uniq_ref_seq','seq_id'), as.character))

annotated_seqs$specimen_cell_subset[annotated_seqs$specimen_cell_subset == 'naïve'] <- 'naive'

annotated_seqs <- annotated_seqs %>% filter(!is.na(n_mutations_partis_nt), !is.na(vgene_mutations_partis_nt))
annotated_seqs <- get_info_from_mouse_id(annotated_seqs)
annotated_seqs <- annotated_seqs %>% dplyr::rename(tissue = specimen_tissue, cell_type = specimen_cell_subset)

# Ignore sequences inferred to be unproductive by partis
annotated_seqs <- annotated_seqs %>% filter(productive_partis)

# What if we exclude naive sequences that are not IGD?
#annotated_seqs <- annotated_seqs %>%
#  filter(!(cell_type == 'naive' & isotype != 'IGD'))


# Gets distribution of the number of mutations by mouse, tissue and cell type
get_distribution_of_mutations <- function(annotated_seqs, n_mutations_variable){
  # n_mutations_variable: the name of the variable with the number of mutations to summarize
  # e.g. 'n_mutations_partis_nt' is the number of nt mutations across the whole sequence,
  # 'vgene_mutations_partis_nt' is the number of mutations in the V gene region only
  
  grouping_vars <- c('mouse_id','day','infection_status','group_controls_pooled',
                     'tissue', 'cell_type', n_mutations_variable)
  
  dist_n_mutations <- annotated_seqs %>%
    group_by(across(grouping_vars)) %>%
    summarise(n_seqs = n()) %>%
    ungroup() %>%
    group_by(across(grouping_vars[grouping_vars!= n_mutations_variable])) %>%
    mutate(compartment_seqs = sum(n_seqs),
           obs_fraction = n_seqs / compartment_seqs) %>%
    ungroup() %>%
    mutate(cell_type = factor(cell_type, levels = c('naive','GC','PC','mem')),
           tissue = factor(tissue, levels = c('LN','spleen','BM')))
  return(dist_n_mutations)
  
}

# Gets distribution of sequence length by mouse, tissue and cell type
get_seq_length_distribution <- function(annotated_seqs, seq_length_variable){
  # seq_length_variable: e.g. 'seq_length_partis' or 'sequenced_bases_in_vgene_region_partis'
  grouping_vars <- c('mouse_id', 'day', 'infection_status', 'group_controls_pooled', 'tissue', 'cell_type',
                     seq_length_variable)
  
  seq_length_dist <- annotated_seqs %>%
    group_by(across(grouping_vars)) %>%
    summarise(n_seqs = n()) %>%
    ungroup() %>%
    group_by(across(grouping_vars[grouping_vars!= seq_length_variable])) %>%
    mutate(compartment_seqs = sum(n_seqs),
           obs_fraction = n_seqs / compartment_seqs) %>%
    ungroup() %>%
    mutate(cell_type = factor(cell_type, levels = c('naive','GC','PC','mem')),
           tissue = factor(tissue, levels = c('LN','spleen','BM')))
  
  return(seq_length_dist)
  
}

# Calculates expected null distribution of mutations given the observed distribution of sequence lengths
get_null_mutation_distribution_given_length_distribution <- function(annotated_seqs, n_mutations_variable, seq_length_variable,
                                                                     estimated_seq_error_rate){
  # Get observed distribution of sequence lengths 
  seq_length_dist <- get_seq_length_distribution(annotated_seqs, seq_length_variable)
  
  # For a range of possible sequence lengths, calculate probability of observing 0,1,2,... mutations
  null_model_base_probs <- generate_mutation_null_model(annotated_seqs, estimated_seq_error_rate, n_mutations_variable = n_mutations_variable,
                                                        seq_length_variable = seq_length_variable)
  
  names(null_model_base_probs)[names(null_model_base_probs) == 'length'] <- seq_length_variable
  
  # Calculate expected null distribution of mutations given the observed distribution of sequence lengths
  null_distribution_given_obs_lengths <- left_join(seq_length_dist  %>%
              select(mouse_id, tissue, cell_type, matches(seq_length_variable), obs_fraction),
            null_model_base_probs) %>%
    group_by(mouse_id, tissue, cell_type, n_mutations) %>%
    # For each number of mutations, calculate null probability as 
    # a weighted average across the observed freq distribution of lengths
    summarise(null_prob = sum(obs_fraction*null_prob)) %>%
    ungroup()
  
  sums <- null_distribution_given_obs_lengths %>% group_by(mouse_id, tissue, cell_type) %>%
    summarise(S = sum(null_prob)) %>% ungroup() %>% select(S) %>% unique() %>% pull(S)
  
  stopifnot(all(abs(sums-1) < 1e-5))
  
  return(null_distribution_given_obs_lengths)
  
}

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

