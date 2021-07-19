library(tidyverse)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

source('gene_frequency_functions.R')

annotated_seq_files_Greiff2017 <- list.files('../results/partis/seq_data_Greiff2017/', pattern = 'annotated_seqs.csv', full.names = T)
# annotated_seq_files_Greiff2017 <- list.files('~/Desktop/v_gene_selection_files/seq_data_Greiff2017/', pattern = 'annotated_seqs.csv', full.names = T)

annotated_seqs_Greiff2017 <- lapply(as.list(annotated_seq_files_Greiff2017),
                                    FUN = function(path){
                                      dataset = str_remove(rev(strsplit(path,'/')[[1]])[1], '_annotated_seqs\\.csv')
                                      dataset = str_extract(dataset, 'untreated_mouse_[0-9]*')
                                      read_csv(path) %>% mutate(mouse_id = dataset) %>%
                                        select(mouse_id, everything())
                                    })

annotated_seqs_Greiff2017 <- bind_rows(annotated_seqs_Greiff2017)

annotated_seqs_Greiff2017 <- annotated_seqs_Greiff2017 %>%
  filter(productive_partis)

# Do the sequences sorted as naive by Greiff et al. 2017 really look naive?

# I.e., are they not too mutated?

distribution_nt_mutations_v_region_Greiff2017 <- get_distribution_of_mutations(annotated_seqs_Greiff2017, n_mutations_variable = "vgene_mutations_partis_nt",
                                                                               disable_grouping = T)
null_distribution_nt_mutations_v_region_Greiff2017 <- get_null_mutation_distribution_given_length_distribution(annotated_seqs_Greiff2017,
                                                                                                               n_mutations_variable = "vgene_mutations_partis_nt",
                                                                                                               seq_length_variable = "sequenced_bases_in_vgene_region_partis",
                                                                                                               estimated_seq_error_rate = estimated_seq_error_rate,
                                                                                                               disable_grouping = T)
distribution_nt_mutations_v_region_Greiff2017 <- left_join(distribution_nt_mutations_v_region_Greiff2017,
                                                           null_distribution_nt_mutations_v_region_Greiff2017 %>%
                                                             dplyr::rename(vgene_mutations_partis_nt = n_mutations))

distribution_nt_mutations_v_region_Greiff2017 %>% filter(vgene_mutations_partis_nt <=2) %>%
  group_by(mouse_id) %>%
  summarise(fraction_seqs_with_two_or_more_mutations = 1 - sum(obs_fraction))

distribution_nt_mutations_v_region_Greiff2017 %>%
  filter(vgene_mutations_partis_nt <= 20) %>% 
  pivot_longer(cols = c('obs_fraction','null_prob'),
               names_to = 'value_type') %>%
  mutate(obs_and_pred_number = value * compartment_seqs) %>%
  mutate(obs_and_pred_number = ifelse(value_type == 'obs_fraction',n_seqs, obs_and_pred_number)) %>%
  ggplot(aes(x = vgene_mutations_partis_nt, y = value)) +
  geom_point(aes(color = value_type)) +
  geom_line(aes(color = value_type)) +
  geom_text(aes(label = round(obs_and_pred_number,2)),
            position = position_dodge(width = 10), size = 2) +
  scale_y_log10() +
  geom_hline(aes(yintercept = 1/compartment_seqs), linetype = 2) +
  xlab('Number of mutations in V gene region') +
  ylab('Frequency') +
  scale_x_continuous(limits = c(0,20)) +
  scale_color_discrete(name = '', 
                       labels = c('Null model','Observed')) +
  theme(legend.position = c(0.75,0.3),
        plot.title = element_text(size = 10)) +
  ggtitle('Omiting sequences with more than 20 mutations to improve clarity ') +
  facet_wrap('mouse_id', scales = 'free')
  
# Do they form large clones?

clone_size_dist_Greiff2017 <- annotated_seqs_Greiff2017 %>%
  group_by(mouse_id, clone_id_partis) %>%
  dplyr::summarise(n_seqs = dplyr::n()) %>%
  group_by(mouse_id, n_seqs) %>%
  dplyr::summarise(n_clones = dplyr::n()) %>%
  ungroup() %>%
  group_by(mouse_id) %>%
  mutate(freq = n_clones/sum(n_clones)) %>%
  ungroup()

clone_size_dist_Greiff2017 %>% 
  filter(n_seqs >= 3) %>%
  group_by(mouse_id) %>%
  summarise(fraction_clones_with_more_than_3_seqs = sum(freq))

clone_size_dist_Greiff2017 %>%
  ggplot(aes(x = n_seqs, y = freq)) +
  geom_point() +
  geom_line() +
  facet_wrap('mouse_id') +
  xlab('Clone size (n. seqs)') +
  ylab('Frequency') +
  background_grid()
  

# ================= How do their V gene frequencies correlate with those from our experiments? =================
naive_freqs_Greiff2017 <- annotated_seqs_Greiff2017 %>%
  group_by(mouse_id, v_segment_partis) %>%
  dplyr::summarise(n_vgene_seqs = dplyr::n()) %>%
  group_by(mouse_id) %>%
  mutate(total_compartment_seqs = sum(n_vgene_seqs),
         vgene_seq_freq = n_vgene_seqs / total_compartment_seqs) %>%
  dplyr::rename(v_gene = v_segment_partis) %>%
  ungroup()
  
# Read pre-computed naive frequencies from our data
load('../results/precomputed_gene_freqs.RData')
#load('~/Desktop/v_gene_selection_files/precomputed_gene_freqs.RData')
naive_freqs_this_study <- naive_freqs %>%
  dplyr::rename(n_vgene_seqs = n_naive_vgene_seqs,
                vgene_seq_freq = naive_vgene_seq_freq,
                total_compartment_seqs = total_mouse_naive_seqs)

get_Greiff2017_mouse_pairs <- function(naive_freqs_Greiff2017){
  unique_pairs <- naive_freqs_Greiff2017 %>% select(mouse_id) %>% unique() %>%
    dplyr::rename(mouse_id_i = mouse_id) %>%
    mutate(mouse_id_j = mouse_id_i) %>%
    complete(mouse_id_i, mouse_id_j) %>%
    rowwise() %>%
    mutate(pair = paste0(sort(c(mouse_id_i, mouse_id_j)), collapse = ';')) %>%
    ungroup() %>%
    filter(mouse_id_i != mouse_id_j) %>%
    select(pair) %>%
    unique() %>% pull(pair)
  
  compartment_sizes <- naive_freqs_Greiff2017 %>%
    select(mouse_id, total_compartment_seqs) %>%
    unique()
  
  
  internal_function <- function(mouse_pair){
    # Parse mouse_pair, determine type of pair (e.g. 'primary', 'primary/control', 'secondary/control')
    mice <- str_split(mouse_pair,';')[[1]]
    
    mouse1_freqs <- naive_freqs_Greiff2017 %>% filter(mouse_id == mice[1]) %>% select(-total_compartment_seqs)
    mouse2_freqs <- naive_freqs_Greiff2017 %>% filter(mouse_id == mice[2]) %>% select(-total_compartment_seqs)
    
    for(var in colnames(mouse1_freqs)){
      if(var != 'v_gene'){
        colnames(mouse1_freqs)[colnames(mouse1_freqs) == var] <- paste0(var,'_i')
        colnames(mouse2_freqs)[colnames(mouse2_freqs) == var] <- paste0(var,'_j')
      }
    }
    
    # Assemble pair tibble
    pair_changes <- full_join(mouse1_freqs, mouse2_freqs) %>%
      mutate(mouse_pair = mouse_pair,
             mouse_id_i = mice[1],
             mouse_id_j = mice[2]) %>%
      select(mouse_pair,everything())  
    #%>% filter(!is.na(vgene_seq_freq_i), !is.na(vgene_seq_freq_j))
    
    return(pair_changes)
  }
  
  paired_gene_freqs <- lapply(unique_pairs, FUN = internal_function)
  paired_gene_freqs <- bind_rows(paired_gene_freqs) 
  
  paired_gene_freqs <- left_join(paired_gene_freqs,
                                 compartment_sizes %>% dplyr::rename_with(.fn = function(x){paste0(x,'_i')},
                                                                          .cols = c(mouse_id, total_compartment_seqs)))
  paired_gene_freqs <- left_join(paired_gene_freqs,
                                 compartment_sizes %>% dplyr::rename_with(.fn = function(x){paste0(x,'_j')},
                                                                          .cols = c(mouse_id, total_compartment_seqs)))
  paired_gene_freqs <- paired_gene_freqs %>%
    mutate(pair_type = 'Greiff2017', cell_type = 'naive', tissue = NA, day_i = NA, day_j = NA, infection_status_i = NA, infection_status_j = NA,
           group_i = NA, group_j = NA, group_controls_pooled_i = NA, group_controls_pooled_j = NA, exp_naive_ratio_i = NA,
           exp_naive_ratio_j = NA, total_mouse_naive_seqs_i = total_compartment_seqs_i, total_mouse_naive_seqs_j = total_compartment_seqs_j,) %>%
    select(mouse_pair, pair_type, matches('mouse_id'), matches('day'), matches('infection_status'), group_i, group_j, 
           matches('group_controls_pooled'), matches('total_compartment_seqs'), matches('naive_seqs'), v_gene, 
           tissue, cell_type, n_vgene_seqs_i, vgene_seq_freq_i, exp_naive_ratio_i, n_vgene_seqs_j, vgene_seq_freq_j, exp_naive_ratio_j)
  
  return(paired_gene_freqs)
  
}

get_cross_dataset_mouse_pairs <- function(naive_freqs_Greiff2017, naive_freqs_this_study){
  base_function <- function(mouse_id_Greiff2017){
    
    total_Greiff2017_mouse_naive_seqs <- unique(naive_freqs_Greiff2017 %>% filter(mouse_id == mouse_id_Greiff2017) %>% pull(total_compartment_seqs))
    stopifnot(length(total_Greiff2017_mouse_naive_seqs) == 1)
    
    cross_dataset_pair <- left_join(naive_freqs_this_study %>% 
                rename_with(function(x){paste0(x,'_i')},
                            any_of(c('mouse_id','day','infection_status','group','group_controls_pooled','total_compartment_seqs'))),
              naive_freqs_Greiff2017 %>% filter(mouse_id == mouse_id_Greiff2017) %>% select(-mouse_id, -total_compartment_seqs),
              by = 'v_gene', suffix = c('_i','_j')) %>%
      mutate(mouse_id_j = mouse_id_Greiff2017,
             total_compartment_seqs_j = total_Greiff2017_mouse_naive_seqs) %>%
      mutate(tissue = NA, cell_type = 'naive', exp_naive_ratio_i = NA, exp_naive_ratio_j = NA,
             total_mouse_naive_seqs_i = total_compartment_seqs_i, total_mouse_naive_seqs_j = total_compartment_seqs_j,
             mouse_pair = paste(mouse_id_i, mouse_id_j, sep = ';'),
             pair_type = 'cross-dataset',
             day_j = NA, group_j = NA, infection_status_j = NA, group_controls_pooled_j = NA) %>%
      select(mouse_pair, pair_type, matches('mouse_id'), matches('day'), matches('infection_status'), group_i, group_j, 
             matches('group_controls_pooled'), matches('total_compartment_seqs'), matches('naive_seqs'), v_gene, 
             tissue, cell_type, n_vgene_seqs_i, vgene_seq_freq_i, exp_naive_ratio_i, n_vgene_seqs_j, vgene_seq_freq_j, exp_naive_ratio_j)
    
    return(cross_dataset_pair)
    
  }
  return(bind_rows(lapply(as.list(unique(naive_freqs_Greiff2017$mouse_id)),
                          FUN = base_function)))
}

pairwise_gene_freqs_Greiff2017 <- get_Greiff2017_mouse_pairs(naive_freqs_Greiff2017)
pairwise_gene_freqs_cross_dataset <- get_cross_dataset_mouse_pairs(naive_freqs_Greiff2017, naive_freqs_this_study)
pairwise_gene_freqs_this_study <- pairwise_gene_freqs %>% 
  filter(cell_type == 'naive') %>%
  mutate(pair_type = 'this study')


pairwise_naive_freqs <- bind_rows(list(pairwise_gene_freqs_this_study,
                                       pairwise_gene_freqs_cross_dataset,
                                       pairwise_gene_freqs_Greiff2017))


pairwise_naive_freq_correlations <- get_pairwise_correlations(
  pairwise_gene_freqs = pairwise_naive_freqs, include_freq_ratios = F)


pairwise_naive_freq_correlations %>%
  filter(pair_type == 'cross-dataset') %>%
  ggplot(aes(x = cor_coef_freqs)) +
  geom_histogram(color = 'white', size = 0.3) +
  geom_density() +
  facet_wrap('mouse_id_j') +
  xlab('Spearman correlation in naive V gene frequencies') +
  ylab('Number mice in our study') +
  xlim(0,1) +
  background_grid() +
  scale_x_continuous()

pairwise_naive_freq_correlations %>%
  ggplot(aes(x = pair_type, y = cor_coef_freqs)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(aes(size = total_compartment_seqs_i + total_compartment_seqs_j),
             alpha = 0.5, shape = 21) +
  scale_size_continuous(name = 'Total number of naive\nsequences in mouse pair') +
  xlab('Type of mouse pair') +
  ylab('Spearman correlation in naive V gene frequencies')


get_vgene_freq_correlation_clustering(
  pairwise_correlations = list(freqs = pairwise_naive_freq_correlations),
  cell_type = 'naive', tissue = NULL, metric = 'freqs'
)




