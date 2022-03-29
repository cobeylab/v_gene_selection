library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggdendro)
library(cowplot)
theme_set(theme_cowplot())
library(scales)
library(viridis)

source('gene_frequency_functions.R')
source('plot_options.R')

min_compartment_size = 100

# These analyses are based on all productive sequences (as opposed to unique sequences only)

# ====== Some functions ======


# Vectorized internal function used by compute_CDR3_similarity
compute_similarity <- function(str1, str2){
  similarity <- mapply(x = str1, y = str2, 
                       FUN = function(x,y){
                         sum(str_split(x,'')[[1]] == str_split(y,'')[[1]]) / nchar(x)
                       })
  names(similarity) <- c()
  
  return(similarity)
}

# As a test, we should get all ones if passing a character vector of identical sequences
all(compute_similarity(c('AAAA','AAAA','AAAA','AAAA'), c('AAAA','AAAA','AAAA','AAAA')) == 1)
# And all 0s if passing entirely different sequences
all(compute_similarity(c('AAAA','BBBB','CCCC','DDDD'), c('XXXX','ZZZZ','YYYY','OOOO')) == 0)


# ====== Analyses =======


clone_info <- read_csv('../processed_data/clone_info.csv')
# clone_info <- read_csv('~/Desktop/v_gene_selection/processed_data/clone_info.csv')

annotated_seqs <- read_csv('../processed_data/annotated_seqs.csv') %>%
  mutate(seq_id = as.character(seq_id),
         clone_id = as.character(clone_id))
 
# annotated_seqs <- read_csv('~/Desktop/v_gene_selection/processed_data/annotated_seqs.csv') %>% mutate(seq_id = as.character(seq_id), clone_id = as.character(clone_id))

annotated_seqs <- annotated_seqs %>%
  filter(productive_partis)

annotated_seqs <- get_info_from_mouse_id(annotated_seqs) %>%
  mutate(group_controls_pooled = factor(group_controls_pooled, levels = group_controls_pooled_factor_levels),
         cdr3_length = nchar(cdr3_seq_partis))

load('../results/precomputed_gene_freqs_all_seqs.RData')
#load('~/Desktop/v_gene_selection/results/precomputed_gene_freqs_all_seqs.RData')

# Focus on LN B cells and naive B cells from all tissues

CDR3_seqs_naive <- annotated_seqs %>%
  filter(group_controls_pooled != 'control',
         cell_type == 'naive') %>%
  select(mouse_id, day, infection_status, group, group_controls_pooled, seq_id, tissue, cell_type, cdr3_length,
         cdr3_seq_partis) %>%
  group_by(mouse_id, day, infection_status, group, group_controls_pooled) %>%
  mutate(total_compartment_seqs = n(),
         tissue = 'all') %>%
  ungroup() 


CDR3_seqs_LN_experienced <- annotated_seqs %>%
  filter(group_controls_pooled != 'control',
         cell_type %in% c('GC','PC','mem'),
         tissue == 'LN') %>%
  select(mouse_id, day, infection_status, group, group_controls_pooled, seq_id, tissue, cell_type, cdr3_length,
         cdr3_seq_partis) %>%
  group_by(mouse_id, day, infection_status, group, group_controls_pooled, tissue, cell_type) %>%
  mutate(total_compartment_seqs = n()) %>%
  ungroup() 

# Plot CDR3 length distributions for each mouse for each cell type in the LN (compared with all-tissue naive)


CDR3_lengths <- bind_rows(CDR3_seqs_naive,
                          CDR3_seqs_LN_experienced) %>%
  mutate(cell_type = case_when(
    cell_type == 'naive' ~ 'Naive cells (all tissues)',
    cell_type == 'GC' ~ 'Lymph node GC cells',
    cell_type == 'PC' ~ 'Lymph node PC cells',
    cell_type == 'mem' ~ 'Lymph node memory cells'
  )) %>%
  mutate(cell_type = factor(cell_type,
                            levels = c('Naive cells (all tissues)', 'Lymph node GC cells',
                                       'Lymph node PC cells', 'Lymph node memory cells')))


CDR3_length_distribution_pl <- CDR3_lengths %>%
  ggplot(aes(x = group_controls_pooled, y = cdr3_length, group = mouse_id, color = group_controls_pooled)) +
  geom_boxplot() +
  facet_wrap('cell_type', nrow = 1) +
  theme(legend.position = 'none',
        axis.text.x = element_text(size = 10)) +
  xlab('Group') +
  ylab('CDR3 amino acid length') +
  scale_x_discrete(labels = function(x){str_replace(x, '-','\n')}) +
  scale_y_continuous(breaks = seq(0,25,5), limits = c(0, NA))
  
save_plot('../figures/all_seqs_freqs/CDR3_length_distribution.pdf',
          CDR3_length_distribution_pl,
          base_height = 5, base_width = 16)

# CDR3 similarity

# For each pair of individuals, take a sample of CDR3 sequences from mouse i and pair it with a sample of sequences from mouse j
sample_sequence_pairs <- function(pairwise_gene_freqs, annotated_seqs, sample_size, selected_tissues, selected_cell_types){
  
  
  # pairwise_gene_freqs object just used to conveniently select pairs of mice
  mouse_pairs <- pairwise_gene_freqs %>%
    filter(day_i == day_j, # Only samples infected mice pair from the same time point
           mouse_id_i != mouse_id_j,
           group_controls_pooled_i != 'control', group_controls_pooled_j != 'control',
           group_controls_pooled_i == group_controls_pooled_j) %>%
    select(mouse_pair, mouse_id_i, mouse_id_j, day_i, group_controls_pooled_i) %>%
    unique() 

   
  # Internal function to sample paired sequences for a single mouse pair
  sample_pair_seqs <- function(annotated_seqs, mouse_pair, selected_tissues, selected_cell_types){
    
    if(all(selected_tissues == 'all')){
      adjust_tissue_column <- T
      selected_tissues <- c('LN','BM','spleen')
    }else{
      adjust_tissue_column <- F
    }
    
    mouse_ids <- str_split(mouse_pair, ';')[[1]]
    mouse_id_i <- mouse_ids[1]
    mouse_id_j <- mouse_ids[2]
    
 
    seqs_mouse_i <- annotated_seqs %>%
      filter(mouse_id == mouse_id_i, tissue %in% selected_tissues, cell_type %in% selected_cell_types) %>%
      select(mouse_id, tissue, cell_type, cdr3_seq_partis, cdr3_length)
    
    seqs_mouse_j <- annotated_seqs %>%
      filter(mouse_id == mouse_id_j, tissue %in% selected_tissues, cell_type %in% selected_cell_types) %>%
      select(mouse_id, tissue, cell_type, cdr3_seq_partis, cdr3_length) 
    
    if(adjust_tissue_column){
      seqs_mouse_i$tissue <- 'all tissues'
      seqs_mouse_j$tissue <- 'all tissues'
    }
    

    # Take a sample of sequences from mouse i for each tissue / cell type combination
    sample_mouse_i <- seqs_mouse_i %>%
      group_by(mouse_id, tissue, cell_type) %>%
      slice_sample(n = sample_size) %>%
      arrange(tissue, cell_type, cdr3_length)
      
    cdr3_length_distribution_in_mouse_i_sample <- sample_mouse_i %>%
      group_by(mouse_id, tissue, cell_type, cdr3_length) %>% count() %>%
      dplyr::rename(mouse_i_sample_size = n) %>%
      ungroup() %>%
      select(tissue, cell_type, cdr3_length, mouse_i_sample_size)
      
    # For each CDR3 length in mouse i sample (for each tissue and cell type), sample the same number of sequences from mouse_j
    # (sample as many as there are, if fewer than the number in the mouse i sample)
    sample_mouse_j <- left_join(seqs_mouse_j, cdr3_length_distribution_in_mouse_i_sample,
                                by = c("tissue", "cell_type", "cdr3_length")) %>%
        arrange(tissue, cell_type, cdr3_length) %>%
        group_by(mouse_id, tissue, cell_type, cdr3_length) %>%
        mutate(sampling_index = sample(1:n())) %>%
        ungroup() %>%
        filter(sampling_index <= mouse_i_sample_size) %>%
        select(mouse_id, tissue, cell_type, cdr3_seq_partis, cdr3_length)
        
    
    paired_sample <- left_join(sample_mouse_i %>%
                                 dplyr::rename(mouse_id_i = mouse_id, cdr3_seq_i = cdr3_seq_partis) %>%
                                 group_by(tissue, cell_type, cdr3_length) %>%
                                 mutate(positional_index = 1:n()) %>%
                                 ungroup(),
                               sample_mouse_j %>%
                                 dplyr::rename(mouse_id_j = mouse_id, cdr3_seq_j = cdr3_seq_partis) %>%
                                 group_by(tissue, cell_type, cdr3_length) %>%
                                 mutate(positional_index = 1:n()),
                               by = c("tissue", "cell_type", "cdr3_length", "positional_index")) %>%
      mutate(mouse_pair = mouse_pair) %>%
      select(mouse_pair, mouse_id_i, mouse_id_j, tissue, cell_type, cdr3_length, everything()) %>%
      select(-positional_index) %>%
      filter(!is.na(cdr3_seq_j))
    
    return(paired_sample)
    
  }
  
  paired_samples <- bind_rows(lapply(as.list(mouse_pairs$mouse_pair), FUN = sample_pair_seqs,
         annotated_seqs = annotated_seqs, selected_tissues = selected_tissues, selected_cell_types = selected_cell_types))
  
  
  paired_samples <- left_join(paired_samples, mouse_pairs) %>%
    dplyr::rename(group_controls_pooled = group_controls_pooled_i, day = day_i) %>%
    select(mouse_pair, matches('mouse_id'), group_controls_pooled, day, everything())
  
  return(paired_samples)
}

# As a control, compare paired CDR3 sequences from naive cells (across all tissues)
paired_naive_CDR3_seqs <- sample_sequence_pairs(pairwise_gene_freqs = pairwise_gene_freqs, annotated_seqs = annotated_seqs,
                                           sample_size = 500, selected_tissues = 'all', selected_cell_types = 'naive')

paired_experienced_LN_CDR3_seqs <- sample_sequence_pairs(pairwise_gene_freqs = pairwise_gene_freqs, annotated_seqs = annotated_seqs,
                                                sample_size = 500, selected_tissues = 'LN', selected_cell_types = c('GC','PC','mem'))

paired_CDR3_seqs <- bind_rows(paired_naive_CDR3_seqs %>% select(-tissue),
                              paired_experienced_LN_CDR3_seqs %>% select(-tissue))

CDR3_similarity <- paired_CDR3_seqs %>%
  mutate(cdr3_similarity = compute_similarity(str1 = cdr3_seq_i, str2 = cdr3_seq_j)) %>%
  group_by(mouse_pair, cell_type) %>%
  mutate(n_seqs_in_comparison = n()) %>%
  ungroup() %>%
  filter(n_seqs_in_comparison >= min_compartment_size) %>% # Remove comparisons with too few sequences
  mutate(cell_type = case_when(
    cell_type == 'naive' ~ 'Naive cells (all tissues)',
    cell_type == 'GC' ~ 'Lymph node GC cells',
    cell_type == 'PC' ~ 'Lymph node PC cells',
    cell_type == 'mem' ~ 'Lymph node memory cells'
  )) %>%
  mutate(cell_type = factor(cell_type,
                            levels = c('Naive cells (all tissues)', 'Lymph node GC cells',
                                       'Lymph node PC cells', 'Lymph node memory cells')))
  

  
CDR3_similarity_plot <- CDR3_similarity %>%
  ggplot(aes(x = group_controls_pooled, y = cdr3_similarity, color = group_controls_pooled)) +
  geom_boxplot() +
  facet_wrap('cell_type', nrow = 1) +
  theme(legend.position = 'none') +
  xlab('Group') +
  ylab('CDR3 similarity between sequences\nfrom different mice') +
  scale_x_discrete(labels = function(x){str_replace(x, '-','\n')})

save_plot('../figures/all_seqs_freqs/CDR3_similarity.pdf',
          CDR3_similarity_plot,
          base_height = 5, base_width = 16)


