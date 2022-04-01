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

# Following Hershberg & Schlomchik 2006, Hershberg & Saini 2015, after Chothia et al. 1998 
amino_acid_classes <- tibble(
  amino_acid = c('F','L','I','M','V','C','W',
                 'Q','R','N','K','D','E',
                 'S','P','T','A','Y','H','G'),
  class = c(rep('hydrophobic', 7),
            rep('hydrophilic', 6),
            rep('neutral', 7))
) %>%
  mutate(code = case_when(class == 'hydrophobic' ~'H',
                          class == 'hydrophilic' ~ 'L',
                          class == 'neutral' ~'N'))

# Converts an amino acid sequence into a string with the biochem class code of each aa (as above)
# (so we can quickly compute biochemical similarity using the function for calculating sequence similarity)
translate_aa_to_biochem_class <- function(seqs){
  chartr(old = paste(amino_acid_classes$amino_acid, collapse = ''),
         new = paste(amino_acid_classes$code, collapse = ''),
         seqs)
}



# These analyses are based on all productive sequences (as opposed to unique sequences only)

# ====== Some functions ======
compute_seq_similarity <- function(str1, str2){
  similarity <- mapply(x = str1, y = str2, 
                       FUN = function(x,y){
                         sum(str_split(x,'')[[1]] == str_split(y,'')[[1]]) / nchar(x)
                       })
  names(similarity) <- c()
  
  return(similarity)
}

# As a test, we should get all ones if passing a character vector of identical sequences
all(compute_seq_similarity(c('AAAA','AAAA','AAAA','AAAA'), c('AAAA','AAAA','AAAA','AAAA')) == 1)
# And all 0s if passing entirely different sequences
all(compute_seq_similarity(c('AAAA','BBBB','CCCC','DDDD'), c('XXXX','ZZZZ','YYYY','OOOO')) == 0)


# For each pair of individuals, take a sample of CDR3 sequences from mouse i and pair it with a sample of sequences from mouse j
# Sequences will be matched by CDR3 length and (optionally) by V allele
sample_sequence_pairs <- function(pairwise_gene_freqs, annotated_seqs, sample_size, selected_tissues, selected_cell_types,
                                  match_v_allele){
  
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
      select(mouse_id, tissue, cell_type, v_gene, cdr3_seq_partis, cdr3_length)
    
    seqs_mouse_j <- annotated_seqs %>%
      filter(mouse_id == mouse_id_j, tissue %in% selected_tissues, cell_type %in% selected_cell_types) %>%
      select(mouse_id, tissue, cell_type, v_gene, cdr3_seq_partis, cdr3_length) 
    
    if(adjust_tissue_column){
      seqs_mouse_i$tissue <- 'all tissues'
      seqs_mouse_j$tissue <- 'all tissues'
    }
    
    
    # Take a sample of sequences from mouse i for each tissue / cell type combination
    sample_mouse_i <- seqs_mouse_i %>%
      group_by(mouse_id, tissue, cell_type) %>%
      slice_sample(n = sample_size) %>%
      arrange(tissue, cell_type, cdr3_length)
    
    
    sample_distribution_grouping_vars <- c('mouse_id', 'tissue', 'cell_type', 'cdr3_length')
    if(match_v_allele){
      sample_distribution_grouping_vars <- c(sample_distribution_grouping_vars, 'v_gene')
    }
    
    cdr3_length_distribution_in_mouse_i_sample <- sample_mouse_i %>%
      group_by(across(sample_distribution_grouping_vars)) %>% count() %>%
      dplyr::rename(mouse_i_sample_size = n) %>%
      ungroup() %>%
      select(tissue, cell_type, cdr3_length, matches('v_gene'), mouse_i_sample_size)
    
    # For each CDR3 length (and, optionally, same v allele) in mouse i sample (for each tissue and cell type),
    # sample the same number of sequences from mouse_j
    # (sample as many as there are, if fewer than the number in the mouse i sample)
    sample_mouse_j <- left_join(seqs_mouse_j, cdr3_length_distribution_in_mouse_i_sample) %>%
      arrange(tissue, cell_type, cdr3_length) %>%
      group_by(across(sample_distribution_grouping_vars)) %>%
      mutate(sampling_index = sample(1:n())) %>%
      ungroup() %>%
      filter(sampling_index <= mouse_i_sample_size) %>%
      select(mouse_id, tissue, cell_type, matches('v_gene'), cdr3_seq_partis, cdr3_length)
    
    if(!match_v_allele){
      sample_mouse_i <- sample_mouse_i %>% select(-v_gene)
      sample_mouse_j <- sample_mouse_j %>% select(-v_gene)
    }
    
    if(nrow(sample_mouse_i) >0 & nrow(sample_mouse_j) > 0){
      paired_sample <- left_join(sample_mouse_i %>%
                                   group_by(across(sample_distribution_grouping_vars)) %>%
                                   mutate(positional_index = 1:n()) %>%
                                   dplyr::rename(mouse_id_i = mouse_id, cdr3_seq_i = cdr3_seq_partis) %>%
                                   ungroup(),
                                 sample_mouse_j %>%
                                   group_by(across(sample_distribution_grouping_vars)) %>%
                                   mutate(positional_index = 1:n()) %>%
                                   dplyr::rename(mouse_id_j = mouse_id, cdr3_seq_j = cdr3_seq_partis) %>%
                                   ungroup()) %>%
        mutate(mouse_pair = mouse_pair) %>%
        select(mouse_pair, mouse_id_i, mouse_id_j, tissue, cell_type, matches('v_gene'), cdr3_length, everything()) %>%
        select(-positional_index) %>%
        filter(!is.na(cdr3_seq_j))
    }else{
      paired_sample <- c()
    }
    
    return(paired_sample)
    
  }
  
  paired_samples <- bind_rows(lapply(as.list(mouse_pairs$mouse_pair), FUN = sample_pair_seqs,
                                     annotated_seqs = annotated_seqs, selected_tissues = selected_tissues, selected_cell_types = selected_cell_types))
  
  
  paired_samples <- left_join(paired_samples, mouse_pairs) %>%
    dplyr::rename(group_controls_pooled = group_controls_pooled_i, day = day_i) %>%
    select(mouse_pair, matches('mouse_id'), group_controls_pooled, day, everything())
  
  return(paired_samples)
}

# Computes sequence and biochemical similarity for paired CDR3 sequence samples produced by sample_sequence_pairs
compute_similarity_matched_samples <- function(matched_samples){
  matched_samples %>%
    mutate(cdr3_seq_similarity = compute_seq_similarity(str1 = cdr3_seq_i, str2 = cdr3_seq_j),
           cdr3_biochem_similarity = compute_seq_similarity(str1 = biochem_cdr3_seq_i, str2 = biochem_cdr3_seq_j)) %>%
    group_by(mouse_pair, cell_type) %>%
    mutate(n_seqs_in_comparison = n()) %>%
    ungroup() %>%
    filter(n_seqs_in_comparison >= min_compartment_size) %>% # Remove comparisons with too few sequences
    mutate(cell_type = case_when(
      cell_type == 'naive' ~ 'Naive cells (all tissues)',
      cell_type == 'GC' ~ 'Lymph node GC cells',
      cell_type == 'PC' ~ 'Lymph node plasma cells',
      cell_type == 'mem' ~ 'Lymph node memory cells'
    )) %>%
    mutate(cell_type = factor(cell_type,
                              levels = c('Naive cells (all tissues)', 'Lymph node GC cells',
                                         'Lymph node plasma cells', 'Lymph node memory cells')))
}

plot_cdr3_similarity <- function(matched_samples_similarity){
  matched_samples_similarity %>%
    pivot_longer(cols = c('cdr3_seq_similarity', 'cdr3_biochem_similarity'), names_to = 'similarity_type',
                 values_to = 'similarity') %>%
    mutate(similarity_type = case_when(
      similarity_type == 'cdr3_seq_similarity' ~ 'sequence similarity',
      similarity_type == 'cdr3_biochem_similarity' ~ 'biochemical similarity'
    )) %>%
    mutate(similarity_type = factor(similarity_type, levels = c('sequence similarity', 'biochemical similarity'))) %>%
    ggplot(aes(x = group_controls_pooled, y = similarity, color = group_controls_pooled)) +
    geom_boxplot() +
    facet_grid(similarity_type~cell_type) +
    theme(legend.position = 'none', axis.text.x = element_text(size = 10)) +
    xlab('Group') +
    ylab('Similarity of length-matched\nCDR3 sequences from different mice') +
    scale_x_discrete(labels = function(x){str_replace(x, '-','\n')})
}


# Get CDR3 sequences for each clone, organize pairs as vectors of sequences, runs vectorized compute_seq_similarity function
# compare_CDR3s_for_all_clone_pairs <- function(cdr3_seqs, mouse_ids){
#   stopifnot(length(unique(nchar(cdr3_seqs))) == 1)
#   
#   seq_vector <- paste(mouse_ids, cdr3_seqs, sep = ';')
#   
#   if(length(seq_vector) > 1){
#     
#     paired_seqs <- as_tibble(t(combn(seq_vector, 2))) %>%
#       separate(V1, into = c('mouse_i', 'seq_i'), sep = ';') %>%
#       separate(V2, into = c('mouse_j', 'seq_j'), sep = ';') %>%
#       mutate(pair_type = case_when(
#         mouse_i == mouse_j ~ 'within mice',
#         mouse_i != mouse_j ~ 'between mice')
#       )
#     
#     similarity <- compute_seq_similarity(str1 = paired_seqs$seq_i, str2 = paired_seqs$seq_j)
#     pair_type <- paired_seqs$pair_type
#   }else{
#     similarity <- NA
#     pair_type <- NA
#   }
#   
#   return(tibble(pair_type = pair_type, similarity = similarity))
#   
# }


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

annotated_seqs <- left_join(annotated_seqs, clone_info %>% select(mouse_id, clone_id, v_gene) %>% mutate(clone_id = as.character(clone_id)))

load('../results/precomputed_gene_freqs_all_seqs.RData')
#load('~/Desktop/v_gene_selection/results/precomputed_gene_freqs_all_seqs.RData')

# Focus on LN B cells and naive B cells from all tissues

CDR3_seqs_naive <- annotated_seqs %>%
  filter(group_controls_pooled != 'control',
         cell_type == 'naive') %>%
  select(mouse_id, day, infection_status, group, group_controls_pooled, clone_id, seq_id, tissue, cell_type, cdr3_length,
         cdr3_seq_partis) %>%
  group_by(mouse_id, day, infection_status, group, group_controls_pooled) %>%
  mutate(total_compartment_seqs = n(),
         tissue = 'all') %>%
  ungroup() 

CDR3_seqs_LN_experienced <- annotated_seqs %>%
  filter(group_controls_pooled != 'control',
         cell_type %in% c('GC','PC','mem'),
         tissue == 'LN') %>%
  select(mouse_id, day, infection_status, group, group_controls_pooled, clone_id, seq_id, tissue, cell_type, cdr3_length,
         cdr3_seq_partis) %>%
  group_by(mouse_id, day, infection_status, group, group_controls_pooled, tissue, cell_type) %>%
  mutate(total_compartment_seqs = n()) %>%
  ungroup() 

# Add v gene info for each sequence (kept separate in a per-clone tibble to reduce size)
CDR3_seqs_naive <- left_join(CDR3_seqs_naive, clone_info %>% select(mouse_id, clone_id, v_gene) %>% mutate(clone_id = as.character(clone_id)))
CDR3_seqs_LN_experienced <- left_join(CDR3_seqs_LN_experienced, clone_info %>% select(mouse_id, clone_id, v_gene) %>% mutate(clone_id = as.character(clone_id)))

CDR3_seqs <- bind_rows(CDR3_seqs_naive,
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


# ======= CDR3 LENGTH ==========

# Plot CDR3 length distributions for each mouse for each cell type in the LN (compared with all-tissue naive)


CDR3_length_distribution_pl <- CDR3_seqs %>%
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

# ======== CDR3 SIMILARITY =======

# Sample length-matched CDR3 sequences from different mice, compute similarity, compare with length-matched CDR3 sequences from naive cells

length_matched_sample_NAIVE <- sample_sequence_pairs(pairwise_gene_freqs = pairwise_gene_freqs, annotated_seqs = annotated_seqs,
                                           sample_size = 500, selected_tissues = 'all', selected_cell_types = 'naive', match_v_allele = F)

length_matched_sample_LNEXP <- sample_sequence_pairs(pairwise_gene_freqs = pairwise_gene_freqs, annotated_seqs = annotated_seqs,
                                                     sample_size = 500, selected_tissues = 'LN', selected_cell_types = c('GC','PC','mem'), match_v_allele = F)

length_matched_CDR3_sample <-  bind_rows(length_matched_sample_NAIVE %>% select(-tissue),
                                    length_matched_sample_LNEXP %>% select(-tissue))  %>%
  mutate(biochem_cdr3_seq_i = translate_aa_to_biochem_class(cdr3_seq_i),
         biochem_cdr3_seq_j = translate_aa_to_biochem_class(cdr3_seq_j))


length_matched_CDR3_similarity <- compute_similarity_matched_samples(matched_samples = length_matched_CDR3_sample)
  
length_matched_CDR3_similarity_plot <- plot_cdr3_similarity(length_matched_CDR3_similarity)

save_plot('../figures/all_seqs_freqs/length_matched_CDR3_similarity.pdf',
          length_matched_CDR3_similarity_plot,
          base_height = 8, base_width = 16)


# Sample length-matched **AND ALLELE-MATCHED** CDR3 sequences from different mice, compute similarity,
# compare with length- and allele-matched CDR3 sequences from naive cells

length_and_allele_matched_sample_NAIVE <- sample_sequence_pairs(pairwise_gene_freqs = pairwise_gene_freqs, annotated_seqs = annotated_seqs,
                                                     sample_size = 500, selected_tissues = 'all', selected_cell_types = 'naive', match_v_allele = T)

length_and_allele_matched_sample_LNEXP <- sample_sequence_pairs(pairwise_gene_freqs = pairwise_gene_freqs, annotated_seqs = annotated_seqs,
                                                     sample_size = 500, selected_tissues = 'LN', selected_cell_types = c('GC','PC','mem'),
                                                     match_v_allele = T)
length_and_allele_matched_sample <- bind_rows(length_and_allele_matched_sample_NAIVE %>% select(-tissue),
                                              length_and_allele_matched_sample_LNEXP %>% select(-tissue))  %>%
  mutate(biochem_cdr3_seq_i = translate_aa_to_biochem_class(cdr3_seq_i),
         biochem_cdr3_seq_j = translate_aa_to_biochem_class(cdr3_seq_j))

length_and_allele_matched_CDR3_similarity <-
  compute_similarity_matched_samples(matched_samples = length_and_allele_matched_sample) 
  
  
# Looks like length- and allele-matched CDR3 sequences become much more similar on day 56 LN plasma cells:
length_and_allele_matched_CDR3_similarity_plot <- plot_cdr3_similarity(length_and_allele_matched_CDR3_similarity)


save_plot('../figures/all_seqs_freqs/length_and_allele_matched_CDR3_similarity.pdf',
          length_and_allele_matched_CDR3_similarity_plot,
          base_height = 5, base_width = 16)

# =Looking closely at pairs of mice from day 56:
length_and_allele_matched_CDR3_similarity_day56_LN_PC <- length_and_allele_matched_CDR3_similarity %>%
  filter(cell_type == 'Lymph node plasma cells', group_controls_pooled == 'secondary-56') %>%
  ggplot(aes(x = mouse_pair, y = cdr3_seq_similarity)) +
  geom_boxplot(outlier.alpha = 0, color = 'blue') +
  geom_point(alpha = 0.5, position = position_jitter(width = 0.1, height = 0)) + 
  xlab("Mouse pair") +
  ylab("Similarity of length- and allele-matched\nCDR3 sequences from lymph node plasma cells")

save_plot('../figures/all_seqs_freqs/length_and_allele_matched_CDR3_similarity_day56_LN_PC.pdf',
          length_and_allele_matched_CDR3_similarity_day56_LN_PC,
          base_height = 5, base_width = 16)

# Looking at sequences driving this pattern (those associated with seq. pairs of over 75% similarity)
cdr3_seqs_driving_similarity <-  length_and_allele_matched_CDR3_similarity %>%
  filter(cell_type == 'Lymph node plasma cells', group_controls_pooled == 'secondary-56') %>%
  filter(cdr3_seq_similarity > 0.75) %>%
  select(v_gene,cdr3_seq_i, cdr3_seq_j) %>%
  pivot_longer(cols = c('cdr3_seq_i', 'cdr3_seq_j'), values_to = 'cdr3_seq') %>%
  group_by(v_gene, cdr3_seq) %>%
  count() %>% arrange(desc(n)) 

# Plot showing what those sequences look like
high_similarity_length_and_allele_matched_seqs_day56_LN_PCs <- cdr3_seqs_driving_similarity %>%
  group_by(cdr3_seq) %>%
  summarise(n = sum(n)) %>%
  mutate(seq_rank = rank(-n, ties.method = 'first')) %>%
  ggplot(aes(x = 1, y = seq_rank, label = cdr3_seq)) +
  geom_text() +
  scale_y_reverse() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  xlab('Sequences driving high similarity in length-\nand allele-matched CDR3 sequences from day 56 LN PCs') + ylab('Rank in paired samples')
save_plot('../figures/all_seqs_freqs/high_similarity_length_and_allele_matched_seqs_day56_LN_PCs.pdf',
          high_similarity_length_and_allele_matched_seqs_day56_LN_PCs,
          base_height = 10, base_width = 6)

# Collectively, how common were those sequences in LN cells from mice at all time points?
combined_freq_of_day56_LN_PC_convergent_CDRs <- CDR3_seqs %>% 
  filter(total_compartment_seqs >= min_compartment_size) %>%
  filter(cdr3_seq_partis %in% unique(cdr3_seqs_driving_similarity$cdr3_seq)) %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, cell_type, total_compartment_seqs) %>%
  summarise(total_n_focal_seqs = n()) %>%
  ungroup() %>%
  mutate(combined_freq_of_focal_seqs = total_n_focal_seqs / total_compartment_seqs) %>%
  ggplot(aes(x = group_controls_pooled, y = combined_freq_of_focal_seqs, color = group_controls_pooled)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(size = 3, position = position_jitter(width = 0.1, height = 0),
             alpha = 0.7) +
  facet_wrap('cell_type', nrow = 1) +
  xlab('Group') +
  ylab('Combined frequency of convergent\nday-56 LN PC CDR3 sequences') +
  theme(legend.position = 'none',
        axis.text.x = element_text(size = 10)) +
  scale_x_discrete(labels = function(x){str_replace(x, '-','\n')})

save_plot('../figures/all_seqs_freqs/combined_freq_of_day56_LN_PC_convergent_CDRs.pdf',
          combined_freq_of_day56_LN_PC_convergent_CDRs,
          base_height = 5, base_width = 16)
  

allele_usage_day56_LN_PC_convergent_CDRs <- CDR3_seqs %>% 
  filter(total_compartment_seqs >= min_compartment_size, group_controls_pooled == 'secondary-56',
         tissue == 'LN') %>%
  filter(cdr3_seq_partis %in% unique(cdr3_seqs_driving_similarity$cdr3_seq)) %>%
  group_by(mouse_id, cell_type, v_gene) %>%
  count() %>%
  group_by(mouse_id, cell_type) %>%
  mutate(rank_gene_in_convergent_seqs = rank(-n, ties.method = 'first'),
         gene_freq_in_convergent_seqs = n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x = rank_gene_in_convergent_seqs, y = n, label = v_gene)) +
  geom_col() +
  facet_grid(cell_type~mouse_id, scales = 'free') +
  geom_text(angle = 45)
  
  group_by(mouse_id, day, infection_status, group_controls_pooled, cell_type, total_compartment_seqs) %>%
  mutate(v_gene_freq_in_convergent_seqs = n/sum(n))

  
# Some other stuff


# CDR3 length distributions by allele
CDR3_length_distributions_by_allele_NAIVE <- CDR3_seqs_naive %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, v_gene, cell_type, total_compartment_seqs) %>%
  summarise(naive_median_cdr3_length = median(cdr3_length),
            naive_cdr3_lowerq = quantile(cdr3_length, 0.25),
            naive_cdr3_upperq = quantile(cdr3_length, 0.75),
            n_naive_seqs = n()) %>%
  ungroup() %>%
  dplyr::rename(total_mouse_naive_seq = total_compartment_seqs) %>%
  select(-cell_type)

CDR3_length_distributions_by_allele_LN_EXPERIENCED <- CDR3_seqs_LN_experienced %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, v_gene, cell_type, total_compartment_seqs) %>%
  summarise(median_cdr3_length = median(cdr3_length),
            cdr3_lowerq = quantile(cdr3_length, 0.25),
            cdr3_upperq = quantile(cdr3_length, 0.75),
            n_seqs = n()) %>%
  ungroup()

CDR3_length_distributions_by_allele <- left_join(CDR3_length_distributions_by_allele_LN_EXPERIENCED,
                                                 CDR3_length_distributions_by_allele_NAIVE) %>%
  mutate(v_gene = factor(v_gene, levels = sort(as.character(unique(v_gene)))))

plot_CDR3_lengths_by_allele <- function(CDR3_length_distributions_by_allele, cell_type,
                                        group_controls_pooled){
  
  CDR3_length_distributions_by_allele %>%
    filter(total_compartment_seqs >= min_compartment_size,
           total_mouse_naive_seq >= min_compartment_size) %>%
    filter(group_controls_pooled == !!group_controls_pooled, cell_type == !!cell_type) %>%
    filter((cdr3_lowerq > naive_cdr3_upperq) | (cdr3_upperq < naive_cdr3_lowerq)) %>%
    ggplot(aes(x = v_gene)) +
    geom_linerange(aes(ymin = cdr3_lowerq, ymax = cdr3_upperq)) +
    geom_linerange(aes(ymin = naive_cdr3_lowerq, ymax = naive_cdr3_upperq),
                   color = 'blue') +
    geom_point(aes(y = median_cdr3_length)) +
    geom_point(aes(y = naive_median_cdr3_length), color = 'blue') +
    #geom_label(aes(y = 18, label = v_gene, angle = 90), size = 2) +
    facet_wrap('mouse_id', ncol = 1) +
    background_grid() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
          plot.title = element_text(hjust = 0.5, size = 14)) +
    xlab('V alleles with distinct CDR3 length\ndistributions relative to the naive repertoire') +
    ylab('CDR3 length') +
    ylim(5, 20)
  
}

CDR3_lengths_by_allele_day8 <- 
  plot_grid(plot_CDR3_lengths_by_allele(CDR3_length_distributions_by_allele, 'PC','primary-8') + xlab('') +
              ggtitle('Day 8 plasma cells'),
            plot_CDR3_lengths_by_allele(CDR3_length_distributions_by_allele, 'GC','primary-8') + ylab('') +
              ggtitle('Day 8 GC cells'),
            plot_CDR3_lengths_by_allele(CDR3_length_distributions_by_allele, 'mem','primary-8') + ylab('') +
              ggtitle('Day 8 memory cells') + xlab(''),
            nrow = 1
            )

save_plot('../figures/all_seqs_freqs/CDR3_lengths_by_allele_day8.pdf',
          CDR3_lengths_by_allele_day8,
          base_height = 10, base_width = 18)


# Enumerate all clones with, say, > 100 sequences.

# For clones with the same V gene, how often do they share the same mutation?

# For clones with the same V gene and a CDR3 of similar length, how similar are their CDR3s?

large_clones <- clone_freqs_by_tissue_and_cell_type %>%
  filter(compartment_tissue == 'LN', compartment_cell_type == 'PC') %>%
  filter(total_seqs_in_compartment >= min_compartment_size,
         n_clone_seqs_in_compartment >= 100) %>%
  filter(compartment_cell_type == 'PC', group_controls_pooled != 'control') %>%
  select(mouse_id, clone_id, v_gene, compartment_tissue, compartment_cell_type, clone_rank_in_compartment, n_clone_seqs_in_compartment, clone_freq_in_compartment) %>%
  arrange(desc(clone_freq_in_compartment))

large_clones <- left_join(large_clones, clone_info %>% select(mouse_id, clone_id,  cdr3_length_partis, clone_naive_cdr3_partis))




large_clones %>% group_by(v_gene, cdr3_length_partis) %>% 
  count() %>%
  ungroup() %>%
  arrange(v_gene, cdr3_length_partis, desc(n))

comparisons <- large_clones %>% group_by(v_gene, cdr3_length_partis) %>%
  summarise(compare_CDR3s_for_all_clone_pairs(cdr3_seqs = clone_naive_cdr3_partis, mouse_ids = mouse_id))

comparisons %>%
  ggplot(aes(x = v_gene, y = similarity)) +
  geom_boxplot(outlier.alpha = 0) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))



  






