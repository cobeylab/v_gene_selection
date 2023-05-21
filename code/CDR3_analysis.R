library(readr)
library(scales)

source('gene_frequency_functions.R')
source('plot_options.R')

theme_set(theme_cowplot())

min_compartment_size = 100

# Loading data
args <- commandArgs(trailingOnly = T)
assignment <- as.character(args[1])

# CDR3 analysis only implemented for these assignments
stopifnot(assignment %in% c('partis', 'partis_ogrdb'))

clone_info <- read_csv(paste0('../processed_data/clone_info_', assignment, '.csv'))

annotated_seqs <- read_csv(paste0('../processed_data/annotated_seqs_', assignment, '.csv'))
  
annotated_seqs <- annotated_seqs %>% mutate(seq_id = as.character(seq_id), clone_id = as.character(clone_id)) %>%
  filter(productive)

annotated_seqs <- get_info_from_mouse_id(annotated_seqs) %>%
  mutate(group_controls_pooled = factor(group_controls_pooled, levels = group_controls_pooled_factor_levels),
         cdr3_length = nchar(cdr3_seq_partis))

annotated_seqs <- left_join(annotated_seqs, clone_info %>% select(mouse_id, clone_id, v_gene) %>% mutate(clone_id = as.character(clone_id)))

# These analyses are based on all productive sequences (as opposed to unique sequences only)
load(paste0('../results/precomputed_gene_freqs_all_seqs_', assignment, '.RData'))


figure_output_dir = paste0('../figures/all_seqs_', assignment, '/exported_ggplot_objects/')
dir.create(figure_output_dir, recursive = T, showWarnings = F)


# ====== Some functions ======

# Defining amino acid classes following Hershberg & Schlomchik 2006, Hershberg & Saini 2015, after Chothia et al. 1998 
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
# Quick test
stopifnot(translate_aa_to_biochem_class('FLQRSP') == 'HHLLNN')


# Computes similarity of two sequences
compute_seq_similarity <- function(str1, str2){
  similarity <- mapply(x = str1, y = str2, 
                       FUN = function(x,y){
                         sum(str_split(x,'')[[1]] == str_split(y,'')[[1]]) / nchar(x)
                       })
  names(similarity) <- c()
  
  return(similarity)
}

# As a test, we should get all ones if passing a character vector of identical sequences
stopifnot(all(compute_seq_similarity(c('AAAA','AAAA','AAAA','AAAA'), c('AAAA','AAAA','AAAA','AAAA')) == 1))
# And all 0s if passing entirely different sequences
stopifnot(all(compute_seq_similarity(c('AAAA','BBBB','CCCC','DDDD'), c('XXXX','ZZZZ','YYYY','OOOO')) == 0))
stopifnot(compute_seq_similarity(c('AACC'),c('AADD')) == 0.5) #0.5 for 50% identical sequences


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
    
    n_compartment_seqs_mouse_i <- seqs_mouse_i %>%
      group_by(mouse_id, tissue, cell_type) %>%
      count() %>%
      pull(n)
    
    if(any(n_compartment_seqs_mouse_i > sample_size)){
      # Take a sample of sequences from mouse i for each tissue / cell type combination
      sample_mouse_i <- seqs_mouse_i %>%
        group_by(mouse_id, tissue, cell_type) %>%
        mutate(n_compartment_seqs = n()) %>%
        filter(n_compartment_seqs >= sample_size) %>%
        slice_sample(n = sample_size) %>%
        arrange(tissue, cell_type, cdr3_length) %>%
        ungroup() %>%
        select(-n_compartment_seqs)
      
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

# For each mouse, find all pairs of naive seqs
# (For CDR3 diversity by V gene analyses)
get_same_mouse_naive_seq_pairs <- function(annotated_seqs){
  mice <- unique(annotated_seqs$mouse_id)
  
  naive_seqs <- annotated_seqs %>% filter(cell_type == 'naive', productive == T)
  
  mouse_v_gene_combinations <- naive_seqs %>% select(mouse_id, v_gene) %>% unique()
  
  internal_function <- function(mouse, v_gene, naive_seqs){
    mouse_vgene_naive_seqs <- naive_seqs %>%
      filter(mouse_id == mouse, v_gene == !!v_gene)
    
    # Create object will all pairs of naive seqs
    naive_seq_pairs <- expand_grid(seq_id_i = mouse_vgene_naive_seqs$seq_id,
                                   seq_id_j = mouse_vgene_naive_seqs$seq_id) %>%
      filter(seq_id_i != seq_id_j)
    
    # Keep only pairs with the same CDR3 length
    naive_seq_pairs <- left_join(naive_seq_pairs,
                                 mouse_vgene_naive_seqs %>% select(seq_id, cdr3_length, cdr3_seq_partis) %>%
                                   dplyr::rename(seq_id_i = seq_id, cdr3_length_i = cdr3_length, cdr3_seq_i = cdr3_seq_partis),
                                 by = 'seq_id_i')
    
    naive_seq_pairs <- left_join(naive_seq_pairs,
                                 mouse_vgene_naive_seqs %>% select(seq_id, cdr3_length, cdr3_seq_partis) %>%
                                   dplyr::rename(seq_id_j = seq_id, cdr3_length_j = cdr3_length, cdr3_seq_j = cdr3_seq_partis),
                                 by = 'seq_id_j')
    
    naive_seq_pairs <- naive_seq_pairs %>%
      filter(cdr3_length_i == cdr3_length_j)
    
    naive_seq_pairs <- naive_seq_pairs %>%
      mutate(biochem_cdr3_seq_i = translate_aa_to_biochem_class(cdr3_seq_i),
             biochem_cdr3_seq_j = translate_aa_to_biochem_class(cdr3_seq_j))
    
    naive_seq_pairs <- naive_seq_pairs %>%
      mutate(mouse_id = mouse, v_gene = v_gene, cdr3_length = cdr3_length_i) %>%
      select(mouse_id, v_gene, cdr3_length, matches('_seq_'))

    return(naive_seq_pairs)
  }
  
  naive_seq_pairs <- mapply(FUN = internal_function,
                            mouse = mouse_v_gene_combinations$mouse_id,
                            v_gene = mouse_v_gene_combinations$v_gene,
                            MoreArgs = list(naive_seqs = naive_seqs),
                            SIMPLIFY = F)
  
  return(bind_rows(naive_seq_pairs))
  
}

# Computes sequence and biochemical similarity for paired CDR3 sequence samples produced by sample_sequence_pairs
compute_similarity_matched_samples <- function(matched_samples){
  if('mouse_pair' %in% names(matched_samples)){
    grouping_vars <- c('mouse_pair', 'cell_type')
  }else{
    stopifnot('mouse_id' %in% names(matched_samples))
    grouping_vars <- 'mouse_id'
  }
  
  matched_samples %>%
    mutate(cdr3_seq_similarity = compute_seq_similarity(str1 = cdr3_seq_i, str2 = cdr3_seq_j),
           cdr3_biochem_similarity = compute_seq_similarity(str1 = biochem_cdr3_seq_i, str2 = biochem_cdr3_seq_j)) %>%
    group_by(across(grouping_vars)) %>%
    mutate(n_seqs_in_comparison = n()) %>%
    ungroup() 
}

# Plots CDR3 sequence and biochemical similarity
plot_cdr3_similarity <- function(matched_samples_similarity){
  
  long_format_tibble <- matched_samples_similarity %>%
    pivot_longer(cols = c('cdr3_seq_similarity', 'cdr3_biochem_similarity'), names_to = 'similarity_type',
                 values_to = 'similarity') %>%
    mutate(similarity_type = case_when(
      similarity_type == 'cdr3_seq_similarity' ~ 'sequence similarity',
      similarity_type == 'cdr3_biochem_similarity' ~ 'biochemical similarity'
    )) %>%
    mutate(similarity_type = factor(similarity_type, levels = c('sequence similarity', 'biochemical similarity'))) %>%
    ungroup() %>%
    mutate(infection_status = str_extract(group_controls_pooled,'[A-Z|a-z]*')) %>%
    mutate(day = as.integer(as.character(day)))
  
  n_comparisons_labels <- long_format_tibble %>% group_by(group_controls_pooled, day, cell_type, similarity_type) %>%
    summarise(similarity = median(similarity), 
              total_seqs_in_comparison = n()) %>%
    mutate(infection_status = str_extract(group_controls_pooled,'[A-Z|a-z]*')) %>%
    mutate(day = as.integer(as.character(day)))

  pl <- long_format_tibble %>%
    mutate(day = as.integer(as.character(day))) %>%
    ggplot(aes(x = day, y = similarity)) +
    geom_boxplot(aes(color = infection_status, group = group_controls_pooled)) +
    #geom_label(data = n_comparisons_labels, aes(label = total_seqs_in_comparison)) +
    facet_grid(similarity_type~cell_type) +
    theme(legend.position = 'top', axis.text.x = element_text(size = 10)) +
    xlab('Days after primary infection') +
    ylab('Similarity of length-matched\nCDR3 sequences from different mice') +
    scale_x_continuous(breaks = unique(n_comparisons_labels$day)) +
    scale_color_manual(values = c('#d95f02','#7570b3'), name = 'Infection')
  return(pl)
}


# ====== Analyses =======

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
    cell_type == 'PC' ~ 'Lymph node plasma cells',
    cell_type == 'mem' ~ 'Lymph node memory cells'
  )) %>%
  mutate(cell_type = factor(cell_type,
                            levels = c('Naive cells (all tissues)', 'Lymph node GC cells',
                                       'Lymph node plasma cells', 'Lymph node memory cells')))


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


# ======== CDR3 SIMILARITY =======

# Sample length-matched CDR3 sequences from different mice, compute similarity, compare with length-matched CDR3 sequences from naive cells
sample_size = 500

length_matched_sample_NAIVE <- sample_sequence_pairs(pairwise_gene_freqs = pairwise_gene_freqs, annotated_seqs = annotated_seqs,
                                           sample_size = sample_size, selected_tissues = 'all', selected_cell_types = 'naive', match_v_allele = F)

length_matched_sample_LNEXP <- sample_sequence_pairs(pairwise_gene_freqs = pairwise_gene_freqs, annotated_seqs = annotated_seqs,
                                                     sample_size = sample_size, selected_tissues = 'LN', selected_cell_types = c('GC','PC','mem'), match_v_allele = F)

length_matched_CDR3_sample <-  bind_rows(length_matched_sample_NAIVE %>% select(-tissue),
                                    length_matched_sample_LNEXP %>% select(-tissue))  %>%
  mutate(biochem_cdr3_seq_i = translate_aa_to_biochem_class(cdr3_seq_i),
         biochem_cdr3_seq_j = translate_aa_to_biochem_class(cdr3_seq_j))


length_matched_CDR3_similarity <- compute_similarity_matched_samples(matched_samples = length_matched_CDR3_sample) %>%
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
  
length_matched_CDR3_similarity_plot <- plot_cdr3_similarity(length_matched_CDR3_similarity)


# Sample length-matched **AND ALLELE-MATCHED** CDR3 sequences from different mice, compute similarity,
# compare with length- and allele-matched CDR3 sequences from naive cells

length_and_allele_matched_sample_NAIVE <- sample_sequence_pairs(pairwise_gene_freqs = pairwise_gene_freqs, annotated_seqs = annotated_seqs,
                                                     sample_size = sample_size, selected_tissues = 'all', selected_cell_types = 'naive', match_v_allele = T)

length_and_allele_matched_sample_LNEXP <- sample_sequence_pairs(pairwise_gene_freqs = pairwise_gene_freqs, annotated_seqs = annotated_seqs,
                                                     sample_size = sample_size, selected_tissues = 'LN', selected_cell_types = c('GC','PC','mem'),
                                                     match_v_allele = T)
length_and_allele_matched_sample <- bind_rows(length_and_allele_matched_sample_NAIVE %>% select(-tissue),
                                              length_and_allele_matched_sample_LNEXP %>% select(-tissue))  %>%
  mutate(biochem_cdr3_seq_i = translate_aa_to_biochem_class(cdr3_seq_i),
         biochem_cdr3_seq_j = translate_aa_to_biochem_class(cdr3_seq_j))

length_and_allele_matched_CDR3_similarity <-
  compute_similarity_matched_samples(matched_samples = length_and_allele_matched_sample)  %>%
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
  
length_and_allele_matched_CDR3_similarity_plot <- plot_cdr3_similarity(length_and_allele_matched_CDR3_similarity)  +
  ylab('Similarity of length- and allele-matched\nCDR3 sequences from different mice')

# ====== SIMILARITY OF CDR3 SEQUENCES BETWEEN NAIVE SEQS FROM THE SAME MOUSE =====

same_mouse_naive_seq_pairs <- get_same_mouse_naive_seq_pairs(annotated_seqs)

CDR3_similarity_NAIVE_same_mouse <- compute_similarity_matched_samples(matched_samples = same_mouse_naive_seq_pairs)

CDR3_similarity_NAIVE_same_mouse <- CDR3_similarity_NAIVE_same_mouse %>%
  select(mouse_id, v_gene, cdr3_seq_similarity, cdr3_biochem_similarity) %>%
  pivot_longer(cols = c('cdr3_seq_similarity', 'cdr3_biochem_similarity'),
               names_to = 'type',
               values_to = 'similarity') %>%
  mutate(dissimilarity = 1 - similarity)

focal_alleles <- c('IGHV14-4*01','IGHV1-82*01','IGHV1-69*01', 'IGHV14-1*01','IGHV14-2*01')

CDR3_similarity_NAIVE_same_mouse_pl <- CDR3_similarity_NAIVE_same_mouse %>%
  mutate(type = case_when(
    type == 'cdr3_biochem_similarity' ~ 'Biochemical dissimilarity',
    type == 'cdr3_seq_similarity' ~ 'Sequence dissimilarity'
  )) %>%
  mutate(focal_allele = v_gene %in% focal_alleles) %>%
  ggplot(aes(x = v_gene, y = dissimilarity, color = focal_allele)) +
  scale_color_manual(values = c('black','red')) +
  geom_boxplot() +
  facet_wrap('type', nrow = 2) +
  theme(legend.position = 'None')

# Export plots 

save(CDR3_length_distribution_pl,
     length_matched_CDR3_similarity_plot,
     length_and_allele_matched_CDR3_similarity_plot,
     CDR3_similarity_NAIVE_same_mouse_pl,
     file = paste0(figure_output_dir, 'CDR3_plots.RData'))
