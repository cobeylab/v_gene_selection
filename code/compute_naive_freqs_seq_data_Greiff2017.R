source('gene_frequency_functions.R')
library(readr)
theme_set(theme_cowplot())

# Read annotated sequences from Greiff et al. 2017

args <- commandArgs(trailingOnly = T)

assignment <- as.character(args[1])

stopifnot(assignment %in% c('partis','partis_ogrdb'))

if(assignment == 'partis'){
  annotated_seq_files_Greiff2017 <- list.files('../results/partis/seq_data_Greiff2017/', pattern = 'annotated_seqs.csv',
                                               full.names = T)
  annotated_seq_files_Greiff2017 <- annotated_seq_files_Greiff2017[str_detect(annotated_seq_files_Greiff2017, 'ogrdb', negate = T)]
  
  # Read pre-computed naive frequencies from our data
  load('../results/precomputed_gene_freqs_all_seqs_partis.RData')
  
  output_file = '../processed_data/naive_freqs_Greiff2017_partis.csv'
  
  
  
}else{
  annotated_seq_files_Greiff2017 <- list.files('../results/partis/seq_data_Greiff2017/', pattern = 'ogrdb_annotated_seqs.csv',
                                               full.names = T)
  
  load('../results/precomputed_gene_freqs_all_seqs_partis_ogrdb.RData')
  
  output_file = '../processed_data/naive_freqs_Greiff2017_partis_ogrdb.csv'
}


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

naive_freqs_Greiff2017 <- annotated_seqs_Greiff2017 %>%
  group_by(mouse_id, v_segment_partis) %>%
  dplyr::summarise(n_vgene_seqs = dplyr::n()) %>%
  group_by(mouse_id) %>%
  mutate(total_compartment_seqs = sum(n_vgene_seqs),
         vgene_seq_freq = n_vgene_seqs / total_compartment_seqs) %>%
  dplyr::rename(v_gene = v_segment_partis) %>%
  ungroup()


# Pre-computed naive frequencies from our data
naive_freqs_this_study <- naive_freqs %>%
  dplyr::rename(n_vgene_seqs = n_naive_vgene_seqs,
                vgene_seq_freq = naive_vgene_seq_freq,
                total_compartment_seqs = total_mouse_naive_seqs)


# ================= How do their V gene frequencies correlate with those from our experiments? =================

# Generates mouse pairs from the Greiff et al. 2017 dataset
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

# Generates pairs with one mouse from each study.
get_cross_dataset_mouse_pairs <- function(naive_freqs_Greiff2017, naive_freqs_this_study){
  base_function <- function(mouse_id_Greiff2017){
    
    total_Greiff2017_mouse_naive_seqs <- unique(naive_freqs_Greiff2017 %>% filter(mouse_id == mouse_id_Greiff2017)
                                                %>% pull(total_compartment_seqs))
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

# Mean and interquartile range of cross-dataset Spearman correlations

pairwise_naive_freq_correlations %>% filter(pair_type == 'cross-dataset') %>%
  summarise(mean_cross_dataset_correlation = mean(cor_coef_freqs),
            lower_quartile = quantile(cor_coef_freqs,0.25),
            upper_quartile = quantile(cor_coef_freqs, 0.75))

# Export Greiff et al. 2017 naive freqs.
write_csv(file = output_file,
          naive_freqs_Greiff2017)

