library(tidyr)
library(dplyr)
library(cowplot)
library(stringr)
library(viridis)
library(reshape2)

# Remove mouse/cell-type combinations with fewer than <min_seqs> sequences and mice with fewer than <min_seqs> naive seqs.
min_seqs <- 1000 # Passed down to other scripts that use the functions defined here.

cell_type_list <- list('GC','PC','mem','experienced')
names(cell_type_list) <- c('GC','PC','mem','experienced')

# Order of mice, for plotting...
mouse_id_factor_levels = expand.grid(id = 1:10, day = c(8,16,24,40,56))%>%
  mutate(factor_levels = paste(day,id, sep = '-')) %>% pull(factor_levels)

# Order of groups, for plotting
group_factor_levels <- c('control-8','primary-8','control-16','primary-16',
                         'control-24','primary-24','control-40','secondary-40',
                         'control-56','secondary-56')

# Order with groups (with controls combined), for plotting
group_controls_pooled_factor_levels <- c('control','primary-8','primary-16',
  'primary-24','secondary-40','secondary-56')


get_info_from_mouse_id <- function(data){
  data %>% mutate(day = str_extract(mouse_id,'[0-9]*[^-]'),
           mouse_id_number = str_replace(str_extract(mouse_id,'-[0-9]*'),'-',''),
           infection_status = assign_infection_status(day, mouse_id_number),
           day = factor(day, levels = sort(unique(as.numeric(day))))) %>%
    mutate(group = paste(infection_status, day, sep = '-')) %>%
    mutate(group_controls_pooled = group) %>%
    mutate(group_controls_pooled = ifelse(grepl('control',group_controls_pooled), 'control', group_controls_pooled)) %>%
    select(-mouse_id_number) %>%
    select(mouse_id, day, infection_status, group, group_controls_pooled, everything())
}


# Determines infection status based on mouse id (treats all post day 40 non-controls as secondary)
assign_infection_status <- function(day, mouse_id_number){
  day <- as.numeric(as.character(day)) # Converting factor to numeric through character
  mouse_id_number <- as.numeric(as.character(mouse_id_number))
  infection_status <- ifelse(mouse_id_number %in% c(5,10), 'control','infected')
  infection_status[infection_status == 'infected'] <- ifelse(
    day[infection_status == 'infected'] < 40, 'primary', 'secondary'
  )
  return(infection_status)
}


# compile_clone_info <- function(clone_info_dir){
#   # Read clone files for all mice 
#   # Combined tibble has one row per cell subpopulation in a clone in tissue in a mouse
#   # (e.g. memory B cells from clone X in the spleen of mouse Y)
#   mouse_clone_files = list.files(clone_info_dir, pattern = '_clones.csv')
#   mouse_clone_files = paste0(clone_info_dir, mouse_clone_files)
#   clone_info <- lapply(mouse_clone_files, FUN = function(x){return(as_tibble(read.csv(x)))})
#   clone_info <- bind_rows(clone_info)
# 
#   #lapply(mouse_clone_files, FUN = function(x){return(as_tibble(read.csv(x)) %>%
#   #                                                     filter(is.na(clone_id)))})
#   #lapply(mouse_clone_files, FUN = function(x){return(  unique(as_tibble(read.csv(x))$mouse_id)  )})
#   
#   # Add variables for observation day, infection status
#   clone_info <- get_info_from_mouse_id(clone_info %>% filter(is.na(clone_id) == F))
#   
#   clone_info <- clone_info %>% mutate(mouse_id = factor(mouse_id, levels = mouse_id_factor_levels)) %>%
#     mutate(cell_type = factor(cell_type, levels = c('naive', 'GC','PC','mem')))
#   return(clone_info)
# }

get_clone_purity <- function(unique_seq_counts){
  clone_purity <- unique_seq_counts %>% 
    # Group each row as representing a naive or non-naive cell subpopulation
    mutate(cell_group = ifelse(cell_type == 'naive', 'naive', 'non_naive')) %>%
    group_by(mouse_id, clone_id, cell_group) %>%
    # Count number of unique productive sequences from naive and non-naive cells in each clone
    summarise(uniq_prod_seqs = sum(uniq_prod_seqs)) %>%
    ungroup() %>%
    pivot_wider(names_from = 'cell_group', values_from = 'uniq_prod_seqs', values_fill = 0) %>%
   # spread(key = "cell_group", value = 'uniq_prod_seqs') %>%
    dplyr::rename(naive_seqs_in_clone = naive, non_naive_seqs_in_clone = non_naive) %>%
    mutate(is_pure_naive = naive_seqs_in_clone > 0 & non_naive_seqs_in_clone ==0,
           is_pure_non_naive = naive_seqs_in_clone == 0 & non_naive_seqs_in_clone > 0) %>%
    mutate(
      clone_purity = case_when(
        is_pure_naive ~ 'pure_naive',
        is_pure_non_naive ~ 'pure_non_naive',
        !(is_pure_naive) & !(is_pure_non_naive) ~ 'mixed'
      )
    ) %>%
    select(-is_pure_naive, -is_pure_non_naive)
  return(clone_purity)
} 

# Retain clone info in original format but only with lines from pure clones
retain_pure_clones <- function(clone_info){
  if(all(c('is_pure_naive','is_pure_non_naive') %in% names(clone_info)) == F){
    clone_info <- get_clone_purity(clone_info)
  }

  pure_clone_info <- clone_info %>%
    filter(is_pure_naive | is_pure_non_naive) %>%
    mutate(clone_type = ifelse(is_pure_naive,'naive','experienced')) %>%
    select(-is_pure_naive, -is_pure_non_naive)
  
  return(pure_clone_info %>% select(mouse_id, day, infection_status, clone_id, clone_type, everything()))
}


# Calculates gene frequencies for each cell type in each mouse (pools across tissues)
calc_gene_freqs <- function(unique_seq_counts, long_format = F, by_tissue = F, tissue_subset = NULL){
  
  grouping_vars <- c('mouse_id', 'v_gene', 'cell_type')
  
  if(by_tissue){
    grouping_vars <- c(grouping_vars, 'tissue')
  }
  
  # All mouse-specific columns other than ID are removed from the grouping for convenience, but they're added back at the end
  mouse_info <- unique_seq_counts %>% select(mouse_id, day, infection_status, group, group_controls_pooled) %>%
    unique()
  
  # Tibble with all genes present in each mice at least once, with rows for that gene
  # in all possible tissue /cell type combinations
  if(by_tissue){
    all_genes_in_mice <- unique_seq_counts %>% select(mouse_id, v_gene, tissue, cell_type) %>%
      complete(nesting(mouse_id, v_gene), tissue, cell_type) %>% unique()
  }else{
    all_genes_in_mice <- unique_seq_counts %>% select(mouse_id, v_gene, cell_type) %>%
      complete(nesting(mouse_id, v_gene), cell_type) %>% unique()
  }
  
  # If a character vector is provided as tissue_subset, calculate frequencies using those tissues only 
  if(is.null(tissue_subset) == F){
    base_data <- unique_seq_counts %>% filter(tissue %in% tissue_subset)
    if(by_tissue){
      all_genes_in_mice <- all_genes_in_mice %>% filter(tissue %in% tissue_subset)
    }
  }else{
    base_data <- unique_seq_counts
  }
  
  # Calculate frequencies in each subcompartment (naive, GC, PC, mem)
  gene_seq_freqs <- base_data %>%
    mutate(mouse_id = factor(mouse_id), v_gene = factor(v_gene), cell_type = factor(cell_type)) %>%
    group_by(across(grouping_vars)) %>%
    summarise(n_vgene_seqs = sum(uniq_prod_seqs, na.rm=T)) %>%
    ungroup() 
   
  # Left join to all_genes_in_mice tibble, so that genes present in mice get an explicit zero value
  # in compartments they're missing from
  gene_seq_freqs <- left_join(all_genes_in_mice, gene_seq_freqs) %>%
    replace_na(list(n_vgene_seqs =  0))
  
  gene_seq_freqs <- gene_seq_freqs %>%
    # Calculate mouse-cell-type (optionally tissue) totals
    group_by(across(grouping_vars[grouping_vars != 'v_gene'])) %>%
    mutate(total_mouse_cell_type_seqs = sum(n_vgene_seqs)) %>%
    ungroup() %>%
    # Exclude compartments that are entirely zero
    filter(total_mouse_cell_type_seqs > 0) %>%
    ungroup()
  
  # Calculate frequencies in the experienced repertoire as a whole, across mem/PC/GC
  gene_seq_freqs_all_exp <- gene_seq_freqs %>%
    filter(cell_type != 'naive') %>%
    group_by(across(grouping_vars[grouping_vars != 'cell_type'])) %>%
    summarise(n_vgene_seqs = sum(n_vgene_seqs)) %>%
    mutate(cell_type = 'experienced') %>%
    ungroup() %>%
    group_by(across(grouping_vars[!(grouping_vars %in% c('v_gene','cell_type'))])) %>%
    mutate(total_mouse_cell_type_seqs = sum(n_vgene_seqs)) %>%
    ungroup() %>%
    select(mouse_id, v_gene, cell_type, everything())
  
  # Tibble with experienced frequencies
  exp_freqs <- bind_rows(gene_seq_freqs %>% filter(cell_type != 'naive'),
                         gene_seq_freqs_all_exp) %>%
    arrange(mouse_id, v_gene, cell_type) %>%
    mutate(vgene_seq_freq = n_vgene_seqs/total_mouse_cell_type_seqs)
  # Add complete mouse_information
  exp_freqs <- left_join(exp_freqs, mouse_info, by = 'mouse_id') %>%
    select(mouse_id, day, infection_status, group, group_controls_pooled, everything())

  # Tibble with naive frequencies
  naive_freqs <- gene_seq_freqs %>% filter(cell_type == 'naive') %>%
    select(-cell_type) %>%
    dplyr::rename(n_naive_vgene_seqs = n_vgene_seqs, total_mouse_naive_seqs = total_mouse_cell_type_seqs) %>%
    mutate(naive_vgene_seq_freq = n_naive_vgene_seqs / total_mouse_naive_seqs)
  naive_freqs <- left_join(naive_freqs, mouse_info, by = 'mouse_id') %>%
    select(mouse_id, day, infection_status, group, group_controls_pooled, everything())
  
  if(long_format == T){
    output <- bind_rows(exp_freqs,
                        naive_freqs %>%
                          rename(n_vgene_seqs = n_naive_vgene_seqs, vgene_seq_freq = naive_vgene_seq_freq,
                                 total_mouse_cell_type_seqs = total_mouse_naive_seqs) %>%
                          mutate(cell_type = 'naive')) %>% arrange(mouse_id,v_gene) %>%
      mutate(cell_type = factor(cell_type, levels = c('naive', 'experienced','GC','mem','PC')))
  }else{
    output <- list(exp_freqs = exp_freqs, naive_freqs = naive_freqs)
  }

  return(output)
}

# Resample experienced frequencies 
# I.e. draw from empirical experienced frequencies in each compartment (cell type), re-estimate frequencies from draw
resample_exp_freqs <- function(exp_freqs){
  resampled_exp_freqs <- exp_freqs %>% group_by(mouse_id, cell_type) %>%
    mutate(resampled_exp_count = rmultinom(1, size = unique(total_mouse_cell_type_seqs), prob = vgene_seq_freq)) %>%
    mutate(n_vgene_seqs = resampled_exp_count, vgene_seq_freq = n_vgene_seqs / total_mouse_cell_type_seqs) %>%
    ungroup() %>%
    select(-resampled_exp_count)
  return(resampled_exp_freqs)
}

# Resample naive frequencies
# I.e. draw from empirical naive frequencies, re-estimate frequencies from draw
resample_naive_freqs <- function(naive_freqs){
  
  # This if statement handles the case where the input is a tibble with both experienced freqs and naive freqs, 
  # (in which naive freqs are repeated for each experienced cell type)
  if('cell_type' %in% names(naive_freqs)){
    stopifnot(length(unique(naive_freqs$cell_type)) == 1)
  }
  
  grouping_vars <- 'mouse_id'
  if('tissue' %in% names(naive_freqs)){
    grouping_vars <- c(grouping_vars, 'tissue')
  }

  resampled_naive_freqs <- naive_freqs %>% group_by(across(grouping_vars)) %>%
    mutate(resampled_naive_count = rmultinom(1, size = unique(total_mouse_naive_seqs), prob = naive_vgene_seq_freq)) %>%
    mutate(n_naive_vgene_seqs = resampled_naive_count, naive_vgene_seq_freq = n_naive_vgene_seqs / total_mouse_naive_seqs) %>%
    ungroup() %>%
    select(-resampled_naive_count)
  
  return(resampled_naive_freqs)
}

adjust_zero_naive_freqs <- function(naive_freqs){
  # Assume genes with zero frequency in naive repertoire have a frequency of 1/total number of naive seqs.
  # Because of the way the input is produced, this only applies to genes present in some other cell subset in the mouse
  
  # This if statement handles the case where the input is a tibble with both experienced freqs and naive freqs, 
  # (in which naive freqs are repeated for each experienced cell type)
  if('cell_type' %in% names(naive_freqs)){
    stopifnot(length(unique(naive_freqs$cell_type)) == 1)
  }
  
  grouping_vars <- 'mouse_id'
  if('tissue' %in% names(naive_freqs)){
    grouping_vars <- c(grouping_vars, 'tissue')
  }
  
  adj_naive_freqs <- naive_freqs %>%
    mutate(n_naive_vgene_seqs = n_naive_vgene_seqs + 1 * (naive_vgene_seq_freq == 0),
           naive_vgene_seq_freq = naive_vgene_seq_freq + 1/total_mouse_naive_seqs * (naive_vgene_seq_freq == 0)) %>%
    # Re-normalize frequencies so they sum to 1
    group_by(across(grouping_vars)) %>%
    mutate(naive_vgene_seq_freq = naive_vgene_seq_freq / sum(naive_vgene_seq_freq)) %>%
    mutate(total_mouse_naive_seqs = sum(n_naive_vgene_seqs)) %>%
    ungroup()
  return(adj_naive_freqs)
}

adjust_zero_exp_freqs <- function(exp_freqs){
  # Assume genes with zero frequency in the experienced repertoire have a frequency of 1/total number of exp seqs in the corresponding compartment
  # (i.e., 'experienced', 'GC', 'mem', 'PC' compartments)
  # Because of the way the input is produced, this only applies to genes present in some other compartment of the mouse (incl. possibly naive)
  adj_exp_freqs <- exp_freqs %>%
    mutate(n_vgene_seqs = n_vgene_seqs + 1 * (vgene_seq_freq == 0),
           vgene_seq_freq = vgene_seq_freq + 1/total_mouse_cell_type_seqs * (vgene_seq_freq == 0)) %>%
    # Re-normalize frequencies so they sum to 1
    group_by(mouse_id, cell_type) %>%
    mutate(vgene_seq_freq = vgene_seq_freq / sum(vgene_seq_freq)) %>%
    mutate(total_mouse_cell_type_seqs = sum(n_vgene_seqs)) %>%
    ungroup()
  
  return(adj_exp_freqs)
}

# Creates matrix with the frequency ratios of all pairs of genes
calculate_frequency_ratio_matrix <- function(gene_freqs){
  # Internal function to construct the matrix for a single cell type for a single mouse
  internal_matrix_function <- function(gene_freqs, cell_type, mouse_id){
    sub_tibble <- gene_freqs %>% filter(cell_type == !!cell_type, mouse_id == !!mouse_id)
    
    # All pairs for a mouse
    long_form_matrix <- t(combn(as.character(unique(sub_tibble$v_gene)),2))
    # This ensures genes in each pair will be sorted alphabetically
    long_form_matrix <- t(apply(long_form_matrix, 1, sort))
    long_form_matrix <- tibble(gene_i = long_form_matrix[,1], gene_j = long_form_matrix[,2])
    
    long_form_matrix <- left_join(long_form_matrix,
                                  sub_tibble %>%
                                    dplyr::rename(gene_i = v_gene, vgene_seq_freq_i = vgene_seq_freq) %>%
                                    select(gene_i, vgene_seq_freq_i), by = 'gene_i') 
    long_form_matrix <- left_join(long_form_matrix,
                                  sub_tibble %>%
                                    dplyr::rename(gene_j = v_gene, vgene_seq_freq_j = vgene_seq_freq) %>%
                                    select(gene_j, vgene_seq_freq_j), by = 'gene_j')
    long_form_matrix <- long_form_matrix %>% 
      mutate(log_freq_ratio = log(vgene_seq_freq_i) - log(vgene_seq_freq_j)) %>%
      mutate(mouse_id = mouse_id, cell_type = cell_type) %>%
      select(mouse_id, cell_type, everything())
    
    return(long_form_matrix)
  }
  
  mouse_cell_type_combinations <- gene_freqs %>%
    select(mouse_id, cell_type) %>% unique()
  
  frequency_ratio_matrix <- mapply(internal_matrix_function,
         cell_type = mouse_cell_type_combinations$cell_type,
         mouse_id = mouse_cell_type_combinations$mouse_id,
         MoreArgs = list(gene_freqs = gene_freqs), SIMPLIFY = F)
  frequency_ratio_matrix <- bind_rows(frequency_ratio_matrix)
  
  # add mouse-specific info
  frequency_ratio_matrix <- left_join(frequency_ratio_matrix, 
            gene_freqs %>% select(mouse_id, day, infection_status, group_controls_pooled) %>%
              unique(),
            by = 'mouse_id') %>%
    select(mouse_id, day, infection_status, group_controls_pooled, everything())
  
  return(frequency_ratio_matrix)
}

# For all pairs of mice, rearrange freqs tibble to show changes in each mouse as different vars.
get_pairwise_freqs <- function(unique_seq_counts, by_tissue = F, adjust_naive_zeros, naive_from_tissue = NULL){
  
  mouse_info <- unique_seq_counts %>%
    select(mouse_id, day, infection_status, group, group_controls_pooled) %>%
    unique()
  
  unique_pairs <- unique_seq_counts %>% select(mouse_id) %>% unique() %>%
    dplyr::rename(mouse_id_i = mouse_id) %>%
    mutate(mouse_id_j = mouse_id_i) %>%
    complete(mouse_id_i, mouse_id_j) %>%
    rowwise() %>%
    mutate(pair = paste0(sort(c(mouse_id_i, mouse_id_j)), collapse = ';')) %>%
    ungroup() %>%
    filter(mouse_id_i != mouse_id_j) %>%
    select(pair) %>%
    unique() %>% pull(pair)
  
  # Calculates naive sequences pooled across a specific set of tissues (by_tissue = F on purpose)
  naive_freqs <- calc_gene_freqs(unique_seq_counts, long_format = T, by_tissue = F, tissue_subset = naive_from_tissue) %>%
    filter(cell_type == 'naive') %>% mutate(exp_naive_ratio = NA)
  
  if(!is.null(naive_from_tissue)){
    naive_freqs <- naive_freqs %>% mutate(tissue = paste(naive_from_tissue, collapse = '+'))
  }else{
    naive_freqs <- naive_freqs %>% mutate(tissue = 'all tissues')
  }
  
  if(adjust_naive_zeros){
    naive_freqs <- adjust_zero_naive_freqs(naive_freqs %>% 
                              dplyr::rename(naive_vgene_seq_freq = vgene_seq_freq,
                                            n_naive_vgene_seqs = n_vgene_seqs,
                                            total_mouse_naive_seqs = total_mouse_cell_type_seqs)) %>%
      dplyr::rename(vgene_seq_freq = naive_vgene_seq_freq,
                    n_vgene_seqs = n_naive_vgene_seqs,
                    total_mouse_cell_type_seqs = total_mouse_naive_seqs)
  }
  
  
  # Calculate experienced frequencies (tissue specific if by_tissue is T)
  exp_freqs <- calc_gene_freqs(unique_seq_counts, long_format = F, by_tissue = by_tissue, tissue_subset = NULL)$exp_freqs
  
  # For experienced cell compartments, calculate rho (ratio between obs and naive freqs)
  exp_freqs <- left_join(exp_freqs,
            naive_freqs %>% select(mouse_id, v_gene, vgene_seq_freq) %>% 
              dplyr::rename(naive_vgene_freq = vgene_seq_freq)) %>%
    mutate(log_exp_naive_ratio = log(vgene_seq_freq) - log(naive_vgene_freq),
           exp_naive_ratio = exp(log_exp_naive_ratio)) %>%
    select(-naive_vgene_freq, -log_exp_naive_ratio)
  
  gene_freqs <- bind_rows(exp_freqs, naive_freqs)
  
  compartment_sizes <- gene_freqs %>%
    select(mouse_id, cell_type, matches('tissue'), total_mouse_cell_type_seqs) %>%
    unique()
  
  internal_function <- function(mouse_pair, gene_freqs){
    # Parse mouse_pair, determine type of pair (e.g. 'primary', 'primary/control', 'secondary/control')
    mice <- str_split(mouse_pair,';')[[1]]
    
    days <- sapply(mice, FUN = function(x){str_split(x,'-')[[1]][1]})
    mouse_id_numbers <- sapply(mice, FUN = function(x){str_split(x,'-')[[1]][2]})
    infection_status = assign_infection_status(days, mouse_id_numbers)
    

    pair_type = ifelse(length(unique(infection_status)) == 1,
                       infection_status,
                       paste(sort(unique(infection_status)),collapse = '/'))
    time_i <- days[1]
    time_j <- days[2]
    pair_time <- paste(sort(as.numeric(days)), collapse=';')
    
    mouse1_freqs <- gene_freqs %>% filter(mouse_id == mice[1]) %>% 
      select(-day, -infection_status, -group, -group_controls_pooled, - total_mouse_cell_type_seqs)
    mouse2_freqs <- gene_freqs %>% filter(mouse_id == mice[2]) %>%
      select(-day, -infection_status, -group, -group_controls_pooled, -total_mouse_cell_type_seqs)
    
    for(var in colnames(mouse1_freqs)){
      if((var %in% c('v_gene','cell_type','tissue')) == F){
        colnames(mouse1_freqs)[colnames(mouse1_freqs) == var] <- paste0(var,'_i')
        colnames(mouse2_freqs)[colnames(mouse2_freqs) == var] <- paste0(var,'_j')
      }
    }
    
    # Assemble pair tibble
    pair_changes <- full_join(mouse1_freqs, mouse2_freqs) %>%
      mutate(mouse_pair = mouse_pair,
             mouse_id_i = mice[1], mouse_id_j = mice[2],
             pair_type = factor(pair_type, levels = c('control','control/primary','primary',
                                                      'control/secondary','secondary','primary/secondary'))) %>%
      select(mouse_pair, pair_type, everything())  
      #%>% filter(!is.na(vgene_seq_freq_i), !is.na(vgene_seq_freq_j))
    
    return(pair_changes)
  }
  # Paired changes in gene frequences
  paired_gene_freqs <- lapply(unique_pairs, FUN = internal_function,
                                     gene_freqs = gene_freqs)
  paired_gene_freqs <- bind_rows(paired_gene_freqs) 
  
  # Add mouse information and compartment sizes (n unique seqs.)
  
  paired_gene_freqs <- left_join(paired_gene_freqs,
                                 mouse_info %>% dplyr::rename_with(.fn = function(x){paste0(x,'_i')}, .cols = everything()))
  paired_gene_freqs <- left_join(paired_gene_freqs,
                                 mouse_info %>% dplyr::rename_with(.fn = function(x){paste0(x,'_j')}, .cols = everything())) 
  paired_gene_freqs <- left_join(paired_gene_freqs,
                                 compartment_sizes %>% dplyr::rename_with(.fn = function(x){paste0(x,'_i')},
                                                                          .cols = c(mouse_id, total_mouse_cell_type_seqs)))
  paired_gene_freqs <- left_join(paired_gene_freqs,
                                 compartment_sizes %>% dplyr::rename_with(.fn = function(x){paste0(x,'_j')},
                                                                          .cols = c(mouse_id, total_mouse_cell_type_seqs))) 
  
  naive_compartment_sizes <- compartment_sizes %>% filter(cell_type == 'naive') %>%
    dplyr::rename(total_mouse_naive_seqs = total_mouse_cell_type_seqs) %>% select(-tissue, -cell_type)
    
  
  paired_gene_freqs <- left_join(paired_gene_freqs,
                                 naive_compartment_sizes %>% dplyr::rename_with(.fn = function(x){paste0(x,'_i')},
                                                                                .cols =  everything()))
  
  paired_gene_freqs <- left_join(paired_gene_freqs,
                                 naive_compartment_sizes %>% dplyr::rename_with(.fn = function(x){paste0(x,'_j')},
                                                                                .cols =  everything())) %>%
    select(mouse_pair, pair_type, matches('id'), matches('day'), matches('infection_status'), matches('group'),
           matches('total_mouse_cell_type_seqs'), matches('total_mouse_naive_seqs'),
           everything())
  
  
  return(paired_gene_freqs)
}

# Calculates pairwise correlations v gene freqs. or v genes experienced-to-naive freq.ratios between pairs of mice.
get_pairwise_correlations <- function(pairwise_gene_freqs, min_genes_in_comparison = 10){
  
  grouping_vars <- c('mouse_pair','pair_type','mouse_id_i','mouse_id_j','day_i','day_j',
                     'cell_type', 'total_mouse_cell_type_seqs_i','total_mouse_cell_type_seqs_j',
                     'total_mouse_naive_seqs_i', 'total_mouse_naive_seqs_j')
  if('tissue' %in% names(pairwise_gene_freqs)){
    grouping_vars <- c(grouping_vars, 'tissue')
  }
  
  pairwise_correlations <- pairwise_gene_freqs %>% 
    filter(mouse_id_i != mouse_id_j) %>%
    group_by(across(grouping_vars)) %>%
    mutate(n_genes_in_freqs_comparison = sum(!is.na(n_vgene_seqs_i) & !is.na(n_vgene_seqs_j)),
           n_genes_in_freq_ratio_comparison = sum(!is.na(exp_naive_ratio_i) & !is.na(exp_naive_ratio_j)))
  
  pairwise_correlations_freqs <- pairwise_correlations %>%
    filter(n_genes_in_freqs_comparison >= min_genes_in_comparison) %>%
    summarise(cor_coef_freqs = cor.test(vgene_seq_freq_i, vgene_seq_freq_j,
                                  method = 'spearman')$estimate) %>%
    ungroup()
  
  pairwise_correlations_freq_ratios <- pairwise_correlations %>%
    filter(n_genes_in_freq_ratio_comparison >= min_genes_in_comparison) %>%
    summarise(cor_coef_freq_ratios = cor.test(exp_naive_ratio_i, exp_naive_ratio_j,
                                  method = 'spearman')$estimate) %>%
    ungroup()
  return(list(freqs = pairwise_correlations_freqs, freq_ratios = pairwise_correlations_freq_ratios))
} 


# Generates datasets with selected genes specified for each mouse in synth_data_input_tibble
# If synth_data_input_tibble is 'neutral', assumes neutrality (i.e. genes in exp. repertoire are sampled based on naive frequencies)
simulate_selection_freq_changes <- function(unique_seq_counts, synth_data_input_tibble, min_seqs, extra_mice = 0, by_tissue = F,
                                            naive_from_tissue = NULL, n_reps = 100){
  
  # If extra_mice = 1, duplicate each mouse before simulation, if = 2, make two copies of each mouse before simulation, etc.
  base_data <- unique_seq_counts
  
  # Generate synthetic data with more mice by replicating existing mice
  if(extra_mice > 0){
    for(i in 1:extra_mice){
      duplicate_data <- base_data
      duplicate_data$mouse_id <- paste0(duplicate_data$mouse_id, '_copy_',i + 1)
      base_data <- bind_rows(base_data, duplicate_data)
    }
  }
  
  # Get gene frequencies
  exp_freqs <- calc_gene_freqs(base_data, long_format = F, by_tissue = by_tissue)$exp_freqs
  
  # Naive frequencies can be across all tissues or from a specific set of tissues
  # (in which case they're used for all tissues if experienced reps are simulated for each tissue)
  
  naive_freqs <- calc_gene_freqs(base_data, long_format = F, by_tissue = F,
                                 tissue_subset = naive_from_tissue)$naive_freqs
  
  # Set genes present in mice but with 0 naive seqs to have 1 naive seq, adjust frequencies accordingly
  naive_freqs <- adjust_zero_naive_freqs(naive_freqs)
  
  # Internal function for generating a single realization
  generate_replicate <- function(naive_freqs, exp_freqs){
    # To account for uncertainty in naive frequencies, sample naive counts, then re-calculate naive frequencies
    naive_freqs <- resample_naive_freqs(naive_freqs)
    
    # Re-adjust naive zeros after resampling
    naive_freqs <- adjust_zero_naive_freqs(naive_freqs)
    
    # Merge naive frequencies with observed total cell counts
    simulated_freqs <- left_join(exp_freqs,
                                 naive_freqs %>% select(mouse_id, v_gene, n_naive_vgene_seqs,naive_vgene_seq_freq, total_mouse_naive_seqs),
                                 by = c('mouse_id','v_gene'))  %>%
      # Remove mouse/cell-type combinations with fewer than <min_seqs> sequences and mice with fewer than <min_seqs> naive seqs
      filter(total_mouse_cell_type_seqs >= min_seqs, total_mouse_naive_seqs >= min_seqs) 
    
    if(is.character(synth_data_input_tibble)){
      stopifnot(synth_data_input_tibble == 'neutral')
      simulated_freqs <- simulated_freqs %>% mutate(fitness = 1) # In neutral scenario, set fitness to arbitrary value
    }else{
      simulated_freqs <- left_join(simulated_freqs, synth_data_input_tibble) 
    }
    
    simulated_freqs <- simulated_freqs %>%
      # Set simulated experienced freq to naive freq * fitness
      mutate(log_sim_exp_freq = log(naive_vgene_seq_freq) + log(fitness)) %>%
      group_by(across(any_of(c('mouse_id','tissue','cell_type')))) %>%
      mutate(sim_exp_freq = exp(log_sim_exp_freq)) %>%
      mutate(sim_exp_freq = sim_exp_freq / sum(sim_exp_freq)) %>%
      ungroup()
    
    stopifnot(all(simulated_freqs$sim_exp_freq <= 1))
    
    simulated_counts <- simulated_freqs  %>%
      group_by(across(any_of(c('mouse_id','tissue','cell_type')))) %>%
      mutate(n_vgene_seqs = rmultinom(1, size = unique(total_mouse_cell_type_seqs), prob = sim_exp_freq)) %>%
      ungroup() %>%
      rename(true_sim_exp_freq = sim_exp_freq) %>%
      mutate(vgene_seq_freq = n_vgene_seqs / total_mouse_cell_type_seqs) %>%
      select(mouse_id, day, infection_status, group_controls_pooled, cell_type, matches('tissue'), v_gene, fitness,
             n_naive_vgene_seqs, total_mouse_naive_seqs, naive_vgene_seq_freq, n_vgene_seqs ,total_mouse_cell_type_seqs,
             true_sim_exp_freq, vgene_seq_freq)
    
    return(simulated_counts)
  }
  
  realizations <- replicate(n_reps, generate_replicate(naive_freqs = naive_freqs, exp_freqs = exp_freqs),
                            simplify = F)  
  
  replicates_tibble <- c()
  for(i in 1:length(realizations)){
    replicates_tibble <- bind_rows(replicates_tibble,
                                   realizations[[i]] %>% mutate(replicate = i) %>%
                                     select(replicate, everything()))
  }
  
  # Check frequencies sum to 1 in each mouse-cell-type-(tissue) combination for all replicates
  checks <- replicates_tibble %>% group_by(across(any_of(c('replicate','mouse_id','tissue','cell_type')))) %>% 
    summarise(S = sum(vgene_seq_freq)) %>% ungroup() %>% select(S) %>% unique() %>% pull(S)
  
  stopifnot(all(abs(checks - 1) < 1e-7))
  
  return(replicates_tibble)
  
}

# Pre-computes the probability of observing a range of mutations for seq. lengths observed in the data, given an estimated sequencing error rate
generate_mutation_null_model <- function(seq_cluster_stats, estimated_seq_error_rate, n_mutations_variable, seq_length_variable){
  # n_mutations_variable: name of column with number of mutations, e.g.
  
  length_set <- unique(seq_cluster_stats[, seq_length_variable]) %>% unlist()
  
  # Find the range of the number of nt mutations observed in the data
  n_mutations_range <- seq(min(seq_cluster_stats[, n_mutations_variable] %>% unlist(), na.rm = T),
                           max(seq_cluster_stats[, n_mutations_variable] %>% unlist(), na.rm = T))
  
  null_model_mutations <- expand_grid(length_set, n_mutations_range) %>%
    dplyr::rename(length = length_set, n_mutations = n_mutations_range) %>%
    group_by(length) %>%
    mutate(null_prob = dbinom(x = n_mutations, size = length, prob = estimated_seq_error_rate))
  
  return(null_model_mutations)
  
}

