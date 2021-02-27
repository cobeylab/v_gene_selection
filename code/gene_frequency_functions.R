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

# Dominant cell type threshold
# A clone has to have at least 5 seqs and a cell type and at least 90% belonging to a cell type to be considered a clone with a dominant cell type.
# Since only 'pure' naive clones are used, naive cells are always the dominant cell type of naive clones.
dominant_ctype_threshold <- 0.9
min_dominant_ctype_seqs <- 5

compile_clone_info <- function(clone_info_dir){
  # Read clone files for all mice 
  # Combined tibble has one row per cell subpopulation in a clone in tissue in a mouse
  # (e.g. memory B cells from clone X in the spleen of mouse Y)
  mouse_clone_files = list.files(clone_info_dir, pattern = '_clones.csv')
  mouse_clone_files = paste0(clone_info_dir, mouse_clone_files)
  clone_info <- lapply(mouse_clone_files, FUN = function(x){return(as_tibble(read.csv(x)))})
  clone_info <- bind_rows(clone_info)

  #lapply(mouse_clone_files, FUN = function(x){return(as_tibble(read.csv(x)) %>%
  #                                                     filter(is.na(clone_id)))})
  #lapply(mouse_clone_files, FUN = function(x){return(  unique(as_tibble(read.csv(x))$mouse_id)  )})
  
  # Add variables for observation day, infection status
  clone_info <- get_info_from_mouse_id(clone_info %>% filter(is.na(clone_id) == F))
  
  clone_info <- clone_info %>% mutate(mouse_id = factor(mouse_id, levels = mouse_id_factor_levels)) %>%
    mutate(cell_type = factor(cell_type, levels = c('naive', 'GC','PC','mem')))
  return(clone_info)
}

get_clone_purity <- function(clone_info){
  clone_purity <- clone_info %>% 
    # Group each row as representing a naive or non-naive cell subpopulation
    mutate(cell_group = ifelse(cell_type == 'naive', 'naive', 'non_naive')) %>%
    group_by(mouse_id, day, infection_status, clone_id, v_gene, d_gene,
             j_gene, cell_group) %>%
    # Count number of unique productive sequences from naive and non-naive cells in each clone
    summarise(uniq_prod_seqs = sum(uniq_prod_seqs)) %>%
    ungroup() %>%
    pivot_wider(names_from = 'cell_group', values_from = 'uniq_prod_seqs') %>%
   # spread(key = "cell_group", value = 'uniq_prod_seqs') %>%
    dplyr::rename(naive_seqs_in_clone = naive, non_naive_seqs_in_clone = non_naive) %>%
    mutate(is_pure_naive = naive_seqs_in_clone > 0 & non_naive_seqs_in_clone ==0,
           is_pure_non_naive = naive_seqs_in_clone == 0 & non_naive_seqs_in_clone > 0)
  stopifnot(nrow(clone_purity) == length(unique(paste(clone_info$mouse_id, clone_info$clone_id, sep = ';'))))
  
  clone_info <- left_join(clone_info,
                          clone_purity %>%
                            select(mouse_id, day, infection_status, clone_id,
                                   is_pure_naive, is_pure_non_naive, naive_seqs_in_clone, non_naive_seqs_in_clone),
                          by = c('mouse_id','day','infection_status','clone_id'))
  
  return(clone_info)
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

# Assigns dominant cell type for each clone
# Dominant cell type > 0.9 uniq prod. seqs in clone, if none, NA.
get_dominant_ctype <- function(pure_clone_info){
  # First, sum unique sequences within cell types across tissues
  dominant_ctype <- pure_clone_info %>% group_by(mouse_id, day, infection_status, clone_id,
                               clone_type, cell_type) %>%
    summarise(uniq_prod_seqs = sum(uniq_prod_seqs)) %>%
    # Calculate total seqs in clone
    ungroup() %>% group_by(mouse_id, day, infection_status, clone_id,
                           clone_type) %>%
    mutate(total_clone_seqs = sum(uniq_prod_seqs)) %>% 
    ungroup() %>%
    mutate(fraction_clone_seqs = uniq_prod_seqs/ total_clone_seqs) %>%
    filter((clone_type == 'naive' & cell_type == 'naive') |
             (fraction_clone_seqs >= dominant_ctype_threshold & total_clone_seqs >= min_dominant_ctype_seqs)) %>%
    mutate(dominant_ctype = cell_type) 
  
  clone_by_ctype_info <- left_join(pure_clone_info, dominant_ctype %>% select(mouse_id, day, infection_status, clone_id,
                                                         clone_type, dominant_ctype),
              by = c('mouse_id', 'day', 'infection_status', 'clone_id', 'clone_type')) %>%
    mutate(dominant_ctype = as.character(dominant_ctype))
  clone_by_ctype_info$dominant_ctype[is.na(clone_by_ctype_info$dominant_ctype)] <- 'none'
  clone_by_ctype_info <- mutate(clone_by_ctype_info, dominant_ctype = factor(dominant_ctype, levels = c("naive", "GC", "PC", "mem",'none')))
  return(clone_by_ctype_info)
}

# Calculates gene frequencies for each cell type in each mouse (pools across tissues)
calc_gene_freqs <- function(clone_info, long_format = F, by_tissue = F){
  
  grouping_vars <- c('mouse_id', 'v_gene', 'cell_type')
  
  if(by_tissue){
    grouping_vars <- c(grouping_vars, 'tissue')
  }
  
  # All mouse-specific columns other than ID are removed from the grouping for convenience, but they're added back at the end
  mouse_info <- clone_info %>% select(mouse_id, day, infection_status, group, group_controls_pooled) %>%
    unique()
  
  # Calculate frequencies in each subcompartment (naive, GC, PC, mem)
  gene_seq_freqs <- clone_info %>%
    mutate(mouse_id = factor(mouse_id), v_gene = factor(v_gene), cell_type = factor(cell_type)) %>%
    group_by(across(grouping_vars)) %>%
    summarise(n_vgene_seqs = sum(uniq_prod_seqs, na.rm=T)) %>%
    ungroup() 
   
  # Complete tibble so that genes present in a mouse are represented by 0s when absent from specific subpopulations
  if(by_tissue){
    gene_seq_freqs <- gene_seq_freqs %>%
      complete(nesting(mouse_id, v_gene), tissue, cell_type,
               fill = list(n_vgene_seqs = 0))
  }else{
    gene_seq_freqs <- gene_seq_freqs %>%
      complete(nesting(mouse_id, v_gene), cell_type,
               fill = list(n_vgene_seqs = 0))
  }
  
  gene_seq_freqs <- gene_seq_freqs %>%
    # Calculate mouse-cell-type (optionally tissue) totals
    group_by(across(grouping_vars[grouping_vars != 'v_gene'])) %>%
    mutate(total_mouse_cell_type_seqs = sum(n_vgene_seqs)) %>%
    ungroup() %>%
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

# Generate pairs of mice to compare
gen_pairs <- function(clone_info){
  mouse_ids <- unique(clone_info %>% pull(mouse_id))
  
  pairs <- expand.grid(ids1 = mouse_ids, ids2 = mouse_ids) %>%
    rowwise() %>%
    mutate(pair = paste(ids1, ids2, sep = ';')) %>%
    pull(pair)
  
  unique_pairs <- c()
  # Order pair labels (to then retain only unique pairs)
  # Label shows smaller times/ids first
  for(p in pairs){
    mouse1 <-  strsplit(p, ';')[[1]][1]
    mouse2 <- strsplit(p, ';')[[1]][2]
    
    unique_pairs <- c(unique_pairs,
                      tibble(mouse = c(mouse1, mouse2)) %>%
                        mutate(mouse_id = mouse) %>%
                        separate(mouse_id, c('day','id')) %>%
                        arrange(as.numeric(day),as.numeric(id)) %>%
                        summarise(pair = paste0(mouse, collapse = ';')) %>% 
                        pull(pair))
  }
  unique_pairs <- unique(unique_pairs)
  
  return(unique_pairs)
}


# For pair of mice, rearrange freqs tibble to show changes in each mouse as different vars.
# Pair is a semi-colon separated string, e.g. "8-1;8-5"
# For now, keeps frequencies in a long format (i.e., with naive frequencies listed as another cell type instead of another column)
get_pairwise_freqs <- function(clone_info, mouse_pairs, by_tissue = F){
  
  gene_freqs <- calc_gene_freqs(clone_info, long_format = T, by_tissue = by_tissue)
  
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
    
    mouse1_freqs <- gene_freqs %>% filter(mouse_id == mice[1])
    mouse2_freqs <- gene_freqs %>% filter(mouse_id == mice[2])
    
    for(var in colnames(mouse1_freqs)){
      if((var %in% c('v_gene','cell_type','tissue')) == F){
        colnames(mouse1_freqs)[colnames(mouse1_freqs) == var] <- paste0(var,'_i')
        colnames(mouse2_freqs)[colnames(mouse2_freqs) == var] <- paste0(var,'_j')
      }
    }
    
    
    # Assemble pair tibble
    pair_changes <- full_join(mouse1_freqs, mouse2_freqs) %>%
      mutate(mouse_pair = mouse_pair,
             mouse_i = factor(mice[1], levels = mouse_id_factor_levels),
             mouse_j = factor(mice[2], levels = mouse_id_factor_levels),
             pair_time = pair_time, time_i = time_i, time_j = time_j,
             pair_type = factor(pair_type, levels = c('control','control/primary','primary',
                                                      'control/secondary','secondary','primary/secondary'))) %>%
      select(mouse_pair, pair_type, matches('infection_status'), everything())  %>%
      filter(!is.na(vgene_seq_freq_i), !is.na(vgene_seq_freq_j))
    return(pair_changes)
  }
  # Paired changes in gene frequences
  paired_gene_freqs <- lapply(mouse_pairs, FUN = internal_function,
                                     gene_freqs = gene_freqs)
  paired_gene_freqs <- bind_rows(paired_gene_freqs) %>%
    mutate(time_i = factor(time_i, levels = c('8','16','24','40','56')))
  return(paired_gene_freqs)
}

get_pairwise_correlations <- function(pairwise_freqs){
  
  grouping_vars <- c('mouse_pair','pair_type','mouse_i','mouse_j','time_i','time_j',
                     'cell_type')
  if('tissue' %in% names(pairwise_freqs)){
    grouping_vars <- c(grouping_vars, 'tissue')
  }
  
  pairwise_correlations <- pairwise_freqs %>% 
    filter(mouse_i != mouse_j) %>%
    group_by(across(grouping_vars)) %>%
    mutate(n_genes_in_comparison = length(unique(v_gene))) %>%
    filter(n_genes_in_comparison >= 10) %>%
    summarise(cor_coef = cor.test(vgene_seq_freq_i, vgene_seq_freq_j,
                                  method = 'spearman')$estimate) %>%
    ungroup()
  
  # Check no comparisons were made between mice sampled in different days
  #stopifnot(all(pairwise_correlations %>% mutate(test = time_i == time_j) %>% pull(test) == T))
  return(pairwise_correlations)
} 


# Generates datasets with selected genes specified for each mouse in synth_data_input_tibble
# If synth_data_input_tibble is 'neutral', assumes neutrality (i.e. genes in exp. repertoire are sampled based on naive frequencies)
simulate_selection_freq_changes <- function(clone_info, synth_data_input_tibble, min_seqs, extra_mice = 0, by_tissue = F,
                                            naive_from_tissue = NULL, n_reps = 100){
  
  # If extra_mice = 1, duplicate each mouse before simulation, if = 2, make two copies of each mouse before simulation, etc.
  base_data <- clone_info
  
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
  
  # Naive frequencies can be across all tissues or from a specific tissue 
  # (in which case they're used for all tissues if experienced reps are simulated for each tissue)
  if(is.null(naive_from_tissue)){
    naive_freqs <- calc_gene_freqs(base_data, long_format = F, by_tissue = F)$naive_freqs
  }else{
    naive_freqs <- calc_gene_freqs(base_data, long_format = F, by_tissue = T)$naive_freqs %>%
      filter(tissue == naive_from_tissue) %>% select(-tissue)
  }
  
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













# OLD FUNCTIONS TO REVIEW OR DELETE:


# Calculates changes in gene frequency from naive pool in terms of the % of unique prod. seqs.
# For now looks only at V genes.
calc_gene_seq_freq_changes <- function(clone_info){
  stopifnot('check' == 'Review this function')
  # Calculate gene usage freqs as % of seqs for each mouse and cell type
  freqs <- calc_gene_freqs(clone_info)
  
  gene_seq_freqs <- freqs$exp_freqs

  naive_freqs <- freqs$naive_freqs

  # Join experienced freqs and naive freqs.
  gene_seq_freq_changes <- left_join(gene_seq_freqs,
                                       naive_freqs %>% 
                                       select(mouse_id, day, infection_status, group_controls_pooled, v_gene,
                                              n_naive_vgene_seqs, total_mouse_naive_seqs, naive_vgene_seq_freq),
                                       by = c("mouse_id", "day", "infection_status","group_controls_pooled","v_gene")) 
  # Remove cases where gene is missing from both naive and experienced repertoire (e.g. in a specific cell type)
  gene_seq_freq_changes <- gene_seq_freq_changes %>%
    filter(naive_vgene_seq_freq > 0 | vgene_seq_freq > 0)
  
  # For genes with naive frequency zero but present in experienced repertoire or vice-versa, if any, set 0 freq to 10-7
  # gene_seq_freq_changes <- gene_seq_freq_changes %>%
  #   rowwise() %>%
  #   mutate(naive_vgene_seq_freq = ifelse(naive_vgene_seq_freq == 0, 1e-7, naive_vgene_seq_freq)) %>%
  #   mutate(vgene_seq_freq = ifelse(vgene_seq_freq == 0, 1e-7, vgene_seq_freq)) %>%
  #   ungroup()
  
  gene_seq_freq_changes <- gene_seq_freq_changes %>%
    mutate(vgene_seq_freq_change = log(vgene_seq_freq) - log(naive_vgene_seq_freq)) %>%
    mutate(vgene_seq_freq_change = exp(vgene_seq_freq_change)) %>%
    group_by(mouse_id, infection_status, group_controls_pooled, day, cell_type) %>%
    # Calculate rank abundance of each gene in naive and experienced repertoires (exp. repertoire by cell type)
    mutate(rank_exp_freq = rank(-vgene_seq_freq, ties.method = 'average')) %>%
    ungroup()
  
  # For the ranking of frequency changes, calculate separate ranks with and without Inf changes caused by zero naive freq.
  rank_change_excluding_zerofreqs <- gene_seq_freq_changes %>%
    filter(naive_vgene_seq_freq > 0) %>%
    group_by(mouse_id, infection_status, group_controls_pooled, day, cell_type) %>%
    mutate(rank_change_adjusted = rank(-vgene_seq_freq_change, ties.method = 'average')) %>%
    ungroup() %>% select(mouse_id, infection_status, day, cell_type, v_gene, rank_change_adjusted)
  
  gene_seq_freq_changes <- left_join(gene_seq_freq_changes, rank_change_excluding_zerofreqs, 
                                 by = c('mouse_id', 'infection_status', 'day', 'cell_type','v_gene')) %>%
    # Now calculate raning of frequency changes including  Inf changes caused by zero naive freq.
    group_by(mouse_id, infection_status, day, group_controls_pooled, cell_type) %>%
    mutate(rank_change = rank(-vgene_seq_freq_change, ties.method = 'average')) %>%
    ungroup()
  
  return(gene_seq_freq_changes)
}



get_changes_consecutive_pairs <- function(paired_gene_freq_changes){
  previous_time_point <- function(time_point){
    tpoints_order <- sort(as.numeric(as.character(unique(paired_gene_freq_changes$time_i))))
    previous_tp <- time_point
    previous_tp[time_point == min(time_point)] <- NA
    for(i in 2:length(tpoints_order)){
      previous_tp[time_point == tpoints_order[i]] <- tpoints_order[i-1]
    }
    return(previous_tp)
  }
  
  consecutive_pairs <- paired_gene_freq_changes %>% 
    filter(time_i != time_j) %>%
    mutate(time_point_before_time_j = previous_time_point(time_j)) %>%
    filter(time_i == time_point_before_time_j)
  return(consecutive_pairs)
}

# Pairwise correlations


get_naive_freq_correlations <- function(paired_gene_freq_changes){
  pairwise_naive_freq_correlations <- paired_gene_freq_changes %>% 
    # Naive frequencies are repeated for each cell type in input, arbitrarily filter to mem
    filter(cell_type == 'mem') %>%
    group_by(mouse_pair, pair_type, mouse_i, mouse_j, time_i, time_j) %>%
    # Remove genes where gene was missing from one or both mice
    filter(!is.na(naive_vgene_seq_freq_i) & !is.na(naive_vgene_seq_freq_j)) %>%
    summarise(cor_coef = cor.test(naive_vgene_seq_freq_i, naive_vgene_seq_freq_j,
                                  method = 'spearman')$estimate) %>%
    ungroup()
  return(pairwise_naive_freq_correlations)
}

plot_pairwise_correlations <- function(pairwise_correlations, day){
  # For plotting, fill half of the matrix with NAs (since it's symmetric)
  for(r in 1:nrow(pairwise_correlations)){
    datarow <- pairwise_correlations[r,]
    mouse_i <- datarow$mouse_i
    mouse_j <- datarow$mouse_j
    coef <- datarow$cor_coef
    datarow <- datarow %>% mutate(mouse_i = !!mouse_j, mouse_j = !!mouse_i,
                                  cor_coef = coef)
    pairwise_correlations <- bind_rows(pairwise_correlations, datarow)
  }
  # Factor levels to cluster mice by infection status
  mouse_ids <- unique(pairwise_correlations$mouse_i)
  grouped_mouse_ids <- c( mouse_id_factor_levels[grepl('-5', mouse_id_factor_levels) |
                                                   grepl('-10', mouse_id_factor_levels)],
                          mouse_id_factor_levels[!grepl('-5', mouse_id_factor_levels) &
                                                  !grepl('-10', mouse_id_factor_levels)])
  grouped_mouse_ids <- grouped_mouse_ids[grepl(paste0(day,'-'), grouped_mouse_ids)]
  
  pl <- pairwise_correlations %>% filter(time_i == day) %>%
    mutate(mouse_i = factor(mouse_i, levels = grouped_mouse_ids),
           mouse_j = factor(mouse_j, levels = grouped_mouse_ids)) %>%
    filter(is.na(cor_coef) == F) %>%
    ggplot(aes(x = mouse_i, y = mouse_j, fill = cor_coef)) +
    geom_tile() +
    scale_fill_viridis(name = "Pearson's\ncoef.",
                       option = 'cividis',na.value = 'white') +
    xlab('Mouse') +
    ylab('Mouse') +
    facet_grid(.~cell_type) +
    theme(axis.text.x = element_text(angle = 90))
  
  return(pl)
}

get_exp_vs_naive_correlations <- function(gene_freq_changes, method = 'spearman'){
  exp_vs_naive_correlations <- gene_freq_changes %>%
    group_by(mouse_id, day, infection_status, cell_type) %>%
    summarise(exp_vs_naive_cor = cor.test(naive_vgene_seq_freq, vgene_seq_freq, method = method)$estimate) %>%
    ungroup()
  return(exp_vs_naive_correlations)
}

get_exp_vs_selection_correlation <- function(gene_freq_changes, method = 'spearman',
                                             remove_inf_selection = T){
  
  if(remove_inf_selection){
    exp_vs_selection_correlations <- gene_freq_changes %>% filter(vgene_seq_freq_change != Inf)
  }else{
    exp_vs_selection_correlations <- gene_freq_changes
  }
  
  exp_vs_selection_correlations <- exp_vs_selection_correlations %>%
    group_by(mouse_id, day, infection_status, cell_type) %>%
    summarise(exp_vs_selection_cor = cor.test(vgene_seq_freq, vgene_seq_freq_change,
                                              method = method)$estimate) %>%
    ungroup()
  return(exp_vs_selection_correlations)
}

# Find genes unique to each mouse.
count_unique_mouse_genes <- function(gene_freq_changes){
  mice_by_gene <- gene_freq_changes %>% group_by(v_gene) %>%
    summarise(n_mice_gene_is_present = length(unique(mouse_id)),
              mice_gene_is_present = paste0(unique(mouse_id), collapse = ';')) %>% 
    ungroup()
  
  mouse_ids <- unique(gene_freq_changes$mouse_id)
  
  unique_genes_tibble <- c()
  for(mid in mouse_ids){
    mouse_unique_genes <- mice_by_gene$v_gene[mice_by_gene$mice_gene_is_present == mid]
    unique_genes_tibble <- bind_rows(unique_genes_tibble,
                                     tibble(mouse_id = mid, mouse_unique_genes = length(mouse_unique_genes)))
  }
  
  return(unique_genes_tibble)
}


