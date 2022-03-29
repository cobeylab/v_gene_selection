library(tidyr)
library(dplyr)
library(cowplot)
library(stringr)
library(viridis)
library(reshape2)
#library(heatmaply)
library(gplots)

source('plot_options.R')

estimated_seq_error_rate <- 0.0018

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
  mouse_info <- data %>% 
    filter(grepl('untreated', mouse_id) == F) %>% # Excludes mouse labels from Greiff et al. 2017 data
    select(mouse_id) %>%
    unique() %>%
    mutate(day = str_extract(mouse_id,'[0-9]*[^-]'),
           mouse_id_number = str_replace(str_extract(mouse_id,'-[0-9]*'),'-',''),
           #day = factor(day, levels = sort(unique(as.numeric(day)))),
           infection_status = assign_infection_status(day, mouse_id_number))
  
  mouse_info$day[mouse_info$mouse_id == '40-7'] <- '8'
  
  mouse_info <- mouse_info %>%
    mutate(group = paste(infection_status, day, sep = '-')) %>%
    mutate(group_controls_pooled = group) %>%
    mutate(group_controls_pooled = ifelse(grepl('control',group_controls_pooled), 'control', group_controls_pooled)) %>%
    select(-mouse_id_number) 
  
  if(any(grepl('untreated', unique(data$mouse_id)))){
    Greiff2017_mouse_info <- tibble(mouse_id = unique(data$mouse_id)[grepl('untreated',unique(data$mouse_id))],
                                    day = NA, infection_status = 'Greiff2017', group = 'Greiff2017',
                                    group_controls_pooled = 'Greiff2017')
    mouse_info <- bind_rows(mouse_info, Greiff2017_mouse_info)
  }
  
 
  return(left_join(data, mouse_info, by = 'mouse_id') %>%
           select(mouse_id, day, infection_status, group, group_controls_pooled, everything()))
}


# Determines infection status based on mouse id (treats all post day 40 non-controls as secondary)
# But treats 40-7 as primary infection on day 8.
assign_infection_status <- function(day, mouse_id_number){
  day <- as.numeric(as.character(day)) # Converting factor to numeric through character
  mouse_id_number <- as.numeric(as.character(mouse_id_number))
  infection_status <- ifelse(mouse_id_number %in% c(5,10), 'control','infected')
  infection_status[infection_status == 'infected'] <- ifelse(
    day[infection_status == 'infected'] < 40, 'primary', 'secondary'
  )
  
  # 40-7 was de facto a primary-8 infected mouse 
  infection_status[day == '40' & mouse_id_number == '7'] <- 'primary'
  
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


get_clone_purity <- function(seq_counts){
  
  # The object will have a "unique_prod_seqs" column if has counts obtained after clustering identical sequences
  names(seq_counts)[names(seq_counts) == 'unique_prod_seqs'] <- 'prod_seqs'
  
  clone_purity <- seq_counts %>% 
    # Group each row as representing a Dump-IgD+B220+ or non-Dump-IgD+B220+ cell subpopulation
    mutate(cell_group = ifelse(cell_type == 'IgD+B220+', 'IgD+B220+', 'non_IgD+B220+')) %>%
    group_by(mouse_id, clone_id, cell_group) %>%
    # Count number of unique productive sequences from naive and non-naive cells in each clone
    dplyr::summarise(prod_seqs = sum(prod_seqs)) %>%
    ungroup() %>%
    pivot_wider(names_from = 'cell_group', values_from = 'prod_seqs', values_fill = 0) %>%
    dplyr::rename(IgD_B220_seqs_in_clone = `IgD+B220+`, non_IgD_B220_seqs_in_clone = `non_IgD+B220+`) %>%
    mutate(
      clone_purity = case_when(
        IgD_B220_seqs_in_clone > 0 & non_IgD_B220_seqs_in_clone == 0 ~ 'pure_IgD+B220+',
        IgD_B220_seqs_in_clone == 0 & non_IgD_B220_seqs_in_clone > 0 ~ 'pure_non_IgD+B220+',
        IgD_B220_seqs_in_clone > 0 & non_IgD_B220_seqs_in_clone > 0 ~ 'mixed'
      )
    )
    
  return(clone_purity)
} 

# Selects sequences sorted as naive based on additional filters (n. mutations, isotype, clone size)
process_IgD_B220_seqs <- function(annotated_seqs, max_clone_unique_IgDB220_seqs,
                                  max_v_gene_mutations){
  
  # Naive sequences are those that satisfy ALL the following:
  # - were sorted as naive (DUMP-IgD+B220+)
  # - Have IgD or IgM as the isotype determined from the read (or isotype is NA)
  # - Are in a clone that only has DUMP-IgD+B220+ cells
  # - Are in a clone with at most max_clone_naive_seqs **unique** seqs, where unique = same mouse, clone, tissue, isotype and sequence
  # - Has max_v_gene_mutations or fewer mutations in the V gene region.
  
  # Sequences sorted as DUMP-IgD+B220+ that do not meet all the other criteria are labeled 'nonnaive_IgD+B220+'
  
  unique_productive_seq_counts <- annotated_seqs %>%
    filter(productive_partis) %>%
    select(mouse_id, clone_id, partis_uniq_ref_seq, tissue, cell_type, isotype) %>%
    unique() %>%
    group_by(mouse_id, clone_id, tissue, cell_type) %>%
    dplyr::summarise(unique_prod_seqs = n()) %>%
    ungroup()
  
  clone_purity <- get_clone_purity(unique_productive_seq_counts) %>%
    dplyr::rename(unique_productive_IgDB220_seqs_in_clone = IgD_B220_seqs_in_clone,
                  unique_productive_nonIgDB220_seqs_in_clone = non_IgD_B220_seqs_in_clone)
  annotated_seqs <- left_join(annotated_seqs, clone_purity)
  
  
  igd_b220_seqs <- annotated_seqs %>% 
    filter(cell_type == 'IgD+B220+')

  igd_b220_seqs <- igd_b220_seqs %>%
    mutate(cell_type = case_when(
      (isotype %in% c('IGM','IGD') | is.na(isotype)) & clone_purity == "pure_IgD+B220+" &
        unique_productive_IgDB220_seqs_in_clone <= max_clone_unique_IgDB220_seqs &
        vgene_mutations_partis_nt <= max_v_gene_mutations ~ 'naive',
      !is.na(clone_purity) ~ 'nonnaive_IgD+B220+',
      T ~ 'unassigned_IgD+B220+' # Will contain IgD+B220+ in clones with only unproductive seqs
    ))
  
  processed_tibble <- bind_rows(annotated_seqs %>% filter(cell_type != 'IgD+B220+'),
            igd_b220_seqs %>% select(all_of(names(annotated_seqs))))
  
  stopifnot(nrow(processed_tibble) == nrow(annotated_seqs))
  
  return(processed_tibble)
  
}

calc_naive_freqs <- function(naive_seq_counts, clone_info){
  
  #clone_info used to get tibble with all genes present in each mouse (including ones potentially not observed in naive rep.)
  names(naive_seq_counts)[names(naive_seq_counts) == 'unique_prod_seqs'] <- 'prod_seqs'
  
  naive_freqs <- naive_seq_counts %>%
    group_by(mouse_id, v_gene) %>%
    dplyr::summarise(n_naive_vgene_seqs = sum(prod_seqs)) %>%
    group_by(mouse_id) %>%
    mutate(naive_vgene_seq_freq = n_naive_vgene_seqs / sum(n_naive_vgene_seqs))
  
  naive_freqs <- left_join(clone_info %>% filter(!is.na(v_gene)) %>%
                             select(mouse_id, v_gene) %>% unique(),
                           naive_freqs) %>%
    replace_na(list(n_naive_vgene_seqs = 0, naive_vgene_seq_freq = 0))
  
  naive_freqs <- naive_freqs %>% 
    group_by(mouse_id) %>%
    mutate(total_mouse_naive_seqs = sum(n_naive_vgene_seqs)) %>%
    ungroup()
  
  return(naive_freqs)
  
}

# Calculates gene frequencies for each cell type in each mouse (pools across tissues)
calc_gene_freqs <- function(exp_seq_counts, naive_seq_counts, clone_info, long_format = F, by_tissue = F){
  names(exp_seq_counts)[names(exp_seq_counts) == 'unique_prod_seqs'] <- 'prod_seqs'
  
  grouping_vars <- c('mouse_id', 'v_gene', 'cell_type')
  
  if(by_tissue){
    grouping_vars <- c(grouping_vars, 'tissue')
  }
  
  # All mouse-specific columns other than ID are removed from the grouping for convenience, but they're added back at the end
  mouse_info <- exp_seq_counts %>% select(mouse_id, day, infection_status, group, group_controls_pooled) %>%
    unique()
  
  # Generate tibble so that all genes present at all in a mouse will be explicitly 
  # represented as zeros in compartments they're missing from
  
  all_genes_in_mice <- clone_info %>% filter(!is.na(v_gene)) %>% select(mouse_id, v_gene) %>% unique()
  tissue_cell_type_combinations <- exp_seq_counts %>% select(tissue, cell_type) %>% unique()
  
  if(by_tissue){
    all_genes_in_mice <- mapply(FUN = function(tis, ctype, tib){tib %>% mutate(tissue = tis, cell_type = ctype) %>%
        select(mouse_id, tissue, cell_type, v_gene)}, 
           tis = tissue_cell_type_combinations$tissue,
           ctype = tissue_cell_type_combinations$cell_type,
        MoreArgs = list(tib = all_genes_in_mice),
        SIMPLIFY = F
        )
  }else{
    all_genes_in_mice <- lapply(as.list(unique(tissue_cell_type_combinations$cell_type)),
           FUN = function(ctype, tib){tib %>% mutate(cell_type = ctype) %>%
               select(mouse_id, cell_type, v_gene)},
           tib = all_genes_in_mice)
    
  }
  all_genes_in_mice <- bind_rows(all_genes_in_mice)
    
  
  # Calculate frequencies in each subcompartment (GC, PC, mem)
  gene_seq_freqs <- exp_seq_counts %>%
    mutate(mouse_id = factor(mouse_id), v_gene = factor(v_gene), cell_type = factor(cell_type)) %>%
    group_by(across(grouping_vars)) %>%
    dplyr::summarise(n_vgene_seqs = sum(prod_seqs, na.rm=T)) %>%
    ungroup() 
   
  # Left join to all_genes_in_mice tibble, so that genes present in mice get an explicit zero value
  # in compartments they're missing from
  gene_seq_freqs <- left_join(all_genes_in_mice, gene_seq_freqs) %>%
    replace_na(list(n_vgene_seqs =  0))
  
  gene_seq_freqs <- gene_seq_freqs %>%
    # Calculate mouse-cell-type (optionally tissue) totals
    group_by(across(grouping_vars[grouping_vars != 'v_gene'])) %>%
    mutate(total_compartment_seqs = sum(n_vgene_seqs)) %>%
    ungroup() %>%
    # Exclude compartments that are entirely zero
    filter(total_compartment_seqs > 0) %>%
    ungroup()
  
  # Calculate frequencies in the experienced repertoire as a whole, across mem/PC/GC
  gene_seq_freqs_all_exp <- gene_seq_freqs %>%
    group_by(across(grouping_vars[grouping_vars != 'cell_type'])) %>%
    dplyr::summarise(n_vgene_seqs = sum(n_vgene_seqs)) %>%
    mutate(cell_type = 'experienced') %>%
    ungroup() %>%
    group_by(across(grouping_vars[!(grouping_vars %in% c('v_gene','cell_type'))])) %>%
    mutate(total_compartment_seqs = sum(n_vgene_seqs)) %>%
    ungroup() %>%
    select(mouse_id, v_gene, cell_type, everything())
  
  # Tibble with experienced frequencies
  exp_freqs <- bind_rows(gene_seq_freqs,
                         gene_seq_freqs_all_exp) %>%
    arrange(mouse_id, v_gene, cell_type) %>%
    mutate(vgene_seq_freq = n_vgene_seqs/total_compartment_seqs)
  
  # Add complete mouse_information
  exp_freqs <- left_join(exp_freqs, mouse_info, by = 'mouse_id') %>%
    select(mouse_id, day, infection_status, group, group_controls_pooled, everything())

  # Naive frequencies calculated by a separate function:
  naive_freqs <- calc_naive_freqs(naive_seq_counts, clone_info)
  
  naive_freqs <- left_join(naive_freqs, mouse_info, by = 'mouse_id') %>%
    select(mouse_id, day, infection_status, group, group_controls_pooled, everything())
  
  if(long_format == T){
    output <- bind_rows(exp_freqs,
                        naive_freqs %>%
                          rename(n_vgene_seqs = n_naive_vgene_seqs, vgene_seq_freq = naive_vgene_seq_freq,
                                 total_compartment_seqs = total_mouse_naive_seqs) %>%
                          mutate(cell_type = 'naive')) %>% arrange(mouse_id, v_gene) %>%
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
    mutate(resampled_exp_count = rmultinom(1, size = unique(total_compartment_seqs), prob = vgene_seq_freq)) %>%
    mutate(n_vgene_seqs = resampled_exp_count, vgene_seq_freq = n_vgene_seqs / total_compartment_seqs) %>%
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
           vgene_seq_freq = vgene_seq_freq + 1/total_compartment_seqs * (vgene_seq_freq == 0)) %>%
    # Re-normalize frequencies so they sum to 1
    group_by(mouse_id, cell_type) %>%
    mutate(vgene_seq_freq = vgene_seq_freq / sum(vgene_seq_freq)) %>%
    mutate(total_compartment_seqs = sum(n_vgene_seqs)) %>%
    ungroup()
  
  return(adj_exp_freqs)
}

# Correlation within each mouse between experienced and naive gene frequencies
get_naive_exp_correlations <- function(gene_freqs, method){
  gene_freqs %>%
    group_by(mouse_id, day, infection_status, group, group_controls_pooled, tissue, cell_type,
             total_compartment_seqs, total_mouse_naive_seqs) %>%
    dplyr::summarise(spearman = cor.test(vgene_seq_freq, naive_vgene_seq_freq, method = 'spearman')$estimate,
                     pearson = cor.test(vgene_seq_freq, naive_vgene_seq_freq, method = 'pearson')$estimate) %>%
    ungroup() %>%
    pivot_longer(cols = c('spearman','pearson'), names_to = 'method', values_to = 'naive_exp_corr')
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
get_pairwise_freqs <- function(gene_freqs, adjust_naive_zeros){
  # Assumes gene_freqs in wide format:
  stopifnot('naive_vgene_seq_freq' %in% names(gene_freqs))
  
  mouse_info <- gene_freqs %>%
    select(mouse_id, day, infection_status, group, group_controls_pooled) %>%
    unique()
  
  unique_pairs <- gene_freqs %>% select(mouse_id) %>% unique() %>%
    dplyr::rename(mouse_id_i = mouse_id) %>%
    mutate(mouse_id_j = mouse_id_i) %>%
    complete(mouse_id_i, mouse_id_j) %>%
    rowwise() %>%
    mutate(pair = paste0(sort(c(mouse_id_i, mouse_id_j)), collapse = ';')) %>%
    ungroup() %>%
    filter(mouse_id_i != mouse_id_j) %>%
    select(pair) %>%
    unique() %>% pull(pair)
  
  # Get naive and experienced frequencies into separate tibbles
  exp_freqs <- gene_freqs %>% select(mouse_id, day, infection_status, group, group_controls_pooled,
                                     v_gene, tissue, cell_type, n_vgene_seqs, total_compartment_seqs,
                                     vgene_seq_freq, obs_rho, mean_sim_rho, lbound_sim_rho, ubound_sim_rho,
                                     deviation_from_naive, v_gene_rank)
  
  naive_freqs <- gene_freqs %>%
    select(mouse_id, day, infection_status, group, group_controls_pooled, v_gene,
           n_naive_vgene_seqs, total_mouse_naive_seqs, naive_vgene_seq_freq) %>%
    unique() %>%
    dplyr::rename(n_vgene_seqs = n_naive_vgene_seqs, total_compartment_seqs = total_mouse_naive_seqs,
                  vgene_seq_freq = naive_vgene_seq_freq) %>%
    mutate(cell_type = 'naive', tissue = 'naive_source_tissue') %>%
    select(mouse_id, day, infection_status, group, group_controls_pooled, v_gene, cell_type, n_vgene_seqs,
           total_compartment_seqs, vgene_seq_freq)
  

  # Put experienced and naive frequencies back together in long format
  gene_freqs <- bind_rows(exp_freqs, naive_freqs)
  
  compartment_sizes <- gene_freqs %>%
    select(mouse_id, cell_type, matches('tissue'), total_compartment_seqs) %>%
    unique()
  
  internal_function <- function(mouse_pair, gene_freqs){
    # Parse mouse_pair, determine type of pair (e.g. 'primary', 'primary/control', 'secondary/control')
    mice <- str_split(mouse_pair,';')[[1]]
    
    pair_info <- mouse_info %>%
      filter(mouse_id %in% mice) %>%
      mutate(mouse_id = factor(mouse_id, levels = mice)) %>%
      arrange(mouse_id)
    
    days <- pair_info$day
    
    infection_status <- pair_info$infection_status
    #mouse_id_numbers <- sapply(mice, FUN = function(x){str_split(x,'-')[[1]][2]})
   
    pair_type = ifelse(length(unique(infection_status)) == 1,
                       infection_status,
                       paste(sort(unique(infection_status)),collapse = '/'))
    time_i <- days[1]
    time_j <- days[2]
    #pair_time <- paste(sort(as.numeric(days)), collapse=';')
    
    mouse1_freqs <- gene_freqs %>% filter(mouse_id == mice[1]) %>% 
      select(-day, -infection_status, -group, -group_controls_pooled, - total_compartment_seqs)
    mouse2_freqs <- gene_freqs %>% filter(mouse_id == mice[2]) %>%
      select(-day, -infection_status, -group, -group_controls_pooled, -total_compartment_seqs)
    
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
                                                                          .cols = c(mouse_id, total_compartment_seqs)))
  paired_gene_freqs <- left_join(paired_gene_freqs,
                                 compartment_sizes %>% dplyr::rename_with(.fn = function(x){paste0(x,'_j')},
                                                                          .cols = c(mouse_id, total_compartment_seqs))) 
  
  naive_compartment_sizes <- compartment_sizes %>% filter(cell_type == 'naive') %>%
    dplyr::rename(total_mouse_naive_seqs = total_compartment_seqs) %>% select(-tissue, -cell_type)
    
  
  paired_gene_freqs <- left_join(paired_gene_freqs,
                                 naive_compartment_sizes %>% dplyr::rename_with(.fn = function(x){paste0(x,'_i')},
                                                                                .cols =  everything()))
  
  paired_gene_freqs <- left_join(paired_gene_freqs,
                                 naive_compartment_sizes %>% dplyr::rename_with(.fn = function(x){paste0(x,'_j')},
                                                                                .cols =  everything())) %>%
    select(mouse_pair, pair_type, matches('id'), matches('day'), matches('infection_status'), matches('group'),
           matches('total_compartment_seqs'), matches('total_mouse_naive_seqs'),
           everything())
  
  # Label each gene in a paper in terms of concordant/discordant direction of change from the naive repertoire
  paired_gene_freqs <- paired_gene_freqs %>%
    mutate(concordance_status = case_when(
      (deviation_from_naive_i == 'positive' & deviation_from_naive_j == 'positive') ~ 'concordant-increasing',
      (deviation_from_naive_i == 'negative' & deviation_from_naive_j == 'negative') ~ 'concordant-decreasing',
      (deviation_from_naive_i == 'neutral' & deviation_from_naive_j == 'neutral') ~ 'concordant-stable',
      (deviation_from_naive_i != deviation_from_naive_j) ~ 'discordant',
      (is.na(deviation_from_naive_i) | is.na(deviation_from_naive_j)) ~ 'NA'))
      
    
  return(paired_gene_freqs)
}

# Calculates pairwise correlations v gene freqs. or v genes experienced-to-naive freq.ratios between pairs of mice.
get_pairwise_correlations <- function(pairwise_gene_freqs, min_genes_in_comparison = 10, include_freq_ratios = T){
  
  grouping_vars <- c('mouse_pair','pair_type','mouse_id_i','mouse_id_j','day_i','day_j',
                     'cell_type', 'total_compartment_seqs_i','total_compartment_seqs_j',
                     'total_mouse_naive_seqs_i', 'total_mouse_naive_seqs_j')
  if('tissue' %in% names(pairwise_gene_freqs)){
    grouping_vars <- c(grouping_vars, 'tissue')
  }
  
  pairwise_correlations <- pairwise_gene_freqs %>% 
    filter(mouse_id_i != mouse_id_j) %>%
    group_by(across(grouping_vars)) %>%
    mutate(n_genes_in_freqs_comparison = sum(!is.na(n_vgene_seqs_i) & !is.na(n_vgene_seqs_j)),
           n_genes_in_freq_ratio_comparison = sum(!is.na(obs_rho_i) & !is.na(obs_rho_j)))
  
  pairwise_correlations_freqs <- pairwise_correlations %>%
    filter(n_genes_in_freqs_comparison >= min_genes_in_comparison) %>%
    dplyr::summarise(spearman = cor.test(vgene_seq_freq_i, vgene_seq_freq_j,
                                  method = 'spearman')$estimate,
                     pearson = cor.test(vgene_seq_freq_i, vgene_seq_freq_j,
                                         method = 'pearson')$estimate) %>%
    ungroup() %>%
    pivot_longer(cols = c('spearman', 'pearson'), names_to = 'method', values_to = 'cor_coef_freqs')
  
  if(include_freq_ratios){
    pairwise_correlations_freq_ratios <- pairwise_correlations %>%
      filter(n_genes_in_freq_ratio_comparison >= min_genes_in_comparison) %>%
      dplyr::summarise(spearman = cor.test(obs_rho_i, obs_rho_j,
                                           method = 'spearman')$estimate,
                       pearson = cor.test(obs_rho_i, obs_rho_j,
                                          method = 'pearson')$estimate) %>%
      ungroup() %>%
      pivot_longer(cols = c('spearman', 'pearson'), names_to = 'method', values_to = 'cor_coef_freq_ratios')
    
    return(list(freqs = pairwise_correlations_freqs, freq_ratios = pairwise_correlations_freq_ratios))
  }else{
    return(pairwise_correlations_freqs)
  }

} 


# Computes "concordance" (for a pair of mice, what % of alleles increase in both, decrease in both, etc.)
compute_deviation_concordance <- function(pairwise_gene_freqs, min_genes_in_comparison = 10){
  grouping_vars <- c('mouse_pair','pair_type','mouse_id_i','mouse_id_j','day_i','day_j',
                     'cell_type', 'total_compartment_seqs_i','total_compartment_seqs_j',
                     'total_mouse_naive_seqs_i', 'total_mouse_naive_seqs_j')
  
  if('tissue' %in% names(pairwise_gene_freqs)){
    grouping_vars <- c(grouping_vars, 'tissue')
  }
  
  deviation_concordance <- pairwise_gene_freqs %>% 
    filter(mouse_id_i != mouse_id_j) %>%
    group_by(across(grouping_vars)) %>%
    mutate(n_genes_in_comparison = sum(!is.na(obs_rho_i) & !is.na(obs_rho_j))) %>%
    filter(n_genes_in_comparison >= min_genes_in_comparison) %>%
    filter(concordance_status != 'NA') %>%
    ungroup() %>%
    group_by(across(c(grouping_vars, 'concordance_status'))) %>%
    dplyr::summarise(n_alleles = n()) %>%
    ungroup() %>%
    group_by(across(grouping_vars)) %>%
    mutate(fraction_alleles = n_alleles / sum(n_alleles)) %>%
    ungroup()
    
  
  return(deviation_concordance)


}


# Generates datasets with selected genes specified for each mouse in synth_data_input_tibble
# If synth_data_input_tibble is 'neutral', assumes neutrality (i.e. genes in exp. repertoire are sampled based on naive frequencies)
simulate_selection_freq_changes <- function(exp_seq_counts, naive_seq_counts, clone_info, synth_data_input_tibble, by_tissue = T, n_reps = 100){
  
  # Get gene frequencies
  gene_freqs <- calc_gene_freqs(exp_seq_counts, naive_seq_counts, clone_info, long_format = F, by_tissue = by_tissue)

  exp_freqs <- gene_freqs$exp_freqs
  naive_freqs <- gene_freqs$naive_freqs
  
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
                                 by = c('mouse_id','v_gene')) 
    
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
      mutate(n_vgene_seqs = rmultinom(1, size = unique(total_compartment_seqs), prob = sim_exp_freq)) %>%
      ungroup() %>%
      rename(true_sim_exp_freq = sim_exp_freq) %>%
      mutate(vgene_seq_freq = n_vgene_seqs / total_compartment_seqs) %>%
      select(mouse_id, day, infection_status, group, group_controls_pooled, cell_type, matches('tissue'), v_gene, fitness,
             n_naive_vgene_seqs, total_mouse_naive_seqs, naive_vgene_seq_freq, n_vgene_seqs ,total_compartment_seqs,
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
    dplyr::summarise(S = sum(vgene_seq_freq)) %>% ungroup() %>% select(S) %>% unique() %>% pull(S)
  
  stopifnot(all(abs(checks - 1) < 1e-7))
  
  return(replicates_tibble)
  
}

# Pre-computes the probability of observing a range of mutations for seq. lengths observed in the data, given an estimated sequencing error rate
generate_mutation_null_model <- function(annotated_seqs, estimated_seq_error_rate, n_mutations_variable, seq_length_variable){
  # n_mutations_variable: name of column with number of mutations, e.g.
  
  length_set <- unique(annotated_seqs[, seq_length_variable]) %>% unlist()
  
  # Find the range of the number of nt mutations observed in the data
  n_mutations_range <- seq(min(annotated_seqs[, n_mutations_variable] %>% unlist(), na.rm = T),
                           max(annotated_seqs[, n_mutations_variable] %>% unlist(), na.rm = T))
  
  null_model_mutations <- expand_grid(length_set, n_mutations_range) %>%
    dplyr::rename(length = length_set, n_mutations = n_mutations_range) %>%
    group_by(length) %>%
    mutate(null_prob = dbinom(x = n_mutations, size = length, prob = estimated_seq_error_rate))
  
  return(null_model_mutations)
  
}


# Gets distribution of the number of mutations by mouse, tissue and cell type
get_distribution_of_mutations <- function(annotated_seqs, n_mutations_variable, disable_grouping = F){
  # n_mutations_variable: the name of the variable with the number of mutations to summarize
  # e.g. 'n_mutations_partis_nt' is the number of nt mutations across the whole sequence,
  # 'vgene_mutations_partis_nt' is the number of mutations in the V gene region only
  
  if(disable_grouping == T){ # computes distribution across the entire tibble (used to analyze separate naive seq datasets)
    grouping_vars <- n_mutations_variable
    if('mouse_id' %in% names(annotated_seqs)){
      grouping_vars <- c('mouse_id', grouping_vars)
    }
    
  }else{
    grouping_vars <- c('mouse_id','day','infection_status','group_controls_pooled',
                       'tissue', 'cell_type', n_mutations_variable)
  }
  
  dist_n_mutations <- annotated_seqs %>%
    group_by(across(grouping_vars)) %>%
    dplyr::summarise(n_seqs = dplyr::n()) %>%
    ungroup() %>%
    group_by(across(grouping_vars[grouping_vars!= n_mutations_variable])) %>%
    mutate(compartment_seqs = sum(n_seqs),
           obs_fraction = n_seqs / compartment_seqs) %>%
    ungroup()
  
  if(!disable_grouping){
    dist_n_mutations <- dist_n_mutations %>%
      mutate(cell_type = factor(cell_type, levels = c('naive','GC','PC','mem')),
             tissue = factor(tissue, levels = c('LN','spleen','BM')))
  }

  return(dist_n_mutations)
  
}

# Gets distribution of sequence length by mouse, tissue and cell type
get_seq_length_distribution <- function(annotated_seqs, seq_length_variable, disable_grouping = F){
  # seq_length_variable: e.g. 'seq_length_partis' or 'sequenced_bases_in_vgene_region_partis'

  if(disable_grouping == T){ # computes distribution across the entire tibble (used to analyze separate naive seq datasets)
    grouping_vars <- seq_length_variable
    if('mouse_id' %in% names(annotated_seqs)){
      grouping_vars <- c('mouse_id', grouping_vars)
    }
  }else{
    grouping_vars <- c('mouse_id', 'day', 'infection_status', 'group_controls_pooled', 'tissue', 'cell_type',
                       seq_length_variable)
  }
  
  seq_length_dist <- annotated_seqs %>%
    group_by(across(grouping_vars)) %>%
    dplyr::summarise(n_seqs = n()) %>%
    ungroup() %>%
    group_by(across(grouping_vars[grouping_vars!= seq_length_variable])) %>%
    mutate(compartment_seqs = sum(n_seqs),
           obs_fraction = n_seqs / compartment_seqs) %>%
    ungroup() 
  
  
  if(!disable_grouping){
    seq_length_dist <- seq_length_dist %>%
      mutate(cell_type = factor(cell_type, levels = c('naive','GC','PC','mem')),
             tissue = factor(tissue, levels = c('LN','spleen','BM')))
    
  }
  
  return(seq_length_dist)
  
}

# Calculates expected null distribution of mutations given the observed distribution of sequence lengths
get_null_mutation_distribution_given_length_distribution <- function(annotated_seqs, n_mutations_variable, seq_length_variable,
                                                                     estimated_seq_error_rate, disable_grouping = F){
  
  if(disable_grouping == T){
    grouping_vars <- 'n_mutations'
    if('mouse_id' %in% names(annotated_seqs)){
      grouping_vars <- c('mouse_id', grouping_vars)
    }
  }else{
    grouping_vars <- c('mouse_id','tissue','cell_type','n_mutations')
  }
  
  
  # Get observed distribution of sequence lengths 
  seq_length_dist <- get_seq_length_distribution(annotated_seqs, seq_length_variable, disable_grouping)
  
  # For a range of possible sequence lengths, calculate probability of observing 0,1,2,... mutations
  null_model_base_probs <- generate_mutation_null_model(annotated_seqs, estimated_seq_error_rate, n_mutations_variable = n_mutations_variable,
                                                        seq_length_variable = seq_length_variable)
  
  names(null_model_base_probs)[names(null_model_base_probs) == 'length'] <- seq_length_variable
  
  # Calculate expected null distribution of mutations given the observed distribution of sequence lengths
  null_distribution_given_obs_lengths <- left_join(seq_length_dist  %>%
                                                     select(matches('mouse_id'), matches('tissue'), matches('cell_type'), matches(seq_length_variable), obs_fraction),
                                                   null_model_base_probs) %>%
    group_by(across(grouping_vars)) %>%
    # For each number of mutations, calculate null probability as 
    # a weighted average across the observed freq distribution of lengths
    dplyr::summarise(null_prob = sum(obs_fraction*null_prob)) %>%
    ungroup()
  
  sums <- null_distribution_given_obs_lengths %>% group_by(across(grouping_vars[grouping_vars != 'n_mutations'])) %>%
    dplyr::summarise(S = sum(null_prob)) %>% ungroup() %>% select(S) %>% unique() %>% pull(S)
  
  stopifnot(all(abs(sums-1) < 1e-5))
  
  return(null_distribution_given_obs_lengths)
  
}

get_clone_size_distribution <- function(seq_counts){
  seq_counts %>%
  group_by(mouse_id, day, infection_status, group, group_controls_pooled, tissue, cell_type, clone_id) %>%
  dplyr::summarise(clone_prod_seqs = sum(uniq_prod_seqs)) %>%
  group_by(mouse_id, day, infection_status, group, group_controls_pooled, tissue, cell_type) %>%
  mutate(clone_freq = clone_prod_seqs / sum(clone_prod_seqs)) 
}

estimate_mut_probs_per_vgene_position <- function(annotated_seqs){
  
  # For a vector with position-specific mutation information for a single mouse / cell type / tissue / v gene,
  # counts how many times each position was mutated
  base_function <- function(annotated_seqs, selected_mouse, selected_tissue, selected_cell_type, selected_v_gene){
    vgene_specific_mutation_list <- annotated_seqs %>% filter(mouse_id == selected_mouse, cell_type == selected_cell_type,
                                             tissue == selected_tissue, v_gene == selected_v_gene) %>%
      pull(vgene_mutations_list_partis_nt)
    
    concatenated_string <- paste0(vgene_specific_mutation_list, collapse = ';')
    
    unique_mutated_positions <- unique(str_extract_all(concatenated_string, '[0-9]+')[[1]])
    
    mutation_counts <- lapply(as.list(unique_mutated_positions),
                              FUN = function(pos){str_count(concatenated_string, pattern = paste0('[A-Z\\*]',pos,'[A-Z\\*]'))})
    
    return(tibble(mouse_id = selected_mouse,
                  tissue = selected_tissue,
                  cell_type = selected_cell_type,
                  v_gene = selected_v_gene,
                  n_vgene_seqs_in_compartment = length(vgene_specific_mutation_list),
                  position = as.integer(unique_mutated_positions),
                  n_mutations = as.integer(mutation_counts)) %>%
             mutate(mutation_freq = n_mutations / length(vgene_specific_mutation_list)) %>%
             arrange(position)
    )
  }
  
  unique_combinations <- annotated_seqs %>%
    select(mouse_id, tissue, cell_type, v_gene) %>%
    unique()

  output <- mapply(FUN = base_function,
         selected_mouse = unique_combinations$mouse_id,
         selected_tissue = unique_combinations$tissue,
         selected_cell_type = unique_combinations$cell_type,
         selected_v_gene = unique_combinations$v_gene,
         MoreArgs = list(annotated_seqs = annotated_seqs), 
         SIMPLIFY = F)
  
  return(bind_rows(output))
  
}

get_mutation_frequencies_within_clones <- function(annotated_seqs, seq_counts, by_tissue_and_cell_type){
  
  grouping_vars <- c('mouse_id', 'clone_id')
  
  # IF this is true, will compute mutation frequencies relative to the number of productive sequences from the 
  # clone in each cell type / tissue combination
  if(by_tissue_and_cell_type){
    grouping_vars <- c(grouping_vars, 'tissue', 'cell_type')
    compartment_totals <- seq_counts
  }else{
    # Otherwise, compute frequencies relative to the clone's total number of productive sequences across tissues and cell types
    compartment_totals <- seq_counts %>%
      group_by(mouse_id, clone_id) %>%
      dplyr::summarise(prod_seqs = sum(prod_seqs)) %>%
      ungroup()
  }
  
  # Count how many times each mutation occurs in each clone
  mutation_counts_within_clones <- annotated_seqs %>%
    filter(productive_partis) %>%
    filter(!is.na(vgene_mutations_list_partis_aa)) %>%
    group_by(across(grouping_vars)) %>%
    dplyr::summarise(mutation = str_split(paste0(vgene_mutations_list_partis_aa, collapse = ';'), ';')[[1]]) %>%
    ungroup() %>%
    group_by(across(c(grouping_vars, 'mutation'))) %>%
    dplyr::summarise(n_seqs_with_mutation = dplyr::n()) %>%
    ungroup()
  
  # Arrange mutations by the position number in which they occurred
  mutation_counts_within_clones <- mutation_counts_within_clones %>%
    mutate(position = as.integer(str_remove_all(mutation, '[A-Z]'))) %>%
    arrange(across(c(grouping_vars, 'position')))
  
  
  # Count what fraction of the clone's sequences each mutation is present in
  mutation_freqs_within_clones <- left_join(mutation_counts_within_clones, compartment_totals) %>%
    dplyr::rename(prod_clone_seqs_in_compartment = prod_seqs) %>%
    mutate(mutation_freq_in_compartment = n_seqs_with_mutation / prod_clone_seqs_in_compartment)
  
  if(by_tissue_and_cell_type){
    mutation_freqs_within_clones <- mutation_freqs_within_clones %>%
      dplyr::rename(compartment_tissue = tissue, compartment_cell_type = cell_type)
  }
  
  return(mutation_freqs_within_clones)

}

# For each clone, finds list of V gene mutations that are present at or above the chosen frequency threshold.
list_clone_mutations_above_threshold <- function(mutation_freqs_within_clones, threshold){
  
  grouping_vars <- c('mouse_id', 'clone_id')
  
  if(all(c('compartment_tissue','compartment_cell_type') %in% names(mutation_freqs_within_clones))){
    grouping_vars <- c(grouping_vars, 'compartment_tissue', 'compartment_cell_type')
  }else{
    stopifnot(!any(c('compartment_tissue','compartment_cell_type') %in% names(mutation_freqs_within_clones)))
  }
  
  mutations_above_threshold <-  mutation_freqs_within_clones %>%
    filter(mutation_freq_in_compartment >= threshold) %>%
    group_by(across(c(grouping_vars, 'prod_clone_seqs_in_compartment'))) %>%
    summarise(mutations_above_threshold = paste0(mutation, collapse = ';'),
              mutation_freqs = paste0(round(mutation_freq_in_compartment,2), collapse = ';')) %>%
    mutate(frequency_threshold = threshold) %>%
    dplyr::rename(n_seqs_in_denominator = prod_clone_seqs_in_compartment) %>%
    ungroup() %>%
    select(mouse_id, clone_id, matches('tissue'), matches('cell_type'), frequency_threshold,
           mutations_above_threshold, n_seqs_in_denominator)
    
  return(mutations_above_threshold)
}

# Clustering based on Spearman correlation of V gene frequencies as an inverse measure of distance.
get_vgene_freq_correlation_clustering <- function(pairwise_correlations, cell_type, tissue, metric,
                                                  min_seqs){
  
  
  cell_type_full_name <- case_when(
    cell_type == 'PC' ~ 'plasma cells',
    cell_type == 'GC' ~ 'germinal center cells',
    cell_type == 'mem' ~ 'memory cells'
  )
  
  
  if(metric == 'freqs'){
    data_subset <- pairwise_correlations$freqs %>%
      rename(cor_coef = cor_coef_freqs)
    plot_title <- paste0('Mouse-pair correlations in ',
                         str_sub(cell_type_full_name,1,-2),
                         '\nV gene frequencies')
  }else{
    stopifnot(metric == 'freq_ratios')
    data_subset <- pairwise_correlations$freq_ratios %>%
      rename(cor_coef = cor_coef_freq_ratios)
    plot_title <- paste0('Mouse-pair correlations in frequency deviations\nbetween ',
           cell_type_full_name,
           ' and the naive repertoire')
    }
  
  data_subset <- data_subset %>%
    filter(total_compartment_seqs_i >= min_seqs, total_compartment_seqs_j >= min_seqs) %>%
    filter(cell_type == !!cell_type, !str_detect(pair_type, 'control')) 
  
  if(!is.null(tissue)){
    data_subset <- data_subset %>%
      filter(tissue == !!tissue)
  }
  
  data_subset <- data_subset %>%
    select(mouse_id_i, mouse_id_j, cor_coef)
  
  # Add a diagonal to the correlation matrix (each mouse has correlation 1 with itself)
  data_subset <- bind_rows(data_subset, 
                           tibble(mouse_id_i = unique(c(data_subset$mouse_id_i, data_subset$mouse_id_j))) %>%
                             mutate(mouse_id_j = mouse_id_i, cor_coef = 1)) %>%
    arrange(mouse_id_i, mouse_id_j)
  
  wide_format_correlations <- data_subset %>%
    pivot_wider(names_from = mouse_id_j, values_from = cor_coef)
  
  correlation_matrix <- as.matrix(wide_format_correlations[colnames(wide_format_correlations) != 'mouse_id_i'])
  rownames(correlation_matrix) <- wide_format_correlations$mouse_id_i
  
  # Fill lower triangle
  for(i in 1:nrow(correlation_matrix)){
    for(j in 1:ncol(correlation_matrix)){
      if(j < i){
        correlation_matrix[i,j] <- correlation_matrix[j,i]
      }
    }
  }
  
  # Convert correlations into distances.
  dist_matrix <- as.dist(1 - correlation_matrix)
  cluster <- hclust(dist_matrix, method = 'complete')
  dendrogram <- as.dendrogram(cluster)
  
  # Within topological constraints of dendrogram, tries to order mice with day 8 on top, then day 16, 
  #etc.
  leaf_weights <- 1/as.integer(str_extract(cluster$labels, '[0-9]+'))
  dendrogram <- reorder(dendrogram,
                        wts = leaf_weights)
  
  annotation <- get_info_from_mouse_id(
    tibble(mouse_id = cluster$labels)
  )
  
  margin_colors <- left_join(annotation, group_controls_pooled_palette)$group_color
  
  rwb <- colorRampPalette(colors = c("blueviolet", "white", "darkorange2"))
  heatmap.2(correlation_matrix,
            trace = 'none',
            #col = rev(brewer.pal(11, name = 'RdBu')),
            col = rwb(100),
            Rowv = dendrogram,
            Colv = rev(dendrogram),
            RowSideColors = margin_colors,
            ColSideColors = margin_colors,
            key.xlab = 'Spearman correlation',
            key.ylab = '',
            denscol = 'black',
            key.title = '',
            main = plot_title,
            xlab = 'Individual mice',
            ylab = 'Individual mice')


  
  #return()
  
  
}

randomize_noncontrol_groups <- function(tibble_with_mouse_id){
  
  # Distributes non-control mice randomly between groups (primary-8, secondary-40, etc.), keeping the original number of mice in each group.
  
  mouse_info <- get_info_from_mouse_id(tibble_with_mouse_id %>% select(mouse_id) %>% unique())
  non_id_vars <- names(mouse_info)[names(mouse_info) != 'mouse_id']
  
  
  noncontrols <- mouse_info %>% filter(group_controls_pooled != 'control')
  controls <- mouse_info %>% filter(group_controls_pooled == 'control')
  
  randomized_noncontrols <- noncontrols %>%
    select(mouse_id, group_controls_pooled) %>%
    mutate(group_controls_pooled = sample(group_controls_pooled, size = length(group_controls_pooled), replace = F)) %>%
    mutate(group = group_controls_pooled,
           day = str_extract(group_controls_pooled,'[0-9]+'),
           infection_status = str_extract(group_controls_pooled, '[a-z]*'))

  
  randomized_mouse_info <- bind_rows(controls, randomized_noncontrols)
  

  randomized_tibble <- left_join(tibble_with_mouse_id %>% select(-any_of(non_id_vars)),
                                 randomized_mouse_info, by = 'mouse_id') %>%
    select(mouse_id, any_of(non_id_vars), everything())
  
  return(randomized_tibble)
  
  
}








