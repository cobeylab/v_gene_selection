library(tidyr)
library(dplyr)
library(cowplot)
library(stringr)
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

# Counts productive sequences by mouse, clone, tissue, cell type. 
get_productive_seq_counts <- function(annotated_seqs, unique_only, assignment){
  
  if(assignment == 'partis'){
    counts <- annotated_seqs %>%
      dplyr::rename(clone_id = clone_id_partis,
                    productive = productive_partis)
  }else{
    if(assignment == 'partis_ogrdb'){
      counts <- annotated_seqs %>%
        dplyr::rename(clone_id = clone_id_partis_ogrdb,
                      productive = productive_partis_ogrdb) %>%
        mutate(partis_uniq_ref_seq = partis_ogrdb_uniq_ref_seq)
    }else{
      stopifnot(assignment == 'igblast')
      
      counts <- annotated_seqs %>%
        dplyr::rename(clone_id = clone_id_igblast,
                      productive = productive_igblast) 
      
    }
    
  }
  
  counts <- counts %>% filter(productive) 
  
  # If computing counts of unique sequences only:
  if(unique_only){
    counts <- counts %>%
      # Partis's assignment of a reference id to multiple identical sequences can be used 
      # when computing unique sequences for clones defined using other methods.
      select(mouse_id, clone_id, partis_uniq_ref_seq, tissue, cell_type, isotype) %>%
      unique() %>%
      group_by(mouse_id, clone_id, tissue, cell_type) %>%
      dplyr::summarise(unique_prod_seqs = n()) %>%
      ungroup()
  # Counts of all productive sequences
  }else{
    counts <- counts %>%
      group_by(mouse_id, clone_id, tissue, cell_type) %>%
      summarise(prod_seqs = n()) %>%
      ungroup()
  }
  
  return(counts)
  
}

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

clone_purity_test <- get_clone_purity(tibble(mouse_id = 'A', clone_id = c(rep(1,4),rep(2,4)),
                                             cell_type = rep(c('IgD+B220+','GC','PC','mem'),2),
                                             prod_seqs = c(1,1,1,1,1,0,0,0)))
stopifnot(clone_purity_test$clone_purity == c('mixed','pure_IgD+B220+'))
rm(clone_purity_test)

# Selects sequences sorted as naive based on additional filters (n. mutations, isotype, clone size)
process_IgD_B220_seqs <- function(annotated_seqs, max_clone_unique_IgDB220_seqs,
                                  max_v_gene_mutations, filters_from){
  
  # Naive sequences are those that satisfy ALL the following:
  # - were sorted as naive (DUMP-IgD+B220+)
  # - Have IgD or IgM as the isotype determined from the read (or isotype is NA)
  # - Are in a clone that only has DUMP-IgD+B220+ cells
  # - Are in a clone with at most max_clone_naive_seqs **unique** seqs, where unique = same mouse, clone, tissue, isotype and sequence
  # - Has max_v_gene_mutations or fewer mutations in the V gene region.
  
  # Sequences sorted as DUMP-IgD+B220+ that do not meet all the other criteria are labeled 'nonnaive_IgD+B220+'
  
  # filters_from specificies the assignment that will be used to check this criteria
  stopifnot(filters_from == 'partis' | filters_from == 'partis_ogrdb')
  
  unique_productive_seq_counts <- get_productive_seq_counts(annotated_seqs, unique_only = T, assignment = filters_from)
    
  clone_purity <- get_clone_purity(unique_productive_seq_counts) %>%
    dplyr::rename(unique_productive_IgDB220_seqs_in_clone = IgD_B220_seqs_in_clone,
                  unique_productive_nonIgDB220_seqs_in_clone = non_IgD_B220_seqs_in_clone)
  
  names(clone_purity)[names(clone_purity) == 'clone_id'] <- paste0('clone_id_', filters_from)
  
  annotated_seqs <- left_join(annotated_seqs, clone_purity)
  
  igd_b220_seqs <- annotated_seqs %>% 
    filter(cell_type == 'IgD+B220+')
  
  if(filters_from == 'partis' | filters_from == 'igblast'){
    igd_b220_seqs <- igd_b220_seqs %>%
      mutate(mutations_filter = vgene_mutations_partis_nt)
  }else{
    stopifnot(filters_from == 'partis_ogrdb')
    igd_b220_seqs <- igd_b220_seqs %>%
      mutate(mutations_filter = vgene_mutations_partis_ogrdb_nt)
  }

  igd_b220_seqs <- igd_b220_seqs %>%
    mutate(cell_type = case_when(
      (isotype %in% c('IGM','IGD') | is.na(isotype)) & clone_purity == "pure_IgD+B220+" &
        unique_productive_IgDB220_seqs_in_clone <= max_clone_unique_IgDB220_seqs &
        mutations_filter <= max_v_gene_mutations ~ 'naive',
      !is.na(clone_purity) ~ 'nonnaive_IgD+B220+',
      T ~ 'unassigned_IgD+B220+' # Will contain IgD+B220+ in clones with only unproductive seqs
    ))
  
  processed_tibble <- bind_rows(annotated_seqs %>% filter(cell_type != 'IgD+B220+'),
            igd_b220_seqs %>% select(all_of(names(annotated_seqs))))
  
  stopifnot(nrow(processed_tibble) == nrow(annotated_seqs))
  
  return(processed_tibble)
  
}

process_IgD_B220_test <- process_IgD_B220_seqs(tibble(mouse_id = 'A',
                                                      clone_id_partis = c(1,1,1,1,2,3,4), 
                                                      partis_uniq_ref_seq = c(1,1,1,2,3,4,5),
                                                      cell_type = c('GC',rep('IgD+B220+',6)),
                                                      isotype = c(rep('IGM',6), 'IGG'),
                                                      vgene_mutations_partis_nt = c(0,0,0,0,0,3,0),
                                                      tissue = 'LN', productive_partis = T),
       max_clone_unique_IgDB220_seqs = 1, max_v_gene_mutations = 2, filters_from = 'partis')
stopifnot(process_IgD_B220_test$cell_type == c('GC', rep('nonnaive_IgD+B220+', 3), 'naive', rep('nonnaive_IgD+B220+', 2)))
rm(process_IgD_B220_test)

# Computes fraction of seqs in each clone that belong to each cell type or to each tissue
get_clone_composition <- function(seq_counts, composition_var){
  
  stopifnot(composition_var %in% names(seq_counts) & composition_var %in% c("tissue","cell_type"))

  # Handles input tibbles with unique prod seq counts instead of total seq counts.
  if('unique_prod_seqs' %in% names(seq_counts)){
    stopifnot(('prod_seqs' %in% names(seq_counts)) == F)
    names(seq_counts)[names(seq_counts) == 'unique_prod_seqs'] <- 'prod_seqs'
    col_prefix <- 'unique_seqs_'
  }else{
    col_prefix <- 'prod_seqs_'
  }

  # Computing clone composition
  clone_composition <- seq_counts %>%
    group_by(across(c('mouse_id', 'clone_id', all_of(composition_var)))) %>%
    # For each clone, sum across cell types within each class of the composition variable
    dplyr::summarise(prod_seqs_in_class = sum(prod_seqs)) %>%
    ungroup() %>%
    group_by(mouse_id, clone_id) %>%
    # Now compute the fraction of sequences in a clone that came from each class
    mutate(total_clone_prod_seqs = sum(prod_seqs_in_class)) %>%
    ungroup() %>%
    pivot_wider(id_cols = any_of(c('mouse_id','clone_id','total_clone_prod_seqs')),
                names_from = any_of(composition_var), values_from = prod_seqs_in_class,
                values_fill = 0, names_prefix = col_prefix) %>%
    # Add rel. frequencies in addition to numbers of each cell type / tissue
    mutate(across(matches(col_prefix), function(x){x/total_clone_prod_seqs},
                  .names = "fraction_{.col}")) 
  
  if(col_prefix == 'unique_seqs_'){
    clone_composition <- clone_composition %>%
      dplyr::rename(total_clone_unique_seqs = total_clone_prod_seqs)
  }

  return(clone_composition)
}


# Annotate clone_info objects with output from get_clone_composition
annotate_clone_info_with_composition <- function(clone_info, seq_counts, unique_seq_counts, assignment){
  clone_tissue_composition_prod_seqs <- get_clone_composition(seq_counts, composition_var = 'tissue')
  clone_cell_type_composition_prod_seqs <- get_clone_composition(seq_counts, composition_var = 'cell_type')
  
  # Total number of productive sequences should be equal for each clone in both tibbles
  stopifnot(all(clone_tissue_composition_prod_seqs$total_clone_prod_seqs ==
                  clone_cell_type_composition_prod_seqs$total_clone_prod_seqs))
  
  clone_info <- left_join(clone_info, clone_tissue_composition_prod_seqs)
  clone_info <- left_join(clone_info, clone_cell_type_composition_prod_seqs %>% select(-total_clone_prod_seqs))
  

  clone_tissue_composition_unique_seqs <- get_clone_composition(unique_seq_counts, composition_var = 'tissue')
  clone_cell_type_composition_unique_seqs <- get_clone_composition(unique_seq_counts, composition_var = 'cell_type')
  
  stopifnot(all(clone_cell_type_composition_unique_seqs$total_clone_unique_seqs ==
                  clone_tissue_composition_unique_seqs$total_clone_unique_seqs))
  
  clone_info <- left_join(clone_info, clone_tissue_composition_unique_seqs)
  clone_info <- left_join(clone_info, clone_cell_type_composition_unique_seqs %>% select(-total_clone_unique_seqs))
  
  # Find frequency of most common tissue or cell type in each clone. 
  clone_info <- clone_info %>%
    rowwise() %>%
    mutate(biggest_tissue_fraction_prod_seqs = max(c(fraction_prod_seqs_BM, fraction_prod_seqs_LN, fraction_prod_seqs_spleen)),
           biggest_tissue_fraction_unique_seqs =  max(c(fraction_unique_seqs_BM, fraction_unique_seqs_LN, fraction_unique_seqs_spleen)),
           biggest_cell_type_fraction_prod_seqs = max(c(fraction_prod_seqs_naive, `fraction_prod_seqs_nonnaive_IgD+B220+`,
                                                        fraction_prod_seqs_GC, fraction_prod_seqs_PC, fraction_prod_seqs_mem)),
           biggest_cell_type_fraction_unique_seqs = max(c(fraction_unique_seqs_naive, `fraction_unique_seqs_nonnaive_IgD+B220+`,
                                                          fraction_unique_seqs_GC, fraction_unique_seqs_PC, fraction_unique_seqs_mem))) %>%
    ungroup()
  
}



# Functions for calculating germline allele frequencies
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

# Combine naive and experienced frequencies into a single tibble in wide format
format_gene_freqs_wide <- function(exp_freqs, naive_freqs){
  left_join(exp_freqs,naive_freqs) %>%
    mutate(tissue = factor(tissue, levels = c('LN','spleen','BM')),
           group_controls_pooled = factor(group_controls_pooled, levels = group_controls_pooled_factor_levels))
}


# Rho: ratio of experienced to naive frequency
compute_rho <- function(gene_freqs){
  gene_freqs %>%
    mutate(log_rho = log(vgene_seq_freq) - log(naive_vgene_seq_freq)) %>%
    mutate(obs_rho = exp(log_rho)) 
}


# Correlation within each mouse between experienced and naive gene frequencies
get_naive_exp_correlations <- function(gene_freqs){
  gene_freqs %>%
    group_by(mouse_id, day, infection_status, group, group_controls_pooled, tissue, cell_type,
             total_compartment_seqs, total_mouse_naive_seqs) %>%
    dplyr::summarise(spearman = cor.test(vgene_seq_freq, naive_vgene_seq_freq, method = 'spearman')$estimate,
                     pearson = cor.test(vgene_seq_freq, naive_vgene_seq_freq, method = 'pearson')$estimate) %>%
    ungroup() %>%
    pivot_longer(cols = c('spearman','pearson'), names_to = 'method', values_to = 'naive_exp_corr')
}


# Generates vector of unique mouse pairs, using gene_freqs as a basis
get_unique_pairs <- function(gene_freqs, within_groups_only){
  
  # Tibble containing group of each mouse.
  mouse_group_assignments <- gene_freqs %>% select(mouse_id, group_controls_pooled) %>% unique()
  
  unique_pairs <- gene_freqs %>% select(mouse_id) %>% unique() %>%
    dplyr::rename(mouse_id_i = mouse_id) %>%
    mutate(mouse_id_j = mouse_id_i) %>%
    complete(mouse_id_i, mouse_id_j) %>%
    rowwise() %>%
    mutate(pair = paste0(sort(c(mouse_id_i, mouse_id_j)), collapse = ';')) %>%
    ungroup() %>%
    filter(mouse_id_i != mouse_id_j) %>%
    unique()
  
  # Add info on mouse groups 
  unique_pairs <- left_join(unique_pairs,  mouse_group_assignments %>%
              dplyr::rename(mouse_id_i = mouse_id, group_controls_pooled_i = group_controls_pooled))
  unique_pairs <- left_join(unique_pairs, mouse_group_assignments %>%
                              dplyr::rename(mouse_id_j = mouse_id, group_controls_pooled_j = group_controls_pooled))
  
  # If true, only return pairs of the same group (e.g. both controls, both primary-8, etc.)
  if(within_groups_only){
    unique_pairs <- unique_pairs %>%
      filter(group_controls_pooled_i == group_controls_pooled_j)
  }
  
  unique_pairs <- unique_pairs %>%
    select(pair) %>%
    unique() %>%
    pull(pair)
  
  return(unique_pairs)
  
}

# For all pairs of mice, rearrange freqs tibble to show changes in each mouse as different vars.
get_pairwise_freqs <- function(gene_freqs, adjust_naive_zeros, within_groups_only){
  # Assumes gene_freqs in wide format:
  stopifnot('naive_vgene_seq_freq' %in% names(gene_freqs))
  
  mouse_info <- gene_freqs %>%
    select(mouse_id, day, infection_status, group, group_controls_pooled) %>%
    unique()
  
  unique_pairs <- get_unique_pairs(gene_freqs = gene_freqs, within_groups_only = within_groups_only)
   

  # Get naive and experienced frequencies into separate tibbles
  exp_freqs <- gene_freqs %>% select(mouse_id, day, infection_status, group, group_controls_pooled,
                                     v_gene, tissue, cell_type, n_vgene_seqs, total_compartment_seqs,
                                     vgene_seq_freq, obs_rho, matches('rho'), matches('deviation_from_naive'),
                                     matches('v_gene_rank'))
  
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
  if('deviation_from_naive_i' %in% names(paired_gene_freqs)){
    # Label each gene in a paper in terms of concordant/discordant direction of change from the naive repertoire
    paired_gene_freqs <- paired_gene_freqs %>%
      mutate(concordance_status = case_when(
        (deviation_from_naive_i == 'positive' & deviation_from_naive_j == 'positive') ~ 'concordant-increasing',
        (deviation_from_naive_i == 'negative' & deviation_from_naive_j == 'negative') ~ 'concordant-decreasing',
        (deviation_from_naive_i == 'neutral' & deviation_from_naive_j == 'neutral') ~ 'concordant-stable',
        (deviation_from_naive_i != deviation_from_naive_j) ~ 'discordant',
        (is.na(deviation_from_naive_i) | is.na(deviation_from_naive_j)) ~ 'NA'))
  }
      
    
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
sample_experienced_from_naive <- function(exp_seq_counts, naive_seq_counts, clone_info, synth_data_input_tibble, by_tissue = T, n_reps = 100){
  
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
    
    simulated_counts <- simulated_freqs  %>%
      group_by(across(any_of(c('mouse_id','tissue','cell_type')))) %>%
      mutate(n_vgene_seqs = rmultinom(1, size = unique(total_compartment_seqs), prob = naive_vgene_seq_freq)) %>%
      ungroup() %>%
      mutate(vgene_seq_freq = n_vgene_seqs / total_compartment_seqs) %>%
      select(mouse_id, day, infection_status, group, group_controls_pooled, cell_type, matches('tissue'), v_gene,
             n_naive_vgene_seqs, total_mouse_naive_seqs, naive_vgene_seq_freq, n_vgene_seqs ,total_compartment_seqs, vgene_seq_freq)
    
    return(simulated_counts)
  }
  
  realizations <- replicate(n_reps, generate_replicate(naive_freqs = naive_freqs, exp_freqs = exp_freqs),
                            simplify = F)  
  
  replicates_tibble <- bind_rows(realizations, .id = 'replicate') %>% select(replicate, everything())
  
  replicates_tibble <- compute_rho(replicates_tibble)
  
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

# Computes frequency of each clone within each compartment
# (tissue, cell type, tissue/cell type combination, or whole mouse)
get_clone_freqs <- function(seq_counts, compartment_vars){
  
  # Handles input with unique instead of total productive seq counts
  if('unique_prod_seqs' %in% names(seq_counts)){
    stopifnot(('prod_seqs' %in% names(seq_counts)) == F)
    seq_counts <- seq_counts %>% rename(prod_seqs = unique_prod_seqs)
  }
  
  grouping_vars <- c('mouse_id','day','infection_status', 'group', 'group_controls_pooled')

  # compartment_vars is null if computing clone frequencies across entire mouse.
  if(!is.null(compartment_vars)){
    stopifnot(compartment_vars %in% c('tissue','cell_type'))
    grouping_vars <- c(grouping_vars, compartment_vars)
  }
  
  clone_freqs <- seq_counts %>%
    group_by(across(any_of(c(grouping_vars, 'clone_id')))) %>%
    # Sum clone sequences after grouping by compartment vars to get n clone seqs in compartment
    summarise(n_clone_seqs_in_compartment = sum(prod_seqs)) %>%
    ungroup() %>%
    # Now compute total seqs in compartment across clones, use as denominator
    group_by(across(any_of(grouping_vars))) %>%
    mutate(total_seqs_in_compartment = sum(n_clone_seqs_in_compartment),
           clone_freq_in_compartment = n_clone_seqs_in_compartment / total_seqs_in_compartment) %>%
    mutate(clone_rank_in_compartment = rank(-clone_freq_in_compartment, ties.method = 'first')) %>%
    ungroup() %>%
    arrange(mouse_id, clone_id)
  
  names(clone_freqs)[names(clone_freqs) %in% compartment_vars] <- paste0('compartment_', compartment_vars)
  return(clone_freqs)

}

estimate_mut_probs_per_vgene_position <- function(annotated_seqs, clone_info_partis, is_ogrdb_run){
  
  # Calculate site-specific mutation frequencies for each V gene
  # Only available for the partis assignments
  # (For these calculations, exclude unproductive sequences and sequences without an assigned V gene)
  
  input_data <- annotated_seqs
  
  if(is_ogrdb_run){
    partis_vars <- names(input_data)[str_detect(names(input_data), 'partis')]
    excluded_vars <- partis_vars[str_detect(partis_vars, 'ogrdb') == F]
    
    input_data <- input_data %>%
      select(-all_of(excluded_vars))
    
    names(input_data) <- str_replace(names(input_data), 'partis_ogrdb','partis')
  }
  
  input_data <- input_data %>% filter(!is.na(n_mutations_partis_nt), !is.na(vgene_mutations_partis_nt), productive_partis == T)
  
  input_data <- input_data %>%
    dplyr::rename(clone_id = clone_id_partis) %>%
    mutate(across(c('clone_id','partis_uniq_ref_seq','seq_id'), as.character))
  
  input_data <- left_join(input_data %>% 
                                mutate(clone_id = as.character(clone_id)),
                              clone_info_partis %>% select(mouse_id, clone_id, v_gene))
  
  
  # For a vector with position-specific mutation information for a single mouse / cell type / tissue / v gene,
  # counts how many times each position was mutated
  base_function <- function(input_data, selected_mouse, selected_tissue, selected_cell_type, selected_v_gene){
    vgene_specific_mutation_list <- input_data %>% filter(mouse_id == selected_mouse, cell_type == selected_cell_type,
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
  
  unique_combinations <- input_data %>%
    select(mouse_id, tissue, cell_type, v_gene) %>%
    unique()

  output <- mapply(FUN = base_function,
         selected_mouse = unique_combinations$mouse_id,
         selected_tissue = unique_combinations$tissue,
         selected_cell_type = unique_combinations$cell_type,
         selected_v_gene = unique_combinations$v_gene,
         MoreArgs = list(input_data = input_data), 
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

# Simple function to count mutations from character vector output by list_clone_mutations_above_threshold
# (adds this count as a column to the input tibble)
count_mutations_above_threshold <- function(tibble_with_mutations_list_vector){
  stopifnot('mutations_above_threshold' %in% names(tibble_with_mutations_list_vector))
  tibble_with_mutations_list_vector %>% 
    # If mutations_above_threshold is not '', there's at least one mutation. Plus 1 for each ';' in mutations_above_threshold
    mutate(n_mutations_above_threshold = (mutations_above_threshold != '') + str_count(mutations_above_threshold,';'))
}

# From a tibble of clone frequencies annotated with n. mutations above threshold, count what fraction of clones has at least target_n_mutations above threshold
get_fraction_of_clones_with_mutations_above_threshold <- function(clone_freqs, target_n_mutations, min_clone_size){
  # If input tibble is by compartment, use total number of clones in the compartment (with at least min_clone_size) as denominator
  compartment_vars <- c('compartment_tissue','compartment_cell_type')
  
  grouping_vars <- c('mouse_id', 'day', 'infection_status', 'group','group_controls_pooled',
                     compartment_vars[compartment_vars %in% names(clone_freqs)])

  return(
    clone_freqs %>%
      filter(n_clone_seqs_in_compartment >= min_clone_size) %>%
      group_by(across(grouping_vars)) %>%
      summarise(n_clones_denominator = n(),
                n_clones_meeting_target_mutations = sum(n_mutations_above_threshold >= target_n_mutations)) %>%
      ungroup() %>%
      mutate(fraction_clones_with_at_least_n_mutations = n_clones_meeting_target_mutations / n_clones_denominator)
  )
}

# From a tibble of clone freqs. annotated with a character vector listing mutations above threshold, count mutations shared by clone pairs
# (pairs of clones in the same compartment and with the same v allele, but from either the same mouse or different mice from the same group)
count_mutations_shared_by_clone_pairs <- function(clone_freqs){
  
  compartment_vars <- c('compartment_tissue','compartment_cell_type')
  
  # NOT grouping by mouse id, so we compare clones with same V allele in different mice. (but from the same cell type/tissue)
  grouping_vars <- c('day', 'infection_status', 'group_controls_pooled', 'v_gene',
                     compartment_vars[compartment_vars %in% names(clone_freqs)])
  
  
  # Internal vectorized function. Compares every element of mutations list column against the others
  comparison_function <- function(mutations_list){
    if(length(mutations_list) > 1){
      paired_lists <- t(combn(mutations_list, 2))
      colnames(paired_lists) <- c('mutations_list_clone_i', 'mutations_list_clone_j')
      paired_lists <- as_tibble(paired_lists)
      
      
      compare_mutation_lists <- function(mutations_list_i, mutations_list_j){
        split_list_i <- str_split(mutations_list_i,';')
        split_list_j <- str_split(mutations_list_j,';')
        
        output <- mapply(list_i = split_list_i, list_j = split_list_j,
                         FUN = function(list_i, list_j){
                           shared_mutations = list_i[list_i %in% list_j]
                           
                           return(length(shared_mutations[!is.na(shared_mutations) & shared_mutations != '']))
                         }
        )
        return(output)
      }
      
      shared_mutations = paired_lists %>% mutate(shared_mutations = compare_mutation_lists(mutations_list_i = mutations_list_clone_i,
                                                                                           mutations_list_j = mutations_list_clone_j)) 
    }else{
      shared_mutations <- tibble(mutations_list_clone_i = c(), mutations_list_clone_j = c(), shared_mutations = c())
    }  
    return(shared_mutations)
  }
  # Quick tests for internal function
  stopifnot(comparison_function(c('A123B;C456D','A123B;C456E'))%>% pull(shared_mutations) == 1)
  stopifnot(comparison_function(c('A123B;C456D','A123B;C456D'))%>% pull(shared_mutations) == 2)
  stopifnot(comparison_function(c('A123B;C456D',''))%>% pull(shared_mutations) == 0)
  stopifnot(comparison_function(c('',''))%>% pull(shared_mutations) == 0)
  stopifnot(comparison_function(c(NA,NA))%>% pull(shared_mutations) == 0)
  
  # Apply internal function after grouping 
  return(
    clone_freqs %>%
      group_by(across(grouping_vars)) %>%
      dplyr::summarise(comparison_function(mutations_list = mutations_above_threshold)) %>%
      ungroup()
  )
}

# Using the output of count_mutations_shared_by_clone_pairs, calculates the fraction of clones (same V allele, compartment, group; same or different mice)
get_fraction_clones_with_shared_mutations <- function(shared_mutations){
  
  compartment_vars <- c('compartment_tissue','compartment_cell_type')
  
  # NOT grouping by mouse id, so we compare clones with same V allele in different mice. (but from the same cell type/tissue)
  grouping_vars <- c('day', 'infection_status', 'group_controls_pooled',
                     compartment_vars[compartment_vars %in% names(shared_mutations)])
  
  shared_mutations %>%
    mutate(any_mutations_shared = shared_mutations > 0) %>%
    group_by(across(c(grouping_vars,'any_mutations_shared'))) %>%
    count() %>%
    ungroup() %>%
    group_by(across(grouping_vars)) %>%
    mutate(total_n_clone_pairs = sum(n),
           probability = n/total_n_clone_pairs) %>%
    filter(any_mutations_shared == T)
}


get_freq_ratio_mutability_correlations <- function(gene_freqs, germline_mutability_by_region_type,
                                                   min_compartment_size, method){
  left_join(gene_freqs, germline_mutability_by_region_type) %>%
    filter(total_compartment_seqs >= min_compartment_size, total_mouse_naive_seqs >= min_compartment_size) %>%
    group_by(mouse_id, day, infection_status, group, group_controls_pooled, tissue, cell_type, total_compartment_seqs) %>%
    summarise(across(matches('mutability'),
                     function(x){cor.test(x, obs_rho, method = method)$estimate})) %>%
    pivot_longer(cols = matches('_mutability'), names_to = 'mutability_metric', values_to = 'correlation')
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


# For each clone, samples a new V gene from the naive repertoire of the corresponding mouse (with replacement)
resample_clone_v_alleles <- function(exp_seq_counts, naive_freqs, adjust_naive_zeros){
  
  if(adjust_naive_zeros){
    naive_freqs <- adjust_zero_naive_freqs(naive_freqs)
  }
  
  clone_v_genes <- exp_seq_counts %>% select(mouse_id, clone_id, v_gene) %>% unique()
  
  # Does the resampling for data, where data is a grouped tibble with each group containing a single mouse
  internal_resampling_function <- function(mouse_id, naive_freqs){
    
    n_clones <- length(mouse_id)
    mouse_id = unique(mouse_id)
    # Ensures there's only a single mouse
    stopifnot(length(mouse_id) == 1)
    
    mouse_specific_freqs <- naive_freqs %>% filter(mouse_id == !!mouse_id) %>%
      select(mouse_id, v_gene, naive_vgene_seq_freq)
    
    new_v_genes <- sample(mouse_specific_freqs$v_gene, size = n_clones,
                                 prob = mouse_specific_freqs$naive_vgene_seq_freq,
                                 replace = T) 
    
    return(new_v_genes)
      
  }
  
  # A built-in test to see if the grouping is done correctly: rename all genes as the id of corresponding mouse
  # Run internal resampling function and check that each mouse only has genes labelled with its name
  internal_function_test <- clone_v_genes %>% group_by(mouse_id) %>%
    mutate(new_v_gene = internal_resampling_function(mouse_id = mouse_id, 
                                                     naive_freqs = naive_freqs %>% mutate(v_gene = mouse_id)))
    
  stopifnot(all(internal_function_test$mouse_id == internal_function_test$new_v_gene))
  

  # run resampling for each mouse in clone_v_genes
  resampled_v_genes <- clone_v_genes %>% group_by(mouse_id) %>%
    mutate(new_v_gene = internal_resampling_function(mouse_id = mouse_id, 
                                                     naive_freqs = naive_freqs))
  

  # Return sequence counts with resampled v genes
  return(
    left_join(exp_seq_counts, resampled_v_genes) %>%
    select(-v_gene) %>%
    dplyr::rename(v_gene = new_v_gene) %>%
    select(mouse_id, clone_id, v_gene, everything())
  )
  
}

# Uses resample_clone_v_alleles to generate gene freqs under a null model where lineages are randomly assigned V alleles based on naive freqs
gene_freqs_random_lineage_alleles <- function(exp_seq_counts, naive_seq_counts, naive_freqs, adjust_naive_zeros){
  
  # Randomize lineage alleles in exp_seq_counts object
  null_model_exp_counts <- resample_clone_v_alleles(exp_seq_counts = exp_seq_counts, 
                                                    naive_freqs = naive_freqs, 
                                                    adjust_naive_zeros = adjust_naive_zeros)
  
  # Calculate allele (gene) frequencies
  null_model_gene_freqs <- calc_gene_freqs(exp_seq_counts = null_model_exp_counts,
                                           naive_seq_counts = naive_seq_counts,
                                           clone_info = clone_info, long_format = F, by_tissue = T)
  
  if(adjust_naive_zeros){
    naive_freqs <- adjust_zero_naive_freqs(null_model_gene_freqs$naive_freqs)
  }else{
    naive_freqs <- null_model_gene_freqs$naive_freqs
  }
  
  # Convert gene frequencies to wide format
  null_model_gene_freqs <- format_gene_freqs_wide(exp_freqs = null_model_gene_freqs$exp_freqs,
                                                  naive_freqs = naive_freqs)
  
  # Add values of rho (ratio of experienced to naive frequencies)
  null_model_gene_freqs <- compute_rho(null_model_gene_freqs)
  
  return(null_model_gene_freqs)

}

# Applies functions for getting pairwise correlations to list with multiple gene freqs. objects (obtained by randomization)
get_pairwise_corrs_for_randomized_datasets <- function(randomized_gene_freqs_list, adjust_naive_zeros, within_groups_only,
                                                       focal_tissue = NULL){
  
  # Master objects containing pw correlation in gene frequencies and experienced-to-naive ratios.
  pw_corr_freqs <- tibble()
  pw_corr_freq_ratios <- tibble()
  
  # For each object in input, convert to pairwise format
  pw_freqs <- lapply(randomized_gene_freqs_list, FUN = get_pairwise_freqs,
                     adjust_naive_zeros = adjust_naive_zeros, within_groups_only = within_groups_only)
  
  if(!is.null(focal_tissue)){
    stopifnot(length(focal_tissue) == 1)
    pw_freqs <- lapply(pw_freqs, FUN = function(x){x %>% filter(tissue == focal_tissue)})
  }
  
  
  # For each paired gene freqs object, compute pairwise correlations, add to master objects
  # More efficient do to loop here because get_pairwise_correlations outputs a list
  for(i in 1:length(pw_freqs)){
    pw_corr <- get_pairwise_correlations(pairwise_gene_freqs = pw_freqs[[i]])
    pw_corr_freqs <- bind_rows(pw_corr_freqs, pw_corr$freqs %>% mutate(replicate = i))
    pw_corr_freq_ratios <- bind_rows(pw_corr_freq_ratios, pw_corr$freq_ratios %>% mutate(replicate = i))
  }
  
  pw_corr_freqs <- pw_corr_freqs %>% select(replicate, everything())
  pw_corr_freq_ratios <- pw_corr_freq_ratios %>% select(replicate, everything())
  
  return(list(freqs = pw_corr_freqs, freq_ratios = pw_corr_freq_ratios))
  
}

export_partis_germline_genes <- function(yaml_object, mouse_id, output_dir, is_ogrdb_run){
  
  run_id <- paste0('partis', ifelse(is_ogrdb_run,'_ogrdb',''))
  
  write.fasta(yaml_object$`germline-info`$seqs$v, names = names(yaml_object$`germline-info`$seqs$v),
              file.out = paste0(output_dir, 'v_genes_', run_id,'_', mouse_id,'.fasta'))
  write.fasta(yaml_object$`germline-info`$seqs$d, names = names(yaml_object$`germline-info`$seqs$d),
              file.out = paste0(output_dir, 'd_genes_', run_id,'_', mouse_id,'.fasta'))
  write.fasta(yaml_object$`germline-info`$seqs$j, names = names(yaml_object$`germline-info`$seqs$j),
              file.out = paste0(output_dir, 'j_genes_', run_id,'_', mouse_id,'.fasta'))
  
  
}



