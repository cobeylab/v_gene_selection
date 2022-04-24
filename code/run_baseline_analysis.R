library(readr)
library(dplyr)
library(tidyr)
library(shazam)
library(stringr)

args <- commandArgs(trailingOnly = T)
mouse_id = as.character(args[1])
nproc = as.integer(args[2])

baseline_test_statistic <- 'focused'
collapse_clones_threshold <- 0.5
focal_tissue = 'LN'
min_clone_size <- 10

# Read clone information to get clone's naive ancestral sequences
clone_info <- read_csv('../processed_data/clone_info.csv') %>% mutate(clone_id = as.character(clone_id))
# clone_info <- read_csv('~/Desktop/v_gene_selection/processed_data/clone_info.csv')  %>% mutate(clone_id = as.character(clone_id))

# Read annotated observed sequences
annotated_seqs <- read_csv('../processed_data/annotated_seqs.csv') %>%
  mutate(seq_id = as.character(seq_id),
         clone_id = as.character(clone_id))

# annotated_seqs <- read_csv('~/Desktop/v_gene_selection/processed_data/annotated_seqs.csv') %>% mutate(seq_id = as.character(seq_id), clone_id = as.character(clone_id))

# Read tsv with FR and CDR3 annotation for clone's naive ancestors
naive_ancestor_annotations <- read_tsv('../results/clone_naive_seqs_igblast.tsv')
# naive_ancestor_annotations <- read_tsv('~/Desktop/v_gene_selection/results/clone_naive_seqs_igblast.tsv')

naive_ancestor_annotations <- naive_ancestor_annotations %>%
  separate(sequence_id, into = c('mouse_id', 'clone_id'), sep = '_')

# Subset sequences to look only at clones from focal tissue (productive sequences only)
focal_seqs <- annotated_seqs %>% filter(tissue == focal_tissue, productive_partis == T, mouse_id == !!mouse_id,
                                        cell_type != 'naive') %>%
            select(mouse_id, clone_id, seq_id, tissue, cell_type, partis_processed_seq)

# Only consider clone/cell type combinations with at least min_clone_size sequences

focal_seqs <- focal_seqs %>%
  group_by(mouse_id, clone_id, tissue, cell_type) %>%
  mutate(n_clone_seqs_in_compartment = n()) %>%
  filter(n_clone_seqs_in_compartment >= min_clone_size) %>%
  ungroup() %>%
  select(-n_clone_seqs_in_compartment)

focal_seqs <- left_join(focal_seqs, 
                        clone_info %>% select(mouse_id, clone_id, clone_naive_seq_nt_partis)) %>%
  filter(nchar(partis_processed_seq) == nchar(clone_naive_seq_nt_partis))

# Use shazam's function to collapse clones into effective (representative) sequences
# This grouping is done by clone/cell type combination

focal_seqs <- focal_seqs %>%
  mutate(baseline_grouping_var = paste(mouse_id, clone_id, tissue, cell_type, sep = ';')) %>%
  select(baseline_grouping_var, everything())

clone_effective_seqs <- collapseClones(focal_seqs, sequenceColumn = "partis_processed_seq", cloneColumn = "baseline_grouping_var",
               germlineColumn = 'clone_naive_seq_nt_partis',
               method="thresholdedFreq", minimumFrequency=collapse_clones_threshold, includeAmbiguous=FALSE,
               breakTiesStochastic = FALSE, nproc = nproc) %>% select(-partis_processed_seq, -clonal_germline)

clone_effective_seqs <- left_join(clone_effective_seqs,
                                  naive_ancestor_annotations %>% select(mouse_id, clone_id, matches('fwr'), matches('cdr'))) %>%
  # In all sequences, remove N-sites introduced by partis at the beginning of germline sequences
  mutate(baseline_input_seq = str_sub(clonal_sequence, fwr1_start),
         baseline_input_naive_seq = str_sub(clone_naive_seq_nt_partis, fwr1_start)) %>%
  # Adjust region positions accordingly by subtracting start of FWR1
  mutate(across(matches(c('start','end')), function(x){x - fwr1_start + 1})) %>%
  select(-baseline_grouping_var)
  
# Test that splitting full seq according to FWR1 positions returns the assigned FWR1:
test_seqs <- clone_effective_seqs$baseline_input_naive_seq
stopifnot(str_sub(test_seqs, clone_effective_seqs$fwr1_start, clone_effective_seqs$fwr1_end) ==
            clone_effective_seqs$fwr1)
stopifnot(str_sub(test_seqs, clone_effective_seqs$cdr2_start, clone_effective_seqs$cdr2_end) ==
            clone_effective_seqs$cdr2)

# Add string column representing each site as either "fwr" or "cdr" to be passed to baseline
clone_effective_seqs <- clone_effective_seqs %>%
  rowwise() %>%
  mutate(regional_definition_string = 
           paste(c(rep('fwr', fwr1_end - fwr1_start + 1),
                   rep('cdr', cdr1_end - cdr1_start + 1),
                   rep('fwr', fwr2_end - fwr2_start + 1),
                   rep('cdr', cdr2_end - cdr2_start + 1),
                   rep('fwr', fwr3_end - fwr3_start + 1),
                   rep('cdr', cdr3_end - cdr3_start + 1)), collapse = ';')) %>%
  ungroup()

# Function for running baseline when each clone has potentially different FR/CDR boundaries
run_baseline <- function(clone_effective_seqs, nproc){
  
    baseline_list <- mapply(
      FUN = function(baseline_input_seq, baseline_input_naive_seq, regional_definition_string){
        regional_definition_factor <- str_split(regional_definition_string, ';')[[1]]
        regional_definition_factor <- factor(regional_definition_factor, levels = c('fwr','cdr'))
        regional_definition <- createRegionDefinition(name = "", boundaries = regional_definition_factor)
        
        input_data <- tibble(baseline_input_seq, baseline_input_naive_seq)
        
        baseline_results <- calcBaseline(input_data, sequenceColumn = 'baseline_input_seq', germlineColumn = 'baseline_input_naive_seq',
                     targetingModel = MK_RS5NF, regionDefinition = regional_definition, testStatistic = baseline_test_statistic,
                     nproc = nproc)
        return(baseline_results)
      },
      baseline_input_seq = clone_effective_seqs$baseline_input_seq, 
      baseline_input_naive_seq = clone_effective_seqs$baseline_input_naive_seq,
      regional_definition_string = clone_effective_seqs$regional_definition_string,
      SIMPLIFY = F
    )
    
    # Combine results for each object in the list to generate 1 master object
    combined_results <- baseline_list[[1]]
    
    for(i in 2:length(baseline_list)){
      next_object <- baseline_list[[i]]
      
      combined_results@db <- bind_rows(combined_results@db, next_object@db)
      combined_results@numbOfSeqs <- rbind(combined_results@numbOfSeqs, next_object@numbOfSeqs)
      combined_results@binomK <- rbind(combined_results@binomK, next_object@binomK)
      combined_results@binomN <- rbind(combined_results@binomN, next_object@binomN)
      combined_results@binomP <- rbind(combined_results@binomP, next_object@binomP)
      combined_results@pdfs$fwr <- rbind(combined_results@pdfs$fwr, next_object@pdfs$fwr)
      combined_results@pdfs$cdr <- rbind(combined_results@pdfs$cdr, next_object@pdfs$cdr)
      combined_results@stats <- bind_rows(combined_results@stats, next_object@stats)
    }
    
    stopifnot(nrow(combined_results@db) == nrow(clone_effective_seqs))
    
    # Add mouse and clone identifiers to db component of results object.
    combined_results@db <- combined_results@db %>% mutate(mouse_id = clone_effective_seqs$mouse_id,
                                                          clone_id = clone_effective_seqs$clone_id,
                                                          tissue = clone_effective_seqs$tissue,
                                                          cell_type = clone_effective_seqs$cell_type) %>%
      select(mouse_id, clone_id, tissue, cell_type, everything())
    
    return(combined_results)

}


baseline_results <- run_baseline(clone_effective_seqs, nproc = nproc)
stopifnot(length(unique(baseline_results@db$tissue)) == 1)

# Convolve results for this mice, separately for each cell type
convolved_baseline_results <- groupBaseline(baseline_results, groupBy = 'cell_type', nproc = nproc)

# Add number of clones per compartment to stats component of results

convolved_baseline_results@stats <- left_join(convolved_baseline_results@stats, 
          clone_effective_seqs %>% group_by(tissue, cell_type) %>% count() %>% dplyr::rename(n_clones_in_compartment = n))

save(convolved_baseline_results,
     file = paste0('../results/baseline_analysis/baseline_analysis_', mouse_id,'.RData'))


