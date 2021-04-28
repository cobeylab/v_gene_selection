library(plyr)
library(dplyr)
library(tidyr)
library(yaml)
library(stringr)
library(seqinr)
library(Biostrings)
library(msa)

# Annotates mouse-specific data files with partis annotations / exports a file with clone information a seq-level file with some annotations

args = commandArgs(trailingOnly = T)
mouse_yaml_file_path = args[1] # e.g. mouse_yaml_file_path = '../results/partis/8-5_partis.yaml'
mouse_data_file_path = args[2] # e.g. mouse_data_file_path = '../processed_data/mouse_specific_data_files/8-5.csv'

# Creates tibble linking duplicated sequences to those chosen as the representative sequence
pair_duplicates <- function(clone){
  duplicate_seqs <- c()
  for(i in 1:length(clone$duplicates)){
    representative_seq = clone$unique_ids[i]
    duplicates <- clone$duplicates[[i]]
    if(length(duplicates) > 0){
      duplicate_seqs <- bind_rows(duplicate_seqs,
                                  tibble(duplicate_seq = duplicates, representative_seq))
    }
  }
  return(duplicate_seqs)
}

count_mutations_from_naive <- function(seq_vector, naive_seq, translate_seqs, is_partis_aa_alignment = F){
  # Will assume seqs are NT seqs unless is_partis_aa_alignment is TRUE.
  # Will translate if translate_seqs is TRUE and is_partis_aa_alignment is FALSE.
  if(translate_seqs){
    naive_seq <- c2s(seqinr::translate(s2c(tolower(naive_seq))))
    seq_vector <- sapply(seq_vector,
                         FUN = function(seq){
                           c2s(seqinr::translate(s2c(tolower(seq))))
                         }, USE.NAMES = F)
    
  }
  if(is_partis_aa_alignment){
    stopifnot(translate_seqs == F) # Sequences already supposed to be AA.
    # If comparing aligned AA sequences output by partis, ignore 'X's and gaps
    count_aa_diffs_ignoring_x <- function(naive_seq, seq){
      naive_seq_char_vector <- s2c(naive_seq)
      seq_char_vector <- s2c(seq)
      
      mutated_sites <- which((naive_seq_char_vector %in% c('X','*') == F) & (seq_char_vector %in% c('X','*') ==F) &
        (naive_seq_char_vector != seq_char_vector))
      return(length(mutated_sites))
    }
    n_mutations <- sapply(seq_vector, FUN = count_aa_diffs_ignoring_x, USE.NAMES = F,
                          naive_seq = naive_seq)
    
    
  }else{
    n_mutations <- sapply(seq_vector,
                          FUN = function(seq){
                            seq_char_vector = s2c(seq)
                            naive_seq_char_vector <- s2c(naive_seq)
                            
                            #Some unproductive sequences will not be successfully aligned to the naive seq. by partis.
                            if(length(naive_seq_char_vector) != length(seq_char_vector)){
                              return(NA)
                            }else{
                              mutated_sites <- which((naive_seq_char_vector != 'N') & (seq_char_vector != 'N') &
                                                       (naive_seq_char_vector != seq_char_vector))
                              return(length(mutated_sites))
                            }
                          })
    
  }
  
  return(n_mutations)
}


# Identifies nt and aa mutations (instead of simply counting them, e.g. K123Q)
identify_vgene_mutations <- function(naive_v_gene_region_seq, v_gene_region_seqs){
  vgene_mutations_aa <- rep('', length(v_gene_region_seqs))
  vgene_mutations_nt <- rep('', length(v_gene_region_seqs))
  
  # Some inferred naive V gene regions will start with 'N' nucleotides. If X is the number of Ns they start with, remove the first X letters from all seqs
  Ns_in_naive_v_region <- str_locate_all(naive_v_gene_region_seq,'N')[[1]]
  if(length(Ns_in_naive_v_region) != 0){
    naive_v_gene_region_seq <- str_sub(naive_v_gene_region_seq, Ns_in_naive_v_region[2] +1)
    v_gene_region_seqs <- str_sub(v_gene_region_seqs, Ns_in_naive_v_region[2] +1)
  }
  
  complete_codon_positions <- nchar(naive_v_gene_region_seq) %/% 3
  
  # Efficiently distributes sequences into a matrix, replacing missing characters in short sequences with 'N'
  character_matrix <- rbind.fill.matrix(lapply(as.list(c(naive_v_gene_region_seq, v_gene_region_seqs)),
                                FUN = function(x){matrix(s2c(x), byrow = T, nrow = 1)}))
  character_matrix[is.na(character_matrix)] <- 'N'
  
  # Retain sequences up to last complete codon
  character_matrix <- character_matrix[, 1:(complete_codon_positions * 3)]
  
  
  # Find all nucleotide posititions where naive sequence isn't 'N' and at least one observed sequence is mutated
  position_is_mutated <- apply(character_matrix, MARGIN = 2, FUN = function(column){
    naive_char <- column[1]
    non_n_chars <- unique(column)[unique(column) != 'N']
    
    return(naive_char != 'N' & length(non_n_chars) > 1)})
  mutated_nt_positions <- (1:ncol(character_matrix))[position_is_mutated]
  
  mutated_codon_positions <- unique(((mutated_nt_positions - 1) %/% 3) + 1)
  
  if(length(mutated_codon_positions) > 0){
    for(codon in mutated_codon_positions){
      
      codon_nt_positions <- seq(codon*3 - 2, codon * 3)
      
      # Record nucleotide mutations
      for(nt_pos in codon_nt_positions){
        naive_nt <- character_matrix[1, nt_pos]
        obs_nts <- character_matrix[-1, nt_pos]
        new_mutations_nt <- ifelse(obs_nts %in% c('N', naive_nt),'', paste0(';',naive_nt, nt_pos, obs_nts))
        vgene_mutations_nt <- paste0(vgene_mutations_nt, new_mutations_nt)
      }
      
      naive_codon <- paste0(character_matrix[1, codon_nt_positions], collapse = '')
      naive_aa <- as.character(Biostrings::translate(DNAString(naive_codon), if.fuzzy.codon = 'X', no.init.codon = T))
      
      # rbind is necessary here because R will annoyingly treat a single row taken from a matrix as a character vector, crashing apply
      obs_codons <- apply(rbind(character_matrix[-1, codon_nt_positions]), MARGIN = 1, FUN = paste0, collapse = '')
      obs_aas <- as.character(Biostrings::translate(DNAStringSet(obs_codons), if.fuzzy.codon = 'X', no.init.codon = T))
      
      # Disregard mutations to X (i.e., sequence did not cover this part of inferred naive V gene seq.)
      new_mutations_aa <- ifelse(obs_aas %in% c('X', naive_aa),'', paste0(';',naive_aa, codon,obs_aas))
      
      vgene_mutations_aa <- paste0(vgene_mutations_aa, new_mutations_aa)
    }
  }
  
  # Remove leading ';'
  vgene_mutations_aa[vgene_mutations_aa != ''] <- str_sub(vgene_mutations_aa[vgene_mutations_aa != ''], 2)
  vgene_mutations_nt[vgene_mutations_nt != ''] <- str_sub(vgene_mutations_nt[vgene_mutations_nt != ''], 2)
  
  return(list(nt = vgene_mutations_nt, aa = vgene_mutations_aa))

}


format_partis_info <- function(yaml_object){
  partis_info <- c()
  for(i in 1:length(yaml_object$events)){
    print(paste0('Processing clone ', i, ' of ', length(yaml_object$events), ' (', round(100*i/length(yaml_object$events)),'%)'))
    clone <- yaml_object$events[[i]]
    
    stopifnot(length(unique(nchar(clone$cdr3_seqs))) == 1)
    
    # Productive status for each unique sequence in clone:
    productive_partis <- (clone$stops == F)&(clone$in_frames)&(clone$mutated_invariants == F)
    
    
    #write.fasta(lapply(as.list(clone$input_seqs), FUN = s2c), 
    #            names = clone$unique_ids, file.out = '~/Desktop/test.fasta')
    
    cdr3_start <- clone$codon_positions$v + 1 # Adds one because partis numbering starts at 0
    cdr3_end <- clone$codon_positions$j + 3 # adds +1 and then +2 to go to end of codon
    naive_cdr3_seq <- str_sub(clone$naive_seq, cdr3_start, cdr3_end)
    naive_cdr3_seq_aa <- c2s(seqinr::translate(s2c(tolower(naive_cdr3_seq))))
    
    v_gene_bounds <- clone$regional_bounds$v
    naive_v_gene_region_seq <- str_sub(clone$naive_seq, v_gene_bounds[1], v_gene_bounds[2])
    v_gene_region_seqs <- str_sub(clone$input_seqs, v_gene_bounds[1], v_gene_bounds[2])
    
    
    # Number of mutations in V gene region
    vgene_mutations_nt <- count_mutations_from_naive(seq_vector = v_gene_region_seqs, naive_seq = naive_v_gene_region_seq,
                                                       translate_seqs = F)
    # Number of bases effectively sequenced in V gene region for each sequence (i.e. excluding "N"'s introduced by partis)
    sequenced_bases_in_vgene_region <- nchar(str_remove_all(v_gene_region_seqs,'N'))
    
    # Precise list of mutations in V gene region for each sequences
    vgene_mutations_list <- identify_vgene_mutations(naive_v_gene_region_seq = naive_v_gene_region_seq,
                             v_gene_region_seqs = v_gene_region_seqs)
    
    vgene_mutations_list_nt <- vgene_mutations_list$nt
    vgene_mutations_list_aa <- vgene_mutations_list$aa
    
    
    stopifnot(length(unique(nchar(clone$cdr3_seqs))) == 1)
    stopifnot(unique(nchar(clone$cdr3_seqs)) == nchar(naive_cdr3_seq))
    
    cdr3_mutations_nt <- count_mutations_from_naive(seq_vector = clone$cdr3_seqs,
                                                    naive_seq = naive_cdr3_seq, translate_seqs = F)
    
    cdr3_mutations_aa <- count_mutations_from_naive(seq_vector = clone$cdr3_seqs,
                                                    naive_seq = naive_cdr3_seq, translate_seqs = T)
    
    n_mutations_partis_nt <- clone$n_mutations
    n_mutations_partis_aa <- count_mutations_from_naive(seq_vector = clone$seqs_aa, naive_seq = clone$naive_seq_aa,
                                                        translate_seqs = F, is_partis_aa_alignment = T)
    
    
    cdr3_seq_partis <- sapply(clone$cdr3_seqs,
                              FUN = function(seq){
                                c2s(seqinr::translate(s2c(tolower(seq))))
                              }, USE.NAMES = F)
    
    # Get sequence info for reference seqs (partis clusters identical seqs. choosing one as the ref)
    ref_seq_info <- tibble(seq = clone$unique_ids, seq_length_partis = nchar(str_remove_all(clone$input_seqs,'N')),
                           productive_partis = productive_partis,
                           n_mutations_partis_nt = n_mutations_partis_nt, 
                           n_mutations_partis_aa = n_mutations_partis_aa,
                           cdr3_length_partis = clone$cdr3_length / 3,
                           cdr3_seq_partis = cdr3_seq_partis,
                           cdr3_mutations_partis_nt = cdr3_mutations_nt,
                           cdr3_mutations_partis_aa = cdr3_mutations_aa,
                           vgene_mutations_partis_nt = vgene_mutations_nt,
                           sequenced_bases_in_vgene_region_partis =sequenced_bases_in_vgene_region,
                           vgene_mutations_list_partis_nt = vgene_mutations_list_nt,
                           vgene_mutations_list_partis_aa = vgene_mutations_list_aa
                           )
    
    # Get productivity of duplicate sequences by referring to the productivity of their corresponding reference seq.
    duplicate_seqs <- pair_duplicates(clone)
    
    if(is.null(duplicate_seqs) == F){
      duplicate_seq_info <- left_join(duplicate_seqs, 
                                          ref_seq_info %>% dplyr::rename(representative_seq = seq),
                                          by = 'representative_seq')
      
      stopifnot(all((duplicate_seq_info$duplicate_seq %in% ref_seq_info$seq) == F))
      
      # Combine both representative and duplicate sequences
      clone_tibble <- bind_rows(ref_seq_info %>% mutate(partis_uniq_ref_seq = seq) %>%
                                  select(seq, partis_uniq_ref_seq, everything()),
                                duplicate_seq_info %>% 
                                  dplyr::rename(seq = duplicate_seq, partis_uniq_ref_seq = representative_seq) %>%
                                  select(seq, partis_uniq_ref_seq, matches('partis')))
    }else{
      clone_tibble <- ref_seq_info %>%
        mutate(partis_uniq_ref_seq = seq) %>% select(seq, partis_uniq_ref_seq, everything())
    }
    
    # Add consensus CDR3 seq from partis
    cdr3_start_aa <- ceiling(cdr3_start/3)
    cdr3_end_aa <- ceiling(cdr3_end/3)
    consensus_cdr3_partis <- str_sub(clone$consensus_seq_aa, start = cdr3_start_aa, end = cdr3_end_aa)
    
    # Add germline genes and clone id
    clone_tibble <- clone_tibble %>%
      mutate(clone_id_partis = i - 1,  # Will keep partis clone numbering with zero indexing.
             v_segment_partis = clone$v_gene,
             d_segment_partis = clone$d_gene,
             j_segment_partis = clone$j_gene,
             v_segment_support_partis = clone$v_per_gene_support[[1]],
             clone_consensus_cdr3_partis = consensus_cdr3_partis,
             clone_naive_cdr3_partis = naive_cdr3_seq_aa,
             clone_naive_seq_nt_partis = clone$naive_seq,
             clone_naive_seq_aa_partis = clone$naive_seq_aa) %>%
      dplyr::rename(trimmed_read_id = seq)
    
        
    # For consistency with the Boyd lab notation, replace TRUE/FALSE with 't'/'f'
    clone_tibble$productive_partis[clone_tibble$productive_partis == T] <- 't'
    clone_tibble$productive_partis[clone_tibble$productive_partis == F] <- 'f'
    
    partis_info <- bind_rows(partis_info,
                             clone_tibble)
  
  }
  return(partis_info)
}


merge_info <- function(yaml_object, mouse_data_file_path){
 
  mouse_raw_data <- as_tibble(read.csv(mouse_data_file_path))
  
  # Get partis information into desired format
  partis_info <- format_partis_info(yaml_object)
    stopifnot(length(unique(partis_info$clone_id_partis)) == length(yaml_object$events))
  
  # Find mouse id 
  mouse_id <- unique(mouse_raw_data$participant_alt_label)
  stopifnot(length(mouse_id) == 1)
  
  # Remove mouse id from read ids in partis info
  partis_info$trimmed_read_id <- str_replace(partis_info$trimmed_read_id,paste0(mouse_id,'_'),'')
  
  # Mouse raw data contains raw data sent from the Boyd lab. Add clone id, V/D/J assignment with those from partis
  merged_data <- left_join(mouse_raw_data %>% mutate(trimmed_read_id = as.character(trimmed_read_id)),
                           partis_info, by = c('trimmed_read_id')) %>%
    dplyr::rename(seq_id = trimmed_read_id)
  
  merged_data <- merged_data %>%
    dplyr::rename(v_segment_igblast = v_segment, d_segment_igblast = d_segment, j_segment_igblast = j_segment,
           clone_id_igblast = igh_igblast_clone_id, productive_igblast = productive,
           v_segment_support_igblast = v_score) 
  
  # Overwrite mouse-specific data file with new file including partis annotation.
  write.csv(merged_data %>% select(-clone_naive_cdr3_partis, -clone_consensus_cdr3_partis, -cdr3_mutations_partis_nt, -cdr3_mutations_partis_aa,
                     -n_mutations_partis_nt, -n_mutations_partis_aa, -clone_naive_seq_nt_partis, -clone_naive_seq_aa_partis,-vgene_mutations_list_partis_nt,
                     -vgene_mutations_list_partis_aa,-vgene_mutations_partis_nt, -sequenced_bases_in_vgene_region_partis),
            mouse_data_file_path, row.names = F)
  
  seq_level_file <- merged_data %>% mutate(mouse_id = mouse_id) %>%
    select(mouse_id, clone_id_partis, seq_id, partis_uniq_ref_seq, specimen_tissue, specimen_cell_subset, isotype, seq_length_partis,
           productive_partis, n_mutations_partis_nt, n_mutations_partis_aa, cdr3_seq_partis, cdr3_mutations_partis_nt, cdr3_mutations_partis_aa,
           vgene_mutations_partis_nt, sequenced_bases_in_vgene_region_partis, vgene_mutations_list_partis_nt, vgene_mutations_list_partis_aa)
  write.csv(seq_level_file, paste0('../processed_data/annotated_seq_files/', mouse_id, '_annotated_seqs.csv'),
            row.names = F)
  
  # Export clone info file
  clone_info_file <- merged_data %>% mutate(mouse_id = mouse_id) %>%
    select(mouse_id, clone_id_partis, v_segment_partis, d_segment_partis, j_segment_partis, cdr3_length_partis,
           clone_consensus_cdr3_partis, clone_naive_cdr3_partis, clone_naive_seq_nt_partis, clone_naive_seq_aa_partis) %>%
    unique()
  
  write.csv(clone_info_file, paste0('../processed_data/clone_info_files/', mouse_id, '_clone_info.csv'),
            row.names = F)
  
}


yaml_object <- read_yaml(mouse_yaml_file_path)
merge_info(yaml_object, mouse_data_file_path)

#Write germline gene sequences (including new alleles)
mouse_id = str_extract(mouse_yaml_file_path,'[0-9]*-[0-9]*')

write.fasta(yaml_object$`germline-info`$seqs$v, names = names(yaml_object$`germline-info`$seqs$v),
            file.out = paste0('../results/partis/partis_germline_genes/v_genes_', mouse_id,'.fasta'))
write.fasta(yaml_object$`germline-info`$seqs$d, names = names(yaml_object$`germline-info`$seqs$d),
            file.out = paste0('../results/partis/partis_germline_genes/d_genes_', mouse_id,'.fasta'))
write.fasta(yaml_object$`germline-info`$seqs$j, names = names(yaml_object$`germline-info`$seqs$j),
            file.out = paste0('../results/partis/partis_germline_genes/j_genes_', mouse_id,'.fasta'))
