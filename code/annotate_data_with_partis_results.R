library(dplyr)
library(tidyr)
library(yaml)
library(stringr)
library(seqinr)
library(Biostrings)

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

count_mutations_from_naive <- function(seq_vector, naive_seq, translate_seqs){
  # Will assume seqs are NT seqs. Will translate if translate_seqs is TRUE.
  if(translate_seqs){
    naive_seq <- c2s(seqinr::translate(s2c(tolower(naive_seq))))
    seq_vector <- sapply(seq_vector,
                         FUN = function(seq){
                           c2s(seqinr::translate(s2c(tolower(seq))))
                         }, USE.NAMES = F)
    
  }
  n_mutations <- sapply(seq_vector,
                        FUN = function(seq){
                          stringDist(c(naive_seq, seq),
                                     method = 'hamming')
                        })
  return(n_mutations)
}




format_partis_info <- function(yaml_object){
  partis_info <- c()
  for(i in 1:length(yaml_object$events)){
    #print(paste('clone', i-1))
    clone <- yaml_object$events[[i]]
    
    cdr3_start <- clone$codon_positions$v + 1 # Adds one because partis numbering starts at 0
    cdr3_end <- clone$codon_positions$j + 3 # adds +1 and then +2 to go to end of codon
    naive_cdr3_seq <- str_sub(clone$naive_seq, cdr3_start, cdr3_end)
    naive_cdr3_seq_aa <- c2s(seqinr::translate(s2c(tolower(naive_cdr3_seq))))
    
    stopifnot(length(unique(nchar(clone$cdr3_seqs))) == 1)
    stopifnot(unique(nchar(clone$cdr3_seqs)) == nchar(naive_cdr3_seq))
    
    cdr3_mutations_nt <- count_mutations_from_naive(seq_vector = clone$cdr3_seqs,
                                                    naive_seq = naive_cdr3_seq, translate_seqs = F)
    
    cdr3_mutations_aa <- count_mutations_from_naive(seq_vector = clone$cdr3_seqs,
                                                    naive_seq = naive_cdr3_seq, translate_seqs = T)
    
    # Get sequence info for reference seqs (partis clusters identical seqs. choosing one as the ref)
    ref_seq_info <- tibble(seq = clone$unique_ids, seq_length_partis = nchar(str_remove_all(clone$input_seqs,'N')),
                           n_mutations_partis = clone$n_mutations, 
                           productive_partis = (clone$stops == F)&(clone$in_frames)&(clone$mutated_invariants == F),
                           cdr3_length_partis = clone$cdr3_length,
                           cdr3_mutations_partis_nt = cdr3_mutations_nt,
                           cdr3_mutations_partis_aa = cdr3_mutations_aa)
    
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
 
    # Add germline genes and clone id
    clone_tibble <- clone_tibble %>%
      mutate(clone_id_partis = i - 1,  # Will keep partis clone numbering with zero indexing.
             v_segment_partis = clone$v_gene,
             d_segment_partis = clone$d_gene,
             j_segment_partis = clone$j_gene,
             v_segment_support_partis = clone$v_per_gene_support[[1]]) %>%
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
  # In the future can also add sequence regions (FRs, CDRs, etc.) from partis. For now, just keeping Igblast.
  merged_data <- left_join(mouse_raw_data %>% mutate(trimmed_read_id = as.character(trimmed_read_id)),
                           partis_info, by = c('trimmed_read_id'))
  
  merged_data <- merged_data %>%
    dplyr::rename(v_segment_igblast = v_segment, d_segment_igblast = d_segment, j_segment_igblast = j_segment,
           clone_id_igblast = igh_igblast_clone_id, productive_igblast = productive,
           v_segment_support_igblast = v_score)
  
  # Overwrite mouse-specific data file with new file including partis annotation.
  write.csv(merged_data, mouse_data_file_path, row.names = F)
  
  # Export the sequence level files
  seq_level_file <- merged_data %>% mutate(mouse_id = mouse_id) %>%
    select(mouse_id, clone_id_partis, trimmed_read_id, partis_uniq_ref_seq, specimen_tissue, specimen_cell_subset, isotype,
                         matches('partis')) 
  write.csv(seq_level_file, paste0('../processed_data/seq_level_files/', mouse_id, '_seq_level_file.csv'),
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
