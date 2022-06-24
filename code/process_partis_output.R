source('gene_frequency_functions.R')
source('partis_output_functions.R')

# Annotates mouse-specific data files with partis annotations / exports a file with clone information a seq-level file with some annotations

args = commandArgs(trailingOnly = T)
mouse_yaml_file_path = args[1] # e.g. mouse_yaml_file_path = '../results/partis/8-5_partis.yaml'
mouse_data_file_path = args[2] # e.g. mouse_data_file_path = '../processed_data/mouse_specific_data_files/8-5.csv'


merge_info <- function(yaml_object, mouse_data_file_path){
 
  # Read mouse sequences
  mouse_raw_data <- as_tibble(read.csv(mouse_data_file_path))
  
  # Get partis information into desired format
  partis_info <- format_partis_info(yaml_object)
  stopifnot(length(unique(partis_info$clone_id_partis)) == length(yaml_object$events))
  
  # Find mouse id 
  mouse_id <- unique(mouse_raw_data$participant_alt_label)
  
  # Check the sequence data was from a single mouse.
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
  
  # Export file with sequence-level annotations
  annotated_seqs <- merged_data %>% mutate(mouse_id = mouse_id) %>%
    select(mouse_id, clone_id_partis, seq_id, partis_uniq_ref_seq, specimen_tissue, specimen_cell_subset, isotype, seq_length_partis,
           productive_partis, n_mutations_partis_nt, n_mutations_partis_aa, cdr3_seq_partis, cdr3_mutations_partis_nt, cdr3_mutations_partis_aa,
           vgene_mutations_partis_nt, sequenced_bases_in_vgene_region_partis, vgene_mutations_list_partis_nt, vgene_mutations_list_partis_aa,
           partis_processed_seq) %>%
    dplyr::rename(tissue = specimen_tissue, cell_type = specimen_cell_subset, clone_id = clone_id_partis) %>%
    mutate(cell_type = as.character(cell_type))
  
  annotated_seqs$cell_type[annotated_seqs$cell_type == 'naïve'] <- 'IgD+B220+'
  
  write.csv(annotated_seqs, paste0('../processed_data/annotated_seq_files/', mouse_id, '_annotated_seqs.csv'),
            row.names = F)
  
  # Export clone info file, including some summary statistics
  clone_info <- merged_data %>% mutate(mouse_id = mouse_id) %>%
    select(mouse_id, clone_id_partis, v_segment_partis, d_segment_partis, j_segment_partis, cdr3_length_partis,
           clone_consensus_cdr3_partis, clone_naive_cdr3_partis, clone_naive_seq_nt_partis, clone_naive_seq_aa_partis) %>%
    unique() %>%
    dplyr::rename(clone_id = clone_id_partis, v_gene = v_segment_partis, j_gene = j_segment_partis,
                  d_gene = d_segment_partis) %>%
    mutate(clone_id = as.character(clone_id))
  
  
  clones_summary <- annotated_seqs %>%
    filter(productive_partis == T) %>%
    group_by(mouse_id, clone_id) %>%
    dplyr::summarise(mean_n_mutations_partis_aa = mean(n_mutations_partis_aa),
              mean_cdr3_mutations_partis_aa = mean(cdr3_mutations_partis_aa),
              median_n_mutations_partis_aa = median(n_mutations_partis_aa),
              median_cdr3_mutations_partis_aa = median(cdr3_mutations_partis_aa),
              max_n_mutations_partis_aa = max(n_mutations_partis_aa),
              max_cdr3_mutations_partis_aa = max(cdr3_mutations_partis_aa),
              min_n_mutations_partis_aa = min(n_mutations_partis_aa),
              min_cdr3_mutations_partis_aa = min(cdr3_mutations_partis_aa)) %>%
    ungroup() %>%
    mutate(clone_id = as.character(clone_id))
  
  clone_info <- left_join(clone_info, clones_summary)
  
  
  write.csv(clone_info, paste0('../processed_data/clone_info_files/', mouse_id, '_clone_info.csv'),
            row.names = F)
  
  # Calculate site-specific mutation frequencies for each V gene.
  # (For these calculations, exclude unproductive sequences and sequences without an assigned V gene)
  
  annotated_seqs <- annotated_seqs %>% filter(!is.na(n_mutations_partis_nt), !is.na(vgene_mutations_partis_nt), productive_partis == T)

  annotated_seqs <- annotated_seqs %>%
    mutate(across(c('clone_id','partis_uniq_ref_seq','seq_id'), as.character))
  
  annotated_seqs <- left_join(annotated_seqs %>% 
                                mutate(clone_id = as.character(clone_id)),
                              clone_info %>% select(mouse_id, clone_id, v_gene))
  
  mut_probs_per_gene_position <- estimate_mut_probs_per_vgene_position(annotated_seqs)
  
  write.csv(mut_probs_per_gene_position, paste0('../results/mutations_per_vgene_base/', mouse_id, '_mutations_per_vgene_base.csv'),
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
