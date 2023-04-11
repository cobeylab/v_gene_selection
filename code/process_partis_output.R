source('gene_frequency_functions.R')
source('partis_output_functions.R')

# Annotates mouse-specific data files with partis annotations / exports a file with clone information a seq-level file with some annotations

args = commandArgs(trailingOnly = T)
yaml_file_path = args[1] # e.g. yaml_file_path = '../results/partis/8-5_partis.yaml'
# Yaml file from the partis run with the OGRDB CB57/6 allele set:
yaml_file_path_ogrdb = args[2] #  yaml_file_path_ogrdb = '../results/partis/8-5_partis_ogrdb.yaml' 
mouse_data_file_path = args[3] # e.g. mouse_data_file_path = '../processed_data/mouse_specific_data_files/8-5.csv'

dir.create('../processed_data/annotated_seq_files/', showWarnings = F)
dir.create('../processed_data/clone_info_files/', showWarnings = F)
dir.create('../results/mutations_per_vgene_base/', recursive = T, showWarnings = F)
dir.create('../results/partis/partis_germline_genes/', showWarnings = F)

merge_info <- function(yaml_object, yaml_object_ogrdb, mouse_data_file_path){
 
  # Read mouse sequences
  mouse_raw_data <- as_tibble(read.csv(mouse_data_file_path))
  
  # Get partis information into desired format
  partis_info <- format_partis_info(yaml_object)
  
  partis_info_ogrdb <- format_partis_info(yaml_object_ogrdb)
  
  stopifnot(length(unique(partis_info$clone_id_partis)) == length(yaml_object$events))
  
  # Find mouse id 
  mouse_id <- unique(mouse_raw_data$participant_alt_label)
  
  # Check the sequence data was from a single mouse.
  stopifnot(length(mouse_id) == 1)
  
  # Remove mouse id from read ids in partis info
  partis_info$trimmed_read_id <- str_replace(partis_info$trimmed_read_id,paste0(mouse_id,'_'),'')
  partis_info_ogrdb$trimmed_read_id <- str_replace(partis_info_ogrdb$trimmed_read_id,paste0(mouse_id,'_'),'')
  
  # Mouse raw data contains raw data sent from the Boyd lab. Add clone id, V/D/J assignment with those from partis
  merged_data <- left_join(mouse_raw_data %>% mutate(trimmed_read_id = as.character(trimmed_read_id)),
                           partis_info, by = 'trimmed_read_id')
  
  # Merge with partis ogrdb info
  names(partis_info_ogrdb) <- str_replace(names(partis_info_ogrdb), 'partis', 'partis_ogrdb')
  
  merged_data <- left_join(merged_data,
            partis_info_ogrdb, by = 'trimmed_read_id') %>%
    dplyr::rename(seq_id = trimmed_read_id)
  
  merged_data <- merged_data %>%
    dplyr::rename(v_segment_igblast = v_segment, d_segment_igblast = d_segment, j_segment_igblast = j_segment,
           clone_id_igblast = igh_igblast_clone_id, productive_igblast = productive,
           v_segment_support_igblast = v_score) %>%
    mutate(v_segment_igblast = str_remove(as.character(v_segment_igblast),'mm'),
           d_segment_igblast = str_remove(as.character(d_segment_igblast),'mm'),
           j_segment_igblast = str_remove(as.character(j_segment_igblast),'mm')) %>%
    mutate(d_segment_igblast = str_remove(d_segment_igblast,"\""))
  
  annotated_seqs_partis_vars <- c('clone_id_partis', 'partis_uniq_ref_seq', 'seq_length_partis',
                                  'productive_partis', 'n_mutations_partis_nt', 'n_mutations_partis_aa', 'cdr3_seq_partis',
                                  'cdr3_mutations_partis_nt', 'cdr3_mutations_partis_aa', 'vgene_mutations_partis_nt',
                                  'sequenced_bases_in_vgene_region_partis', 'vgene_mutations_list_partis_nt',
                                  'vgene_mutations_list_partis_aa')
  
  annotated_seqs_partis_vars <- c(annotated_seqs_partis_vars,
                                  str_replace(annotated_seqs_partis_vars,'partis', 'partis_ogrdb'))
  
  # Export file with sequence-level annotations
  annotated_seqs <- merged_data %>% mutate(mouse_id = mouse_id) %>%
    select(mouse_id, seq_id, specimen_tissue, specimen_cell_subset, isotype,
           all_of(annotated_seqs_partis_vars),
           productive_igblast, clone_id_igblast) %>%
    dplyr::rename(tissue = specimen_tissue, cell_type = specimen_cell_subset) %>%
    mutate(cell_type = as.character(cell_type))
  
  annotated_seqs$cell_type[annotated_seqs$cell_type == 'naïve'] <- 'IgD+B220+'
  
  write.csv(annotated_seqs, paste0('../processed_data/annotated_seq_files/', mouse_id, '_annotated_seqs.csv'),
            row.names = F)
  
  # Export clone info files, including some summary statistics
  clone_info_partis <- get_partis_clone_info(merged_data = merged_data, is_ogrdb_run = F)
  clone_info_partis_ogrdb <- get_partis_clone_info(merged_data = merged_data, is_ogrdb_run = T)
  
  write.csv(clone_info_partis, paste0('../processed_data/clone_info_files/', mouse_id, '_clone_info_partis.csv'),
            row.names = F)
  write.csv(clone_info_partis_ogrdb, paste0('../processed_data/clone_info_files/', mouse_id, '_clone_info_partis_ogrdb.csv'),
            row.names = F)
  
  # Export clone info file based on the IgBLAST assignment by the Boyd lab.
  clone_info_igblast <- get_igblast_clone_info(merged_data)
  
  write.csv(clone_info_igblast, paste0('../processed_data/clone_info_files/', mouse_id, '_clone_info_igblast.csv'),
            row.names = F)
  
  # Calculate site-specific mutation frequencies for each V gene
  # Only available for the partis assignments 
  mut_probs_per_gene_position_partis <- estimate_mut_probs_per_vgene_position(annotated_seqs, is_ogrdb_run = F)
  
  write.csv(mut_probs_per_gene_position_partis,
            paste0('../results/mutations_per_vgene_base/', mouse_id, '_mutations_per_vgene_base_partis.csv'),
            row.names = F)
  
  mut_probs_per_gene_position_partis_ogrdb <- estimate_mut_probs_per_vgene_position(annotated_seqs, is_ogrdb_run = T)
  write.csv(mut_probs_per_gene_position_partis_ogrdb,
            paste0('../results/mutations_per_vgene_base/', mouse_id, '_mutations_per_vgene_base_partis_ogrdb.csv'))
}

yaml_object <- read_yaml(yaml_file_path)
yaml_object_ogrdb <- read_yaml(yaml_file_path_ogrdb)
merge_info(yaml_object, yaml_object_ogrdb, mouse_data_file_path)

#Write germline gene sequences (including new alleles)
mouse_id = str_extract(yaml_file_path,'[0-9]*-[0-9]*')

germline_genes_dir <- '../results/partis/partis_germline_genes/'

export_partis_germline_genes(yaml_object = yaml_object, mouse_id = mouse_id,
                             output_dir = germline_genes_dir, is_ogrdb_run = F)

export_partis_germline_genes(yaml_object = yaml_object_ogrdb, mouse_id = mouse_id,
                             output_dir = germline_genes_dir, is_ogrdb_run = T)
