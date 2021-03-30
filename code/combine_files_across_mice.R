library(dplyr)
library(tidyr)
library(readr)

annotated_seq_files <- list.files('../processed_data/annotated_seq_files/', pattern = 'csv', full.names = T)
clone_info_files <- list.files('../processed_data/clone_info_files/', pattern = 'csv', full.names = T)
unique_seq_cluster_files <-  list.files('../processed_data/unique_seq_clusters_files/', pattern = 'csv', full.names = T)


annotated_seqs <- lapply(annotated_seq_files, read_csv)
annotated_seqs <- bind_rows(annotated_seqs)
write_csv(annotated_seqs, '../processed_data/annotated_seqs.csv')

clone_info <- lapply(clone_info_files, read_csv)
clone_info <- bind_rows(clone_info)
write_csv(clone_info, '../processed_data/clone_info.csv')


unique_seq_clusters <- lapply(unique_seq_cluster_files, read_csv)

unique_seq_clusters <- lapply(unique_seq_clusters, FUN = function(x){x %>% mutate(clone_id = as.character(clone_id),
                                                                                  cluster_ref_seq = as.character(cluster_ref_seq),
                                                                                  seq_id = as.character(seq_id),
                                                                                  isotype = as.character(isotype))}) 

unique_seq_clusters <- bind_rows(unique_seq_clusters) %>%
  filter(!is.na(seq_id)) #  Handles incomplete lines for files still being filled

write_csv(unique_seq_clusters, '../processed_data/unique_seq_clusters.csv')

# In addition to file with detailed sequence clustering, export file with unique sequence counts by mouse, clone, cell type, tissue

unique_seq_counts <- unique_seq_clusters %>%
  filter(!is.na(isotype)) %>% # Filter sequences with missing isotype
  group_by(mouse_id, clone_id, tissue, cell_type, isotype) %>%
  summarise(uniq_prod_seqs = length(unique(cluster_ref_seq))) %>%
  # Sum across isotypes
  group_by(mouse_id, clone_id, tissue, cell_type) %>% 
  summarise(uniq_prod_seqs = sum(uniq_prod_seqs)) %>%
  ungroup()

write_csv(unique_seq_counts, '../processed_data/unique_seq_counts.csv')


# File with total productive and unique productive sequences per mouse
total_prod_seqs_per_mouse <- unique_seq_clusters %>%
  group_by(mouse_id) %>%
  summarise(total_prod_seqs = n()) %>%
  ungroup()

unique_prod_seqs_per_mouse <- unique_seq_counts %>%
  group_by(mouse_id) %>%
  summarise(uniq_prod_seqs = sum(uniq_prod_seqs))

prod_seqs_per_mouse <- left_join(total_prod_seqs_per_mouse, unique_prod_seqs_per_mouse) 

prod_seqs_per_mouse %>% summarise(across(c('total_prod_seqs', 'uniq_prod_seqs'), list(mean = mean, sum = sum)))

prod_seqs_per_mouse %>%
  write_csv('../processed_data/prod_seqs_per_mouse.csv')
  


