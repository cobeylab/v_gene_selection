library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(Biostrings)
source('gene_frequency_functions.R')

annotated_seq_files <- list.files('../processed_data/annotated_seq_files/', pattern = 'csv', full.names = T)
# annotated_seq_files <- list.files('~/Desktop/annotated_seq_files/', pattern = 'csv', full.names = T)
clone_info_files <- list.files('../processed_data/clone_info_files/', pattern = 'csv', full.names = T)
mutations_per_vgene_base_files <- list.files('../results/mutations_per_vgene_base/', pattern = 'csv', full.names = T)
germline_v_genes_files <- list.files('../results/partis/partis_germline_genes/', pattern = 'v_genes', full.names = T)

print('Combining annotated sequences')
annotated_seqs <- lapply(annotated_seq_files, read_csv)
annotated_seqs <- bind_rows(annotated_seqs)
annotated_seqs$cell_type[annotated_seqs$cell_type == 'naive'] <- 'IgD+B220+'

# Process Igd+B220+ cells to decide which are naive, which are likely non-naive
annotated_seqs <- process_IgD_B220_seqs(annotated_seqs,
                                        max_clone_unique_IgDB220_seqs = 1,
                                        max_v_gene_mutations = 2)


write_csv(annotated_seqs, '../processed_data/annotated_seqs.csv')

print('Combining clone-level info')
clone_info <- lapply(clone_info_files, read_csv)
clone_info <- bind_rows(clone_info)
write_csv(clone_info, '../processed_data/clone_info.csv')

print('Combining per-base V gene mutations')
mutations_per_vgene_base <- lapply(mutations_per_vgene_base_files, read_csv)
mutations_per_vgene_base <- bind_rows(mutations_per_vgene_base)
write_csv(mutations_per_vgene_base, '../results/mutations_per_vgene_base.csv')

germline_v_genes <- lapply(germline_v_genes_files,
                           FUN = function(path){
                             mouse_id = str_remove(rev(str_split(path,'//')[[1]])[1], c('v_genes_'))
                             mouse_id = str_remove(mouse_id, '\\.fasta')
                             seq_object <- readDNAStringSet(path, format = 'fasta')

                             return(tibble(mouse_id = mouse_id, 
                                           v_gene = names(seq_object),
                                           v_gene_seq = as.character(seq_object)))
                             }
                           )

germline_v_genes <- bind_rows(germline_v_genes)

# Check that each v gene (including new alleles identified in different mice) is associated with a unique sequence
uniq_seqs_per_v_gene <- germline_v_genes %>% group_by(v_gene) %>% summarise(n_uniq_seqs = length(unique(v_gene_seq))) %>%
  ungroup() %>% select(n_uniq_seqs) %>% unique() %>% pull(n_uniq_seqs)
stopifnot(length(uniq_seqs_per_v_gene) == 1 & uniq_seqs_per_v_gene == 1)


germline_v_genes <- germline_v_genes %>% select(v_gene, v_gene_seq) %>% unique()
write_csv(germline_v_genes, '../results/germline_genes.csv')

# ====== Sequences counts for experienced cells, by mouse, cell type, tissue, V gene
seq_counts <- annotated_seqs %>%
  filter(productive_partis) %>%
  group_by(mouse_id, clone_id, tissue, cell_type) %>%
  summarise(prod_seqs = n()) %>%
  ungroup()
write_csv(seq_counts, '../processed_data/seq_counts.csv')

# Same structure, but after grouping sequences from the same mouse, clone, cell type, tissue and isotype that are identical
unique_seq_counts <- annotated_seqs %>%
  filter(productive_partis) %>%
  select(mouse_id, clone_id, partis_uniq_ref_seq, tissue, cell_type, isotype) %>%
  unique() %>%
  group_by(mouse_id, clone_id, tissue, cell_type) %>%
  summarise(unique_prod_seqs = n()) %>%
  ungroup()
write_csv(unique_seq_counts, '../processed_data/unique_seq_counts.csv')

# From the now-discontinued clustering analysis...
#unique_seq_cluster_files <-  list.files('../processed_data/unique_seq_clusters_files/', pattern = 'csv', full.names = T)



#unique_seq_clusters <- lapply(unique_seq_cluster_files, read_csv)
#unique_seq_clusters <- lapply(unique_seq_clusters, FUN = function(x){x %>% mutate(clone_id = as.character(clone_id),
#                                                                                  cluster_ref_seq = as.character(cluster_ref_seq),
#                                                                                  seq_id = as.character(seq_id),
#                                                                                  isotype = as.character(isotype))}) 

#unique_seq_clusters <- bind_rows(unique_seq_clusters) %>%
#  filter(!is.na(seq_id)) #  Handles incomplete lines for files still being filled

#write_csv(unique_seq_clusters, '../processed_data/unique_seq_clusters.csv')

# In addition to file with detailed sequence clustering, export file with unique sequence counts by mouse, clone, cell type, tissue

#unique_seq_counts <- unique_seq_clusters %>%
#  filter(!is.na(isotype)) %>% # Filter sequences with missing isotype
#  group_by(mouse_id, clone_id, tissue, cell_type, isotype) %>%
#  summarise(uniq_prod_seqs = length(unique(cluster_ref_seq))) %>%
  # Sum across isotypes
#  group_by(mouse_id, clone_id, tissue, cell_type) %>% 
#  summarise(uniq_prod_seqs = sum(uniq_prod_seqs)) %>%
#  ungroup()

#write_csv(unique_seq_counts, '../processed_data/unique_seq_counts.csv')

# File with total productive and unique productive sequences per mouse
#total_prod_seqs_per_mouse <- unique_seq_clusters %>%
#  group_by(mouse_id) %>%
#  summarise(total_prod_seqs = n()) %>%
#  ungroup()

#unique_prod_seqs_per_mouse <- unique_seq_counts %>%
#  group_by(mouse_id) %>%
#  summarise(uniq_prod_seqs = sum(uniq_prod_seqs))

#prod_seqs_per_mouse <- left_join(total_prod_seqs_per_mouse, unique_prod_seqs_per_mouse) 

#prod_seqs_per_mouse %>% summarise(across(c('total_prod_seqs', 'uniq_prod_seqs'), list(mean = mean, sum = sum)))

#prod_seqs_per_mouse %>%
#  write_csv('../processed_data/prod_seqs_per_mouse.csv')