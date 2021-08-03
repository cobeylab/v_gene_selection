library(dplyr)
library(purrr)
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

# ========= Annotate clone_info with clone tissue composition and clone cell type composition =======
# Numbers of productive sequences (as opposed to only unique seqs) in each tissue and cell type in each clone

clone_tissue_composition_prod_seqs <- seq_counts %>%
  group_by(mouse_id, clone_id, tissue) %>%
  # For each clone, sum across cell types within each tissue
  dplyr::summarise(prod_seqs_in_tissue = sum(prod_seqs)) %>%
  group_by(mouse_id, clone_id) %>%
  # Now compute the fraction of sequences in a clone that came from each tissue
  mutate(total_clone_prod_seqs = sum(prod_seqs_in_tissue)) %>%
  ungroup() %>%
  pivot_wider(id_cols = any_of(c('mouse_id','clone_id','total_clone_prod_seqs')),
              names_from = tissue, values_from = prod_seqs_in_tissue,
              values_fill = 0, names_prefix = 'prod_seqs_') 

clone_cell_type_composition_prod_seqs <- seq_counts %>%
  group_by(mouse_id, clone_id, cell_type) %>%
  # For each clone, sum across tissue within each cell type
  dplyr::summarise(prod_seqs_in_cell_type = sum(prod_seqs)) %>%
  group_by(mouse_id, clone_id) %>%
  # Now compute the fraction of sequences in a clone that came from each cell_type
  mutate(total_clone_prod_seqs = sum(prod_seqs_in_cell_type)) %>%
  ungroup() %>%
  pivot_wider(id_cols = any_of(c('mouse_id','clone_id','total_clone_prod_seqs')),
              names_from = cell_type, values_from = prod_seqs_in_cell_type,
              values_fill = 0, names_prefix = 'prod_seqs_') 

# Total number of productive sequences should be equal for each clone in both tibbles
stopifnot(all(clone_tissue_composition_prod_seqs$total_clone_prod_seqs == clone_tissue_composition_prod_seqs$total_clone_prod_seqs))


clone_info <- left_join(clone_info, clone_tissue_composition_prod_seqs)
clone_info <- left_join(clone_info, clone_cell_type_composition_prod_seqs %>% select(-total_clone_prod_seqs))


# Numbers of *unique* productive sequences 
# (identical sequences from different tissues, cell types or isotypes are counted as separate unique sequences)
clone_tissue_composition_unique_seqs <- unique_seq_counts %>%
  group_by(mouse_id, clone_id, tissue) %>%
  dplyr::summarise(unique_seqs_in_tissue = sum(unique_prod_seqs)) %>%
  group_by(mouse_id, clone_id) %>%
  mutate(total_clone_unique_seqs = sum(unique_seqs_in_tissue)) %>%
  ungroup() %>%
  pivot_wider(id_cols = any_of(c('mouse_id','clone_id','total_clone_unique_seqs')),
              names_from = tissue, values_from = unique_seqs_in_tissue,
              values_fill = 0, names_prefix = 'unique_seqs_') 

clone_cell_type_composition_unique_seqs <- unique_seq_counts %>%
  group_by(mouse_id, clone_id, cell_type) %>%
  dplyr::summarise(unique_seqs_in_cell_type = sum(unique_prod_seqs)) %>%
  group_by(mouse_id, clone_id) %>%
  # Now compute the fraction of sequences in a clone that came from each cell_type
  mutate(total_clone_unique_seqs = sum(unique_seqs_in_cell_type)) %>%
  ungroup() %>%
  pivot_wider(id_cols = any_of(c('mouse_id','clone_id','total_clone_unique_seqs')),
              names_from = cell_type, values_from = unique_seqs_in_cell_type,
              values_fill = 0, names_prefix = 'unique_seqs_') 

stopifnot(all(clone_cell_type_composition_unique_seqs$total_clone_unique_seqs ==
                clone_tissue_composition_unique_seqs$total_clone_unique_seqs))

clone_info <- left_join(clone_info, clone_tissue_composition_unique_seqs)
clone_info <- left_join(clone_info, clone_cell_type_composition_unique_seqs %>% select(-total_clone_unique_seqs))


# Compute Simpson diversity of tissues and cell types within each clone
# Both with frequencies relative to all productive sequences and with freqs relative to unique seqs only.

clone_info <- clone_info %>%
  # Compute fractions in each compartment (tissue or cell type)
  mutate(across(matches('prod_seqs_'), function(x){x/total_clone_prod_seqs},
                .names = "fraction_{.col}")) %>%
  mutate(across(matches('unique_seqs_'), function(x){x/total_clone_unique_seqs},
                .names = "fraction_{.col}")) %>%
  mutate(across(matches('fraction_'), function(x){x^2},
                .names = 'squared_{.col}'))

tissues <- unique(seq_counts$tissue)
cell_types <- unique(seq_counts$cell_type)

clone_info <- clone_info %>%
  mutate(
    simpson_diversity_tissues_prod_seqs = 1 - purrr::reduce(
      clone_info %>% select(any_of(paste0('squared_fraction_prod_seqs_', tissues))), `+`),
    simpson_diversity_tissues_unique_seqs = 1 - purrr::reduce(
      clone_info %>% select(any_of(paste0('squared_fraction_unique_seqs_', tissues))), `+`),
    simpson_diversity_cell_types_prod_seqs = 1 - purrr::reduce(
      clone_info %>% select(any_of(paste0('squared_fraction_prod_seqs_', cell_types))), `+`),
    simpson_diversity_cell_types_unique_seqs = 1 - purrr::reduce(
      clone_info %>% select(any_of(paste0('squared_fraction_unique_seqs_', cell_types))), `+`)
    ) %>%
  select(-matches('squared_fraction'))

clone_info <- clone_info %>%
  rowwise() %>%
  mutate(biggest_tissue_fraction_prod_seqs = max(c(fraction_prod_seqs_BM, fraction_prod_seqs_LN, fraction_prod_seqs_spleen)),
         biggest_tissue_fraction_unique_seqs =  max(c(fraction_unique_seqs_BM, fraction_unique_seqs_LN, fraction_unique_seqs_spleen)),
         biggest_cell_type_fraction_prod_seqs = max(c(fraction_prod_seqs_naive, `fraction_prod_seqs_nonnaive_IgD+B220+`,
                                                      fraction_prod_seqs_GC, fraction_prod_seqs_PC, fraction_prod_seqs_mem)),
         biggest_cell_type_fraction_unique_seqs = max(c(fraction_unique_seqs_naive, `fraction_unique_seqs_nonnaive_IgD+B220+`,
                                                        fraction_unique_seqs_GC, fraction_unique_seqs_PC, fraction_unique_seqs_mem))) %>%
  ungroup()
    

write_csv(clone_info, '../processed_data/clone_info.csv')





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