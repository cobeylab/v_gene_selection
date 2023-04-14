library(purrr)
library(readr)
library(Biostrings)
library(seqinr)
source('gene_frequency_functions.R')

annotated_seq_files <- list.files('../processed_data/annotated_seq_files/', pattern = 'csv', full.names = T)

partis_clone_info_files <- list.files('../processed_data/clone_info_files/', pattern = 'partis\\.csv', full.names = T)
partis_ogrdb_clone_info_files <- list.files('../processed_data/clone_info_files/', pattern = 'partis_ogrdb\\.csv', full.names = T)
igblast_clone_info_files <-  list.files('../processed_data/clone_info_files/', pattern = 'igblast\\.csv', full.names = T)

mutations_per_vgene_base_partis_files <- list.files('../results/mutations_per_vgene_base/', pattern = 'partis\\.csv', full.names = T)
mutations_per_vgene_base_partis_ogrdb_files <- list.files('../results/mutations_per_vgene_base/', pattern = 'partis_ogrdb\\.csv', full.names = T)
germline_v_genes_partis_files <- list.files('../results/partis/partis_germline_genes/', pattern = 'v_genes_partis_[0-9]', full.names = T)
germline_v_genes_partis_ogrdb_files <- list.files('../results/partis/partis_germline_genes/', pattern = 'v_genes_partis_ogrdb', full.names = T)

# =========================================== Some auxillary functions =================================
read_partis_germline_genes <- function(path){
  mouse_id = str_extract(path,'[0-9]+-[0-9]+')
  seq_object <- readDNAStringSet(path, format = 'fasta')
  
  return(tibble(mouse_id = mouse_id, 
                v_gene = names(seq_object),
                v_gene_seq = as.character(seq_object)))
}

process_partis_germline_genes <- function(germline_v_genes_partis_files){
  germline_v_genes_partis <- bind_rows(lapply(germline_v_genes_partis_files, FUN = read_partis_germline_genes))
  
  # Check that each v gene (including new alleles identified in different mice) is associated with a unique sequence
  uniq_seqs_per_v_gene <- germline_v_genes_partis %>% group_by(v_gene) %>% summarise(n_uniq_seqs = length(unique(v_gene_seq))) %>%
    ungroup() %>% select(n_uniq_seqs) %>% unique() %>% pull(n_uniq_seqs)
  stopifnot(length(uniq_seqs_per_v_gene) == 1 & uniq_seqs_per_v_gene == 1)
  
  germline_v_genes_partis <- germline_v_genes_partis %>% select(v_gene, v_gene_seq) %>% unique()
  
  return(germline_v_genes_partis)
}

# ================================================ Processing ======================================

print('Combining annotated sequences')
annotated_seqs <- lapply(annotated_seq_files, read_csv)
annotated_seqs <- bind_rows(annotated_seqs)
annotated_seqs$cell_type[annotated_seqs$cell_type == 'naive'] <- 'IgD+B220+'

# Process Igd+B220+ cells to decide which are naive, which are likely non-naive
# (Because some criteria depend on clone size, this call is assignment-specific)


annotated_seqs_partis <- process_IgD_B220_seqs(annotated_seqs,
                                        max_clone_unique_IgDB220_seqs = 1,
                                        max_v_gene_mutations = 2,
                                        filters_from = 'partis')

write_csv(annotated_seqs_partis, '../processed_data/annotated_seqs_partis.csv')


annotated_seqs_partis_ogrdb <- process_IgD_B220_seqs(annotated_seqs,
                                                     max_clone_unique_IgDB220_seqs = 1,
                                                     max_v_gene_mutations = 2,
                                                     filters_from = 'partis_ogrdb')
write_csv(annotated_seqs_partis_ogrdb, '../processed_data/annotated_seqs_partis_ogrdb.csv')

rm(annotated_seqs)

print('Combining clone-level info')
clone_info_partis <- bind_rows(lapply(partis_clone_info_files, read_csv))
clone_info_partis_ogrdb <- bind_rows(lapply(partis_ogrdb_clone_info_files, read_csv))
clone_info_igblast <- bind_rows(lapply(igblast_clone_info_files, read_csv))

# Mutations per v gene base (only implemented for the partis assignments)
print('Combining per-base V gene mutations')
mutations_per_vgene_base_partis <- bind_rows(lapply(mutations_per_vgene_base_partis_files, read_csv))
write_csv(mutations_per_vgene_base_partis, '../results/mutations_per_vgene_base_partis.csv')

mutations_per_vgene_base_partis_ogrdb <- bind_rows(lapply(mutations_per_vgene_base_partis_ogrdb_files, read_csv))
write_csv(mutations_per_vgene_base_partis_ogrdb, '../results/mutations_per_vgene_base_partis_ogrdb.csv')


# Export sequences of genes identified by partis

germline_v_genes_partis <- process_partis_germline_genes(germline_v_genes_partis_files)
write_csv(germline_v_genes_partis, '../results/germline_genes_partis.csv')

germline_v_genes_partis_ogrdb <- process_partis_germline_genes(germline_v_genes_partis_ogrdb_files)
write_csv(germline_v_genes_partis_ogrdb, '../results/germline_genes_partis_ogrdb.csv')

# ====== Sequences counts for experienced cells, by mouse, cell type, tissue, V gene
seq_counts_partis <- get_productive_seq_counts(annotated_seqs_partis, unique_only = F, assignment = 'partis')
# Ig blast assignment will use the partis assignment-based filters for naive sequences
seq_counts_igblast <- get_productive_seq_counts(annotated_seqs_partis, unique_only = F, assignment = 'igblast')
seq_counts_partis_ogrdb <- get_productive_seq_counts(annotated_seqs_partis_ogrdb, unique_only = F, assignment = 'partis_ogrdb')

write_csv(seq_counts_partis, '../processed_data/seq_counts_partis.csv')
write_csv(seq_counts_igblast, '../processed_data/seq_counts_igblast.csv')
write_csv(seq_counts_partis_ogrdb, '../processed_data/seq_counts_partis_ogrdb.csv')

# Same structure, but after grouping sequences from the same mouse, clone, cell type, tissue and isotype that are identical

unique_seq_counts_partis <- get_productive_seq_counts(annotated_seqs_partis, unique_only = T, assignment = 'partis')
unique_seq_counts_igblast <- get_productive_seq_counts(annotated_seqs_partis, unique_only = T, assignment = 'igblast')
unique_seq_counts_partis_ogrdb <- get_productive_seq_counts(annotated_seqs_partis_ogrdb, unique_only = T, assignment = 'partis_ogrdb')

write_csv(unique_seq_counts_partis, '../processed_data/unique_seq_counts_partis.csv')
write_csv(unique_seq_counts_igblast, '../processed_data/unique_seq_counts_igblast.csv')
write_csv(unique_seq_counts_partis_ogrdb, '../processed_data/unique_seq_counts_partis_ogrdb.csv')

# ========= Annotate clone_info with clone tissue composition and clone cell type composition =======

# Numbers of productive sequences (as opposed to only unique seqs) in each tissue and cell type in each clone
clone_info_partis <- annotate_clone_info_with_composition(clone_info = clone_info_partis,
                                                          seq_counts = seq_counts_partis,
                                                          unique_seq_counts = unique_seq_counts_partis)

clone_info_igblast <- annotate_clone_info_with_composition(clone_info = clone_info_igblast,
                                                          seq_counts = seq_counts_igblast,
                                                          unique_seq_counts = unique_seq_counts_igblast)

clone_info_partis_ogrdb <- annotate_clone_info_with_composition(clone_info = clone_info_partis_ogrdb,
                                                          seq_counts = seq_counts_partis_ogrdb,
                                                          unique_seq_counts = unique_seq_counts_partis_ogrdb)

write_csv(clone_info_partis, '../processed_data/clone_info_partis.csv')
write_csv(clone_info_igblast, '../processed_data/clone_info_igblast.csv')
write_csv(clone_info_partis_ogrdb, '../processed_data/clone_info_partis_ogrdb.csv')
 
