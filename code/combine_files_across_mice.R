library(purrr)
library(readr)
library(Biostrings)
library(seqinr)
source('gene_frequency_functions.R')

annotated_seq_files <- list.files('../processed_data/annotated_seq_files/', pattern = 'csv', full.names = T)
partis_clone_info_files <- list.files('../processed_data/clone_info_files/', pattern = 'partis\\.csv', full.names = T)
igblast_clone_info_files <-  list.files('../processed_data/clone_info_files/', pattern = 'igblast\\.csv', full.names = T)
mutations_per_vgene_base_files <- list.files('../results/mutations_per_vgene_base/', pattern = 'csv', full.names = T)
germline_v_genes_files <- list.files('../results/partis/partis_germline_genes/', pattern = 'v_genes', full.names = T)

print('Combining annotated sequences')
annotated_seqs <- lapply(annotated_seq_files, read_csv)
annotated_seqs <- bind_rows(annotated_seqs)
annotated_seqs$cell_type[annotated_seqs$cell_type == 'naive'] <- 'IgD+B220+'

# Process Igd+B220+ cells to decide which are naive, which are likely non-naive
# (All criteria based on the Partis assignment)
annotated_seqs <- process_IgD_B220_seqs(annotated_seqs,
                                        max_clone_unique_IgDB220_seqs = 1,
                                        max_v_gene_mutations = 2,
                                        filters_from = 'partis')

write_csv(annotated_seqs, '../processed_data/annotated_seqs.csv')

print('Combining clone-level info')
clone_info_partis <- bind_rows(lapply(partis_clone_info_files, read_csv))
clone_info_igblast <- bind_rows(lapply(igblast_clone_info_files, read_csv))


# Mutations per v gene base (only implemented for the partis assignment)
print('Combining per-base V gene mutations')
mutations_per_vgene_base <- lapply(mutations_per_vgene_base_files, read_csv)
mutations_per_vgene_base <- bind_rows(mutations_per_vgene_base)
write_csv(mutations_per_vgene_base, '../results/mutations_per_vgene_base.csv')

# Export sequences of genes identified by partis
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
write_csv(germline_v_genes, '../results/germline_genes_partis.csv')

# ====== Sequences counts for experienced cells, by mouse, cell type, tissue, V gene
seq_counts_partis <- get_productive_seq_counts(annotated_seqs, unique_only = F, assignment = 'partis')
seq_counts_igblast <- get_productive_seq_counts(annotated_seqs, unique_only = F, assignment = 'igblast')

write_csv(seq_counts_partis, '../processed_data/seq_counts_partis.csv')
write_csv(seq_counts_igblast, '../processed_data/seq_counts_igblast.csv')


# Same structure, but after grouping sequences from the same mouse, clone, cell type, tissue and isotype that are identical

unique_seq_counts_partis <- get_productive_seq_counts(annotated_seqs, unique_only = T, assignment = 'partis')
unique_seq_counts_igblast <- get_productive_seq_counts(annotated_seqs, unique_only = T, assignment = 'igblast')

write_csv(unique_seq_counts_partis, '../processed_data/unique_seq_counts_partis.csv')
write_csv(unique_seq_counts_igblast, '../processed_data/unique_seq_counts_igblast.csv')

# ========= Annotate clone_info with clone tissue composition and clone cell type composition =======

# Numbers of productive sequences (as opposed to only unique seqs) in each tissue and cell type in each clone
clone_info_partis <- annotate_clone_info_with_composition(clone_info = clone_info_partis,
                                                          seq_counts = seq_counts_partis,
                                                          unique_seq_counts = unique_seq_counts_partis)

clone_info_igblast <- annotate_clone_info_with_composition(clone_info = clone_info_igblast,
                                                          seq_counts = seq_counts_igblast,
                                                          unique_seq_counts = unique_seq_counts_igblast)
    

write_csv(clone_info_partis, '../processed_data/clone_info_partis.csv')
write_csv(clone_info_igblast, '../processed_data/clone_info_igblast.csv')

# Export naive sequences of each clone to a fasta file
#write.fasta(names = paste(clone_info_partis$mouse_id, clone_info_partis$clone_id, sep = '_'),
#            sequences = as.list(clone_info_partis$clone_naive_seq_nt_partis),
#            file.out = '../results/partis_clone_naive_seqs.fasta')

