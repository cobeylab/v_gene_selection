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

# OGRDB dataset, for renaming V genes in partis-ogrdb assignment back to their IMGT names
ogrdb_data <- read_csv('../data/OGRDB_C57BL6_genes.csv') %>%
  # Remove space from gene labels
  mutate(Label = str_remove(Label, '\\s')) %>%
  dplyr::rename(ogrdb_name = Label,
                imgt_name = `IMGT Name`)

# Function for renaming ogrdb genes following the traditional nomenclature
convert_ogrdb_to_imgt <- function(input_tibble, ogrdb_data){
  # Removing the '*x' suffix
  input_tibble <- input_tibble %>% mutate(v_gene = str_remove(v_gene,'\\*x'))
  
  # Check if all genes in input tibble are listed in the ogrdb data
  input_tibble_labels <- unique(input_tibble$v_gene)
  stopifnot(all(input_tibble_labels[!is.na(input_tibble_labels)] %in% ogrdb_data$ogrdb_name))
  
  relabelled_tibble <- left_join(input_tibble,
            ogrdb_data %>% dplyr::rename(v_gene = ogrdb_name) %>% select(v_gene, imgt_name)) %>%
    mutate(new_label = ifelse(is.na(imgt_name), v_gene, imgt_name)) %>%
    mutate(v_gene = new_label) %>%
    select(-new_label, imgt_name) %>%
    select(matches('mouse_id'), matches('clone_id'), matches('v_gene'), everything())
  
  stopifnot(nrow(relabelled_tibble) == nrow(input_tibble))
  
  return(relabelled_tibble)

}

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

# Final annotated seqs depend on filters for naive seqs, which are assignment specific
# Hence assignment-specific annotated_seqs objects

# For the partis assignments, standardize variable names
partis_var_names <- names(annotated_seqs)[str_detect(names(annotated_seqs),'partis')]

annotated_seqs_partis <- annotated_seqs %>%
  select(-all_of(partis_var_names[str_detect(partis_var_names,'ogrdb')])) %>%
  select(-matches('igblast')) %>%
  dplyr::rename(clone_id = clone_id_partis, productive = productive_partis)

annotated_seqs_partis_ogrdb <- annotated_seqs %>%
  select(-all_of(partis_var_names[str_detect(partis_var_names,'ogrdb') == F])) %>%
  select(-matches('igblast'))

names(annotated_seqs_partis_ogrdb) <- str_replace(names(annotated_seqs_partis_ogrdb), 'partis_ogrdb','partis')

annotated_seqs_partis_ogrdb <- annotated_seqs_partis_ogrdb %>%
  dplyr::rename(clone_id = clone_id_partis, productive = productive_partis)


# Process Igd+B220+ cells to decide which are naive, which are likely non-naive
max_clone_unique_IgDB220_seqs = 1
max_v_gene_mutations = 2

annotated_seqs_partis <- process_IgD_B220_seqs(annotated_seqs_partis,
                                        max_clone_unique_IgDB220_seqs = max_clone_unique_IgDB220_seqs,
                                        max_v_gene_mutations = max_v_gene_mutations)

annotated_seqs_partis_ogrdb <- process_IgD_B220_seqs(annotated_seqs_partis_ogrdb,
                                                     max_clone_unique_IgDB220_seqs = max_clone_unique_IgDB220_seqs,
                                                     max_v_gene_mutations = max_v_gene_mutations)

# For the igblast assignment, use the naive filters from the partis assignment
annotated_seqs_igblast <- annotated_seqs %>%
  select(-matches('partis')) %>%
  dplyr::rename(clone_id = clone_id_igblast, productive = productive_igblast) %>%
  select(-cell_type)

annotated_seqs_igblast <- left_join(annotated_seqs_igblast, 
                                    annotated_seqs_partis %>% select(mouse_id, seq_id, cell_type, partis_uniq_ref_seq),
                                    by = c('mouse_id','seq_id')) %>%
  select(mouse_id, seq_id, partis_uniq_ref_seq, tissue, cell_type, everything())

rm(annotated_seqs)

write_csv(annotated_seqs_partis, '../processed_data/annotated_seqs_partis.csv')
write_csv(annotated_seqs_partis_ogrdb, '../processed_data/annotated_seqs_partis_ogrdb.csv')
write_csv(annotated_seqs_igblast, '../processed_data/annotated_seqs_igblast.csv')


print('Combining clone-level info')
clone_info_partis <- bind_rows(lapply(partis_clone_info_files, read_csv))

clone_info_partis_ogrdb <- bind_rows(lapply(partis_ogrdb_clone_info_files, read_csv))
clone_info_partis_ogrdb <- convert_ogrdb_to_imgt(clone_info_partis_ogrdb, ogrdb_data)

clone_info_igblast <- bind_rows(lapply(igblast_clone_info_files, read_csv))

# Mutations per v gene base (only implemented for the partis assignments)
print('Combining per-base V gene mutations')
mutations_per_vgene_base_partis <- bind_rows(lapply(mutations_per_vgene_base_partis_files, read_csv))
write_csv(mutations_per_vgene_base_partis, '../results/mutations_per_vgene_base_partis.csv')

mutations_per_vgene_base_partis_ogrdb <- bind_rows(lapply(mutations_per_vgene_base_partis_ogrdb_files, read_csv))
mutations_per_vgene_base_partis_ogrdb <- convert_ogrdb_to_imgt(mutations_per_vgene_base_partis_ogrdb, ogrdb_data)

write_csv(mutations_per_vgene_base_partis_ogrdb, '../results/mutations_per_vgene_base_partis_ogrdb.csv')


# Export sequences of genes identified by partis

germline_v_genes_partis <- process_partis_germline_genes(germline_v_genes_partis_files)
write_csv(germline_v_genes_partis, '../results/germline_genes_partis.csv')

germline_v_genes_partis_ogrdb <- process_partis_germline_genes(germline_v_genes_partis_ogrdb_files)
germline_v_genes_partis_ogrdb <- convert_ogrdb_to_imgt(germline_v_genes_partis_ogrdb, ogrdb_data)

write_csv(germline_v_genes_partis_ogrdb, '../results/germline_genes_partis_ogrdb.csv')

# ====== Productive sequences counts for experienced cells, by mouse, cell type, tissue, V gene
seq_counts_partis <- get_seq_counts(annotated_seqs_partis, unique_only = F, productive_only = T)
# Ig blast assignment will use the partis assignment-based filters for naive sequences
seq_counts_igblast <- get_seq_counts(annotated_seqs_igblast, unique_only = F, productive_only = T)
seq_counts_partis_ogrdb <- get_seq_counts(annotated_seqs_partis_ogrdb, unique_only = F, productive_only = T)

write_csv(seq_counts_partis, '../processed_data/seq_counts_partis.csv')
write_csv(seq_counts_igblast, '../processed_data/seq_counts_igblast.csv')
write_csv(seq_counts_partis_ogrdb, '../processed_data/seq_counts_partis_ogrdb.csv')

# Same structure, but after grouping sequences from the same mouse, clone, cell type, tissue and isotype that are identical

unique_seq_counts_partis <- get_seq_counts(annotated_seqs_partis, unique_only = T, productive_only = T)
unique_seq_counts_igblast <- get_seq_counts(annotated_seqs_igblast, unique_only = T, productive_only = T)
unique_seq_counts_partis_ogrdb <- get_seq_counts(annotated_seqs_partis_ogrdb, unique_only = T, productive_only = T)

write_csv(unique_seq_counts_partis, '../processed_data/unique_seq_counts_partis.csv')
write_csv(unique_seq_counts_igblast, '../processed_data/unique_seq_counts_igblast.csv')
write_csv(unique_seq_counts_partis_ogrdb, '../processed_data/unique_seq_counts_partis_ogrdb.csv')

# ========= Annotate clone_info with clone tissue composition and clone cell type composition =======
#                                     (based on productive sequences)

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
 
# ========== Export sequence counts including unproductive sequences
seq_counts_partis_incl_unprod <- get_seq_counts(annotated_seqs_partis, unique_only = F, productive_only = F)
seq_counts_partis_ogrdb_incl_unprod <- get_seq_counts(annotated_seqs_partis_ogrdb, unique_only = F, productive_only = F)
seq_counts_igblast_incl_unprod <- get_seq_counts(annotated_seqs_igblast, unique_only = F, productive_only = F)

# As a test, for clones with productive sequences, the n of seqs when unproductive sequences are included has to be equal to or greater
# than the number of prod seqs:
seq_count_test <- left_join(seq_counts_partis, seq_counts_partis_incl_unprod) %>%
  mutate(test = seqs >= prod_seqs) %>%
  select(test) %>% unique() %>% pull(test)
stopifnot(seq_count_test)

write_csv(seq_counts_partis_incl_unprod, '../processed_data/seq_counts_partis_incl_unprod.csv')
write_csv(seq_counts_partis_ogrdb_incl_unprod, '../processed_data/seq_counts_partis_ogrdb_incl_unprod.csv')
write_csv(seq_counts_igblast_incl_unprod, '../processed_data/seq_counts_igblast_incl_unprod.csv')




