library(dplyr)
library(purrr)
library(tidyr)
library(readr)
library(stringr)
library(Biostrings)
library(seqinr)
source('gene_frequency_functions.R')

annotated_seq_files <- list.files('../processed_data/annotated_seq_files/', pattern = 'csv', full.names = T)
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
seq_counts <- get_productive_seq_counts(annotated_seqs, unique_only = F)
write_csv(seq_counts, '../processed_data/seq_counts.csv')


# Same structure, but after grouping sequences from the same mouse, clone, cell type, tissue and isotype that are identical
unique_seq_counts <- get_productive_seq_counts(annotated_seqs, unique_only = T)
write_csv(unique_seq_counts, '../processed_data/unique_seq_counts.csv')

# ========= Annotate clone_info with clone tissue composition and clone cell type composition =======
# Numbers of productive sequences (as opposed to only unique seqs) in each tissue and cell type in each clone

clone_tissue_composition_prod_seqs <- get_clone_composition(seq_counts, composition_var = 'tissue')
clone_cell_type_composition_prod_seqs <- get_clone_composition(seq_counts, composition_var = 'cell_type')
  
# Total number of productive sequences should be equal for each clone in both tibbles
stopifnot(all(clone_tissue_composition_prod_seqs$total_clone_prod_seqs == clone_tissue_composition_prod_seqs$total_clone_prod_seqs))

clone_info <- left_join(clone_info, clone_tissue_composition_prod_seqs)
clone_info <- left_join(clone_info, clone_cell_type_composition_prod_seqs %>% select(-total_clone_prod_seqs))


# Numbers of *unique* productive sequences 

clone_tissue_composition_unique_seqs <- get_clone_composition(unique_seq_counts, composition_var = 'tissue')
clone_cell_type_composition_unique_seqs <- get_clone_composition(unique_seq_counts, composition_var = 'cell_type')
  
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

# Export naive sequences of each clone to a fasta file
write.fasta(names = paste(clone_info$mouse_id, clone_info$clone_id, sep = '_'),
            sequences = as.list(clone_info$clone_naive_seq_nt_partis),
            file.out = '../results/clone_naive_seqs.fasta')

