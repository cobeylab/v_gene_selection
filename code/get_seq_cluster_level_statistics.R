library(dplyr)
library(tidyr)
library(readr)
library(Biostrings)

source('gene_frequency_functions.R')

# Seq information (isotype, n mutations in CDRs, etc.)
annotated_seqs <- read_csv('../processed_data/annotated_seqs.csv')

annotated_seqs <- annotated_seqs %>%
  mutate(across(c('clone_id_partis','partis_uniq_ref_seq','seq_id'), as.character))

annotated_seqs$specimen_cell_subset[annotated_seqs$specimen_cell_subset == 'naïve'] <- 'naive'

# Sequence clustering information to estimate unique sequences
unique_seq_clusters <- read_csv('../processed_data/unique_seq_clusters.csv')
unique_seq_clusters <- unique_seq_clusters %>%
  mutate(across(c('clone_id','cluster_ref_seq','seq_id'), as.character)) %>%
  filter(!is.na(isotype))

combined_tibble <- left_join(unique_seq_clusters,
                             annotated_seqs %>% dplyr::rename(clone_id = clone_id_partis,
                                                              cell_type = specimen_cell_subset,
                                                              tissue = specimen_tissue)) 

# Average the number of mutations within each cluster of closely related seqs (so each cluster is counted once)
combined_tibble <- combined_tibble %>%
  group_by(mouse_id, clone_id, tissue, cell_type, isotype, cluster_ref_seq) %>%
  summarise(n_mutations_partis_aa = mean(n_mutations_partis_aa),
            n_mutations_partis_nt = mean(n_mutations_partis_nt),
            seq_length_partis = mean(seq_length_partis),
            cdr3_mutations_partis_aa = mean(cdr3_mutations_partis_aa),
            cdr3_mutations_partis_nt = mean(cdr3_mutations_partis_nt),
            vgene_mutations_partis_nt = mean(vgene_mutations_partis_nt),
            sequenced_bases_in_vgene_region_partis = mean(sequenced_bases_in_vgene_region_partis),
            cluster_consensus_CDR3 = consensusString(AAStringSet(cdr3_seq_partis))) %>%
  ungroup()

write_csv(combined_tibble, '../processed_data/seq_cluster_stats.csv')
