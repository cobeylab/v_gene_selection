# For each clone, calculates summary statistics of mutation across sequences
library(dplyr)
library(readr)
library(tidyr)
library(Biostrings)

# Seq information (isotype, n mutations in CDRs, etc.)
annotated_seqs <- read_csv('../processed_data/annotated_seqs.csv')

annotated_seqs <- annotated_seqs %>%
  mutate(across(c('clone_id_partis','partis_uniq_ref_seq','seq_id'), as.character))

annotated_seqs$specimen_cell_subset[annotated_seqs$specimen_cell_subset == 'naïve'] <- 'naive'

# Sequence clustering information to estimate unique sequences
unique_seq_clusters <- read_csv('../processed_data/unique_seq_clusters.csv')
unique_seq_clusters <- unique_seq_clusters %>%
  mutate(across(c('clone_id','cluster_ref_seq','seq_id'), as.character))


combined_tibble <- left_join(unique_seq_clusters,
          annotated_seqs %>% dplyr::rename(clone_id = clone_id_partis,
                                           cell_type = specimen_cell_subset,
                                           tissue = specimen_tissue)) 



clones_summary <- combined_tibble %>%
  group_by(mouse_id, clone_id, cluster_ref_seq) %>%
  # First, average within each cluster of closely related seqs (so that each cluster is counted once for clone stats)
  summarise(n_mutations_partis_aa = mean(n_mutations_partis_aa),
            cdr3_mutations_partis_aa = mean(cdr3_mutations_partis_aa),
            productive_consensus_CDR3 = consensusString(AAStringSet(cdr3_seq_partis))) %>%
  ungroup()

# Then average for each clone.
final_summary <- clones_summary %>%
  group_by(mouse_id, clone_id) %>%
  summarise(mean_n_mutations_partis_aa = mean(n_mutations_partis_aa),
            mean_cdr3_mutations_partis_aa = mean(cdr3_mutations_partis_aa),
            median_n_mutations_partis_aa = median(n_mutations_partis_aa),
            median_cdr3_mutations_partis_aa = median(cdr3_mutations_partis_aa),
            max_n_mutations_partis_aa = max(n_mutations_partis_aa),
            max_cdr3_mutations_partis_aa = max(cdr3_mutations_partis_aa),
            min_n_mutations_partis_aa = min(n_mutations_partis_aa),
            min_cdr3_mutations_partis_aa = min(cdr3_mutations_partis_aa),
            productive_consensus_CDR3 = consensusString(AAStringSet(productive_consensus_CDR3))) %>%
  ungroup
  
write_csv(final_summary, '../results/clone_summary_statistics.csv')
