# For each clone, calculates summary statistics of mutation across sequences
library(dplyr)
library(readr)
library(tidyr)
library(DECIPHER)

annotated_seqs <- read_csv('../processed_data/annotated_seqs.csv')

clones_summary <- annotated_seqs %>%
  dplyr::rename(clone_id = clone_id_partis) %>%
  filter(productive_partis) %>%
  group_by(mouse_id, clone_id) %>%
  summarise(mean_n_nt_mutations = mean(n_mutations_partis_nt),
            mean_n_aa_mutations = mean(n_mutations_partis_aa),
            mean_n_nt_mutations_CDR3 = mean(cdr3_mutations_partis_nt),
            mean_n_aa_mutations_CDR3 = mean(cdr3_mutations_partis_aa),
            max_n_nt_mutations = max(n_mutations_partis_nt),
            max_n_aa_mutations = max(n_mutations_partis_aa),
            max_n_nt_mutations_CDR3 = max(cdr3_mutations_partis_nt),
            max_n_aa_mutations_CDR3 = max(cdr3_mutations_partis_aa),
            productive_consensus_CDR3 = consensusString(AAStringSet(cdr3_seq_partis))
            )

write_csv(clones_summary, '../results/clone_summary_statistics.csv')