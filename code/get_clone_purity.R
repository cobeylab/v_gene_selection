# Using the unique_seq_counts file, annotates the clone_info file with whether clone is pure (all naive or all non-naive)
library(readr)
library(dplyr)
library(tidyr)
source('gene_frequency_functions.R')

clone_info <- read_csv('../processed_data/clone_info.csv')
unique_seq_counts <- read_csv('../processed_data/unique_seq_counts.csv')

clone_purity <- get_clone_purity(unique_seq_counts) %>%
  dplyr::rename(unique_naive_prod_seqs = naive_seqs_in_clone,
                unique_nonnaive_prod_seqs = non_naive_seqs_in_clone,
                clone_id_partis = clone_id) %>%
  mutate(unique_prod_seqs = unique_naive_prod_seqs + unique_nonnaive_prod_seqs)


clone_info <- left_join(clone_info %>% 
            select(-matches(c('unique_naive_prod_seqs','unique_nonnaive_prod_seqs',
                              'unique_prod_seqs','clone_purity'))),
          clone_purity)

write_csv(clone_info, '../processed_data/clone_info.csv')

clone_purity %>%
  group_by(clone_purity) %>%
  count() %>%
  ungroup() %>%
  mutate(fraction = n/sum(n)) %>%
  dplyr::rename(n_clones_with_prod_seqs = n) %>%
  write_csv('../processed_data/clone_purity.csv')

