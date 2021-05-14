# Using the unique_seq_counts file, annotates the clone_info file with whether clone is pure (all naive or all non-naive)
library(readr)
library(dplyr)
library(tidyr)
source('gene_frequency_functions.R')

seq_counts <- read_csv('../processed_data/seq_counts_unclustered.csv')

clone_info <- read_csv('../processed_data/clone_info.csv')
clone_purity <- get_clone_purity(seq_counts) %>%
  dplyr::rename(naive_prod_seqs = naive_seqs_in_clone,
                nonnaive_prod_seqs = non_naive_seqs_in_clone) %>%
  mutate(prod_seqs = naive_prod_seqs + nonnaive_prod_seqs) %>%
  mutate(fraction_naive = naive_prod_seqs / prod_seqs)
  

clone_info <- left_join(clone_info, clone_purity) %>%
  select(mouse_id, clone_id, prod_seqs, clone_purity, fraction_naive, everything())
write_csv(clone_info, '../processed_data/clone_info.csv')

# Fraction of clones that are pure naive, pure non-naive, or mixed
clone_purity %>%
  group_by(clone_purity) %>%
  count() %>%
  ungroup() %>%
  mutate(fraction = n/sum(n)) %>%
  dplyr::rename(n_clones_with_prod_seqs = n) %>%
  write_csv('../processed_data/clone_purity.csv')

# Fraction of clones WITH NAIVE SEQs that are pure naive
clone_purity %>%
  filter(naive_prod_seqs > 0) %>%
  group_by(clone_purity) %>%
  count() %>%
  ungroup() %>%
  mutate(fraction = n / sum(n))

# Fraction of clones with NON-NAIVE SEQS that are pure-non-naive
clone_purity %>%
  filter(nonnaive_prod_seqs > 0) %>%
  group_by(clone_purity) %>%
  count() %>%
  ungroup() %>%
  mutate(fraction = n / sum(n))


