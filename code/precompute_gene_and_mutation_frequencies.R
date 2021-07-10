library(readr)
library(dplyr)
library(tidyr)

source('gene_frequency_functions.R')

# Annotated sequences
annotated_seqs <- read_csv('../processed_data/annotated_seqs.csv')
# annotated_seqs <- read_csv('~/Desktop/v_gene_selection_files/annotated_seqs.csv')

# Basic info for each clone (germline genes, CDR lenght, naive CDR seq)
clone_info <- read_csv('../processed_data/clone_info.csv')
# clone_info <- read_csv('~/Desktop/v_gene_selection_files/clone_info.csv')

seq_counts <- read_csv('../processed_data/seq_counts.csv')
# seq_counts <- read_csv('~/Desktop/v_gene_selection_files/seq_counts.csv')

seq_counts <- get_info_from_mouse_id(seq_counts)

seq_counts <- left_join(seq_counts, clone_info %>% select(mouse_id, clone_id, v_gene)) %>%
  select(mouse_id, clone_id, v_gene, everything())


# ======= Calculate V gene frequencies =========

naive_seq_counts <- seq_counts %>% filter(cell_type == 'naive')
exp_seq_counts <- seq_counts %>% filter(cell_type != 'naive')

gene_freqs <- calc_gene_freqs(exp_seq_counts = exp_seq_counts,
                              naive_seq_counts = naive_seq_counts,
                              clone_info = clone_info, long_format = F, by_tissue = T)

naive_freqs <- gene_freqs$naive_freqs
exp_freqs <- gene_freqs$exp_freqs

# Combine naive and experienced frequencies into a single tibble
gene_freqs <- left_join(exp_freqs,naive_freqs) %>%
  mutate(tissue = factor(tissue, levels = c('LN','spleen','BM')),
         group_controls_pooled = factor(group_controls_pooled, levels = group_controls_pooled_factor_levels))

# Gene frequencies with adjusted naive zeroes. 
# i.e. genes present in mouse but with 0 obs seqs in naive rep. are assigned 1 sequence, with freqs adjusted accordingly
gene_freqs_adj_naive_zeros <- left_join(exp_freqs, adjust_zero_naive_freqs(naive_freqs)) %>%
  mutate(tissue = factor(tissue, levels = c('LN','spleen','BM')),
         group_controls_pooled = factor(group_controls_pooled,
                                        levels = group_controls_pooled_factor_levels))

# Generate realizations under neutral model
neutral_realizations <- simulate_selection_freq_changes(exp_seq_counts = exp_seq_counts,
                                                        naive_seq_counts = naive_seq_counts,
                                                        clone_info = clone_info,
                                                        synth_data_input_tibble = 'neutral',
                                                        by_tissue = T, 
                                                        n_reps = 100)

# =========== PAIRWISE CORRELATIONS IN GENE FREQUENCIES BETWEEN MICE ==========
pairwise_gene_freqs <- get_pairwise_freqs(gene_freqs, adjust_naive_zeros = T)
pairwise_correlations <- get_pairwise_correlations(pairwise_gene_freqs)


neutral_pairwise_correlations_freqs <- lapply(neutral_realizations %>% group_by(replicate) %>% group_split(),
                                              FUN = function(x){
                                                get_pairwise_correlations(get_pairwise_freqs(x, adjust_naive_zeros = T))$freqs
                                              })

neutral_pairwise_correlations_freq_ratios <- lapply(neutral_realizations %>% group_by(replicate) %>% group_split(),
                                                    FUN = function(x){
                                                      get_pairwise_correlations(get_pairwise_freqs(x, adjust_naive_zeros = T))$freq_ratios
                                                    })

for(i in 1:length(neutral_pairwise_correlations_freqs)){
  neutral_pairwise_correlations_freqs[[i]] <- neutral_pairwise_correlations_freqs[[i]] %>%
    mutate(replicate = i) %>% select(replicate, everything())
  neutral_pairwise_correlations_freq_ratios[[i]] <- neutral_pairwise_correlations_freq_ratios[[i]] %>%
    mutate(replicate = i) %>% select(replicate, everything())
}

neutral_pairwise_correlations_freqs <- bind_rows(neutral_pairwise_correlations_freqs)
neutral_pairwise_correlations_freq_ratios <- bind_rows(neutral_pairwise_correlations_freq_ratios)

neutral_pairwise_correlations <- list(freqs = neutral_pairwise_correlations_freqs,
                                      freq_ratios = neutral_pairwise_correlations_freq_ratios)

# =========== V-GENE REGION MUTATION FREQUENCIES WITHIN CLONES ==========

mutation_freqs_within_clones <- get_mutation_frequencies_within_clones(annotated_seqs, seq_counts, by_tissue_and_cell_type = F)
mutation_freqs_within_clones_by_tissue_and_cell_type <- get_mutation_frequencies_within_clones(annotated_seqs, seq_counts, by_tissue_and_cell_type = T)



# =========== EXPORT RData ==========
save(naive_seq_counts, exp_seq_counts, gene_freqs, naive_freqs, exp_freqs, gene_freqs_adj_naive_zeros,
          neutral_realizations, pairwise_gene_freqs,
          pairwise_correlations,
          neutral_pairwise_correlations,
          mutation_freqs_within_clones,
          mutation_freqs_within_clones_by_tissue_and_cell_type,
     file = '../results/precomputed_gene_freqs.RData')


