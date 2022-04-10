library(readr)
library(dplyr)
library(tidyr)

source('gene_frequency_functions.R')


args <- commandArgs(trailingOnly = T)

frequency_type <- as.character(args[1])
use_Greiff2017_naive_freqs <- as.logical(args[2])
collapse_novel_alleles <- as.logical(args[3])


if(is.na(use_Greiff2017_naive_freqs)){
  use_Greiff2017_naive_freqs <- F
}
if(is.na(collapse_novel_alleles)){
  collapse_novel_alleles <- F
}


if(frequency_type == 'all_seqs'){
  seq_counts <- read_csv('../processed_data/seq_counts.csv')
  # seq_counts <- read_csv('~/Desktop/v_gene_selection/processed_data/seq_counts.csv')
}else{
  stopifnot(frequency_type == 'unique_seqs')
  stopifnot(!use_Greiff2017_naive_freqs)
  seq_counts <- read_csv('../processed_data/unique_seq_counts.csv')
}

output_file <- paste0('../results/precomputed_gene_freqs_', frequency_type, '.RData')


# Annotated sequences
annotated_seqs <- read_csv('../processed_data/annotated_seqs.csv')
# annotated_seqs <- read_csv('~/Desktop/v_gene_selection/processed_data/annotated_seqs.csv')

# Basic info for each clone (germline genes, CDR length, naive CDR seq)
clone_info <- read_csv('../processed_data/clone_info.csv')
# clone_info <- read_csv('~/Desktop/v_gene_selection/processed_data/clone_info.csv')

  
seq_counts <- get_info_from_mouse_id(seq_counts)

seq_counts <- left_join(seq_counts, clone_info %>% select(mouse_id, clone_id, v_gene)) %>%
  select(mouse_id, clone_id, v_gene, everything())

# If collapse_novel_alleles is true, assigns novel alleles to their reference allele
if(collapse_novel_alleles){
  stopifnot(use_Greiff2017_naive_freqs == F)
  seq_counts <- seq_counts %>% mutate(v_gene = str_remove(v_gene, '\\+.*'))
  clone_info <- clone_info %>% mutate(v_gene = str_remove(v_gene, '\\+.*'))
  output_file <- paste0('../results/precomputed_gene_freqs_', frequency_type, '_collapsed_novel_alleles.RData') 
}


# ======= Calculate V gene frequencies =========

naive_seq_counts <- seq_counts %>% filter(cell_type == 'naive')
exp_seq_counts <- seq_counts %>% filter(cell_type != 'naive')

gene_freqs <- calc_gene_freqs(exp_seq_counts = exp_seq_counts,
                              naive_seq_counts = naive_seq_counts,
                              clone_info = clone_info, long_format = F, by_tissue = T)

naive_freqs <- gene_freqs$naive_freqs
exp_freqs <- gene_freqs$exp_freqs

if(use_Greiff2017_naive_freqs){
  naive_freqs_greiff2017 <- read_csv('../processed_data/naive_freqs_Greiff2017.csv')
  # Compute average naive frequency of each V gene across mice in the naive dataset
  naive_freqs_greiff2017 <- naive_freqs_greiff2017 %>%
    group_by(v_gene) %>%
    summarise(average_greiff2017_naive_freq = mean(vgene_seq_freq)) %>%
    ungroup() %>%
    # Re-normalize so frequencies sum to 1
    mutate(average_greiff2017_naive_freq = average_greiff2017_naive_freq/sum(average_greiff2017_naive_freq))
  
  naive_freqs <- left_join(naive_freqs, 
            naive_freqs_greiff2017) %>%
    mutate(naive_vgene_seq_freq = case_when(
      is.na(average_greiff2017_naive_freq) ~ 0,
      T ~ average_greiff2017_naive_freq
    )) %>%
    select(-average_greiff2017_naive_freq) %>%
    group_by(mouse_id) %>%
    mutate(naive_vgene_seq_freq = naive_vgene_seq_freq / sum(naive_vgene_seq_freq)) %>%
    mutate(n_naive_vgene_seqs = round(naive_vgene_seq_freq * total_mouse_naive_seqs)) %>%
    ungroup()
 
  output_file <- paste0('../results/precomputed_gene_freqs_', frequency_type, '_Greiff2017_naive_freqs.RData') 
}

if(collapse_novel_alleles){
  
}


# Combine naive and experienced frequencies into a single tibble
gene_freqs <- left_join(exp_freqs,naive_freqs) %>%
  mutate(tissue = factor(tissue, levels = c('LN','spleen','BM')),
         group_controls_pooled = factor(group_controls_pooled, levels = group_controls_pooled_factor_levels))

# Gene frequencies with adjusted naive zeros. 
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

# Rho: ratio of experienced to naive frequency
obs_rhos <- gene_freqs_adj_naive_zeros %>%
  mutate(log_rho = log(vgene_seq_freq) - log(naive_vgene_seq_freq)) %>%
  mutate(obs_rho = exp(log_rho)) %>%
  rename(obs_n_vgene_seqs = n_vgene_seqs) %>%
  select(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type, v_gene,
         obs_n_vgene_seqs,obs_rho)

neutral_rhos <- neutral_realizations %>%
  mutate(log_rho = log(vgene_seq_freq) - log(naive_vgene_seq_freq)) %>%
  mutate(rho = exp(log_rho)) %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type, v_gene) %>%
  summarise(mean_sim_rho = mean(rho),
            lbound_sim_rho = quantile(rho, 0.025),
            ubound_sim_rho = quantile(rho, 0.975)) %>%
  ungroup()

deviation_from_naive <- left_join(obs_rhos, neutral_rhos) %>%
  mutate(deviation_from_naive = case_when(
    obs_rho > ubound_sim_rho ~ 'positive',
    obs_rho < lbound_sim_rho ~ 'negative',
    (obs_rho >= lbound_sim_rho) & (obs_rho <= ubound_sim_rho) ~ 'neutral'
  )) %>%
  mutate(group_controls_pooled = factor(group_controls_pooled, levels = group_controls_pooled_factor_levels)) %>%
  # deviation from naive is NA for mice lacking naive sequences in the tissues in naive_from_tissue
  filter(!is.na(deviation_from_naive))

gene_freqs <- left_join(gene_freqs, 
                        deviation_from_naive %>%
                          select(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type, v_gene, matches('rho'), deviation_from_naive))

gene_freqs <- gene_freqs %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type) %>%
  mutate(v_gene_rank = rank(-vgene_seq_freq, ties.method = 'first')) %>%
  ungroup()



# =========== PAIRWISE CORRELATIONS IN GENE FREQUENCIES BETWEEN MICE ==========

# Pairwise correlations in observed data
pairwise_gene_freqs <- get_pairwise_freqs(gene_freqs, adjust_naive_zeros = T)
pairwise_correlations <- get_pairwise_correlations(pairwise_gene_freqs)
deviation_concordance <- compute_deviation_concordance(pairwise_gene_freqs)


# Pairwise correlations with randomized noncontrol groups
gene_freqs_randomized_noncontrol_groups <- replicate(500, randomize_noncontrol_groups(gene_freqs),
                                   simplify = F)

pairwise_correlations_randomized_noncontrol_groups_FREQS <- tibble()
pairwise_correlations_randomized_noncontrol_groups_FREQ_RATIOS <- tibble()

# 
for(i in 1:500){
  gene_freqs_randomized_noncontrol_groups <- randomize_noncontrol_groups(gene_freqs)
  pairwise_gene_freqs_randomized_noncontrol_groups <- get_pairwise_freqs(gene_freqs_randomized_noncontrol_groups, 
                                                                                adjust_naive_zeros = T)
  pw_corr <- get_pairwise_correlations(pairwise_gene_freqs_randomized_noncontrol_groups)
  
  pairwise_correlations_randomized_noncontrol_groups_FREQS <- bind_rows(
    pairwise_correlations_randomized_noncontrol_groups_FREQS,
    pw_corr$freqs %>% mutate(replicate = i)
  )
  
  pairwise_correlations_randomized_noncontrol_groups_FREQ_RATIOS <- bind_rows(
    pairwise_correlations_randomized_noncontrol_groups_FREQ_RATIOS,
    pw_corr$freq_ratios %>% mutate(replicate = i)
  )
  
}


pairwise_correlations_randomized_noncontrol_groups <- list(freqs = pairwise_correlations_randomized_noncontrol_groups_FREQS %>%
                                                             select(replicate, everything()),
                                                           freq_ratios = pairwise_correlations_randomized_noncontrol_groups_FREQ_RATIOS %>%
                                                             select(replicate, everything()))


# =========== CLONE FREQUENCIES RELATIVE TO THE TOTAL IN EACH TISSUE / CELL TYPE COMBINATION ==========

if(frequency_type == 'unique_seqs'){
  seq_counts <- seq_counts %>%
    rename(prod_seqs = unique_prod_seqs)
}

clone_freqs_by_tissue_and_cell_type <- seq_counts %>%
  dplyr::rename(n_clone_seqs_in_compartment = prod_seqs) %>%
  #group_by(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type, clone_id, v_gene) %>%
  #summarise(n_clone_seqs_in_compartment = sum(prod_seqs)) %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type) %>%
  mutate(total_seqs_in_compartment = sum(n_clone_seqs_in_compartment),
         clone_freq_in_compartment = n_clone_seqs_in_compartment / total_seqs_in_compartment) %>%
  mutate(clone_rank_in_compartment = rank(-clone_freq_in_compartment, ties.method = 'first')) %>%
  ungroup() %>%
  mutate(group_controls_pooled = factor(group_controls_pooled, levels = group_controls_pooled_factor_levels)) %>%
  dplyr::rename(compartment_tissue = tissue, compartment_cell_type = cell_type)

clone_freqs_by_tissue <- seq_counts %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, tissue, clone_id, v_gene) %>%
  summarise(n_clone_seqs_in_compartment = sum(prod_seqs)) %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, tissue) %>%
  mutate(total_seqs_in_compartment = sum(n_clone_seqs_in_compartment),
         clone_freq_in_compartment = n_clone_seqs_in_compartment / total_seqs_in_compartment) %>%
  mutate(clone_rank_in_compartment = rank(-clone_freq_in_compartment, ties.method = 'first')) %>%
  ungroup() %>%
  mutate(group_controls_pooled = factor(group_controls_pooled, levels = group_controls_pooled_factor_levels)) %>%
  dplyr::rename(compartment_tissue = tissue)

# Clone frequencies across the entire mouse
clone_freqs <- clone_freqs_by_tissue %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, clone_id) %>%
  summarise(n_clone_seqs_in_compartment = sum(n_clone_seqs_in_compartment)) %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled) %>%
  mutate(total_seqs_in_compartment = sum(n_clone_seqs_in_compartment)) %>%
  mutate(compartment_tissue = 'all') %>%
  ungroup()

# Add general clone information to object with clone size distributions
clone_freqs_by_tissue_and_cell_type <- left_join(clone_freqs_by_tissue_and_cell_type, clone_info)


# Add information on whether the V gene used by each clone has evidence of gene-level selection
clone_freqs_by_tissue_and_cell_type <- left_join(clone_freqs_by_tissue_and_cell_type,
                                                 deviation_from_naive %>%
                                                   select(mouse_id, day, infection_status, group_controls_pooled,
                                                          tissue, cell_type, v_gene, matches('rho'), deviation_from_naive) %>%
                                                   dplyr::rename(compartment_tissue = tissue, compartment_cell_type = cell_type))

clone_freqs_by_tissue_and_cell_type <- left_join(clone_freqs_by_tissue_and_cell_type,
                                                 gene_freqs %>% select(mouse_id, total_mouse_naive_seqs) %>% unique()) %>%
  mutate(compartment_tissue = factor(compartment_tissue, levels = c('LN','spleen','BM')),
         compartment_cell_type = factor(compartment_cell_type, levels = c('naive', 'nonnaive_IgD+B220+','GC','PC','mem')))

# =========== V-GENE REGION MUTATION FREQUENCIES WITHIN CLONES ==========

if(frequency_type == 'all_seqs'){
  mutation_freqs_within_clones <- get_mutation_frequencies_within_clones(annotated_seqs, seq_counts, by_tissue_and_cell_type = F)
  mutation_freqs_within_clones_by_tissue_and_cell_type <- get_mutation_frequencies_within_clones(annotated_seqs, seq_counts, by_tissue_and_cell_type = T)
  
  # For each clone, include a list of mutations above a certain frequency threshold
  # (For now, frequency is calculated relative to the number of productive sequences from a clone in a particular cell type and tissue combination)
  mutations_above_threshold <- list_clone_mutations_above_threshold(mutation_freqs_within_clones_by_tissue_and_cell_type,
                                                                    threshold = 0.5)
  
  
  clone_freqs_by_tissue_and_cell_type <- left_join(clone_freqs_by_tissue_and_cell_type,
                                                   mutations_above_threshold %>% select(mouse_id, clone_id, compartment_tissue,
                                                                                        compartment_cell_type,
                                                                                        mutations_above_threshold)) %>%
    mutate(mutations_above_threshold = ifelse(is.na(mutations_above_threshold), '', mutations_above_threshold))
  
  
}else{
  mutation_freqs_within_clones <- 'Not defined when using unique seqs'
  mutation_freqs_within_clones_by_tissue_and_cell_type <- 'Not defined when using unique seqs'
}





# =========== EXPORT RData ==========
save(naive_seq_counts, exp_seq_counts, gene_freqs, naive_freqs, exp_freqs, gene_freqs_adj_naive_zeros,
     clone_freqs_by_tissue_and_cell_type, 
     clone_freqs_by_tissue,
     clone_freqs,
     neutral_realizations,
     deviation_from_naive,
     pairwise_gene_freqs,
     pairwise_correlations,
     #neutral_pairwise_correlations,
     pairwise_correlations_randomized_noncontrol_groups,
     mutation_freqs_within_clones,
     mutation_freqs_within_clones_by_tissue_and_cell_type,
     file = output_file)


