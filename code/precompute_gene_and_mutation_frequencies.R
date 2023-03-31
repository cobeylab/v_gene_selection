library(readr)
source('gene_frequency_functions.R')

high_frequency_mutation_threshold <- 0.5 # Threshold for calling high-frequency mutations.

args <- commandArgs(trailingOnly = T)

# Frequencies based on either "all_seqs" (all productive sequences) or "unique_seqs" (unique productive sequences)
frequency_type <- as.character(args[1]) 
use_Greiff2017_naive_freqs <- as.logical(args[2]) # Use alternative dataset to estimate naive freqs?
collapse_novel_alleles <- as.logical(args[3]) # Ignore novel alleles identified by partis?
assignment <- as.character(args[4]) # VDJ/Clone assignment method ('partis', 'igblast')

if(is.na(use_Greiff2017_naive_freqs)){
  use_Greiff2017_naive_freqs <- F
}
if(is.na(collapse_novel_alleles)){
  collapse_novel_alleles <- F
}
if(is.na(assignment)){
  assignment <- 'partis'
}

n_null_model_realizations <- 500

if(assignment == 'partis'){
  # Import basic info for each clone (germline genes, CDR length, naive CDR seq)
  clone_info <- read_csv('../processed_data/clone_info_partis.csv')
  
  # Read sequences counts by mouse/clone/cell type
  if(frequency_type == 'all_seqs'){
    seq_counts <- read_csv('../processed_data/seq_counts_partis.csv')
  }else{
    stopifnot(frequency_type == 'unique_seqs')
    stopifnot(!use_Greiff2017_naive_freqs)
    seq_counts <- read_csv('../processed_data/unique_seq_counts_partis.csv')
  }
  
  # Path to exported RData object.
  output_file <- paste0('../results/precomputed_gene_freqs_', frequency_type, '.RData')
  
  # If collapse_novel_alleles is true, assigns novel alleles to their reference allele
  if(collapse_novel_alleles){
    stopifnot(use_Greiff2017_naive_freqs == F)
    seq_counts <- seq_counts %>% mutate(v_gene = str_remove(v_gene, '\\+.*'))
    clone_info <- clone_info %>% mutate(v_gene = str_remove(v_gene, '\\+.*'))
    output_file <- paste0('../results/precomputed_gene_freqs_', frequency_type, '_collapsed_novel_alleles.RData') 
  }
}else{
  stopifnot(assignment == 'igblast')
  stopifnot(frequency_type == 'all_seqs')
  stopifnot(!use_Greiff2017_naive_freqs)
  stopifnot(!collapse_novel_alleles | is.na(collapse_novel_alleles))

  
  clone_info <- read_csv('../processed_data/clone_info_igblast.csv')
  seq_counts <- read_csv('../processed_data/unique_seq_counts_igblast.csv')
  
  output_file <- paste0('../results/precomputed_gene_freqs_', frequency_type, '_igblast_assignment.RData')
  
}

# Import annotated sequences
annotated_seqs <- read_csv('../processed_data/annotated_seqs.csv')

# Add mouse information to sequence counts object, then left join clone info
seq_counts <- get_info_from_mouse_id(seq_counts)

seq_counts <- left_join(seq_counts, clone_info %>% select(mouse_id, clone_id, v_gene)) %>%
  select(mouse_id, clone_id, v_gene, everything())


# =========== CALCULATING V GENE FREQUENCIES ==========
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

# Combine naive and experienced frequencies into a single tibble
gene_freqs <- format_gene_freqs_wide(exp_freqs = exp_freqs, naive_freqs = naive_freqs)

# Gene frequencies with adjusted naive zeros. 
# i.e. genes present in mouse but with 0 obs seqs in naive rep. are assigned 1 sequence, with freqs adjusted accordingly
gene_freqs_adj_naive_zeros <- format_gene_freqs_wide(exp_freqs = exp_freqs, naive_freqs = adjust_zero_naive_freqs(naive_freqs))

# Rho: ratio of experienced to naive frequency
obs_rhos <- compute_rho(gene_freqs_adj_naive_zeros) %>%
  rename(obs_n_vgene_seqs = n_vgene_seqs) %>%
  select(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type, v_gene,
         obs_n_vgene_seqs,obs_rho)
  
# Generate allele freqs under null model where experienced frequencies are sampled from naive freqs.
neutral_realizations <- sample_experienced_from_naive(exp_seq_counts = exp_seq_counts,
                                                        naive_seq_counts = naive_seq_counts,
                                                        clone_info = clone_info,
                                                        synth_data_input_tibble = 'neutral',
                                                        by_tissue = T, 
                                                        n_reps = n_null_model_realizations)

# Compare observed vs neutral experienced-to-naive ratios
neutral_rhos <- neutral_realizations %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type, v_gene) %>%
  summarise(mean_sim_rho = mean(obs_rho),
            lbound_sim_rho = quantile(obs_rho, 0.025),
            ubound_sim_rho = quantile(obs_rho, 0.975)) %>%
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

# Add comparison with null model to gene_freqs object.
gene_freqs <- left_join(gene_freqs, 
                        deviation_from_naive %>%
                          select(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type, v_gene, matches('rho'), deviation_from_naive))

gene_freqs <- gene_freqs %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type) %>%
  mutate(v_gene_rank = rank(-vgene_seq_freq, ties.method = 'first')) %>%
  ungroup()


# =========== PAIRWISE CORRELATIONS IN GENE FREQUENCIES BETWEEN MICE ==========

# Pairwise correlations in observed data
pairwise_gene_freqs <- get_pairwise_freqs(gene_freqs, adjust_naive_zeros = T, within_groups_only =F)
pairwise_correlations <- get_pairwise_correlations(pairwise_gene_freqs)
deviation_concordance <- compute_deviation_concordance(pairwise_gene_freqs)


# Pairwise correlations under null model where lineage sizes are kept fixed but lineage alleles rdmly assigned based on naive freqs
null_model_random_lineage_alleles <- replicate(n_null_model_realizations,
                                               gene_freqs_random_lineage_alleles(exp_seq_counts = exp_seq_counts,
                                                                                 naive_seq_counts = naive_seq_counts,
                                                                                 naive_freqs = naive_freqs,
                                                                                 adjust_naive_zeros = T),
                                               simplify = F)

pairwise_correlations_randomized_lineage_V_alleles <- get_pairwise_corrs_for_randomized_datasets(
  randomized_gene_freqs_list = null_model_random_lineage_alleles, adjust_naive_zeros = T, within_groups_only = T, focal_tissue = 'LN')


# Pairwise correlations with randomized noncontrol groups
null_model_randomized_noncontrol_groups <- replicate(n_null_model_realizations,
                                                     randomize_noncontrol_groups(gene_freqs %>% filter(tissue == 'LN')),
                                                     simplify = F)

pairwise_correlations_randomized_noncontrol_groups <- get_pairwise_corrs_for_randomized_datasets(
  randomized_gene_freqs_list = null_model_randomized_noncontrol_groups, adjust_naive_zeros = T, within_groups_only = T
)


# =========== CLONE FREQUENCIES RELATIVE TO THE TOTAL IN EACH COMPARTMENT ==========
# Clone frequencies by tissue 
clone_freqs_by_tissue <- get_clone_freqs(seq_counts, compartment_vars = c('tissue')) %>%
  mutate(group_controls_pooled = factor(group_controls_pooled, levels = group_controls_pooled_factor_levels)) 

# Clone frequencies by tissue and cell type combination
clone_freqs_by_tissue_and_cell_type <- get_clone_freqs(seq_counts, compartment_vars = c('tissue','cell_type')) %>%
  mutate(group_controls_pooled = factor(group_controls_pooled, levels = group_controls_pooled_factor_levels)) 
  
# Add general clone information to objects clone freqs objects
clone_freqs_by_tissue <- left_join(clone_freqs_by_tissue, clone_info)
clone_freqs_by_tissue_and_cell_type <- left_join(clone_freqs_by_tissue_and_cell_type, clone_info)


# Add information on V gene used by each clone
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

if(frequency_type == 'all_seqs' & assignment == 'partis'){
  mutation_freqs_within_clones <- get_mutation_frequencies_within_clones(annotated_seqs, seq_counts, by_tissue_and_cell_type = F)
  mutation_freqs_within_clones_by_tissue_and_cell_type <- get_mutation_frequencies_within_clones(annotated_seqs, seq_counts, by_tissue_and_cell_type = T)
  
  # For each clone, include a list of mutations above a certain frequency threshold
  # (For now, frequency is calculated relative to the number of productive sequences from a clone in a particular cell type and tissue combination)
  mutations_above_threshold <- list_clone_mutations_above_threshold(mutation_freqs_within_clones_by_tissue_and_cell_type,
                                                                    threshold = high_frequency_mutation_threshold)
  
  
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
     clone_freqs_by_tissue,
     clone_freqs_by_tissue_and_cell_type, 
     neutral_realizations,
     deviation_from_naive,
     pairwise_gene_freqs,
     pairwise_correlations,
     pairwise_correlations_randomized_noncontrol_groups,
     pairwise_correlations_randomized_lineage_V_alleles,
     mutation_freqs_within_clones,
     mutation_freqs_within_clones_by_tissue_and_cell_type,
     file = output_file)


