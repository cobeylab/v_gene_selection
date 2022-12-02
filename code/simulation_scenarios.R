library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(readr)
theme_set(theme_cowplot())

source('gene_frequency_functions.R')

selected_allele_eligibility_threshold <- 40 # Only alleles occurring in at least these many mice can be under positive allele-level selection
selected_allele_naive_freq_interval <- c(0.02, 0.03) # Only alleles with a mean naive freq. across mice of 2-3% can be chosen.

min_naive_seqs <- 1000 # Only use mice with at least 1000 naive seqs as a base for simulations

# To use realistic naive frequencies, import precomputed gene frequencies object

load('../results/precomputed_gene_freqs_all_seqs.RData')


# Remove mice with fewer than min_naive_seqs
obs_naive_freqs <- naive_freqs %>% 
  filter(total_mouse_naive_seqs >= min_naive_seqs)

# Adjust naive zeros.
obs_naive_freqs <- adjust_zero_naive_freqs(obs_naive_freqs)

# ===== Functions for creating scenario input files =====

# For testing, takes obs_naive_freqs and assigns one allele in each mouse (same across mice) a naive frequency of dominant_allele_naive_freq
# (then assign uniform frequencies for the rest of the alleles)
set_dominant_allele <- function(obs_naive_freqs, dominant_allele_naive_freq){
  total_n_mice <- length(unique(obs_naive_freqs$mouse_id))
  
  shared_alleles <- obs_naive_freqs %>%
    group_by(v_gene) %>%
    summarise(n_mice_allele_occurs = length(unique(mouse_id))) %>%
    filter(n_mice_allele_occurs == total_n_mice) %>%
    pull(v_gene)
  
  dominant_allele_in_naive_rep <- sample(shared_alleles, size = 1)
  
  changed_naive_freqs <- obs_naive_freqs %>%
    # Sets frequency of selected dominant allele in naive repertoire to dominant_allele_naive_freq, others to NA
    mutate(n_naive_vgene_seqs = NA,
           naive_vgene_seq_freq = ifelse(v_gene == dominant_allele_in_naive_rep, dominant_allele_naive_freq, NA)) %>%
    group_by(mouse_id, day, infection_status, group, group_controls_pooled, total_mouse_naive_seqs) %>%
    mutate(n_alleles_in_mice = length(unique(v_gene))) %>%
    ungroup() %>%
    # Assign uniform frequencies to non-dominant alleles
    mutate(naive_vgene_seq_freq = case_when(
      is.na(naive_vgene_seq_freq) ~ (1 - dominant_allele_naive_freq)/(n_alleles_in_mice-1),
      !is.na(naive_vgene_seq_freq) ~ naive_vgene_seq_freq
    )) 
  
  return(changed_naive_freqs)

}

generate_allele_info <- function(obs_naive_freqs, n_high_avg_alleles, n_high_mutability_alleles,
                                 selected_allele_eligibility_threshold, selected_allele_naive_freq_interval,
                                 dominant_allele_naive_freq = NULL){
  # Assigns affinity distributions to alleles. For each allele, the distribution is the same in all individuals where it occurs 
  # Uses empirical allele sets and naive frequencies. 
  # The number of alleles with "high-average" is input, and so is the increase in mean affinity associated with using them
  # Only alleles that occur in at least n = selected_allele_eligibility_threshold mice go can go in those categories
  # All alleles not in the "high-average" category go in the "low-average category".
  # Same logic with highly mutable alleles (mutability is 1 for low mutability alleles, gamma for highly mutable alleles)
  
  # if dominant_allele_naive_freq is specified, assign that naive frequency to a randomly chosen allele shared by mice
  if(!is.null(dominant_allele_naive_freq)){
    # We'll use this option just for the neutral case where all alleles are equivalent
    stopifnot(n_high_avg_alleles == 0)
    stopifnot(n_high_mutability_alleles == 0)
    
    obs_naive_freqs <- set_dominant_allele(obs_naive_freqs, dominant_allele_naive_freq = dominant_allele_naive_freq)
    
  }
  
  n_alleles_in_data <- length(obs_naive_freqs %>% select(v_gene) %>% unique() %>% pull(v_gene))
  
  #First, assign all alleles in all mice to the low_avg, low mutability category
  allele_info <- obs_naive_freqs %>%
    select(mouse_id, v_gene, naive_vgene_seq_freq) %>%
    dplyr::rename(naive_freq = naive_vgene_seq_freq) %>%
    mutate(allele_type_affinity = 'low_avg',
           allele_type_mutability = 'low_mut')
  
  # Pick alleles that meet criteria for potential assignment to high avg or high mut. class. 
  candidate_selected_alleles <- obs_naive_freqs %>%
    group_by(mouse_id) %>%
    mutate(rank_in_naive_rep = rank(naive_vgene_seq_freq, ties.method = 'average')) %>%
    ungroup() %>%
    group_by(v_gene) %>%
    summarise(n_mice = length(unique(mouse_id)),
              mean_naive_freq = mean(naive_vgene_seq_freq)) %>%
    filter(n_mice >= selected_allele_eligibility_threshold,
           mean_naive_freq >= selected_allele_naive_freq_interval[1], mean_naive_freq <= selected_allele_naive_freq_interval[2]) %>%
    pull(v_gene)
  
  if(n_high_avg_alleles > 0){
    if(n_high_avg_alleles > length(candidate_selected_alleles)){
      stop("Selected number of high-affinity alleles is greater than number of alleles meeting selection criteria")
    }
    
    high_avg_alleles <- sample(candidate_selected_alleles,size = n_high_avg_alleles, replace = F)
    allele_info$allele_type_affinity[allele_info$v_gene %in% high_avg_alleles] <- 'high_avg'
  }
  
  if(n_high_mutability_alleles > 0){
    if(n_high_mutability_alleles > length(candidate_selected_alleles)){
      stop("Selected number of high-mutability alleles is greater than number of alleles meeting selection criteria")
    }
    high_mutability_alleles <- sample(candidate_selected_alleles, size = n_high_mutability_alleles, replace = F)
    allele_info$allele_type_mutability[allele_info$v_gene %in% high_mutability_alleles] <- 'high_mut'
  }
  
  # Replace mouse and gene ids with arbitrary integer ids
  individual_integer_ids <- obs_naive_freqs %>% select(mouse_id) %>% unique() %>% mutate(individual = 1:n())
  
  allele_integer_ids <- obs_naive_freqs %>% select(v_gene) %>% unique() %>% mutate(allele = 1:n()) %>%
    mutate(allele = paste0('V', allele))
  
  allele_info <- left_join(allele_info, individual_integer_ids) %>% select(-mouse_id)
  allele_info <- left_join(allele_info, allele_integer_ids) %>% select(-v_gene) %>%
    select(individual, allele, allele_type_affinity, allele_type_mutability, naive_freq)
  

  return(allele_info)
  

}


create_scenario <- function(scenario_directory, obs_naive_freqs, selected_allele_eligibility_threshold,
                            selected_allele_naive_freq_interval, n_high_avg_alleles, s, sigma_r, n_high_mutability_alleles,
                            gamma, K, I_total, t_imm, mu_max, delta, mutation_rate, beta, tmax, uniform_naive_freqs,
                            dominant_allele_naive_freq = NULL){
  
  # Create base directory for the scenario.
  dir.create(scenario_directory, showWarnings = F, recursive = T)
  raw_simulations_dir <- paste0(scenario_directory, 'raw_simulation_files/')
  dir.create(raw_simulations_dir, showWarnings = F)
  
  # If multiple parameter combinations sampled within scenario, they all share same allele info (naive freqs, selection)
  allele_info <- generate_allele_info(obs_naive_freqs = obs_naive_freqs,
                                      n_high_avg_alleles = n_high_avg_alleles,
                                      n_high_mutability_alleles = n_high_mutability_alleles,
                                      selected_allele_eligibility_threshold = selected_allele_eligibility_threshold,
                                      selected_allele_naive_freq_interval = selected_allele_naive_freq_interval,
                                      dominant_allele_naive_freq = dominant_allele_naive_freq)
  
  if(uniform_naive_freqs){
    stopifnot(is.null(dominant_allele_naive_freq))
    allele_info <- allele_info %>% 
      group_by(individual) %>%
      mutate(naive_freq = 1 /n())
  }
  
  # Export allele information for simulations
  write_csv(allele_info, paste0(scenario_directory, 'allele_info.csv'))
  

  # Identify parameters with more than one value
  par_combinations <- expand_grid(K = K, I_total = I_total, t_imm = t_imm, mu_max = mu_max, delta = delta,
                                  s = s, sigma_r = sigma_r, gamma = gamma, mutation_rate = mutation_rate,
                                  beta = beta, tmax = tmax, uniform_naive_freqs = uniform_naive_freqs)
  
    
  n_par_values <- par_combinations %>%
    summarise(across(everything(),.fns = function(x){length(unique(x))})) %>% unlist()
  variable_pars <- names(n_par_values)[n_par_values > 1]
  
  
  for(i in 1:nrow(par_combinations)){
    if(length(variable_pars) > 0){
      
      par_values_label <- paste(paste(variable_pars, par_combinations[i, variable_pars] %>% unlist(), sep = '_'),
                                collapse = '_')
    }else{
      par_values_label <- 'single_par_combination'
    }
      par_combination_subfolder <- paste0(raw_simulations_dir, par_values_label)
      dir.create(par_combination_subfolder, showWarnings = F)
      
      # Export model parameters for this specific parameter combination
      write_csv(par_combinations[i,],
                paste0(par_combination_subfolder, '/model_parameters.csv'))
  }
  
}


# ============================ NEUTRAL SCENARIO ===================================
create_scenario(scenario_directory = '../results/simulations/neutral_scenario/',
                obs_naive_freqs = obs_naive_freqs,
                selected_allele_eligibility_threshold = selected_allele_eligibility_threshold,
                selected_allele_naive_freq_interval = selected_allele_naive_freq_interval,
                n_high_avg_alleles = 0,
                s = 0,
                sigma_r = 1,
                n_high_mutability_alleles = 0,
                gamma = 1,
                K = 2000,
                I_total = c(50,100,200),
                t_imm = 6,
                mu_max = 3,
                delta = 0.2,
                mutation_rate = c(0, 0.01, 0.05),
                beta = c(1,2,3,4),
                tmax = 50,
                uniform_naive_freqs = F)

# ============================ NEUTRAL SCENARIO WITH UNIFORM NAIVE FREQS ===================================
create_scenario(scenario_directory = '../results/simulations/neutral_uniform_freqs_scenario/',
                obs_naive_freqs = obs_naive_freqs,
                selected_allele_eligibility_threshold = selected_allele_eligibility_threshold,
                selected_allele_naive_freq_interval = selected_allele_naive_freq_interval,
                n_high_avg_alleles = 0,
                s = 0,
                sigma_r = 1,
                n_high_mutability_alleles = 0,
                gamma = 1,
                K = 2000,
                I_total = c(50,100,200),
                t_imm = 6,
                mu_max = 3,
                delta = 0.2,
                mutation_rate = c(0, 0.01, 0.05),
                beta = c(1,2,3,4),
                tmax = 50,
                uniform_naive_freqs = T)

# ============================ HIGH AFFINITY SCENARIO ===================================
# FINAL SCENARIO SHOWING DIFFERENCES IN ALLELE'S AFFINITIES
create_scenario(scenario_directory = '../results/simulations/high_affinity_scenario/',
                obs_naive_freqs = obs_naive_freqs,
                selected_allele_eligibility_threshold = selected_allele_eligibility_threshold,
                selected_allele_naive_freq_interval = selected_allele_naive_freq_interval,
                n_high_avg_alleles = 5,
                s = c(0.5, 1, 1.5, 2),
                sigma_r = 1,
                n_high_mutability_alleles = 0,
                gamma = 1,
                K = 2000,
                I_total = 200,
                t_imm = 6,
                mu_max = 3,
                delta = 0.2,
                mutation_rate = c(0, 0.01, 0.05),
                beta = c(1,2,3,4),
                tmax = 50,
                uniform_naive_freqs = F)

# ============================ HIGH MUTATION SCENARIO ===================================
# some alleles have higher mutation rate by a factor gamma
create_scenario(scenario_directory = '../results/simulations/high_mutation_scenario/',
                obs_naive_freqs = obs_naive_freqs,
                selected_allele_eligibility_threshold = selected_allele_eligibility_threshold,
                selected_allele_naive_freq_interval = selected_allele_naive_freq_interval,
                n_high_avg_alleles = 0,
                s = 1,
                sigma_r = 1,
                n_high_mutability_alleles = 5,
                gamma = c(1.5,2,4,6),
                K = 2000,
                I_total = 200,
                t_imm = 6,
                mu_max = 3,
                delta = 0.2,
                mutation_rate = c(0, 0.01, 0.05),
                beta = c(1,2,3,4),
                tmax = 50,
                uniform_naive_freqs = F)


# Test scenarios

# === TEST SCENARIO: neutral with a highly dominant allele in naive repertoire (0.9)
create_scenario(scenario_directory = '../results/simulations/test_scenario_neutral_highly_dominant/',
                obs_naive_freqs = obs_naive_freqs,
                selected_allele_eligibility_threshold = selected_allele_eligibility_threshold,
                selected_allele_naive_freq_interval = selected_allele_naive_freq_interval,
                n_high_avg_alleles = 0,
                s = 0,
                sigma_r = 1,
                n_high_mutability_alleles = 0,
                gamma = 1,
                K = 2000,
                I_total = c(50,200),
                t_imm = 6,
                mu_max = 3,
                delta = 0.2,
                mutation_rate = c(0, 0.05),
                beta = c(1,4),
                tmax = 50,
                uniform_naive_freqs = F,
                dominant_allele_naive_freq = 0.999)

# === TEST SCENARIO: neutral with uniform frequencies in the naive repertoire
create_scenario(scenario_directory = '../results/simulations/test_scenario_neutral_uniform/',
                obs_naive_freqs = obs_naive_freqs,
                selected_allele_eligibility_threshold = selected_allele_eligibility_threshold,
                selected_allele_naive_freq_interval = selected_allele_naive_freq_interval,
                n_high_avg_alleles = 0,
                s = 0,
                sigma_r = 1,
                n_high_mutability_alleles = 0,
                gamma = 1,
                K = 2000,
                I_total = c(50,200),
                t_imm = 6,
                mu_max = 3,
                delta = 0.2,
                mutation_rate = c(0, 0.05),
                beta = c(1,4),
                tmax = 50,
                uniform_naive_freqs = T,
                dominant_allele_naive_freq = NULL)





