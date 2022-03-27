library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(readr)
theme_set(theme_cowplot())

source('gene_frequency_functions.R')

selected_allele_eligibility_threshold <- 40 # Only alleles occurring in at least these many mice can be under positive allele-level selection

min_naive_seqs <- 1000 # Only use mice with at least 1000 naive seqs as a base for simulations

# Scenarios are specified by alleles' affinity distributions and allele frequencies and 
# by germinal center parameters (see simulation_functions.R)
# Each class of allele has a mean alpha and a mean beta parameter.
# Specific alphas and betas for each allele are obtained by adding a normally distributed error around 0 to their class mean values
# Affinities for B cells using each allele are then distributed according to those specific alphas and betas.


# To use realistic naive frequencies, import precomputed gene frequencies object

load('../results/precomputed_gene_freqs_all_seqs.RData')
#load('~/Desktop/v_gene_selection/results/precomputed_gene_freqs_all_seqs.RData')
obs_naive_freqs <- naive_freqs %>%
  filter(total_mouse_naive_seqs >= min_naive_seqs)

obs_naive_freqs <- adjust_zero_naive_freqs(obs_naive_freqs)


generate_allele_info <- function(obs_naive_freqs, n_high_avg_alleles, n_high_mutability_alleles,
                                 selected_allele_eligibility_threshold){
  # Assigns affinity distributions to alleles. For each allele, the distribution is the same in all individuals where it occurs 
  # Uses empirical allele sets and naive frequencies. 
  # The number of alleles with "high-average" is input, and so is the increase in mean affinity associated with using them
  # Only alleles that occur in at least n = selected_allele_eligibility_threshold mice go can go in those categories
  # All alleles not in the "high-average" category go in the "low-average category".
  # Same logic with highly mutable alleles (mutability is 1 for low mutability alleles, gamma for highly mutable alleles)
  
  
  n_alleles_in_data <- length(obs_naive_freqs %>% select(v_gene) %>% unique() %>% pull(v_gene))
  
  
  #First, assign all alleles in all mice to the low_avg, low mutability category
  allele_info <- obs_naive_freqs %>%
    select(mouse_id, v_gene, naive_vgene_seq_freq) %>%
    dplyr::rename(naive_freq = naive_vgene_seq_freq) %>%
    mutate(allele_type_affinity = 'low_avg',
           allele_type_mutability = 'low_mut')
  
  candidate_selected_alleles <- obs_naive_freqs %>% group_by(v_gene) %>% filter(naive_vgene_seq_freq > 0) %>%
    summarise(n_mice = length(unique(mouse_id))) %>%
    filter(n_mice >= selected_allele_eligibility_threshold) %>% pull(v_gene)
  
  if(n_high_avg_alleles > 0){
    high_avg_alleles <- sample(candidate_selected_alleles,size = n_high_avg_alleles, replace = F)
    allele_info$allele_type_affinity[allele_info$v_gene %in% high_avg_alleles] <- 'high_avg'
  }
  
  if(n_high_mutability_alleles > 0){
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
                            n_high_avg_alleles, s, sigma_r, n_high_mutability_alleles, gamma, K, I_total, t_imm, mu_max, delta,
                            mutation_rate, beta, tmax, uniform_naive_freqs, fix_initial_affinities){
  
  # Create base directory for the scenario.
  dir.create(scenario_directory, showWarnings = F)
  raw_simulations_dir <- paste0(scenario_directory, 'raw_simulation_files/')
  dir.create(raw_simulations_dir, showWarnings = F)
  
  # If multiple parameter combinations sampled within scenario, they all share same allele info (naive freqs, selection)
  allele_info <- generate_allele_info(obs_naive_freqs = obs_naive_freqs,
                                      n_high_avg_alleles = n_high_avg_alleles,
                                      n_high_mutability_alleles = n_high_mutability_alleles,
                                      selected_allele_eligibility_threshold = selected_allele_eligibility_threshold)
  
  if(uniform_naive_freqs){
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



# ============================ SCENARIO 0 ===================================
# Truly neutral scenario, where all naive B cells have exactly the same affinity regardless of V gene
# This is achieved by making sigma r = 0. Since the effect of mutations is sigma_r times beta, 
# we set sigma_r = 0.0001 (effectively no variation in naive B cell affinity) and beta = 20000 to allow mutations to have an effect with sd=2

# create_scenario(scenario_directory = '../results/simulations/scenario_0/',
#                 obs_naive_freqs = obs_naive_freqs,
#                 selected_allele_eligibility_threshold = selected_allele_eligibility_threshold,
#                 n_high_avg_alleles = 0,
#                 s = 0,
#                 sigma_r = 0.0001,
#                 n_high_mutability_alleles = 0,
#                 gamma = 1,
#                 K = 2000,
#                 I_total = c(50,100,200),
#                 t_imm = 6,
#                 mu_max = 3,
#                 delta = 0.2,
#                 mutation_rate = c(0,0.01, 0.05),
#                 beta = 20000,
#                 tmax = 50,
#                 uniform_naive_freqs = F)


# ============================ SCENARIO 1 ===================================
# create_scenario(scenario_directory = '../results/simulations/scenario_1/',
#                 obs_naive_freqs = obs_naive_freqs,
#                 selected_allele_eligibility_threshold = selected_allele_eligibility_threshold,
#                 n_high_avg_alleles = 0,
#                 s = 0,
#                 sigma_r = c(0.5,1,1.5,2),
#                 n_high_mutability_alleles = 0,
#                 gamma = 1,
#                 K = 2000,
#                 I_total = 100,
#                 t_imm = 6,
#                 mu_max = 3,
#                 delta = 0.2,
#                 mutation_rate = 0.01,
#                 beta = c(0,0.5,1,2),
#                 tmax = 50,
#                 uniform_naive_freqs = F)




# ============================ SCENARIO 2 ===================================
# create_scenario(scenario_directory = '../results/simulations/scenario_2/',
#                 obs_naive_freqs = obs_naive_freqs,
#                 selected_allele_eligibility_threshold = selected_allele_eligibility_threshold,
#                 n_high_avg_alleles = 5,
#                 s = c(0, 0.5, 1, 2),
#                 sigma_r = 1,
#                 n_high_mutability_alleles = 0,
#                 gamma = 1,
#                 K = 2000,
#                 I_total = 100,
#                 t_imm = 6,
#                 mu_max = 3,
#                 delta = 0.2,
#                 mutation_rate = c(0, 0.01,0.05),
#                 beta = 1,
#                 tmax = 50,
#                 uniform_naive_freqs = F)


# ============================ SCENARIO 2b ===================================
# Like 2, but with much higher effect of mutation relative to sigma_r (beta = 4)
# create_scenario(scenario_directory = '../results/simulations/scenario_2b/',
#                 obs_naive_freqs = obs_naive_freqs,
#                 selected_allele_eligibility_threshold = selected_allele_eligibility_threshold,
#                 n_high_avg_alleles = 5,
#                 s = c(0, 0.5, 1, 2),
#                 sigma_r = 1,
#                 n_high_mutability_alleles = 0,
#                 gamma = 1,
#                 K = 2000,
#                 I_total = 100,
#                 t_imm = 6,
#                 mu_max = 3,
#                 delta = 0.2,
#                 mutation_rate = c(0, 0.01,0.05),
#                 beta = 4,
#                 tmax = 50,
#                 uniform_naive_freqs = F)

# ============================ SCENARIO 2C ===================================
# Like 2B, but with more GCs (currently n GCs is chosen upon job submission)
# create_scenario(scenario_directory = '../results/simulations/scenario_2c/',
#                 obs_naive_freqs = obs_naive_freqs,
#                 selected_allele_eligibility_threshold = selected_allele_eligibility_threshold,
#                 n_high_avg_alleles = 5,
#                 s = c(0, 0.5, 1, 2),
#                 sigma_r = 1,
#                 n_high_mutability_alleles = 0,
#                 gamma = 1,
#                 K = 2000,
#                 I_total = 100,
#                 t_imm = 6,
#                 mu_max = 3,
#                 delta = 0.2,
#                 mutation_rate = c(0, 0.01,0.05),
#                 beta = 4,
#                 tmax = 50,
#                 uniform_naive_freqs = F)

# ============================ SCENARIO 2d ===================================
# Like 2,b but with I_total = 200 and different ranges for s and the mutation rate
# create_scenario(scenario_directory = '../results/simulations/scenario_2d/',
#                 obs_naive_freqs = obs_naive_freqs,
#                 selected_allele_eligibility_threshold = selected_allele_eligibility_threshold,
#                 n_high_avg_alleles = 5,
#                 s = c(0.5, 1, 1.5, 2),
#                 sigma_r = 1,
#                 n_high_mutability_alleles = 0,
#                 gamma = 1,
#                 K = 2000,
#                 I_total = 200,
#                 t_imm = 6,
#                 mu_max = 3,
#                 delta = 0.2,
#                 mutation_rate = c(0.01,0.05,0.10),
#                 beta = 4,
#                 tmax = 50,
#                 uniform_naive_freqs = F)


# ============================ HIGH AFFINITY SCENARIO ===================================
# FINAL SCENARIO SHOWING DIFFERENCES IN ALLELE'S AFFINITIES
# create_scenario(scenario_directory = '../results/simulations/high_affinity_scenario/',
#                 obs_naive_freqs = obs_naive_freqs,
#                 selected_allele_eligibility_threshold = selected_allele_eligibility_threshold,
#                 n_high_avg_alleles = 5,
#                 s = c(0.5, 1, 1.5, 2),
#                 sigma_r = 1,
#                 n_high_mutability_alleles = 0,
#                 gamma = 1,
#                 K = 2000,
#                 I_total = 200,
#                 t_imm = 6,
#                 mu_max = 3,
#                 delta = 0.2,
#                 mutation_rate = c(0, 0.01, 0.05),
#                 beta = c(1,2,3,4),
#                 tmax = 50,
#                 uniform_naive_freqs = F)

# ============================ NEUTRAL SCENARIO ===================================
# FINAL NEUTRAL SCENARIO
# create_scenario(scenario_directory = '../results/simulations/neutral_scenario/',
#                 obs_naive_freqs = obs_naive_freqs,
#                 selected_allele_eligibility_threshold = selected_allele_eligibility_threshold,
#                 n_high_avg_alleles = 0,
#                 s = 0,
#                 sigma_r = 1,
#                 n_high_mutability_alleles = 0,
#                 gamma = 1,
#                 K = 2000,
#                 I_total = c(50,100,200),
#                 t_imm = 6,
#                 mu_max = 3,
#                 delta = 0.2,
#                 mutation_rate = c(0, 0.01, 0.05),
#                 beta = c(1,2,3,4),
#                 tmax = 50,
#                 uniform_naive_freqs = F)

# ============================ HIGH MUTABILITY SCENARIO ===================================
create_scenario(scenario_directory = '../results/simulations/high_mutability_scenario/',
                obs_naive_freqs = obs_naive_freqs,
                selected_allele_eligibility_threshold = selected_allele_eligibility_threshold,
                n_high_avg_alleles = 0,
                s = 1,
                sigma_r = 1,
                n_high_mutability_alleles = 5,
                gamma = c(1.5,2,4,6), # WILL ADD 1.1, 1.5, 3
                K = 2000,
                I_total = 200,
                t_imm = 6,
                mu_max = 3,
                delta = 0.2,
                mutation_rate = c(0, 0.01, 0.05),
                beta = c(1,2,3,4),
                tmax = 50,
                uniform_naive_freqs = F)


