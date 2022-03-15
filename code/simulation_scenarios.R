library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library('DirichletReg')
library(readr)
theme_set(theme_cowplot())

source('simulation_functions.R')
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


generate_allele_info <- function(obs_naive_freqs, n_high_avg_alleles, n_long_tail_alleles,
                                 selected_allele_eligibility_threshold){
  
  # Assigns simulated affinity distributions to alleles. For each allele, the distribution is the same in all individuals where it occurs 
  # Uses empirical allele sets and naive frequencies. 
  # The number of alleles with "high-average" or "long-tail" affinity distributions is input
  # Only alleles that occur in at least n = selected_allele_eligibility_threshold mice go can go in those categories
  # All alleles not in the "high-average" or "long-tail" categories go in the "low-average category".
  
  n_alleles_in_data <- length(obs_naive_freqs %>% select(v_gene) %>% unique() %>% pull(v_gene))
  stopifnot(n_high_avg_alleles + n_long_tail_alleles < n_alleles_in_data)
  
  
  # Object specifying gamma distribution parameter for each allele category
  allele_types <-   tibble(allele_type = c('low_avg', 'high_avg', 'long_tail'),
                           alpha = c(6, 10, 1), # FOR NOW THESE ARE FIXED
                           beta = c(2, 2, 0.25), # FOR NOW THESE ARE FIXED
                           #alpha = c(8, 10, 2), 
                           #beta = c(2, 2, 0.5)
                           expected_affinity = alpha / beta
                           )
  
  
  #First, assign all alleles in all mice to the low_avg category
  allele_info <- obs_naive_freqs %>%
    select(mouse_id, v_gene, naive_vgene_seq_freq) %>%
    dplyr::rename(naive_freq = naive_vgene_seq_freq) %>%
    mutate(allele_type = 'low_avg')
  
  
  if(n_high_avg_alleles > 0 | n_long_tail_alleles > 0){
    candidate_selected_alleles <- obs_naive_freqs %>% group_by(v_gene) %>% filter(naive_vgene_seq_freq > 0) %>%
      summarise(n_mice = length(unique(mouse_id))) %>%
      filter(n_mice >= selected_allele_eligibility_threshold) %>% pull(v_gene)
    
    high_avg_alleles <- sample(candidate_selected_alleles,size = n_high_avg_alleles, replace = F)
    long_tail_alleles <- sample(candidate_selected_alleles[!(candidate_selected_alleles %in% high_avg_alleles)],
                               size = n_long_tail_alleles, replace = F)
    
    allele_info$allele_type[allele_info$v_gene %in% high_avg_alleles] <- 'high_avg'
    allele_info$allele_type[allele_info$v_gene %in% long_tail_alleles] <- 'long_tail'
    

  }

  # Replace mouse and gene ids with arbitrary integer ids
  
  individual_integer_ids <- obs_naive_freqs %>% select(mouse_id) %>% unique() %>% mutate(individual = 1:n())
  
  allele_integer_ids <- obs_naive_freqs %>% select(v_gene) %>% unique() %>% mutate(allele = 1:n()) %>%
    mutate(allele = paste0('V', allele))
  
  allele_info <- left_join(allele_info, individual_integer_ids) %>% select(-mouse_id)
  allele_info <- left_join(allele_info, allele_integer_ids) %>% select(-v_gene) %>%
    select(individual, allele, allele_type, naive_freq)
  
  # Add gamma-distribution parameters to the allele_info tibble
  allele_info <- left_join(allele_info, 
                           allele_types)
  return(allele_info)
  
}


# Old function using Dirichlet draws to generate naive frequencies (now using empirical frequencies)
# generate_naive_freqs <- function(allele_info){
#   # Allele info is the output of generate_affinity_distributions
#   # In addition to a distribution of affinities, each gene is assigned a naive frequency
#   # Naive frequencies sampled from a Dirichlet distribution using average freqs. for each rank as alpha parameters
#   # For development, using a Dirichlet with Poisson alphas + 1
#   alphas <- rpois(nrow(allele_info),1) + 1
#   naive_freqs <- rdirichlet(nrow(allele_info), alphas)[1,]
#   #hist(naive_freqs)
#   
#   allele_info <- allele_info %>%
#     mutate(naive_freq = naive_freqs)
#   return(allele_info)
# }


create_scenario <- function(scenario_directory, obs_naive_freqs, selected_allele_eligibility_threshold, 
                            n_high_avg_alleles, n_long_tail_alleles, K, I_total, t_imm, mu_max, delta,
                            mutation_rate, mutation_sd, tmax,
                            uniform_naive_freqs, fixed_initial_affinities){
  
  # Create base directory for the scenario.
  dir.create(scenario_directory, showWarnings = F)
  raw_simulations_dir <- paste0(scenario_directory, 'raw_simulation_files/')
  dir.create(raw_simulations_dir, showWarnings = F)
  
  # If multiple parameter combinations sampled within scenario, they all share same allele info (naive freqs, selection)
  allele_info <- generate_allele_info(obs_naive_freqs = obs_naive_freqs,
                                      n_high_avg_alleles = n_high_avg_alleles,
                                      n_long_tail_alleles = n_long_tail_alleles, 
                                      selected_allele_eligibility_threshold = selected_allele_eligibility_threshold)
  
  if(uniform_naive_freqs){
    allele_info <- allele_info %>% 
      group_by(individual) %>%
      mutate(naive_freq = 1 /n())
  }
  
  # Export allele information for simulations
  write_csv(allele_info, paste0(scenario_directory, 'allele_info.csv'))
  
  # Showing what the affinity distributions look like
  save_plot(paste0(scenario_directory, 'expected_affinities_by_allele_type.pdf'),
            plot_expected_affinity_boxplots(allele_info = allele_info),
            base_height = 5, base_width = 6)
  save_plot(paste0(scenario_directory, 'affinity_distributions.pdf'),
            plot_affinity_distributions(allele_info = allele_info),
            base_height = 5, base_width = 6)
  #plot_affinity_distributions(allele_info = allele_info, log_density = T)
  
  
  # Identify parameters with more than one value
  par_combinations <- expand_grid(K = K, I_total = I_total, t_imm = t_imm, mu_max = mu_max, delta = delta,
                      mutation_rate = mutation_rate, mutation_sd = mutation_sd, tmax = tmax,
                      uniform_naive_freqs = uniform_naive_freqs,
                      fixed_initial_affinities = fixed_initial_affinities)
  
    
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



# ============================ NEUTRAL SCENARIO 1 ===================================
create_scenario(scenario_directory = '../results/simulations/neutral_scenario_1/',
                obs_naive_freqs = obs_naive_freqs,
                selected_allele_eligibility_threshold = selected_allele_eligibility_threshold,
                n_high_avg_alleles = 0,
                n_long_tail_alleles = 0,
                K = 2000, 
                I_total = c(50,100,200), 
                t_imm = 6,
                mu_max = 3,
                delta = 0.2, 
                mutation_rate = c(0,0.01, 0.05),
                mutation_sd = 5, 
                tmax = 50,
                uniform_naive_freqs = F,
                fixed_initial_affinities = T)

