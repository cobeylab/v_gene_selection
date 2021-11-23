library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library('DirichletReg')
library(readr)
theme_set(theme_cowplot())

source('simulation_functions.R')

# Scenarios are specified by alleles' affinity distributions and allele frequencies and 
# by germinal center parameters (see simulation_functions.R)
# Each class of allele has a mean alpha and a mean beta parameter.
# Specific alphas and betas for each allele are obtained by adding a normally distributed error around 0 to their class mean values
# Affinities for B cells using each allele are then distributed according to those specific alphas and betas.

generate_affinity_specs <- function(n_low_avg_alleles, n_high_avg_alleles, n_long_tail_alleles,
                                    alpha_sd, beta_sd){
  tibble(allele_type = c('low_avg', 'high_avg', 'long_tail'),
         expected_alpha = c(1.5, 10, 2), # FOR NOW THESE ARE FIXED
         expected_beta = c(1, 2, 0.5), # FOR NOW THESE ARE FIXED
         n_alleles = c(n_low_avg_alleles, n_high_avg_alleles, n_long_tail_alleles),
         sd_alpha = 0,
         sd_beta = 0)
  
}


generate_affinity_distributions <- function(affinity_specs){
  affinity_specs %>%
    uncount(n_alleles) %>%
    mutate(allele = paste0('V',1:n())) %>%
    select(allele, everything()) %>%
    mutate(alpha = expected_alpha + rnorm(n(), mean = 0, sd = sd_alpha),
           beta = expected_beta + rnorm(n(), mean = 0, sd = sd_beta)) %>%
    select(-expected_alpha, -expected_beta, -sd_alpha, -sd_beta) %>%
    mutate(expected_affinity = alpha / beta)
}

generate_naive_freqs <- function(allele_info){
  # Allele info is the output of generate_affinity_distributions
  # In addition to a distribution of affinities, each gene is assigned a naive frequency
  # Naive frequencies sampled from a Dirichlet distribution using average freqs. for each rank as alpha parameters
  # For development, using a Dirichlet with Poisson alphas + 1
  alphas <- rpois(nrow(allele_info),1) + 1
  naive_freqs <- rdirichlet(nrow(allele_info), alphas)[1,]
  #hist(naive_freqs)
  
  allele_info <- allele_info %>%
    mutate(naive_freq = naive_freqs)
  return(allele_info)
}

create_scenario <- function(scenario_directory, n_low_avg_alleles, n_high_avg_alleles, n_long_tail_alleles,
                            alpha_sd, beta_sd, nGCs, K, mu, lambda_max, mutation_rate, mutation_sd,
                            tmax, uniform_naive_freqs){

  affinity_specs <- generate_affinity_specs(n_low_avg_alleles = n_low_avg_alleles,
                                            n_high_avg_alleles = n_high_avg_alleles,
                                            n_long_tail_alleles = n_long_tail_alleles,
                                            alpha_sd = alpha_sd, beta_sd = beta_sd)
  
  allele_info <- generate_affinity_distributions(affinity_specs = affinity_specs)
  
  if(uniform_naive_freqs){
    allele_info <- allele_info %>% mutate(naive_freq = 1 / length(unique(allele_info$allele)))
  }else{
    allele_info <- generate_naive_freqs(allele_info = allele_info)
  }
  #hist(allele_info$naive_freq)
  
  GC_parameters <- tibble(
    nGCs = nGCs, # Number of germinal centers in an individual
    K = K, # carrying capacity of germinal centers
    mu = mu, # expected number of newly recruited clones arriving at germinal centers per time step
    lambda_max = lambda_max, # expected reproductive rate per B cell in an empty germinal center
    mutation_rate = mutation_rate, # mutation probability per B cell per time step
    mutation_sd = mutation_sd, # standard deviation for normal distribution of mutational effects (mean 0)
    tmax = tmax, # Number of time steps observed,
    uniform_naive_freqs = uniform_naive_freqs
  )
  
  # Showing what the affinity distributions look like
  save_plot(paste0(scenario_directory, 'expected_affinities_by_allele_type.pdf'),
            plot_expected_affinity_boxplots(allele_info = allele_info),
            base_height = 5, base_width = 6)
  save_plot(paste0(scenario_directory, 'affinity_distributions.pdf'),
            plot_affinity_distributions(allele_info = allele_info),
            base_height = 5, base_width = 6)
  #plot_affinity_distributions(allele_info = allele_info, log_density = T)
  
  
  # Export allele information for simulations
  write_csv(allele_info, 
            paste0(scenario_directory, 'allele_info.csv'))
  
  # Export GC parameters
  write_csv(GC_parameters,
            paste0(scenario_directory, 'GC_parameters.csv'))
  
}

# Neutral scenarios where all alleles have IDENTICAL affinity distributions
# ============================ NEUTRAL SCENARIO 1 ===================================
create_scenario(scenario_directory = '../results/simulations/neutral_scenario_1/',
                n_low_avg_alleles = 80,
                n_high_avg_alleles = 0,
                n_long_tail_alleles = 0,
                alpha_sd = 0, beta_sd = 0,
                nGCs = 100,
                K = 10000, 
                mu = 5, 
                lambda_max = 1.5, 
                mutation_rate = 0.01,
                mutation_sd = 0.5, 
                tmax = 200,
                uniform_naive_freqs = F)

# ============================ NEUTRAL SCENARIO 2 ===================================
# Like neutral scenario 1, but with uniform naive frequencies in all individuals
create_scenario(scenario_directory = '../results/simulations/neutral_scenario_2/',
                n_low_avg_alleles = 80,
                n_high_avg_alleles = 0,
                n_long_tail_alleles = 0,
                alpha_sd = 0, beta_sd = 0,
                nGCs = 100,
                K = 10000, 
                mu = 5, 
                lambda_max = 1.5, 
                mutation_rate = 0.01,
                mutation_sd = 0.5, 
                tmax = 200,
                uniform_naive_freqs = T)

# ============================ NEUTRAL SCENARIO 3 ===================================
# Like neutral scenario 1, but with slower immigration
create_scenario(scenario_directory = '../results/simulations/neutral_scenario_3/',
                n_low_avg_alleles = 80,
                n_high_avg_alleles = 0,
                n_long_tail_alleles = 0,
                alpha_sd = 0, beta_sd = 0,
                nGCs = 100,
                K = 10000, 
                mu = 1, 
                lambda_max = 1.5, 
                mutation_rate = 0.01,
                mutation_sd = 0.5, 
                tmax = 200,
                uniform_naive_freqs = F)

# ============================ NEUTRAL SCENARIO 4 ===================================
# Like scenario 1, but fewer GCs 
create_scenario(scenario_directory = '../results/simulations/neutral_scenario_4/',
                n_low_avg_alleles = 80,
                n_high_avg_alleles = 0,
                n_long_tail_alleles = 0,
                alpha_sd = 0, beta_sd = 0,
                nGCs = 30,
                K = 10000, 
                mu = 5, 
                lambda_max = 1.5, 
                mutation_rate = 0.01,
                mutation_sd = 0.5, 
                tmax = 200, 
                uniform_naive_freqs = F)

# ============================ NEUTRAL SCENARIO 5 ===================================
# Like scenario 1, but with a thousand GCs 
create_scenario(scenario_directory = '../results/simulations/neutral_scenario_5/',
                n_low_avg_alleles = 80,
                n_high_avg_alleles = 0,
                n_long_tail_alleles = 0,
                alpha_sd = 0, beta_sd = 0,
                nGCs = 1000,
                K = 10000, 
                mu = 5, 
                lambda_max = 1.5, 
                mutation_rate = 0.01,
                mutation_sd = 0.5, 
                tmax = 200, 
                uniform_naive_freqs = F)




# ============================ NON-NEUTRAL SCENARIO 1 ===================================
# Like neutral scenario 1, but with 70 low average and 10 high average alleles
create_scenario(scenario_directory = '../results/simulations/non_neutral_scenario_1/',
                n_low_avg_alleles = 70,
                n_high_avg_alleles = 10,
                n_long_tail_alleles = 0,
                alpha_sd = 0, beta_sd = 0,
                nGCs = 100,
                K = 10000, 
                mu = 5, 
                lambda_max = 1.5, 
                mutation_rate = 0.01,
                mutation_sd = 0.5, 
                tmax = 200,
                uniform_naive_freqs = F)






