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
         expected_alpha = c(6, 10, 1), # FOR NOW THESE ARE FIXED
         expected_beta = c(2, 2, 0.25), # FOR NOW THESE ARE FIXED
         #expected_alpha = c(8, 10, 2), 
         #expected_beta = c(2, 2, 0.5), 
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
                            alpha_sd, beta_sd, nGCs, K, mu, theta, lambda_max, mutation_rate, mutation_sd,
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
    theta = theta,
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

# Neutral scenarios where all alleles have IDENTICAL affinity distributions (default no migration)
# ============================ NEUTRAL SCENARIO 1 ===================================
create_scenario(scenario_directory = '../results/simulations/neutral_scenario_1/',
                n_low_avg_alleles = 80,
                n_high_avg_alleles = 0,
                n_long_tail_alleles = 0,
                alpha_sd = 0, beta_sd = 0,
                nGCs = 10,
                K = 1000, 
                mu = 10, 
                theta = 0,
                lambda_max = 1.5, 
                mutation_rate = 0.01,
                mutation_sd = 1, 
                tmax = 200,
                uniform_naive_freqs = F)

# ============================ NEUTRAL SCENARIO 2 ===================================
# Like neutral scenario 1, slower max growth and WITHOUT MUTATIONS
create_scenario(scenario_directory = '../results/simulations/neutral_scenario_2/',
                n_low_avg_alleles = 80,
                n_high_avg_alleles = 0,
                n_long_tail_alleles = 0,
                alpha_sd = 0, beta_sd = 0,
                nGCs = 10,
                K = 1000, 
                mu = 10, 
                theta = 0,
                lambda_max = 1.05, 
                mutation_rate = 0,
                mutation_sd = 0, 
                tmax = 200,
                uniform_naive_freqs = F)

# ============================ NEUTRAL SCENARIO 3 ===================================
# like neutral scenario 1, but with uniform naive freqs.
create_scenario(scenario_directory = '../results/simulations/neutral_scenario_3/',
                n_low_avg_alleles = 80,
                n_high_avg_alleles = 0,
                n_long_tail_alleles = 0,
                alpha_sd = 0, beta_sd = 0,
                nGCs = 10,
                K = 1000, 
                mu = 10, 
                theta = 0,
                lambda_max = 1.5, 
                mutation_rate = 0.01,
                mutation_sd = 0.33, 
                tmax = 200,
                uniform_naive_freqs = T)

# ============================ NEUTRAL SCENARIO 4 ===================================
# Like neutral scenario 1, but with a ton of seeding (high mu)
create_scenario(scenario_directory = '../results/simulations/neutral_scenario_4/',
                n_low_avg_alleles = 80,
                n_high_avg_alleles = 0,
                n_long_tail_alleles = 0,
                alpha_sd = 0, beta_sd = 0,
                nGCs = 10,
                K = 1000, 
                mu = 100, 
                theta = 0,
                lambda_max = 1.5, 
                mutation_rate = 0.01,
                mutation_sd = 1, 
                tmax = 200,
                uniform_naive_freqs = F)



# ============================ NEUTRAL SCENARIO 5 ===================================
# Like neutral scenario 1, but with uniform naive frequencies in all individuals




# ============================ NON-NEUTRAL SCENARIO 1 ===================================
# Like neutral scenario 1, but with 79 low average and 1 high average alleles, also uniform naive freqs.
create_scenario(scenario_directory = '../results/simulations/non_neutral_scenario_1/',
                n_low_avg_alleles = 79,
                n_high_avg_alleles = 1,
                n_long_tail_alleles = 0,
                alpha_sd = 0, beta_sd = 0,
                nGCs = 10,
                K = 1000, 
                mu = 10, 
                theta = 0,
                lambda_max = 1.5, 
                mutation_rate = 0.01,
                mutation_sd = 1, 
                tmax = 200,
                uniform_naive_freqs = T)

# ============================ NON-NEUTRAL SCENARIO 2 ===================================
# Like non-neutral scenario 1, but with bigger mutation step, slower dynamics
create_scenario(scenario_directory = '../results/simulations/non_neutral_scenario_2/',
                n_low_avg_alleles = 79,
                n_high_avg_alleles = 1,
                n_long_tail_alleles = 0,
                alpha_sd = 0, beta_sd = 0,
                nGCs = 10,
                K = 1000, 
                mu = 10, 
                theta = 0,
                lambda_max = 1.1, 
                mutation_rate = 0.01,
                mutation_sd = 4, 
                tmax = 200,
                uniform_naive_freqs = T)

# ============================ NON-NEUTRAL SCENARIO 3 ===================================
# Like non-neutral scenario 2, but with 5 high affinity alleles and higher K
create_scenario(scenario_directory = '../results/simulations/non_neutral_scenario_3/',
                n_low_avg_alleles = 75,
                n_high_avg_alleles = 5,
                n_long_tail_alleles = 0,
                alpha_sd = 0, beta_sd = 0,
                nGCs = 10,
                K = 5000, 
                mu = 10, 
                theta = 0,
                lambda_max = 1.1, 
                mutation_rate = 0.01,
                mutation_sd = 4, 
                tmax = 200,
                uniform_naive_freqs = T)

# ============================ NON-NEUTRAL SCENARIO 4 ===================================
# Like non-neutral scenario 3, but with higher K
create_scenario(scenario_directory = '../results/simulations/non_neutral_scenario_4/',
                n_low_avg_alleles = 75,
                n_high_avg_alleles = 5,
                n_long_tail_alleles = 0,
                alpha_sd = 0, beta_sd = 0,
                nGCs = 10,
                K = 10000, 
                mu = 10, 
                theta = 0,
                lambda_max = 1.1, 
                mutation_rate = 0.01,
                mutation_sd = 4, 
                tmax = 200,
                uniform_naive_freqs = T)

# ============================ NON-NEUTRAL SCENARIO 5 ===================================
# Like non-neutral scenario 3, but with migration between GCs
create_scenario(scenario_directory = '../results/simulations/non_neutral_scenario_5/',
                n_low_avg_alleles = 75,
                n_high_avg_alleles = 5,
                n_long_tail_alleles = 0,
                alpha_sd = 0, beta_sd = 0,
                nGCs = 10,
                K = 5000, 
                mu = 10, 
                theta = 0.0001,
                lambda_max = 1.1, 
                mutation_rate = 0.01,
                mutation_sd = 4, 
                tmax = 200,
                uniform_naive_freqs = T)

# ============================ NON-NEUTRAL SCENARIO 6 ===================================
# Like non-neutral scenario 5, but with non-uniform naive frequencies
create_scenario(scenario_directory = '../results/simulations/non_neutral_scenario_6/',
                n_low_avg_alleles = 75,
                n_high_avg_alleles = 5,
                n_long_tail_alleles = 0,
                alpha_sd = 0, beta_sd = 0,
                nGCs = 10,
                K = 5000, 
                mu = 10, 
                theta = 0.0001,
                lambda_max = 1.1, 
                mutation_rate = 0.01,
                mutation_sd = 4, 
                tmax = 200,
                uniform_naive_freqs = F)

# ============================ NON-NEUTRAL SCENARIO 7 ===================================
# Like non-neutral scenario 3, but with ....???
# create_scenario(scenario_directory = '../results/simulations/non_neutral_scenario_7/',
#                 n_low_avg_alleles = 75,
#                 n_high_avg_alleles = 5,
#                 n_long_tail_alleles = 0,
#                 alpha_sd = 0, beta_sd = 0,
#                 nGCs = 10,
#                 K = 5000, 
#                 mu = 10, 
#                 theta = 0,
#                 lambda_max = 1.1, 
#                 mutation_rate = 0.01,
#                 mutation_sd = 4, 
#                 tmax = 200,
#                 uniform_naive_freqs = T)










