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


# First, neutral scenarios where all alleles have IDENTICAL affinity distributions

# ============================ NEUTRAL SCENARIO 1 ===================================

scenario_directory <- '../results/simulations/neutral_scenario_1/'

n_low_avg_alleles <- 80 # Alleles with low affinity overall
n_high_avg_alleles <- 0 # Alleles with high average affinity
n_long_tail_alleles <- 0 # Alleles with intermediate average but long tails
sd_alpha <- 0
sd_beta <- 0

GC_parameters_neutral <- tibble(
  nGCs = 30, # Number of germinal centers in an individual
  K = 10000, # carrying capacity of germinal centers
  mu = 5, # expected number of newly recruited clones arriving at germinal centers per time step
  lambda_max = 1.5, # expected reproductive rate per B cell in an empty germinal center
  mutation_rate = 0.01, # mutation probability per B cell per time step
  mutation_sd = 0.05, # standard deviation for normal distribution of mutational effects (mean 0)
  tmax = 200 # Number of time steps observed
)

neutral_affinity_specs <- tibble(allele_type = c('low_avg', 'high_avg', 'long_tail'),
                           expected_alpha = c(1.5, 10, 2),
                           expected_beta = c(1, 2, 0.5),
                           n_alleles = c(n_low_avg_alleles, n_high_avg_alleles, n_long_tail_alleles),
                           sd_alpha = 0,
                           sd_beta = 0)

allele_info_neutral <- generate_affinity_distributions(affinity_specs = neutral_affinity_specs)
allele_info_neutral <- generate_naive_freqs(allele_info = allele_info_neutral)
hist(allele_info_neutral$naive_freq)


# Showing what the affinity distributions look like
save_plot(paste0(scenario_directory, 'expected_affinities_by_allele_type.pdf'),
          plot_expected_affinity_boxplots(allele_info = allele_info_neutral),
          base_height = 5, base_width = 6)
save_plot(paste0(scenario_directory, 'affinity_distributions.pdf'),
          plot_affinity_distributions(allele_info = allele_info_neutral),
          base_height = 5, base_width = 6)
#plot_affinity_distributions(allele_info = allele_info_neutral, log_density = T)



# Export allele information for simulations
write_csv(allele_info_neutral, 
          paste0(scenario_directory, 'allele_info.csv'))

# Export GC parameters
write_csv(GC_parameters_neutral,
          paste0(scenario_directory, 'GC_parameters.csv'))




n_low_avg_alleles <- 1 # Alleles with low affinity overall
n_high_avg_alleles <- 1 # Alleles with high average affinity
n_long_tail_alleles <- 1 # Alleles with intermediate average but long tails



nonneutral_scenario <- tibble(allele_type = c('low_avg', 'high_avg', 'long_tail'),
         expected_alpha = c(1.5, 10, 2),
         expected_beta = c(1, 2, 0.5),
         n_alleles = c(n_low_avg_alleles, n_high_avg_alleles, n_long_tail_alleles),
         sd_alpha = 0,
         sd_beta = 0)

allele_info_nonneutral <- generate_affinity_distributions(scenario = nonneutral_scenario)

# Showing what the distributions look like
plot_expected_affinity_boxplots(allele_info = allele_info_nonneutral)
plot_affinity_distributions(allele_info = allele_info_nonneutral)
plot_affinity_distributions(allele_info = allele_info_nonneutral, log_density = T)







