library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library('DirichletReg')
library(readr)
theme_set(theme_cowplot())

# Scenarios are specified by the parameters controlling gamma-distributed affinities for 3 classes of alleles
# Each class of allele has a mean alpha and a mean beta parameter.
# Specific alphas and betas for each allele are obtained by adding a normally distributed error around 0 to their class mean values
# Affinities for B cells using each allele are then distributed according to those specific alphas and betas.

generate_affinity_distributions <- function(scenario){
  scenario %>%
    uncount(n_alleles) %>%
    mutate(allele = paste0('V',1:n())) %>%
    select(allele, everything()) %>%
    mutate(alpha = expected_alpha + rnorm(n(), mean = 0, sd = sd_alpha),
           beta = expected_beta + rnorm(n(), mean = 0, sd = sd_beta)) %>%
    select(-expected_alpha, -expected_beta, -sd_alpha, -sd_beta) %>%
    mutate(expected_affinity = alpha / beta)
}

plot_expected_affinity_boxplots <- function(allele_info){
  allele_info %>%
    ggplot(aes(x = allele_type, y = expected_affinity)) +
    geom_boxplot(outlier.alpha = 0) +
    geom_point() +
    xlab('Allele type') +
    ylab('Expected affinity')
}

plot_affinity_distributions <- function(allele_info, log_density = F){
  densities_for_plotting <- full_join(tibble(x = seq(0,15,0.1)),
                                      allele_info, by = character()) %>%
    select(allele, allele_type, alpha, beta, x) %>%
    mutate(density = dgamma(x = x, shape = alpha, rate = beta)) 
  
  pl <- densities_for_plotting %>%
    ggplot(aes(x = x, y = density, group = allele)) +
    geom_line(aes(color = allele_type)) +
    theme(legend.position = 'top') +
    xlab('Affinity') +
    ylab('Density')
  scale_color_discrete(name = 'Allele type')
  
  if(log_density){
    pl <- pl + scale_y_log10()
  }
  return(pl)
}


# First, a neutral scenario all alleles have IDENTICAL affinity distributions
n_low_avg_alleles <- 0 # Alleles with low affinity overall
n_high_avg_alleles <- 0 # Alleles with high average affinity
n_long_tail_alleles <- 20 # Alleles with intermediate average but long tails
sd_alpha <- 0
sd_beta <- 0


neutral_scenario <- tibble(allele_type = c('low_avg', 'high_avg', 'long_tail'),
                           expected_alpha = c(1.5, 10, 2),
                           expected_beta = c(1, 2, 0.5),
                           n_alleles = c(n_low_avg_alleles, n_high_avg_alleles, n_long_tail_alleles),
                           sd_alpha = 0,
                           sd_beta = 0)

allele_info_neutral <- generate_affinity_distributions(scenario = neutral_scenario)

# Showing what the distributions look like
plot_expected_affinity_boxplots(allele_info = allele_info_neutral)
plot_affinity_distributions(allele_info = allele_info_neutral)
plot_affinity_distributions(allele_info = allele_info_neutral, log_density = T)

# In addition to a distribution of affinities, each gene is assigned a naive frequency
# Naive frequencies sampled from a Dirichlet distribution using average freqs. for each rank as alpha parameters
# For development, using a Dirichlet with Poisson alphas + 1
alphas <- rpois(nrow(allele_info_neutral),1) + 1
naive_freqs <- rdirichlet(nrow(allele_info_neutral), alphas)[1,]
hist(naive_freqs)

allele_info_neutral <- allele_info_neutral %>%
  mutate(naive_freq = naive_freqs)

# Export allele information for simulations
write_csv(allele_info_neutral, 
          '../results/simulations/neutral_scenario/replicate_1/allele_info.csv')




#n_low_avg_alleles <- 12 # Alleles with low affinity overall
#n_high_avg_alleles <- 6 # Alleles with high average affinity
#n_long_tail_alleles <- 2 # Alleles with intermediate average but long tails



#scenario <- tibble(allele_type = c('low_avg', 'high_avg', 'long_tail'),
#         expected_alpha = c(1.5, 10, 2),
#         expected_beta = c(1, 2, 0.5),
#         n_alleles = c(n_low_avg_alleles, n_high_avg_alleles, n_long_tail_alleles),
#         sd_alpha = 0,
#         sd_beta = 0)







