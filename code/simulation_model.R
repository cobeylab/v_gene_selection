library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library('DirichletReg')
theme_set(theme_cowplot())

# Functions for simulating the recruitment, competition and evolution of B cell lineages in GCs
# While keeping track of V allele usage and clone's identities

K <- 1000 # carrying capacity of germinal centers
mu <- 1 # expected number of newly recruited clones arriving at germinal centers per time step
lambda_max <- 2 # expected reproductive rate per B cell in an empty germinal center

get_pop_lambda <- function(N, lambda_max, K){
  alpha = log(lambda_max)/K
  return(lambda_max*exp(-alpha*N))
}
# Quick tests
stopifnot(abs(get_pop_lambda(N = K, lambda_max = lambda_max, K = K) - 1) <=1e-7) # lambda = 1 at N = K
stopifnot(get_pop_lambda(N = 0, lambda_max = lambda_max, K = K) == lambda_max) # lambda = lambda_max at N = 0
stopifnot(get_pop_lambda(N = K*1.1, lambda_max = lambda_max, K = K) < 1) # lambda <1 if N>K


poisson_mean_n_immigrants <- 1 # mean number of (Poisson distributed) number of clones clones arriving at GCs in each time step.
mutation_rate <- 0.1
mutation_sd <- 0.2


# Three types of alleles, with gamma distributions for their possible affinites
# Alleles with low affinity overall
n_low_alleles <- 10
low_affinity_alphas <- 1.5 + rnorm(n_low_alleles,0, 0)
low_affinity_betas <- 1 + rnorm(n_low_alleles,0,0)

# Alleles with high average affinity
n_high_avg_alleles <- 5 
high_avg_alphas <- 10 + rnorm(n_high_avg_alleles,0, 0)
high_avg_betas <-  2 + rnorm(n_high_avg_alleles,0,0)

# Alleles with intermediate average but long tails
n_long_tail_alleles <- 2
long_tail_alphas <- 2 + rnorm(n_long_tail_alleles, 0, 0)
long_tail_betas <- 0.5 + rnorm(n_long_tail_alleles,0, 0)

total_n_allelles <- n_low_alleles + n_high_avg_alleles + n_long_tail_alleles

allele_info <- 
  tibble(allele = paste0('V', 1:total_n_allelles),
         allele_type = c(rep('low', n_low_alleles),
                         rep('high average', n_high_avg_alleles),
                         rep('long-tail', n_long_tail_alleles)),
         alpha = c(low_affinity_alphas, high_avg_alphas, long_tail_alphas),
         beta = c(low_affinity_betas, high_avg_betas, long_tail_betas),
         ) %>%
  mutate(expected_affinity = alpha / beta)

allele_info %>%
  ggplot(aes(x = allele_type, y = expected_affinity)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point()
  
  
# Showing what the distributions look like
densities_for_plotting <- full_join(tibble(x = seq(0,15,0.1)),
                                    allele_info, by = character()) %>%
  select(allele, allele_type, alpha, beta, x) %>%
  mutate(density = dgamma(x = x, shape = alpha, rate = beta)) 

densities_for_plotting %>%
  ggplot(aes(x = x, y = density, group = allele)) +
  geom_line(aes(color = allele_type)) +
  theme(legend.position = 'none')

densities_for_plotting %>%
  ggplot(aes(x = x, y = density, group = allele)) +
  geom_line(aes(color = allele_type)) +
  theme(legend.position = 'none') +
  scale_y_log10()

# In addition to a distribution of affinities, each gene is assigned a naive frequency
# Naive frequencies sampled from a Dirichlet distribution using average freqs. for each rank as alpha parameters


# For development, using a Dirichlet with Poisson alphas + 1
alphas <- rpois(total_n_allelles,1) + 1
naive_freqs <- rdirichlet(total_n_allelles, alphas)[1,]
hist(naive_freqs)

allele_info <- allele_info %>%
  mutate(naive_freq = naive_freqs)

# Creates a tibble with one cell per row, representing clones newly arrived at germinal center
# The expected number of new clones (each arriving as a single cell) is mu
# Alleles contribute different numbers of clones to mu depending on their naive freqs. and expected affinities
sample_immigrants <- function(allele_info, mu, clone_numbering_start){
  
  average_naive_affinity <- allele_info %>%
    summarise(A = sum(expected_affinity * naive_freqs)) %>% pull(A)
  
  # Compute number of clones arriving by allele (this is set up so, on average, mu clones will arrive)
  arrivals_by_allele <- allele_info %>%
    rowwise() %>%
    mutate(new_clones = rpois(n = 1, lambda = naive_freq*expected_affinity*mu/average_naive_affinity)) %>%
    ungroup() %>%
    filter(new_clones > 0)
  
  if(nrow(arrivals_by_allele) > 0){
    
    immigrants_tibble <- arrivals_by_allele %>%
      uncount(new_clones) %>%
      rowwise() %>%
      mutate(affinity = rgamma(n = 1, shape = alpha, rate = beta)) %>%
      ungroup() %>%
      select(allele, affinity)
    
    clone_ids <- seq(from = clone_numbering_start, by = 1, length.out = nrow(immigrants_tibble))
    immigrants_tibble$clone_id = clone_ids
    
    immigrants_tibble <- immigrants_tibble %>% select(clone_id, everything())
    
  }else{
    immigrants_tibble <- tibble()
  }
  
  return(immigrants_tibble)
}
# For a test, the command below should return a value close to 100:
#mean(replicate(100, nrow(sample_immigrants(allele_info, mu = 100, 1)), simplify = T))


# Generates GC tibble at time t+1 given state at time t-1 (GC_t) and immigration, growth and mutation parameters
get_GC_tplus1 <- function(allele_info, GC_t, mu, lambda_max, K, mutation_rate, mutation_sd, clone_numbering_start){
  
  # At beginning of time step, potential immigrants arrive
  new_immigrants <- sample_immigrants(allele_info = allele_info, mu = mu, clone_numbering_start = clone_numbering_start)
  
  # The pop size used to determine the expected growth rate per B cell (lambda) is given by the number of
  # cells already in the GC plus the number of new arrivals (each new clone arrives as a single cell)
  GC_pop_size = nrow(GC_t) + nrow(new_immigrants)
  
  if(GC_pop_size > 0){
    GC_t_with_new_immigrants <- bind_rows(GC_t, new_immigrants)
    lambda <- get_pop_lambda(N = GC_pop_size, lambda_max = lambda_max, K = K)  
    
    average_affinity_in_GC <- mean(GC_t_with_new_immigrants$affinity)
    
    # Tibble with one cell per row, corresponding to state of GC at t+1
    GC_tplus1 <- GC_t_with_new_immigrants %>%
      select(-t) %>%
      mutate(descendants = rpois(n = n(), lambda = lambda * affinity/average_affinity_in_GC)) %>%
      uncount(descendants)
    
    # Introduce mutations if at least one cell left descendants
    if(nrow(GC_tplus1) >0){
      GC_tplus1 <- GC_tplus1 %>%
        mutate(has_mutation = rbinom(n = n(), size = 1,prob = mutation_rate),
               mutation_effect = rnorm(n = n(), mean = 0, sd = mutation_sd),
               affinity = affinity + mutation_effect*has_mutation) %>%
        select(-has_mutation, -mutation_effect) %>%
        mutate(affinity = ifelse(affinity < 0, 0, affinity))
    }else{
      GC_tplus1 <- tibble()
    }
    
  }else{
    # This covers the case when the entire GC went extinct and no immigrants arrived at this timestep.
    GC_tplus1 <- tibble()
  }
  
  return(GC_tplus1)

}

# As a test, this mean should be close to the value given by the line below it.
# I.e., a single cell should produce on average get_pop_lambda(1,lambda_max, K)
# in one time step if mu = 0 (no immigrants arrive at beginning of time step)
# mean(replicate(100,
#          nrow(get_GC_tplus1(allele_info, GC_t = tibble(t = 1, clone_id = 1,
#                                         allele = 'A',
#                                         affinity = 1),
#              mu = 0,lambda_max, K, mutation_rate, mutation_sd, clone_numbering_start)),
#          simplify = T))
# get_pop_lambda(1,lambda_max, K)

simulate_GC_dynamics <- function(allele_info, lambda_max, K, mu, mutation_rate, mutation_sd, tmax){
  
  # To initialize GC, sample immigrants until there's at least one clone arriving
  time = 1
  GC_tibble <- tibble()
  while(nrow(GC_tibble) == 0){
    GC_tibble <- sample_immigrants(allele_info = allele_info, mu = mu,
                                   clone_numbering_start = 1) 
  }
  GC_tibble <- GC_tibble %>%
    mutate(t = time) %>% select(t, everything())
  
  # Maximum numeric clone idalready  used (so new clones will always have different ids)
  max_used_clone_id <- max(GC_tibble$clone_id)
  
  
  while(time < tmax){
    
    # Part of tibble corresponding to current time step
    GC_t <- GC_tibble %>% filter(t == time)
    
    # Generate tibble (1 row per cell) for next time step
    
    GC_tplus1 <- get_GC_tplus1(allele_info = allele_info, GC_t = GC_t, mu = mu, lambda_max = lambda_max,
                               K = K, mutation_rate = mutation_rate, mutation_sd = mutation_sd,
                               clone_numbering_start = max_used_clone_id +1)
    
    if(nrow(GC_tplus1) > 0){
      GC_tplus1 <- GC_tplus1 %>%
        mutate(t = time + 1) %>%
        select(t, everything())
    }
    
    max_used_clone_id <- max(max_used_clone_id, GC_tplus1$clone_id)
    
    # Store new time step in master tibble
    GC_tibble <- bind_rows(GC_tibble, GC_tplus1)
    time <- time + 1
  }
  return(GC_tibble)

}

# Quick plots for visual tests:
quick_plotting_function <- function(GC_tibble){
  left_join(GC_tibble, allele_info %>% select(allele, allele_type)) %>%
    group_by(t,clone_id, allele, allele_type) %>% count() %>% 
    ggplot(aes(x = t, y = n, group = clone_id)) +
    geom_line(aes(color = allele_type)) +
    theme(legend.position = 'top')
}

# If lambda_max is < 1, GC populations should die out
# quick_plotting_function(
#   simulate_GC_dynamics(allele_info, lambda_max = 0.9, K, mu, mutation_rate, mutation_sd, tmax)
# ) + ylim(0, K) + geom_hline(aes(yintercept = K), linetype = 2)

# Otherwise the total GC population size should remain around K
# quick_plotting_function(
#   simulate_GC_dynamics(allele_info, lambda_max = 2, K, mu, mutation_rate, mutation_sd, tmax)
# ) + geom_hline(aes(yintercept = K), linetype = 2)

# Super high lambda_max should lead to crazy oscillations due to populations overshooting past K
# quick_plotting_function(
#   simulate_GC_dynamics(allele_info, lambda_max = 25, K, mu, mutation_rate, mutation_sd, tmax)
# ) + geom_hline(aes(yintercept = K), linetype = 2)

GC_tibble <- simulate_GC_dynamics(allele_info, lambda_max, K, mu, mutation_rate, mutation_sd, tmax)



#simulate_repertoire <- function(N_GCs, allele_info, lambda_max, K, mu, mutation_rate, mutation_sd, tmax)


