library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(readr)
theme_set(theme_cowplot())

args <- commandArgs(trailingOnly = T)

allele_info_file_path <- args[1] # allele_info_file_path = '../results/simulations/neutral_scenario/replicate_1/allele_info.csv'
individual_id <- args[2]

allele_info <- read_csv(allele_info_file_path)
output_directory <- paste0(dirname(allele_info_file_path),'/')
output_file_path <- paste0(output_directory,'repertoire_counts_individual_', individual_id, '.csv')


# Functions for simulating the recruitment, competition and evolution of B cell lineages in GCs
# While keeping track of V allele usage and clone's identities

nGCs <- 100 # Number of germinal centers in an individual
tmax <- 500 # Number of timesteps observed
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

mutation_rate <- 0.01
mutation_sd <- 1

# Creates a tibble with one cell per row, representing clones newly arrived at germinal center
# The expected number of new clones (each arriving as a single cell) is mu
# Alleles contribute different numbers of clones to mu depending on their naive freqs. and expected affinities
sample_immigrants <- function(allele_info, mu, clone_numbering_start){
  
  average_naive_affinity <- allele_info %>%
    summarise(A = sum(expected_affinity * naive_freq)) %>% pull(A)
  
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

# As a test, the values produced by the next two lines should be close
# I.e., a single cell should produce on average get_pop_lambda(1,lambda_max, K)
# in one time step if mu = 0 (no immigrants arrive at beginning of time step)
# mean(replicate(100,
#          nrow(get_GC_tplus1(allele_info, GC_t = tibble(t = 1, clone_id = 1,
#                                         allele = 'A',
#                                         affinity = 1),
#              mu = 0,lambda_max, K, mutation_rate, mutation_sd, clone_numbering_start)),
#          simplify = T))
# get_pop_lambda(1,lambda_max, K)


#Combines previous functions to simulate GC dynamics with recruitment, growth, evolution and competition:
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
#quick_plotting_function(
#  simulate_GC_dynamics(allele_info, lambda_max = 0.9, K, mu, mutation_rate, mutation_sd, tmax = 100)
#) + ylim(0, K) + geom_hline(aes(yintercept = K), linetype = 2)

# Otherwise the total GC population size should remain around K
#quick_plotting_function(
# simulate_GC_dynamics(allele_info, lambda_max = 2, K, mu, mutation_rate, mutation_sd, tmax = 100)
#) + geom_hline(aes(yintercept = K), linetype = 2)

# Super high lambda_max should lead to crazy oscillations due to populations overshooting past K
# quick_plotting_function(
#   simulate_GC_dynamics(allele_info, lambda_max = 25, K, mu, mutation_rate, mutation_sd, tmax = 100)
# ) + geom_hline(aes(yintercept = K), linetype = 2)

# For a single realization of a germinal center, compute allele counts
extract_allele_counts <- function(GC_tibble){
  GC_tibble %>%
    group_by(t, allele) %>%
    count() %>%
    ungroup()
}

# Runs realizations for multiple germinal centers, returns combined allele counts across GCs for each time point
simulate_repertoire_allele_counts <- function(nGCs, allele_info, lambda_max, K, mu, mutation_rate, mutation_sd, tmax){
  allele_counts <- replicate(nGCs,
            extract_allele_counts(
              simulate_GC_dynamics(allele_info, lambda_max, K, mu, mutation_rate, mutation_sd, tmax)
            ),
            simplify = F)
  
  # Aggregate counts across GCs for each time point
  allele_counts <- bind_rows(allele_counts, .id = 'GC') %>%
    group_by(t, allele) %>%
    summarise(n = sum(n)) %>%
    ungroup()
  return(allele_counts)
    
}


repertoire_allele_counts <- simulate_repertoire_allele_counts(nGCs = nGCs,
                                                              allele_info = allele_info,
                                                              lambda_max = lambda_max,
                                                              K = K,
                                                              mu = mu,
                                                              mutation_rate = mutation_rate,
                                                              mutation_sd = mutation_sd,
                                                              tmax = tmax)

write_csv(repertoire_allele_counts, file = output_file_path)

#system.time(
#  simulate_repertoire_allele_counts(nGCs = 10, allele_info, lambda_max, K, mu, mutation_rate, mutation_sd, tmax)
#)  


