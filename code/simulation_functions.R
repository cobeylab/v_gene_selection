# Functions for simulating the recruitment, competition and evolution of B cell lineages in GCs
# While keeping track of V allele usage and clone's identities
# See simulation_scenarios.R for

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

# For each sampling of clones to seed GCs, sample based on affinity 
# from a sample of this size taken from the naive repertoire with replacement
recruitment_pool_size <- 1000


# Paths for easy access during development:
# allele_info <- read_csv('../results/simulations/neutral_scenario_1/allele_info.csv')
# GC_parameters <-  read_csv('../results/simulations/neutral_scenario_1/GC_parameters.csv')

# nGCs: Number of germinal centers in an individual
#tmax: Number of time steps observed
# K: carrying capacity of germinal centers
# mu: expected number of newly recruited clones arriving at germinal centers per time step
# theta: Probability that a B cell migrates to another germinal center per time step
# lambda_max: expected reproductive rate per B cell in an empty germinal center
# mutation_rate: mutation probability per B cell per time step
# mutation_sd: standard deviation for normal distribution of mutational effects (mean 0)
# allele_info: tibble with alleles' affinity distributions and naive frequencies (see simulation_scenarios.R)

# (Load pre-generated files with allele_info and the other parameters to run the tests here.)


# Describes expected growth rate per B cell (lambda) as a function of carrying capacity K and total GC pop size (N)
get_pop_lambda <- function(N, lambda_max, K){
  alpha = log(lambda_max)/K
  return(lambda_max*exp(-alpha*N))
}
# Quick tests
stopifnot(abs(get_pop_lambda(N = 100, lambda_max = 2, K = 100) - 1) <=1e-7) # lambda = 1 at N = K
stopifnot(get_pop_lambda(N = 0, lambda_max = 2, K = 100) == 2) # lambda = lambda_max at N = 0
stopifnot(get_pop_lambda(N = 200, lambda_max = 2, K = 100) < 1) # lambda <1 if N>K


# Creates a tibble with one cell per row, representing clones newly arrived at germinal center
# The expected number of new clones (each arriving as a single cell) is mu
# Alleles contribute different numbers of clones to mu depending on their naive freqs. and expected affinities
recruit_naive_clones_single_GC <- function(allele_info, mu, clone_numbering_start){
  
  n_arrivals <- rpois(n = 1, lambda = mu)

  if(n_arrivals > 0){
    recruitment_pool <- sample(allele_info$allele, size = recruitment_pool_size,
                               prob = allele_info$naive_freq, replace = T)
    
    recruitment_pool <- left_join(tibble(allele = recruitment_pool),
                                  allele_info %>% select(allele, alpha, beta, expected_affinity)) %>%
      rowwise() %>%
      mutate(affinity = rgamma(n = 1, shape = alpha, rate = beta)) %>%
      ungroup() %>%
      mutate(sum_affinity = sum(affinity),
             log_normalized_affinity = log(affinity) - log(sum_affinity),
             normalized_affinity = exp(log_normalized_affinity))
    
    recruited_cell_numbers <- sample(1:nrow(recruitment_pool), size = n_arrivals, replace = F,
                                     prob = recruitment_pool$normalized_affinity)
    
    # Sample n_arrivals from the recruitment pool based on each cell's normalized affinity
    immigrants_tibble <- recruitment_pool %>%
      mutate(cell_number = 1:n()) %>%
      filter(cell_number %in% recruited_cell_numbers) %>%
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
#mean(replicate(100, nrow(recruit_naive_clones_single_GC(allele_info, mu = 100, 1)), simplify = T))

# Applies the previous function to a collection of GCs, with their specificied clone numbering starts 
recruit_naive_clones <- function(allele_info, mu, clone_numbering_start){
  new_recruits <- mapply(
    FUN = function(GC, clone_numbering_start, allele_info, mu){
      recruit_naive_clones_single_GC(allele_info = allele_info, mu = mu, clone_numbering_start = clone_numbering_start) %>%
        mutate(GC = GC) %>% select(GC, everything())},
    GC = clone_numbering_start$GC,
    clone_numbering_start = clone_numbering_start$clone_numbering_start,
    MoreArgs = list(allele_info = allele_info, mu = mu),
    SIMPLIFY = F
  )
  new_recruits <- bind_rows(new_recruits)
  
  # Clones are labelled based on a number and the GC where they first appeared (e.g. GC1_1 and GC2_1 are diff.)
  if(nrow(new_recruits) >0 ){
    new_recruits <- new_recruits %>%
      mutate(clone_id = paste0('GC',GC,'_' ,clone_id))
  }
  return(new_recruits)
}

# Generates GC tibble at time t+1 given state at time t-1 (GC_t) and immigration, growth and mutation parameters
get_GC_tplus1 <- function(allele_info, GC_t, mu, theta, lambda_max, K, mutation_rate, mutation_sd, clone_numbering_start){
  
  # theta gives the true migration probability (i.e. prob that a cell will go to *ANOTHER* GC)
  # For convenience, we define theta_adj, which gives the probability a cell will be randomly assigned to a GC which may be its current one. 
  # theta_adj is chosen such that the probability of going to ANOTHER gc is the chosen value theta
  number_of_GCs <- length(clone_numbering_start$GC)
  
  if(theta != 0){
    stopifnot(number_of_GCs > 1) # Stops if a non-zero migration rate is specified but there's only one GC.
    theta_adj <- theta * number_of_GCs / (number_of_GCs - 1)
    # i.e, theta_adj * (1 - 1/nGCs) = theta
  }else{
    theta_adj <- theta # Only happens if theta = 0
  }
  
  
  # At beginning of time step, new clones from the naive repertoire potentially arrive in each GC
  # clone_numbering_start is used to pass the numbering of clones in each GC

  new_recruits <- recruit_naive_clones(allele_info = allele_info, mu = mu, clone_numbering_start = clone_numbering_start)
  
  # Skip if GCs are currently empty and no new recruits arrived. Otherwise proceed.
  if(nrow(GC_t) >0 | nrow(new_recruits) >0){
    
    # After potential new recruits arrive, there's a chance for cells to migrate between GCs:
    # (New recruits won't migrate to another GC in the same time step, they're not in this tibble yet)
    if(nrow(GC_t) > 0){
      GC_t_after_cross_GC_migration <- GC_t %>%
        mutate(cell_migrates = (rbinom(n(), size = 1, prob = theta_adj)) == 1,
               new_GC = sample(unique(GC), size = n(), replace = T)) %>%
        mutate(GC = ifelse(cell_migrates, new_GC, GC)) %>%
        select(-new_GC, -cell_migrates)
    }else{
      GC_t_after_cross_GC_migration <- GC_t
    }
      
    GC_t_with_new_recruits <- bind_rows(GC_t_after_cross_GC_migration, new_recruits)
      
    
    # The pop size used to determine the expected growth rate per B cell (lambda) is given by the number of
    # cells already in the GC (including migrants from other GCs)
    # plus the number of new arrivals (each new clone arrives as a single cell)
    
    lambdas <- GC_t_with_new_recruits %>%
      arrange(GC) %>%
      group_by(GC) %>% 
      summarise(N_in_GC = n(),
             average_affinity_in_GC = mean(affinity)) %>%
      rowwise() %>%
      mutate(lambda_in_GC = get_pop_lambda(N = N_in_GC, lambda_max = lambda_max, K = K)) %>%
      ungroup()
    
    # Tibble with one cell per row, corresponding to state of GC at t+1:
    
    GC_tplus1 <-  left_join(GC_t_with_new_recruits, lambdas, by = 'GC') %>%
      arrange(GC) %>%
      select(-any_of(c('t'))) %>%
      mutate(descendants = rpois(n = n(), lambda = lambda_in_GC * affinity/average_affinity_in_GC)) %>%
      select(-average_affinity_in_GC, -lambda_in_GC, -N_in_GC) %>%
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


#Combines previous functions to simulate GC dynamics with recruitment, growth, evolution and competition:
simulate_GC_dynamics <- function(nGCs, allele_info, lambda_max, K, mu, theta, mutation_rate, mutation_sd, tmax){
  
  time = 0

  # Internal function used to keep track of clone numbering as newly recruited clones arrive
  update_clone_numbering_start <- function(clone_numbering_start,  GC_tibble){
    
    if(nrow(GC_tibble) == 0){
      new_numbering_start <- clone_numbering_start 
    }else{
      max_id_n_in_each_GC <- GC_tibble %>%
        select(-GC) %>% # Numbering considers GCs where clones first arrive, not their current GC w possible migration
        separate(clone_id, sep = '_', into = c('GC','number_within_GC')) %>%
        mutate(GC = as.integer(str_remove(GC,'GC')),
               number_within_GC = as.integer(number_within_GC)) %>%
        group_by(GC) %>%
        summarise(max_n_within_GC = max(number_within_GC))
      # Each GC has a counter identifying clones that start there.
      # If there's been GC1_1 and GC1_2 and GC1_3 in GC1, the next clone arriving from the naive repertoire
      # will be GC1_4, even if, say, GC1_3 went extinct in that GC.
      # Clones will retain their labels based on original GC even when they migrate between GCs.
      
      new_numbering_start <- left_join(clone_numbering_start %>% rename(current_numbering_start = clone_numbering_start),
                max_id_n_in_each_GC) %>%
        replace_na(list(max_n_within_GC = 0)) %>%
        mutate(
          clone_numbering_start = case_when(
            max_n_within_GC >= current_numbering_start ~ max_n_within_GC + 1,
            T ~ current_numbering_start
          )
        ) %>%
        select(GC, clone_numbering_start)
      
    }
    return(new_numbering_start)

  }
  
  clone_numbering_start <- tibble(GC = 1:nGCs, clone_numbering_start = 1)
  
  # Master tibble keeping results for each time step
  GC_tibble <- tibble()
  
  # Initialize GC state as an empty tibble 
  GC_t <- tibble()

  while(time < tmax){
    time = time + 1
    
    # Generate tibble (1 row per cell) for next time step
    
    GC_tplus1 <- get_GC_tplus1(allele_info = allele_info, GC_t = GC_t, mu = mu, theta = theta, lambda_max = lambda_max,
                               K = K, mutation_rate = mutation_rate, mutation_sd = mutation_sd,
                               clone_numbering_start = clone_numbering_start)
    
    if(nrow(GC_tplus1) > 0){
      GC_tplus1 <- GC_tplus1 %>%
        mutate(t = time) %>%
        select(t, everything())
    }

    # Store new time step in master tibble every 10 steps
    if(time == 1 | (time%%10 == 0)){
      GC_tibble <- bind_rows(GC_tibble, GC_tplus1)
    }

    # Update starting id number for clones in each GC
    clone_numbering_start <- update_clone_numbering_start(clone_numbering_start, GC_tplus1)
    
    # Set GC_tplus1 as the new GC_t
    GC_t <- GC_tplus1
  }
  return(GC_tibble)
  
}

# Quick plots for visual tests:
quick_plotting_function <- function(GC_tibble){
  left_join(GC_tibble, allele_info %>% select(allele, allele_type)) %>%
    group_by(t, GC, clone_id, allele, allele_type) %>% count() %>% 
    ggplot(aes(x = t, y = n, group = clone_id)) +
    geom_line(aes(color = allele_type)) +
    theme(legend.position = 'top') +
    facet_wrap('GC')
}

# If lambda_max is < 1, GC populations should die out
# quick_plotting_function(
#   simulate_GC_dynamics(nGCs = 2, allele_info, lambda_max = 0.9, K, mu, theta, mutation_rate, mutation_sd,
#                        tmax = 100)) + ylim(0, K) + geom_hline(aes(yintercept = K), linetype = 2)

# Otherwise the total GC population size should remain around K
# quick_plotting_function(
#  simulate_GC_dynamics(nGCs = 2, allele_info, lambda_max = 2, K, mu, theta, mutation_rate, mutation_sd, tmax = 100)
# ) + geom_hline(aes(yintercept = K), linetype = 2)

# Super high lambda_max should lead to crazy oscillations due to populations overshooting past K
# quick_plotting_function(
#   simulate_GC_dynamics(nGCs = 2, allele_info, lambda_max = 25, K, mu, theta, mutation_rate, mutation_sd, tmax = 100)
# ) + geom_hline(aes(yintercept = K), linetype = 2)


# Runs realizations for multiple germinal centers in an individual, returns combined allele counts across GCs for each time point
# Output:
# - allele counts across multiple GCs within an individual
# - Clone and allele diversity statistics for each GC
# - Also exports first GC trajectory in detailed format (abundance of each clone for each time point)

run_simulation <- function(nGCs, allele_info, lambda_max, K, mu, theta, mutation_rate, mutation_sd, tmax){
  
  
  simulation <- simulate_GC_dynamics(nGCs = nGCs, allele_info = allele_info, lambda_max = lambda_max,
                                     K = K, mu = mu, theta = theta, mutation_rate = mutation_rate,
                                     mutation_sd = mutation_sd, tmax = tmax)
  
  # Repertoire allele counts (aggregated across individual GCs)
  allele_counts <- simulation %>%
    group_by(t, allele) %>%
    count() %>%
    ungroup()
  
  # Export some statistics (but not full population structure) for each GC at each time step
  clone_diversity_by_GC <- simulation %>%
    group_by(t, GC, clone_id) %>%
    count() %>%
    group_by(t, GC) %>%
    mutate(clone_freq = n / sum(n)) %>%
    summarise(n_clones = length(unique(clone_id)),
              total_GC_pop = sum(n),
              fraction_biggest_clone = max(clone_freq),
              clone_diversity = 1 - sum(clone_freq^2)) %>%
    ungroup()
  
  allele_diversity_by_GC <- simulation %>%
    group_by(t, GC, allele) %>%
    count() %>%
    group_by(t, GC) %>%
    mutate(allele_freq = n / sum(n)) %>%
    summarise(n_alleles = length(unique(allele)),
              fraction_most_common_allele = max(allele_freq),
              allele_diversity = 1 - sum(allele_freq^2)) %>%
    ungroup()
  
  GC_statistics <- left_join(clone_diversity_by_GC, allele_diversity_by_GC, by = c('t','GC'))
  
  # Example GC
  # Export the detailed trajectory (size of each clone at each time step) for an example GC
  example_GC <- sample(1:nGCs, size = 1, replace = F)
  example_GC_trajectory <- simulation %>% filter(GC == example_GC)
  
  return(list(allele_counts = allele_counts, GC_statistics = GC_statistics,
              example_GC_trajectory = example_GC_trajectory))
  
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
    geom_line(aes(color = allele_type), size = 1.5) +
    theme(legend.position = 'top') +
    xlab('Affinity') +
    ylab('Density')
  scale_color_discrete(name = 'Allele type')
  
  if(log_density){
    pl <- pl + scale_y_log10()
  }
  return(pl)
}


