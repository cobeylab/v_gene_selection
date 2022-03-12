# Functions for simulating the recruitment, competition and evolution of B cell lineages in GCs
# While keeping track of V allele usage and clone's identities
# (Alternative using Gillespie algorithm)
# See simulation_scenarios.R for

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

# For each sampling of clones to seed GCs, sample based on affinity 
# from a sample of this size taken from the naive repertoire with replacement
recruitment_pool_size <- 1000

# Example allele_info object for testing
test_allele_info <- read_csv('test_allele_info.csv')


#tmax: Maximum simulation time.
# K: carrying capacity of germinal centers
# I_total: expected total number of clones that will seed GCs
# t_imm: by at which seeding of GCs will stop.
# mu_max: expected reproductive rate per B cell (day^-1) in an empty germinal center
# delta: per-capita death rate per day
# mutation_rate: mutation probability per B cell per division
# mutation_sd: standard deviation for normal distribution of mutational effects (mean 0)
# allele_info: tibble with alleles' affinity distributions and naive frequencies (see simulation_scenarios.R)
# fixed_initial_affinities: if TRUE, all naive B cells have exactly the expected affinity or their V allele.

# Describes expected growth rate per B cell (mu) as a function of carrying capacity K and total GC pop size (N)
get_pop_mu <- function(N, mu_max, alpha){
  return(mu_max*exp(-alpha*N))
}

# Given lambda_imm, mu_max and K, calculate lambda so that expected population growth is zero with N = K
get_alpha <- function(mu_max, delta, K){
  
  alpha <- (log(mu_max) - log(delta)) / K
  mu_at_K <- get_pop_mu(N = K, mu_max = mu_max, alpha = alpha)
  # Expected population growth is zero at N = K if this equality is true:
  stopifnot(abs(mu_at_K * K - delta*K) < 1e-7)
  names(alpha) = ''
  return(alpha)
}

# Creates a tibble with one row, representing a clone newly arrived at the germinal center
recruit_naive_clone <- function(allele_info, fixed_initial_affinities){

    recruitment_pool <- sample(allele_info$allele, size = recruitment_pool_size,
                               prob = allele_info$naive_freq, replace = T)
    
    recruitment_pool <- left_join(tibble(allele = recruitment_pool),
                                  allele_info %>% select(allele, alpha, beta, expected_affinity))
    
    # Sample affinities of naive B cells
    if(fixed_initial_affinities){
      # If fixed_initial_affinities is TRUE, affinity is exactly the allele expectation
      recruitment_pool <- recruitment_pool %>% mutate(affinity = expected_affinity)
    }else{
      recruitment_pool <- recruitment_pool %>%
        rowwise() %>%
        #Else it's sampled from a Gamma distribution
        mutate(affinity = rgamma(n = 1, shape = alpha, rate = beta)) %>%
        ungroup() 
    }

    
    # Sample one arriving clone from the recruitment pool based on each cell's normalized affinity    
    recruitment_pool <- recruitment_pool %>%
      mutate(sum_affinity = sum(affinity),
             log_normalized_affinity = log(affinity) - log(sum_affinity),
             normalized_affinity = exp(log_normalized_affinity))
    
    recruited_cell_number <- sample(1:nrow(recruitment_pool), size = 1, replace = F,
                                     prob = recruitment_pool$normalized_affinity)
    
    immigrant_tibble <- recruitment_pool %>%
      mutate(cell_number = 1:n()) %>%
      filter(cell_number == recruited_cell_number) %>%
      select(allele, affinity)
    
  
  return(immigrant_tibble)
}

# Get immigration rate:
# Immigration rate declines linearly until reaching 0 on day t_imm.
# We choose t_mm and the total number of immigrants by t_imm, I_total.
# Then compute coefficients accordingly
get_immigration_parameters <- function(I_total, t_imm){
  slope = -2*I_total/(t_imm ^ 2)
  intercept = 2 * I_total / t_imm 
  return(list(slope = slope, intercept = intercept))
}


simulate_GC_dynamics <- function(K, I_total, t_imm, mu_max, delta, mutation_rate, mutation_sd, allele_info, tmax,
                                 observation_times, initial_GC_state = NULL, fixed_initial_affinities){
  
  # Compute alpha (how fast division rate decreases with increasing pop size in GC.)
  alpha <- get_alpha(mu_max = mu_max, delta = delta, K = K)
  
  # Immigration parameters
  imm_pars <- get_immigration_parameters(I_total = I_total, t_imm = t_imm)
  
  time <- 0
  
  clone_numbering_start <- 1
  
  if(!is.null(initial_GC_state)){
    GC_t <- initial_GC_state
    # Master tibble keeping results for each time step
    GC_tibble <- initial_GC_state
    
  }else{
    #If no initial state provided, initialize GC as an empty tibble
    GC_t <- tibble()
    # Master tibble keeping results for each time step
    GC_tibble <- tibble()
    
  }
  
  while(time < tmax){
    
    # Compute immigration rate (immigration stops at time == t_imm)
    if(time <= t_imm){
      lambda_imm = imm_pars$intercept + time * imm_pars$slope
    }else{
      lambda_imm <- 0
    } 
    
    # Current total population size in GC (all clones)
    current_GC_pop <- nrow(GC_t)
  
    # Compute rates of division and death
    lambda_div <- current_GC_pop * get_pop_mu(N = current_GC_pop, mu_max = mu_max, alpha = alpha)
    lambda_death <- current_GC_pop * delta
    
    # Compute total rate of events lambda
    lambda <- lambda_imm + lambda_div + lambda_death
    
    # If lambda is zero (GC population extinct after the end of immigration), stop simulation
    if(lambda == 0){
      time <- tmax
    # Otherwise simulate next event
    }else{
      # Sample time to next event:
      time_to_next_event <- rexp(n = 1, rate = lambda)
      
      # Decide if it's time to export a snapshot
      # For instance, if the next event is the first event after time t=10, export current state as the state at t=10.
      snapshot_interval_current <- findInterval(time, observation_times)
      snapshot_interval_next <- findInterval(time + time_to_next_event, observation_times)
      
      if(snapshot_interval_current != snapshot_interval_next){
        export_snapshot <- T
        snapshot_time <- observation_times[snapshot_interval_next]
        GC_tibble <- bind_rows(GC_tibble, GC_t %>% mutate(t = snapshot_time) %>% select(t, everything()))
        
      }
      
      # Move to the time of the next event
      time <- time + time_to_next_event
      
      # Sample what the next event is
      next_event <- sample(c('immigration', 'division', 'death'), size = 1, 
                           prob = c(lambda_imm, lambda_div, lambda_death))
      
      # Implement next event
      if(next_event == 'immigration'){
        GC_next_t <- bind_rows(GC_t,
                               recruit_naive_clone(allele_info, fixed_initial_affinities) %>%
                                 mutate(clone_id = clone_numbering_start)
        ) %>%
          mutate(t = time)
        clone_numbering_start <- clone_numbering_start +1
      }else{
        if(next_event == 'division'){
          # Probability of cell being sampled to divide is proportional to its affinity
          dividing_cell <- sample(1:nrow(GC_t), size = 1, prob = GC_t$affinity, replace = F)
          # (replace = F is irrelevant, but just to be safe in case we start mutating more than 1 cell at same time)
          
          # Creates new row for the daughter cell, which mutates with some probability
          GC_next_t <- bind_rows(GC_t,
                                 GC_t %>%
                                   slice(dividing_cell) %>%
                                   mutate(cell_mutates = rbinom(n = 1, size = 1, prob = mutation_rate),
                                          affinity = affinity + rnorm(1, mean = 0, sd = mutation_sd * cell_mutates)) %>%
                                   # If mutation takes affinity below zero, set it to zero.
                                   mutate(affinity = ifelse(affinity < 0, 0, affinity)) %>%
                                   select(-cell_mutates)) %>%
            mutate(t = time)
          
        }else{
          stopifnot(next_event == 'death')
          # All cells have the same probability of dying
          dying_cell <- sample(1:nrow(GC_t), size = 1, replace = F) 
          # (replace = F is irrelevant, but just to be safe in case we start more than 1 cells to die at same time)
          GC_next_t <- GC_t %>%
            slice(-dying_cell) %>%
            mutate(t = time)
        }
      }
      # Redefine current state
      GC_t <- GC_next_t
      
    }
    
  }
  
  return(GC_tibble)
}

master_simulation_function <- function(K, I_total, t_imm, mu_max, delta, mutation_rate, mutation_sd, allele_info, tmax,
                           initial_GC_state = NULL, fixed_initial_affinities){
  
  observation_times <- c(1, seq(5, tmax, 5))
  
  simulation <- simulate_GC_dynamics(K = K, I_total = I_total, t_imm = t_imm, mu_max = mu_max,
                                     delta = delta, mutation_rate = mutation_rate,
                                     mutation_sd = mutation_sd, allele_info = allele_info,
                                     tmax = tmax, observation_times = observation_times,
                                     initial_GC_state = initial_GC_state,
                                     fixed_initial_affinities = fixed_initial_affinities)
  
  # Return clone counts per GC over time, together with clone statistics (mean, median and SD affinity)
  clone_stats_per_GC <- simulation %>%
    group_by(t, clone_id, allele) %>%
    summarise(clone_size = n(),
              mean_affinity = mean(affinity),
              median_affinity = median(affinity),
              sd_affinity = sd(affinity)) %>%
    group_by(t) %>%
    mutate(clone_freq = clone_size / sum(clone_size)) %>%
    ungroup()
  
  return(clone_stats_per_GC)
}


# ----- Tests ------


# Quick plots for visual tests:
quick_plotting_function <- function(GC_tibble, allele_info){
  left_join(GC_tibble, allele_info %>% select(allele, allele_type)) %>%
    group_by(t, clone_id, allele, allele_type) %>% count() %>%
    ggplot(aes(x = t, y = n, group = clone_id)) +
    geom_line(aes(color = allele_type)) +
    theme(legend.position = 'top')
}

# Some visual tests:
# If mu_max is < delta, clones should consistently die out
# (this will stop simulations before tmax)
  # quick_plotting_function(
  #  simulate_GC_dynamics(K = 500, I_total = 100, t_imm = 10, mu_max = 0.25, delta = 0.5, mutation_rate = 0.01, mutation_sd = 0.1,
  #                       allele_info = test_allele_info, tmax = 100, observation_times = c(seq(5,100,5)),
  #                       fixed_initial_affinities = F),
  #  allele_info = test_allele_info)

# Otherwise the total GC population size should remain around K
 # simulate_GC_dynamics(K = 300, I_total = 100, t_imm = 6, mu_max = 2, delta = 0.5, mutation_rate = 0.01, mutation_sd = 0.1,
 #                       allele_info = test_allele_info, tmax = 50, observation_times = c(1,seq(5,50,5)),
 #                      fixed_initial_affinities = F) %>%
 #   group_by(t) %>%
 #   count() %>%
 #   ggplot(aes(x = t, y = n)) +
 #   geom_line() +
 #   scale_y_continuous(limits = c(0, NA)) +
 #   geom_hline(yintercept = 300, linetype = 2) +
 #   ylab('Total population in germinal center')
 
 
# An arbitrary initial state for testing purposes. Two clones with the same initial abundance but different affinities
# test_initial_GC <- tibble(allele = c(rep('V1',20), rep('V2', 20)),
#                          affinity = c(rep(1,20), rep(1.2, 20)),
#                          clone_id = c(rep(1,20), rep(2, 20)),
#                          t = 0)

# Without immigration and mutation, clone 2 should consistently win out:
# bind_rows(replicate(n = 5, simulate_GC_dynamics(K = 100, I_total = 0, t_imm = 1, mu_max = 2, delta = 0.5,
#                                                      mutation_rate = 0, mutation_sd = 0,
#                                                      allele_info = test_allele_info, tmax = 50,
#                                                      observation_times = c(seq(5,50,5)),
#                                                      initial_GC_state = test_initial_GC,
#                                                      fixed_initial_affinities = F),
#           simplify = F), .id = 'replicate') %>%
#   group_by(clone_id, replicate, t) %>%
#   count() %>%
#   ggplot(aes(x = t, y = n)) +
#   geom_line(aes(group = clone_id, color = factor(clone_id))) +
#   scale_y_continuous(limits = c(0, NA)) +
#   facet_wrap('replicate') +
#   theme(legend.position = 'top')

# With mutation, clone 1 will sometimes win out by overcoming clone 2's initial advantage

# bind_rows(replicate(n = 5, simulate_GC_dynamics(K = 100, I_total=0, t_imm = 1, mu_max = 2, delta = 0.5,
#                                                 mutation_rate = 0.05, mutation_sd = 4,
#                                                 allele_info = test_allele_info, tmax = 50,
#                                                 observation_times = c(seq(5,50,5)),
#                                                 initial_GC_state = test_initial_GC,
#                                                 fixed_initial_affinities = F),
#                     simplify = F), .id = 'replicate') %>%
#   group_by(clone_id, replicate, t) %>%
#   count() %>%
#   ggplot(aes(x = t, y = n)) +
#   geom_line(aes(group = clone_id, color = factor(clone_id))) +
#   scale_y_continuous(limits = c(0, NA)) +
#   facet_wrap('replicate') +
#   theme(legend.position = 'top')

 

