# Functions for simulating the recruitment, competition and evolution of B cell lineages in GCs
# While keeping track of V allele usage and clone's identities
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(cowplot)
library(truncnorm)
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
# mutation_sd: standard deviation for normal distribution of mutational effects (mean 0) (given by beta * sigma_r)
# allele_info: tibble with alleles' affinity distributions and naive frequencies (see simulation_scenarios.R)
# baseline_mean: baseline mean affinity (i.e., average for 'low affinity' alleles.)
# s: increase in mean affinity associated with high-affinity alleles, relative to sigma_r
# sigma_r: baseline variation in naive affinity due to VDJ recombination
# gamma relative mutability of high mutability alleles (compared to 1, for regular alleles)


# Describes expected growth rate per B cell (mu) as a function of carrying capacity K and total GC pop size (N)
get_pop_mu <- function(N, mu_max, alpha){
  return(mu_max*exp(-alpha*N))
}

# Given mu_max and K, calculate lambda so that expected population growth is zero with N = K
get_alpha <- function(mu_max, delta, K){
  
  alpha <- (log(mu_max) - log(delta)) / K
  mu_at_K <- get_pop_mu(N = K, mu_max = mu_max, alpha = alpha)
  # Expected population growth is zero at N = K if this equality is true:
  stopifnot(abs(mu_at_K * K - delta*K) < 1e-7)
  names(alpha) = ''
  return(alpha)
}

# Creates a tibble with one row, representing a clone newly arrived at the germinal center
recruit_naive_clone <- function(allele_info){

    recruitment_pool <- sample(allele_info$allele, size = recruitment_pool_size,
                               prob = allele_info$naive_freq, replace = T)
    
    recruitment_pool <- left_join(tibble(allele = recruitment_pool),
                                  allele_info %>% select(allele, mean_affinity, sd_affinity, relative_mutability))
  
    recruitment_pool <- recruitment_pool %>%
      rowwise() %>%
      #Else it's sampled from a Truncated normal distribution
      mutate(affinity = max(0,rtruncnorm(n = 1, a = 0, b = Inf, mean = mean_affinity, sd = sd_affinity))) %>%
      ungroup() 
    
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
      select(allele, affinity, relative_mutability)
    
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
                                 observation_times, initial_GC_state = NULL){
  
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
                               recruit_naive_clone(allele_info) %>%
                                 mutate(clone_id = clone_numbering_start)
        ) %>%
          mutate(t = time)
        clone_numbering_start <- clone_numbering_start +1
      }else{
        if(next_event == 'division'){
          # Probability of cell being sampled to divide is proportional to its affinity
          dividing_cell <- sample(1:nrow(GC_t), size = 1, prob = GC_t$affinity, replace = F)
          stopifnot(length(dividing_cell) == 1)
          # (replace = F is irrelevant, but just to be safe in case we start mutating more than 1 cell at same time)
          
          # Creates new row for the daughter cell, which mutates with some probability
          daughter_cell <- GC_t[dividing_cell,]
          
          cell_mutates <- rbinom(n = 1, size = 1, prob = mutation_rate*daughter_cell$relative_mutability)
          if(cell_mutates == T){
            # If mutation takes affinity below zero, set it to very small positive value.
            new_affinity <- max(1e-10, daughter_cell$affinity + rnorm(1, mean = 0, sd = mutation_sd))
            daughter_cell$affinity <- new_affinity
          }
           
          GC_next_t <- bind_rows(GC_t, daughter_cell) %>% 
            mutate(t = time)
          
        }else{
          stopifnot(next_event == 'death')
          
          # All cells have the same probability of dying
          # (replace = F is irrelevant, but just to be safe in case we start more than 1 cells to die at same time)
          dying_cell <- sample(1:nrow(GC_t), size = 1, replace = F) 
          stopifnot(length(dying_cell) == 1)
          
          GC_next_t <- GC_t[-dying_cell,] %>%
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
                           initial_GC_state = NULL){
  
  observation_times <- c(1, seq(5, tmax, 5))
  
  simulation <- simulate_GC_dynamics(K = K, I_total = I_total, t_imm = t_imm, mu_max = mu_max,
                                     delta = delta, mutation_rate = mutation_rate,
                                     mutation_sd = mutation_sd, allele_info = allele_info,
                                     tmax = tmax, observation_times = observation_times,
                                     initial_GC_state = initial_GC_state)
  
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

# Adds info on affinity distributions and mutability to tibble of allele info
assign_allele_properties <- function(allele_info, baseline_mean, s, sigma_r, gamma){
  
  allele_types_affinity <- tibble(allele_type_affinity = c('low_avg', 'high_avg'),
                                  mean_affinity = c(baseline_mean, baseline_mean + s*sigma_r),
                                  sd_affinity = c(sigma_r, sigma_r))
  
  allele_types_mutability <- tibble(allele_type_mutability = c('low_mut', 'high_mut'),
                                    relative_mutability = c(1, gamma))
  
  allele_info <- left_join(allele_info, allele_types_affinity)
  
  allele_info <- left_join(allele_info, allele_types_mutability)
  
  return(allele_info)

}


# ========= Functions for summarizing simulations ===============
# simulations: clone counts per time point per GC per individual
# variable_pars: if simulations object includes simulations with multiple par. combinations, does grouping by variable parameters

compute_allele_freqs_per_GC <- function(simulations, variable_pars){
  simulations %>%
    group_by(across(c(any_of(variable_pars),'individual', 'base_individual', 't', 'GC', 'allele'))) %>%
    summarise(n_cells = sum(clone_size)) %>%
    ungroup() %>%
    group_by(across(c(any_of(variable_pars),'individual', 'base_individual','t', 'GC'))) %>%
    mutate(allele_freq = n_cells / sum(n_cells)) %>%
    ungroup()
}

compute_allele_diversity_per_GC <- function(allele_freqs_by_GC, variable_pars){
  allele_freqs_by_GC %>%
    group_by(across(c(any_of(variable_pars),'individual', 'base_individual', 't', 'GC'))) %>%
    summarise(n_alleles = length(unique(allele)),
              fraction_most_common_allele = max(allele_freq),
              allele_diversity = 1 - sum(allele_freq^2)) %>%
    ungroup()
}

compute_clone_diversity_per_GC <- function(simulations, variable_pars){
  simulations %>%
    group_by(across(c(any_of(variable_pars),'individual', 'base_individual', 't', 'GC'))) %>%
    summarise(n_clones = length(unique(clone_id)),
              total_GC_pop = sum(clone_size),
              fraction_biggest_clone = max(clone_freq),
              clone_diversity = 1 - sum(clone_freq^2)) %>%
    ungroup()
}

sample_random_GCs <- function(allele_freqs_by_GC, sample_size){
  GC_indices <- unique(allele_freqs_by_GC$GC)
  stopifnot(sample_size <= length(GC_indices))
  
  sampled_indices <- sample(GC_indices, size = sample_size)
  
  subsampled_tibble <- allele_freqs_by_GC %>%
    filter(GC %in% sampled_indices) %>%
    mutate(nGCs = sample_size) %>%
    select(nGCs, everything())
  return(subsampled_tibble)
}

# Computes repertoire-wide allele freqs per individual per time from GC-specific freqs.
compute_repertoire_allele_freqs <- function(allele_freqs_by_GC, nGC_values, allele_info, variable_pars){
  
  # For each nGCs value, create a subsample of allele_freqs_by_GC
  extended_freqs_tibble <- lapply(as.list(nGC_values), FUN = sample_random_GCs,
         allele_freqs_by_GC = allele_freqs_by_GC)
  
  extended_freqs_tibble <- bind_rows(extended_freqs_tibble)
  
  repertoire_allele_freqs <- extended_freqs_tibble %>%
    group_by(across(c(any_of(variable_pars),'nGCs','individual', 'base_individual', 't', 'allele'))) %>%
    summarise(n_cells = sum(n_cells)) %>%
    ungroup() %>%
    group_by(across(c(any_of(variable_pars), 'nGCs', 'individual', 'base_individual', 't'))) %>%
    mutate(experienced_freq = n_cells / sum(n_cells)) %>%
    ungroup()
  
  # Adjust so zeros are explicitly represented
  repertoire_allele_freqs <- complete_repertoire_allele_freqs(repertoire_allele_freqs = repertoire_allele_freqs,
                                                              allele_info = allele_info, variable_pars = variable_pars)
  
  # Add total number of cells, rank alleles
  repertoire_allele_freqs <- repertoire_allele_freqs %>%
    group_by(across(c(any_of(variable_pars), 'nGCs','t', 'individual'))) %>%
    mutate(total_time_point_cells = sum(n_cells),
           allele_rank = rank(-experienced_freq, ties.method = 'average')) %>%
    ungroup()
  
  # Add allele affinities /types, add naive allele frequencies and compute experienced-to-naive ratios,
  repertoire_allele_freqs <- left_join(repertoire_allele_freqs, allele_info %>%
                                         select(base_individual, allele, naive_freq, allele_type_affinity, allele_type_mutability)) %>%
    mutate(freq_ratio_log = log(experienced_freq) - log(naive_freq),
           freq_ratio = exp(freq_ratio_log)) %>% select(-freq_ratio_log)
  
  # Order alleles as a factor
  allele_order <- paste0('V',sort(as.integer(str_remove(unique(repertoire_allele_freqs$allele), 'V'))))
  repertoire_allele_freqs <- repertoire_allele_freqs %>%
    mutate(allele = factor(allele, levels = allele_order))
  
  # Check allele frequencies sum to 1
  freq_sums <- repertoire_allele_freqs %>%
    group_by(across(c(any_of(variable_pars), 'nGCs','t','individual','base_individual'))) %>%
    summarise(freq_sum = sum(experienced_freq)) %>% ungroup() %>% select(freq_sum) %>%
    unique() %>% pull(freq_sum)
  
  stopifnot(all(abs(freq_sums-1) < 1e-6))
  
  return(repertoire_allele_freqs)

}

# Used by compute_repertoire_allele_freqs to have zeros explicitly represented 
# (uses allele_info to find the full set of alleles present in the naive repertoire for each individual)
complete_repertoire_allele_freqs <- function(repertoire_allele_freqs, allele_info, variable_pars){
  
  # All alleles of an individual explicitly represented at all time points, for all nGCs
  complete_scaffold <-  expand_grid(individual = unique(repertoire_allele_freqs$individual),
                                              t = unique(repertoire_allele_freqs$t),
                                              nGCs = unique(repertoire_allele_freqs$nGCs))
  complete_scaffold <- left_join(complete_scaffold, repertoire_allele_freqs %>% select(individual, base_individual) %>% unique())
  
  
  
  complete_scaffold <- left_join(complete_scaffold,
                                 allele_info %>% select(base_individual, allele))
  
  if(length(variable_pars) > 0){
    complete_scaffold <- expand_grid(repertoire_allele_freqs %>% select(any_of(variable_pars)) %>% unique(),
                                     complete_scaffold)
  }
  
  completed_tibble <- left_join(complete_scaffold, repertoire_allele_freqs) %>%
    replace_na(list(n_cells = 0, experienced_freq = 0))
  
  return(completed_tibble)
}

compute_repertoire_allele_diversity <- function(repertoire_allele_freqs, variable_pars){
  repertoire_allele_freqs %>%
    group_by(across(c(any_of(variable_pars), 'nGCs', 't', 'individual', 'base_individual', 'total_time_point_cells'))) %>%
    summarise(repertoire_allele_diversity = 1 - sum(experienced_freq^2),
              n_alleles_in_experienced_repertoire = sum(experienced_freq>0)) %>%
    ungroup()
}

# Not worth trying to use get_pairwise_freqs function written for obs data (too many other variables/groupings) 
# Best to write a function specific for simulations even though it will look similar
get_pairwise_sim_freqs <- function(repertoire_allele_freqs, variable_pars){
  unique_pairs <- repertoire_allele_freqs %>% select(individual) %>% unique() %>%
    dplyr::rename(ind_i = individual) %>%
    mutate(ind_j = ind_i) %>%
    complete(ind_i, ind_j) %>%
    rowwise() %>%
    mutate(pair = paste0(sort(c(ind_i, ind_j)), collapse = ';')) %>%
    ungroup() %>%
    filter(ind_i != ind_j) %>%
    select(pair) %>%
    unique() %>% pull(pair)
  
  # Will add these back later to fill full-join missing values in rows where a gene is missing from one individual
  total_time_point_cells <- repertoire_allele_freqs %>% select(any_of(variable_pars),
                                                               individual, nGCs, t, total_time_point_cells) %>% unique()
  
  internal_function <- function(pair, repertoire_allele_freqs){
    
    ind_specific_vars <- c('individual','n_cells', 'total_time_point_cells', 'experienced_freq', 'naive_freq', 'freq_ratio',
                           'allele_rank')
    
    individual_ids <- str_split(pair,';')[[1]]
    
    ind_i_values <- repertoire_allele_freqs %>% filter(individual == individual_ids[1])  %>%
      rename_with(.cols = any_of(ind_specific_vars), .fn = function(x){paste0(x,'_i')})
    ind_j_values <- repertoire_allele_freqs %>% filter(individual == individual_ids[2])  %>%
      rename_with(.cols = any_of(ind_specific_vars), .fn = function(x){paste0(x,'_j')})
    
    pair_values <- full_join(ind_i_values %>% select(-total_time_point_cells_i),
                             ind_j_values %>% select(-total_time_point_cells_j)) %>%
      mutate(pair = pair) %>%
      select(any_of(variable_pars), nGCs, t, pair, allele, matches('_i'), matches('_j')) %>%
      mutate(individual_i = individual_ids[1], individual_j = individual_ids[2])
    
    return(pair_values)
  }
  
  paired_tibble <- bind_rows(lapply(as.list(unique_pairs), FUN = internal_function, repertoire_allele_freqs = repertoire_allele_freqs)) %>%
    mutate(individual_i = as.integer(individual_i), individual_j = as.integer(individual_j))
  
  # Add back total cells for each individual at each time point.
  paired_tibble <- left_join(paired_tibble,
                             total_time_point_cells %>% rename(individual_i = individual, total_time_point_cells_i = total_time_point_cells))
  paired_tibble <- left_join(paired_tibble,
                             total_time_point_cells %>% rename(individual_j = individual, total_time_point_cells_j = total_time_point_cells)) %>%
    select(any_of(variable_pars), nGCs, t, pair, allele, matches('_i'), matches('_j')) 
  
  # keep only alleles present in the germline set of both individuals. For plotting purposes, compute Spearman ranks
  # (i.e., with smaller values assigned smaller ranks)
  paired_tibble <- paired_tibble %>%
    filter(!is.na(naive_freq_i), !is.na(naive_freq_j)) %>%
    group_by(across(c(any_of(variable_pars), nGCs, 'pair', 't'))) %>%
    mutate(spearman_rank_freq_i = rank(experienced_freq_i),
           spearman_rank_freq_j = rank(experienced_freq_j),
           spearman_rank_freq_ratio_i = rank(freq_ratio_i),
           spearman_rank_freq_ratio_j = rank(freq_ratio_j)) %>%
    ungroup()
  
  return(paired_tibble)
  
}

compute_pairwise_correlations <- function(pairwise_allele_freqs, variable_pars){
  bind_rows(pairwise_allele_freqs %>% mutate(method = 'pearson'),
            pairwise_allele_freqs %>% mutate(method = 'spearman')) %>%
    group_by(across(c(any_of(variable_pars), 'pair', 'nGCs','t', 'method'))) %>%
    mutate(n_points_in_pairwise_comparison = n()) %>%
    filter(n_points_in_pairwise_comparison >= 3) %>%
    summarise(freq_correlation = cor.test(experienced_freq_i, experienced_freq_j, method = unique(method))$estimate,
              freq_ratio_correlation = cor.test(freq_ratio_i, freq_ratio_j, method = unique(method))$estimate) %>%
    ungroup()
}


# Summarise statistics across individual (or pairs of individuals) per time point
# Median, 1st and 3rd quartiles

summarise_across_individuals <- function(data, vars_to_summarise){
  data %>%
    summarise(across(all_of(vars_to_summarise),
                     list(median = median, lowerq = ~ quantile(.x, 0.25, na.rm = F),
                          upperq = ~ quantile(.x, 0.75, na.rm = F))))
}

count_increases_and_decreases <- function(repertoire_allele_freqs, variable_pars){
  repertoire_allele_freqs %>% 
    mutate(direction_of_change = case_when(
      freq_ratio > 1 ~ 'increasing',
      freq_ratio == 1 ~ 'stable',
      freq_ratio < 1 ~ 'decreasing'
    )) %>%
    group_by(across(c(any_of(variable_pars), 'nGCs' ,'t','allele', 'allele_type_affinity', 'allele_type_mutability', 'direction_of_change'))) %>%
    summarise(n_individuals = n()) %>%
    ungroup() %>%
    pivot_wider(names_from = direction_of_change, values_from = n_individuals, values_fill = 0) %>%
    mutate(net_direction = increasing - decreasing,
           fraction_increasing = increasing / (increasing + decreasing),
           fraction_decreasing = decreasing / (increasing + decreasing))
    
}

find_variable_parameters <- function(model_parameters){
  n_par_values <- model_parameters %>% summarise(across(everything(), function(x){length(unique(x))})) %>% unlist()
  variable_pars <- names(n_par_values)[n_par_values >1]
  
  return(variable_pars)
}

# ----- Tests ------


# Quick plots for visual tests:
quick_plotting_function <- function(GC_tibble, allele_info){
  left_join(GC_tibble, allele_info %>% select(allele, allele_type_affinity)) %>%
    group_by(t, clone_id, allele, allele_type_affinity) %>% count() %>%
    ggplot(aes(x = t, y = n, group = clone_id)) +
    geom_line(aes(color = allele_type_affinity)) +
    theme(legend.position = 'top')
}


test_allele_info <- assign_allele_properties(allele_info = test_allele_info,
                                             baseline_mean = 1, s = 0, sigma_r = 1, gamma = 1)

# Some visual tests:
# If mu_max is < delta, clones should consistently die out
# (this will stop simulations before tmax)
   # quick_plotting_function(
   #   GC_tibble = simulate_GC_dynamics(K = 500, I_total = 100, t_imm = 10, mu_max = 0.25, delta = 0.5, mutation_rate = 0.01, mutation_sd = 0.1,
   #                        allele_info = test_allele_info, tmax = 100, observation_times = c(seq(5,100,5))),
   #   allele_info = test_allele_info)

# Otherwise the total GC population size should remain around K
 # simulate_GC_dynamics(K = 300, I_total = 100, t_imm = 6, mu_max = 2, delta = 0.5, mutation_rate = 0.01, mutation_sd = 0.1,
 #                       allele_info = test_allele_info, tmax = 50, observation_times = c(1,seq(5,50,5))) %>%
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
#                           relative_mutability = c(rep(1,20), rep(1, 20)),
#                          clone_id = c(rep(1,20), rep(2, 20)),
#                          t = 0)

# Without immigration and mutation, clone 2 should consistently win out:
# bind_rows(replicate(n = 5, simulate_GC_dynamics(K = 100, I_total = 0, t_imm = 1, mu_max = 2, delta = 0.5,
#                                                      mutation_rate = 0, mutation_sd = 0,
#                                                      allele_info = test_allele_info, tmax = 50,
#                                                      observation_times = c(seq(5,50,5)),
#                                                      initial_GC_state = test_initial_GC),
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
#                                                 initial_GC_state = test_initial_GC),
#                     simplify = F), .id = 'replicate') %>%
#   group_by(clone_id, replicate, t) %>%
#   count() %>%
#   ggplot(aes(x = t, y = n)) +
#   geom_line(aes(group = clone_id, color = factor(clone_id))) +
#   scale_y_continuous(limits = c(0, NA)) +
#   facet_wrap('replicate') +
#   theme(legend.position = 'top')

 

