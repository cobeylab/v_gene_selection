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
lambda_max <- 8 # expected reproductive rate per B cell in an empty germinal center

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

# Initializes GCs with a poisson distributed number of clones, each seeded by 1 cell
# Mean of distribution = base_n_seed_clones
sample_immigrants <- function(allele_info, mu, clone_numbering_start){
  
  average_naive_affinity <- allele_info %>%
    summarise(A = sum(expected_affinity * naive_freqs)) %>% pull(A)
  
  # Compute number of clones arriving by allele (this is set up so, on average, mu clones will arrive)
  arrivals_by_allele <- allele_info %>%
    rowwise() %>%
    mutate(new_clones = rpois(n = 1, lambda = naive_freq*expected_affinity*mu/average_naive_affinity)) %>%
    filter(new_clones > 0)
  
  if(nrow(arrivals_by_allele) > 0){
    
    immigrants_tibble <- arrivals_by_allele %>%
      uncount(new_clones) %>%
      rowwise() %>%
      mutate(affinity = rgamma(n = 1, shape = alpha, rate = beta)) %>%
      select(allele, affinity)
    
    clone_ids <- seq(from = clone_numbering_start, by = 1, length.out = nrow(immigrants_tibble))
    immigrants_tibble$clone_id = clone_ids
    
    immigrants_tibble <- immigrants_tibble %>% select(clone_id, everything())
    
    
  }else{
    immigrants_tibble <- NULL
  }
  
  return(immigrants_tibble)
}
# For a test, the command below should return a value close to 100:
#mean(replicate(100, nrow(sample_immigrants(allele_info, mu = 100, 1)), simplify = T))


# Samples population at time t given population at time t-1 (GC_tibble) and total size expected at Nt under logistic backbone
draw_GC_tplus1 <- function(allele_info, Ntplus1, GC_t, poisson_mean_n_immigrants, mutation_rate, mutation_sd){
  
  current_t <- unique(GC_t$t)
  stopifnot(length(current_t) == 1)
  
  
  # At beginning of time step, potential immigrants arrive
  new_immigrants <- sample_immigrants(allele_info = allele_info, poisson_mean_n_immigrants = poisson_mean_n_immigrants)
  
  # New immigrants will compete with current occupants to be parents to the next generation
  GC_t_with_new_immigrants <- bind_rows(GC_t, new_immigrants) %>%
    mutate(cell_id = 1:n())
  
  sampling_probs <- GC_t_with_new_immigrants$affinity/sum(GC_t_with_new_immigrants$affinity)
  
  parent_cells <- sample(GC_t_with_new_immigrants$cell_id, size = Ntplus1, replace = T, prob = sampling_probs)
  
  GC_tplus1 <- left_join(tibble(parent_id = parent_cells),
                         GC_t_with_new_immigrants %>% rename(parent_id = cell_id), by = 'parent_id') %>%
    mutate(t = current_t + 1) %>%
    select(t, everything()) %>%
    select(-cell_id)
  
}


simulate_GC_dynamics <- function(allele_info, lambda_max, K, mu, mutation_rate, mutation_sd, tmax){
  
  # To initialize GC, sample immigrants from poisson_mean_n_immigrants until there's at least one clone arriving
  time = 1
  GC_tibble <- NULL
  while(is.null(GC_tibble)){
    GC_tibble <- sample_immigrants(allele_info = allele_info, mu = mu,
                                   clone_numbering_start = 1) %>%
      mutate(t = time) %>% select(t, everything())
  }
  
  # Generate logistic backbone with N0 set by initial number of cells
  logistic_backbone <- create_logistic_backbone(r = r, K = K, N0 = nrow(GC_tibble))
  time_N_reaches_K <- min(which(logistic_backbone == K))
  
  unique_clone_ids <- unique(GC_tibble$clone_id)
  
  while(time < timax){
    if(time < time_N_reaches_K){
      Ntplus1 = logistic_backbone[time + 1]
    }else{
      Ntplus1 = K
    }
    # Part of tibble corresponding to current time step
    GC_t <- GC_tibble %>% filter(t == time)
    
    
    next_tibble
  }
  
  
  
  
  
}





