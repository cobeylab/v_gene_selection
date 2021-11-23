library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(cowplot)
library(stringr)
theme_set(theme_cowplot())
source('simulation_functions.R')

# Directory containing simulation results and allele_info.csv file with allele properties
args <- commandArgs(trailingOnly = T) 
results_directory <- args[1] # results_directory <- '../results/simulations/neutral_scenario_1/'

allele_info <- read_csv(paste0(results_directory, 'allele_info.csv'))
GC_parameters <- read_csv(paste0(results_directory, 'GC_parameters.csv'))
  

example_GC_files <- list.files(results_directory, pattern = 'example_GC', full.names = T)[1]
GC_statistics_files <- list.files(results_directory, pattern = 'GC_statistics', full.names = T)
repertoire_allele_counts_files <- list.files(results_directory, pattern = 'repertoire_counts',
                                          full.names = T)

results <- lapply(list(example_GCs = example_GC_files, GC_statistics = GC_statistics_files,
                       allele_counts = repertoire_allele_counts_files),
       FUN = function(file_paths){
         bind_rows(lapply(as.list(file_paths), FUN = read_csv), .id = 'individual') %>%
           select(individual, everything())
       })

example_GCs <- results$example_GCs
allele_counts <- results$allele_counts
GC_statistics <- results$GC_statistics

# Completes allele_counts tibble so that zeros are explicitly represented in allele counts in the experienced repertoire
# (uses allele_info to find the full set of alleles present in the naive repertoire)
complete_allele_counts <- function(allele_counts, allele_info){
  complete_scaffold <- as_tibble(expand.grid(t = unique(allele_counts$t),
                                             individual = unique(allele_counts$individual),
                                             allele = unique(allele_info$allele))) %>%
    select(t, allele, individual) %>%
    arrange(t, allele, individual)
  
  left_join(complete_scaffold, allele_counts) %>%
    replace_na(list(n = 0))
  
}

allele_counts <- complete_allele_counts(allele_counts = allele_counts, allele_info = allele_info)

# Compute experienced frequencies and total cells in the repertoire per time point
allele_counts <- allele_counts %>%
  group_by(t, individual) %>%
  mutate(total_time_point_cells = sum(n),
         experienced_freq = n/total_time_point_cells) %>%
  ungroup()

# Add allele affinities /types, add naive allele frequencies and compute experienced-to-naive ratios,
allele_counts <- left_join(allele_counts, allele_info %>%
                             select(allele, allele_type, naive_freq, expected_affinity)) %>%
  mutate(freq_ratio_log = log(experienced_freq) - log(naive_freq),
         freq_ratio = exp(freq_ratio_log)) %>% select(-freq_ratio_log)


repertoire_allele_diversity <- allele_counts %>%
  group_by(t, individual, total_time_point_cells) %>%
  summarise(repertoire_allele_diversity = 1 - sum(experienced_freq^2),
            n_alleles_in_experienced_repertoire = sum(experienced_freq>0)) %>%
  ungroup()



# Not worth trying to use get_pairwise_freqs function written for obs data (too many other variables/groupings) 
# Best to write a function specific for these simulations even though it will look similar
get_pairwise_values <- function(allele_counts){
  unique_pairs <- allele_counts %>% select(individual) %>% unique() %>%
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
  total_time_point_cells <- allele_counts %>% select(individual, t, total_time_point_cells) %>% unique()
  
  
  internal_function <- function(pair, allele_counts){
    
    ind_specific_vars <- c('individual','n', 'total_time_point_cells', 'experienced_freq', 'naive_freq', 'freq_ratio')
    
    individual_ids <- str_split(pair,';')[[1]]
    
    ind_i_values <- allele_counts %>% filter(individual == individual_ids[1])  %>%
      rename_with(.cols = any_of(ind_specific_vars), .fn = function(x){paste0(x,'_i')})
    ind_j_values <- allele_counts %>% filter(individual == individual_ids[2])  %>%
      rename_with(.cols = any_of(ind_specific_vars), .fn = function(x){paste0(x,'_j')})
    
    pair_values <- full_join(ind_i_values %>% select(-total_time_point_cells_i),
                             ind_j_values %>% select(-total_time_point_cells_j)) %>%
      mutate(pair = pair) %>%
      select(pair, t, allele, matches('_i'), matches('_j')) %>%
      mutate(individual_i = individual_ids[1], individual_j = individual_ids[2])
  
    return(pair_values)
  }
  
  paired_tibble <- bind_rows(lapply(as.list(unique_pairs), FUN = internal_function, allele_counts = allele_counts))
  
  # Add back total cells for each individual at each time point.
  paired_tibble <- left_join(paired_tibble,
                             total_time_point_cells %>% rename(individual_i = individual, total_time_point_cells_i = total_time_point_cells))
  paired_tibble <- left_join(paired_tibble,
                             total_time_point_cells %>% rename(individual_j = individual, total_time_point_cells_j = total_time_point_cells)) %>%
    select(pair, t, allele, matches('_i'), matches('_j')) 
  
  # Fill NAs with zeros
  paired_tibble <- paired_tibble %>%
   mutate(across(matches(c('n_','freq_','ratio')), function(x){replace_na(x,0)}) )
  
  return(paired_tibble)
    
}

pairwise_allele_freqs <- get_pairwise_values(allele_counts)

pairwise_correlations <- pairwise_allele_freqs %>%
  #filter(experienced_freq_i >0, experienced_freq_j > 0) %>%
  group_by(pair, t) %>%
  summarise(freq_correlation = cor.test(experienced_freq_i, experienced_freq_j, method = 'spearman')$estimate,
            freq_ratio_correlation = cor.test(freq_ratio_i, freq_ratio_j, method = 'spearman')$estimate) %>%
  ungroup()

mean_pairwise_correlations <- pairwise_correlations %>%
  group_by(t) %>%
  summarise(freq_correlation_lowerq = quantile(freq_correlation, 0.25, na.rm = T),
            freq_correlation_upperq = quantile(freq_correlation,0.75, na.rm = T),
            freq_ratio_correlation_lowerq = quantile(freq_ratio_correlation, 0.25, na.rm = T),
            freq_ratio_correlation_upperq = quantile(freq_ratio_correlation,0.75, na.rm = T),
            freq_correlation = mean(freq_correlation),
            freq_ratio_correlation = mean(freq_ratio_correlation),
            ) %>%
  ungroup()
  
pairwise_correlations %>%
  ggplot(aes(x = t, y = freq_correlation)) +
  geom_line(aes(group = pair), alpha = 0.1) +
  geom_linerange(data = mean_pairwise_correlations, 
                 aes(ymin = freq_correlation_lowerq, ymax = freq_correlation_upperq),
                 color = 'red', alpha = 0.5) +
  geom_line(data = mean_pairwise_correlations, color = 'red', size = 1.5) +
  geom_hline(yintercept = 0, linetype =2) +
  xlab('Time') +
  ylab('Pairwise correlation in allele frequencies') +
  ylim(-1,1)
  

pairwise_correlations %>%
  ggplot(aes(x = t, y = freq_ratio_correlation)) +
  geom_line(aes(group = pair), alpha = 0.2) +
  geom_linerange(data = mean_pairwise_correlations, 
                 aes(ymin = freq_ratio_correlation_lowerq, ymax = freq_ratio_correlation_upperq),
                 color = 'red', alpha = 0.5) +
  geom_line(data = mean_pairwise_correlations, color = 'red', size = 1.5)  +
  geom_hline(yintercept = 0, linetype =2) +
  xlab('Time') +
  ylab('Pairwise correlation in experienced-to-naive ratios') +
  ylim(-1,1)


repertoire_allele_diversity %>%
  ggplot(aes(x = t, y = n_alleles_in_experienced_repertoire)) +
  geom_line(aes(group = individual), alpha = 0.5) +
  geom_line(data = repertoire_allele_diversity %>%
              group_by(t) %>% 
              summarise(n_alleles_in_experienced_repertoire = mean(n_alleles_in_experienced_repertoire)),
            color = 'red', size = 1.5) +
  xlab('Time') + ylab('Number of V alleles in the experienced repertoire')


allele_counts %>%
  filter(t == max(t)) %>%
  group_by(allele, naive_freq, allele_type) %>%
  summarise(n_inds_allele_present = sum(experienced_freq >0)) %>%
  ggplot(aes(x = naive_freq, y = n_inds_allele_present)) +
  geom_point(aes(color = allele_type)) +
  xlab('Naive frequency') +
  ylab('N. individuals where the allele is\npresent at the end of simulation') +
  theme(legend.position = 'top')


mean_clones_and_alleles_per_GC <- GC_statistics %>%
  group_by(individual, t) %>%
  summarise(mean_clones_per_GC = mean(n_clones),
            mean_alleles_per_GC = mean(n_alleles)) %>%
  ungroup()

mean_clones_and_alleles_per_GC %>%
  ggplot(aes(x = t, y = mean_alleles_per_GC)) +
  geom_line(aes(group = individual), alpha = 0.5) +
  xlab('Time') +
  ylab('Mean number of alleles per GC')

mean_clones_and_alleles_per_GC %>%
  ggplot(aes(x = t, y = mean_clones_per_GC)) +
  geom_line(aes(group = individual), alpha = 0.5) +
  xlab('Time') +
  ylab('Mean number of clones per GC')


allele_counts %>%
  filter(t == max(t)) %>%
  mutate(change_in_frequency = case_when(
    freq_ratio ==1 ~ 'none',
    freq_ratio > 1 ~ 'increase',
    freq_ratio < 1 ~ 'decrease'
  )) %>%
  group_by(allele, naive_freq, change_in_frequency) %>%
  count() %>%
  group_by(allele, naive_freq) %>%
  mutate(prob = n/sum(n)) %>%
  pivot_wider(names_from = change_in_frequency, values_from = prob) %>%
  ggplot(aes(x = naive_freq, y = decrease)) +
  geom_point() +
  geom_smooth() +
  xlab('Naive frequency') +
  ylab('Fraction of individuals with frequency decrease') +
  ylim(0,1)

allele_counts %>%
  filter(t == max(t)) %>%
  ggplot(aes(x = naive_freq, y = experienced_freq)) +
  geom_point(alpha = 0.5, aes(color = allele_type)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  xlab('Naive frequency') +
  ylab('Experienced frequency') +
  theme(legend.position = 'top')

allele_counts %>%
  filter(t == max(t)) %>%
  ggplot(aes(naive_freq, freq_ratio)) +
  geom_point(alpha = 0.2) +
  geom_hline(yintercept = 1, linetype = 2)



allele_counts %>%
  filter(t == max(t)) %>%
  ggplot(aes(allele, freq_ratio)) +
  geom_point(alpha = 0.2) 









allele_counts %>%
  group_by(individual, t, allele_type) %>%
  summarise(combined_freq = sum(experienced_freq)) %>%
  ungroup() %>%
  ggplot(aes(x = t, y = combined_freq)) +
  geom_point(aes(color = allele_type)) +
  scale_x_log10()



allele_counts %>%
  filter(individual == 2) %>%
  ggplot(aes(x = t, y = freq_ratio)) +
  geom_line(aes(group = allele, color = allele_type)) +
  geom_hline(yintercept = 1)
  







example_GC <- example_GCs %>% filter(individual == 1)


quick_plotting_function(example_GC)  +
  xlab('Time') +
  ylab('Number of cells') +
  theme(legend.position = 'none')


quick_plotting_function(example_GC) +
  xlab('Time') +
  ylab('Number of cells') +
  theme(legend.position = 'none') +
  #xlim(0,100) +
  scale_y_log10()

example_GC_statistics <- example_GC  %>%
  group_by(t) %>%
  summarise(n_clones = length(unique(clone_id)),
            n_alleles = length(unique(allele)),
            mean_affinity = mean(affinity))

example_GC_statistics %>% ggplot(aes(x = t, y = mean_affinity)) +
  geom_line() +
  xlab('Time') +
  ylab('Average affinity within GC')

example_GC_statistics %>% ggplot(aes(x = t, y = n_clones)) + geom_line() +
  xlab('Time') +
  ylab('Number of clones in GC')

GC_statistics %>% ggplot(aes(x = t, y = n_alleles)) + geom_line()  +
  xlab('Time') +
  ylab('Number of V alleles in GC')


test_function <- function(){
  x <- c(rep(10,30),rep(0,40))
  x <- sample(x, size = length(x), replace  =F)
  y <- sample(x, size = length(x), replace  =F)
  cor.test(x,y, method = 'spearman')$estimate
}
hist(replicate(1000,test_function(), simplify = T), breaks = 10)





