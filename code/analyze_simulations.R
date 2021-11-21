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

### ADD A CHECK TO SEE IF NAIVE FREQS WERE RANDOMIZED IN EACH INDIVIDUAL




example_GC_files <- list.files(results_directory, pattern = 'example_GC', full.names = T)[1]

example_GCs <- lapply(as.list(example_GC_files), FUN = read_csv)
example_GCs <- bind_rows(example_GCs, .id = 'individual') %>%
  select(individual, everything())

individual_repertoire_files <- list.files(results_directory, pattern = 'repertoire_counts',
                                         full.names = T)

allele_counts <- lapply(as.list(individual_repertoire_files), FUN = read_csv)
allele_counts <- bind_rows(allele_counts, .id = 'individual') %>%
  select(individual, everything())


allele_counts <- allele_counts %>%
  group_by(t, individual) %>%
  mutate(total_time_point_cells = sum(n),
         experienced_freq = n/total_time_point_cells) %>%
  ungroup()

# Add naive allele frequencies and allele expected affinitie
allele_counts <- left_join(allele_counts, allele_info %>% select(allele, allele_type, naive_freq, expected_affinity)) %>%
  mutate(freq_ratio_log = log(experienced_freq) - log(naive_freq),
         freq_ratio = exp(freq_ratio_log)) %>% select(-freq_ratio_log)

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
                 color = 'red') +
  geom_line(data = mean_pairwise_correlations, color = 'black', size = 1.5) +
  geom_hline(yintercept = 0, linetype =2) +
  xlab('time') +
  ylab('Pairwise correlation in allele frequencies') +
  ylim(-1,1)

pairwise_correlations %>%
  ggplot(aes(x = t, y = freq_ratio_correlation)) +
  geom_line(aes(group = pair), alpha = 0.2) +
  geom_linerange(data = mean_pairwise_correlations, 
                 aes(ymin = freq_ratio_correlation_lowerq, ymax = freq_ratio_correlation_upperq),
                 color = 'red') +
  geom_line(data = mean_pairwise_correlations, color = 'black', size = 1.5)  +
  geom_hline(yintercept = 0, linetype =2) +
  xlab('time') +
  ylab('Pairwise correlation in experienced-to-naive ratios') +
  ylim(-1,1)





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
  







example_GC <- example_GCs %>% filter(individual == 5)


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

GC_statistics <- example_GC  %>%
  group_by(t) %>%
  summarise(n_clones = length(unique(clone_id)),
            n_alleles = length(unique(allele)),
            mean_affinity = mean(affinity))

GC_statistics %>% ggplot(aes(x = t, y = mean_affinity)) +
  geom_line() +
  xlab('Time') +
  ylab('Average affinity within GC')

GC_statistics %>% ggplot(aes(x = t, y = n_clones)) + geom_line() +
  xlab('Time') +
  ylab('Number of clones in GC')

GC_statistics %>% ggplot(aes(x = t, y = n_alleles)) + geom_line()  +
  xlab('Time') +
  ylab('Number of V alleles in GC')







