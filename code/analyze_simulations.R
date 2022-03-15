library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(cowplot)
library(stringr)
library(viridis)
theme_set(theme_cowplot())
source('simulation_functions.R')

# Directory containing simulation results and allele_info.csv file with allele properties
args <- commandArgs(trailingOnly = T) 
results_directory <- args[1] # results_directory <- '../results/simulations/neutral_scenario_1/'

print(results_directory)

allele_info <- read_csv(paste0(results_directory, 'allele_info.csv'))


parameter_files <- list.files(paste0(results_directory, 'raw_simulation_files'), recursive = T, pattern = 'model_parameters',
                              full.names = T)
model_parameters <- bind_rows(lapply(as.list(parameter_files), FUN = read_csv))


model_parameters <- read_csv(paste0(results_directory, 'model_parameters.csv'))
simulations <- read_csv(paste0(results_directory, 'combined_simulations.csv'))

# Plot annotation with parameters (have to redo now that looking at >1 parameter combination)
#parameter_annotation <- paste(paste(names(model_parameters) , model_parameters[1,], sep = ' = '), collapse = ' ; ')
parameter_annotation <- ''

# Figure out which parameters vary:
n_par_values <- model_parameters %>% summarise(across(everything(), function(x){length(unique(x))})) %>% unlist()
variable_pars <- names(n_par_values)[n_par_values >1]


# Compute allele frequencies over time and metrics of diversity for each GC in each individual
allele_freqs_by_GC <- simulations %>%
  group_by(across(c(any_of(variable_pars),'individual', 't', 'GC', 'allele'))) %>%
  summarise(n_cells = sum(clone_size)) %>%
  group_by(across(c(any_of(variable_pars),'individual', 't', 'GC'))) %>%
  mutate(allele_freq = n_cells / sum(n_cells)) %>%
  ungroup()

allele_diversity_by_GC <-  allele_freqs_by_GC %>%
  group_by(across(c(any_of(variable_pars),'individual', 't', 'GC'))) %>%
  summarise(n_alleles = length(unique(allele)),
            fraction_most_common_allele = max(allele_freq),
            allele_diversity = 1 - sum(allele_freq^2)) %>%
  ungroup()

clone_diversity_by_GC <- simulations %>%
  group_by(across(c(any_of(variable_pars),'individual', 't', 'GC'))) %>%
  summarise(n_clones = length(unique(clone_id)),
            total_GC_pop = sum(clone_size),
            fraction_biggest_clone = max(clone_freq),
            clone_diversity = 1 - sum(clone_freq^2)) %>%
  ungroup()

GC_statistics <- left_join(clone_diversity_by_GC, allele_diversity_by_GC, by = c(variable_pars, 'individual','t','GC'))

mean_GC_stats_per_time_per_individual <- GC_statistics %>%
  group_by(across(c(any_of(variable_pars),'individual', 't'))) %>%
  summarise(across(-any_of('GC'), mean))

mean_GC_stats_per_time <- mean_GC_stats_per_time_per_individual %>%
  group_by(across(c(any_of(variable_pars),'t'))) %>%
  summarise(across(-any_of('individual'), mean))
  
# Repertoire-wide allele freqs (aggregated across individual GCs)
repertoire_allele_freqs <- allele_freqs_by_GC %>%
  group_by(across(c(any_of(variable_pars),'individual', 't', 'allele'))) %>%
  summarise(n_cells = sum(n_cells)) %>%
  group_by(across(c(any_of(variable_pars),'individual', 't'))) %>%
  mutate(experienced_freq = n_cells / sum(n_cells)) %>%
  ungroup()


# Completes repertoire_allele_freqs tibble so that zeros are explicitly represented in allele counts in the experienced repertoire
# (uses allele_info to find the full set of alleles present in the naive repertoire for each individual)

complete_repertoire_allele_freqs <- function(repertoire_allele_freqs, allele_info, variable_pars){
  
  # All alleles of an individual explicitly represented at all time points
  complete_scaffold <- left_join(expand_grid(individual = unique(repertoire_allele_freqs$individual),
                                             t = unique(repertoire_allele_freqs$t)),
                                 allele_info %>% select(individual, allele))
  
  if(length(variable_pars) > 0){
    complete_scaffold <- expand_grid(model_parameters[,variable_pars], complete_scaffold)
  }

  completed_tibble <- left_join(complete_scaffold, repertoire_allele_freqs) %>%
    replace_na(list(n_cells = 0, experienced_freq = 0))
  
  return(completed_tibble)
}


repertoire_allele_freqs <- complete_repertoire_allele_freqs(repertoire_allele_freqs = repertoire_allele_freqs,
                                                            allele_info = allele_info, variable_pars = variable_pars)

# Compute experienced frequencies and total cells in the repertoire per time point
repertoire_allele_freqs <- repertoire_allele_freqs %>%
  group_by(across(c(any_of(variable_pars), 't', 'individual'))) %>%
  mutate(total_time_point_cells = sum(n_cells),
         allele_rank = rank(-experienced_freq, ties.method = 'average')) %>%
  ungroup()

# Add allele affinities /types, add naive allele frequencies and compute experienced-to-naive ratios,
repertoire_allele_freqs <- left_join(repertoire_allele_freqs, allele_info %>%
                             select(individual, allele, allele_type, naive_freq, expected_affinity)) %>%
  mutate(freq_ratio_log = log(experienced_freq) - log(naive_freq),
         freq_ratio = exp(freq_ratio_log)) %>% select(-freq_ratio_log)

# Order alleles as a factor
allele_order <- paste0('V',sort(as.integer(str_remove(unique(repertoire_allele_freqs$allele), 'V'))))
repertoire_allele_freqs <- repertoire_allele_freqs %>%
  mutate(allele = factor(allele, levels = allele_order))

repertoire_allele_diversity <- repertoire_allele_freqs %>%
  group_by(across(c(any_of(variable_pars), 't', 'individual', 'total_time_point_cells'))) %>%
  summarise(repertoire_allele_diversity = 1 - sum(experienced_freq^2),
            n_alleles_in_experienced_repertoire = sum(experienced_freq>0)) %>%
  ungroup()



# Not worth trying to use get_pairwise_freqs function written for obs data (too many other variables/groupings) 
# Best to write a function specific for these simulations even though it will look similar
get_pairwise_values <- function(repertoire_allele_freqs, variable_pars){
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
                                                               individual, t, total_time_point_cells) %>% unique()
  
  
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
      select(any_of(variable_pars), pair, t, allele, matches('_i'), matches('_j')) %>%
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
    select(any_of(variable_pars), pair, t, allele, matches('_i'), matches('_j')) 
  
  # keep only alleles present in the germline set of both individuals. For plotting purposes, compute Spearman ranks
  # (i.e., with smaller values assigned smaller ranks)
  paired_tibble <- paired_tibble %>%
    filter(!is.na(naive_freq_i), !is.na(naive_freq_j)) %>%
    group_by(across(c(any_of(variable_pars), 'pair', 't'))) %>%
    mutate(spearman_rank_freq_i = rank(experienced_freq_i),
           spearman_rank_freq_j = rank(experienced_freq_j),
           spearman_rank_freq_ratio_i = rank(freq_ratio_i),
           spearman_rank_freq_ratio_j = rank(freq_ratio_j)) %>%
    ungroup()
  
  return(paired_tibble)
    
}

pairwise_allele_freqs <- get_pairwise_values(repertoire_allele_freqs, variable_pars = variable_pars)


pairwise_correlations <- bind_rows(pairwise_allele_freqs %>% mutate(method = 'pearson'),
                                   pairwise_allele_freqs %>% mutate(method = 'spearman')) %>%
  group_by(across(c(any_of(variable_pars), 'pair', 't', 'method'))) %>%
  mutate(n_points_in_pairwise_comparison = n()) %>%
  filter(n_points_in_pairwise_comparison >= 3) %>%
  summarise(freq_correlation = cor.test(experienced_freq_i, experienced_freq_j, method = unique(method))$estimate,
            freq_ratio_correlation = cor.test(freq_ratio_i, freq_ratio_j, method = unique(method))$estimate) %>%
  ungroup()

# Now making plots

# Some plots have alleles' naive frequencies on x axis.
# If using uniform naive frequencies, plot the alleles themselves as the x axis
if(any(model_parameters$uniform_naive_freqs)){
  x_axis_var <- 'allele'
  xlabel <- 'Allele'
  col_width <- 0.7 # Bar width for directionality plot
}else{
  x_axis_var <- 'naive_freq'
  xlabel <- 'Naive frequency'
  col_width <- 0.001 # Bar width for directionality plot
  # (For some reason different on a continuous vs. discrete axis)
}

median_pairwise_correlations <- pairwise_correlations %>%
  group_by(across(c(any_of(variable_pars), 't', 'method'))) %>%
  summarise(freq_correlation_lowerq = quantile(freq_correlation, 0.25, na.rm = T),
            freq_correlation_upperq = quantile(freq_correlation,0.75, na.rm = T),
            freq_ratio_correlation_lowerq = quantile(freq_ratio_correlation, 0.25, na.rm = T),
            freq_ratio_correlation_upperq = quantile(freq_ratio_correlation,0.75, na.rm = T),
            freq_correlation = median(freq_correlation),
            freq_ratio_correlation = median(freq_ratio_correlation),
            ) %>%
  ungroup()

median_pairwise_correlations %>%
  filter(method == 'pearson') %>%
  ggplot(aes(x = t, y = freq_correlation, color = factor(mutation_rate), group = mutation_rate)) +
  geom_line(size = 1.5) +
  geom_point(size = 3) +
  geom_linerange(aes(ymin = freq_correlation_lowerq, ymax = freq_correlation_upperq),
                 alpha = 0.5, size = 1.2) +
  facet_wrap('mu_max') +
  xlab('Time') +
  ylab('Pairwise correlation in allele frequencies') +
  ylim(-1,1) +
  geom_hline(yintercept = 0, linetype =2) +
  theme(legend.position = 'top')
  
  
pairwise_corr_freqs <- pairwise_correlations  %>%
  ggplot(aes(x = t, y = freq_correlation)) +
  geom_line(aes(group = pair), alpha = 0.1) +
  geom_linerange(data = median_pairwise_correlations, 
                 aes(ymin = freq_correlation_lowerq, ymax = freq_correlation_upperq),
                 color = 'red', alpha = 0.5) +
  geom_line(data = median_pairwise_correlations, color = 'red', size = 1.5) +
  geom_hline(yintercept = 0, linetype =2) +
  xlab('Time') +
  ylab('Pairwise correlation in allele frequencies') +
  ylim(-1,1) +
  facet_wrap('method')


pairwise_corr_freq_ratios <- pairwise_correlations %>%
  ggplot(aes(x = t, y = freq_ratio_correlation)) +
  geom_line(aes(group = pair), alpha = 0.2) +
  geom_linerange(data = median_pairwise_correlations, 
                 aes(ymin = freq_ratio_correlation_lowerq, ymax = freq_ratio_correlation_upperq),
                 color = 'red', alpha = 0.5) +
  geom_line(data = median_pairwise_correlations, color = 'red', size = 1.5)  +
  geom_hline(yintercept = 0, linetype =2) +
  xlab('Time') +
  ylab('Pairwise correlation in experienced-to-naive ratios') +
  ylim(-1,1) +
  facet_wrap('method')


n_alleles_in_rep <- repertoire_allele_diversity %>%
  ggplot(aes(x = t, y = n_alleles_in_experienced_repertoire)) +
  geom_line(aes(group = individual), alpha = 0.5) +
  geom_line(data = repertoire_allele_diversity %>%
              group_by(t) %>% 
              summarise(n_alleles_in_experienced_repertoire = mean(n_alleles_in_experienced_repertoire)),
            color = 'red', size = 1.5) +
  xlab('Time') + ylab('Number of V alleles in the experienced repertoire')

repertoire_allele_diversity_pl <- repertoire_allele_diversity %>%
  ggplot(aes(x = t, y = repertoire_allele_diversity)) +
  geom_line(aes(group = individual), alpha = 0.5) +
  geom_line(data = repertoire_allele_diversity %>%
              group_by(t) %>% 
              summarise(repertoire_allele_diversity = mean(repertoire_allele_diversity)),
            color = 'red', size = 1.5) +
  xlab('Time') + ylab('Repertoire allele diversity')

naive_freq_vs_n_inds_present <- repertoire_allele_freqs %>%
  filter(t == max(t)) %>%
  group_by(allele, naive_freq, allele_type) %>%
  summarise(n_inds_allele_present = sum(experienced_freq >0)) %>%
  ggplot(aes_string(x = x_axis_var, y = 'n_inds_allele_present')) +
  geom_point(aes(color = allele_type)) +
  xlab(xlabel) +
  ylab('N. individuals with allele \npresent at the end of simulation') +
  #theme(legend.position = c(0.7,0.1)) +
  theme(legend.position = 'top')  +
  scale_color_discrete(name = 'Allele type')


naive_vs_exp_freq <- repertoire_allele_freqs %>%
  filter(t == max(t)) %>%
  ggplot(aes_string(x = x_axis_var, y = 'experienced_freq')) +
  geom_point(alpha = 0.5, aes(color = allele_type)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  xlab(xlabel) +
  ylab('Experienced frequency\nat the end of simulation') +
  theme(legend.position = 'top')  +
  scale_color_discrete(name = 'Allele type')

mean_GC_total_pop <-  mean_GC_stats_per_time_per_individual  %>%
  ggplot(aes(x = t, y = total_GC_pop)) +
  geom_point(aes(group = individual), alpha = 0.5) +
  geom_line(data = mean_GC_stats_per_time,
            color = 'red', size = 1.5) +
  xlab('Time') +
  ylab('Mean total GC population') +
  geom_hline(yintercept = model_parameters$K, linetype = 2)

mean_freq_dominant_clone <- mean_GC_stats_per_time_per_individual  %>%
  ggplot(aes(x = t, y = fraction_biggest_clone)) +
  geom_point(aes(group = individual), alpha = 0.5) +
  geom_line(data = mean_GC_stats_per_time,
            color = 'red', size = 1.5) +
  geom_vline(xintercept = 16, linetype = 2) +
  geom_hline(yintercept = 0.4, linetype = 2) +
  xlab('Time') +
  ylab('Mean frequency of biggest clone within GCs')

mean_freq_dominant_allele  <- mean_GC_stats_per_time_per_individual  %>%
  ggplot(aes(x = t, y = fraction_most_common_allele)) +
  geom_point(aes(group = individual), alpha = 0.5) +
  geom_line(data = mean_GC_stats_per_time,
            color = 'red', size = 1.5) +
  xlab('Time') +
  ylab('Mean frequency of most common allele within GCs')

mean_alleles_per_gc <- mean_GC_stats_per_time_per_individual %>%
  ggplot(aes(x = t, y = n_alleles)) +
  geom_line(aes(group = individual), alpha = 0.5) +
  geom_line(data = mean_GC_stats_per_time, color = 'red', size = 1.5) +
  xlab('Time') +
  ylab('Mean number of alleles per GC')

mean_clones_per_gc <- mean_GC_stats_per_time_per_individual %>%
  ggplot(aes(x = t, y = n_clones)) +
  geom_line(aes(group = individual), alpha = 0.5) +
  geom_line(data = mean_GC_stats_per_time, color = 'red', size = 1.5) +
  xlab('Time') +
  ylab('Mean number of clones per GC')


directionality_per_allele <- repertoire_allele_freqs %>%
  filter(t %in% c(10,50,max(t))) %>%
  mutate(direction_of_change = case_when(
    freq_ratio > 1 ~ 'increase',
    freq_ratio == 1 ~ 'stable',
    freq_ratio < 1 ~ 'decrease'
  )) %>%
  group_by(t, allele, allele_type, naive_freq, direction_of_change) %>%
  summarise(n_individuals = n()) %>%
  ungroup() %>%
  mutate(n_individuals = ifelse(direction_of_change == 'decrease',-n_individuals, n_individuals)) %>%
  ggplot(aes_string(x = x_axis_var, y = 'n_individuals', fill = 'allele_type')) +
  geom_col(width = col_width) +
  scale_y_continuous(limits = c(-length(unique(repertoire_allele_freqs$individual)),
                                length(unique(repertoire_allele_freqs$individual)))) +
  geom_hline(yintercept = 0, linetype = 2, size = 1.5) +
  facet_wrap('t') +
  ylab('Number of individuals\n decreasing | increasing') +
  xlab(xlabel) +
  theme(panel.border = element_rect(color = 'black'),
        legend.position = 'top') +
  background_grid()

  
pairwise_corr_row <- plot_grid(pairwise_corr_freqs, pairwise_corr_freq_ratios, nrow = 1)
repertoire_stats_row <- plot_grid(n_alleles_in_rep, repertoire_allele_diversity_pl, nrow = 1)
GC_stats_row <- plot_grid(mean_GC_total_pop, mean_freq_dominant_clone,
                          mean_freq_dominant_allele, nrow = 1)
naive_freq_scatterplots_row <- plot_grid(naive_vs_exp_freq,naive_freq_vs_n_inds_present, nrow = 1)


main_panel <- plot_grid(pairwise_corr_row,
                        repertoire_stats_row,
                        GC_stats_row,
                        naive_freq_scatterplots_row,
                        directionality_per_allele,
                        nrow = 5)
  

save_plot(paste0(results_directory, basename(results_directory), '_main_panel.pdf'),
          plot_grid(ggdraw() + 
                      draw_label(
                        parameter_annotation,
                        fontface = 'bold',
                        x = 0,
                        hjust = -0.05,
                        size = 12
                      ),
                    main_panel,
                    nrow = 2,
                    rel_heights = c(1, 20)),
          base_height = 25, base_width = 13)


# Results for Pearson correlation specifically
pearson_pairwise_freqs <- pairwise_correlations  %>%
  filter(method == 'pearson') %>%
  ggplot(aes(x = t, y = freq_correlation)) +
  geom_line(aes(group = pair), alpha = 0.1) +
  geom_linerange(data = median_pairwise_correlations %>% filter(method == 'pearson'), 
                 aes(ymin = freq_correlation_lowerq, ymax = freq_correlation_upperq),
                 color = 'red', alpha = 0.8, size = 1) +
  geom_line(data = median_pairwise_correlations %>%
              filter(method == 'pearson'), color = 'red', size = 1.5) +
  geom_hline(yintercept = 0, linetype =2) +
  xlab('Time') +
  ylab('Pairwise correlation in\nallele frequencies') +
  ylim(-1,1) 

pearson_pairwise_freq_ratios <- pairwise_correlations %>%
  filter(method == 'pearson') %>%
  ggplot(aes(x = t, y = freq_ratio_correlation)) +
  geom_line(aes(group = pair), alpha = 0.2) +
  geom_linerange(data = median_pairwise_correlations %>%
                   filter(method == 'pearson'), 
                 aes(ymin = freq_ratio_correlation_lowerq, ymax = freq_ratio_correlation_upperq),
                 color = 'red', alpha = 0.8, size = 1) +
  geom_line(data = median_pairwise_correlations %>%
              filter(method == 'pearson'), color = 'red', size = 1.5)  +
  geom_hline(yintercept = 0, linetype =2) +
  xlab('Time') +
  ylab('Pairwise correlation in\nexperienced-to-naive ratios') +
  ylim(-1,1)

save_plot(paste0(results_directory,'pearson_pairwise_cors.pdf'),
          plot_grid(pearson_pairwise_freqs, 
                    pearson_pairwise_freq_ratios,
                    nrow = 2),
          base_width = 5, base_height = 8)
  
# Arrow plots
# Arrow plots
individual_sample <- sample(unique(repertoire_allele_freqs$individual), size = 10, replace = F)

arrow_plot <- repertoire_allele_freqs %>%
  filter(t %in% c(10,50, max(t)), individual %in% individual_sample) %>%
  filter(allele_rank <= 20) %>%
  mutate(label_position = ifelse(experienced_freq > naive_freq, 1.05*experienced_freq, 1.05*naive_freq)) %>%
  
  ggplot(aes(x = allele_rank)) +
  geom_segment(aes(xend = allele_rank, y = naive_freq, yend = experienced_freq, color = allele_type),
               arrow = arrow(ends = 'last', length = unit(10, "pt"), type = 'closed')) +
  facet_grid(t~individual, scales = 'free') +
  geom_text(aes(label = allele,
                x = allele_rank, y = label_position), angle = 20, size = 3, alpha = 0.8) +
  scale_x_continuous(expand = c(0.15,0)) +
  xlab('Top 20 alleles in each individual') +
  ylab('Frequency') +
  theme(legend.position = 'top') +
  background_grid()

save_plot(paste0(results_directory, basename(results_directory), '_arrow_plot.pdf'),
          arrow_plot,
          base_height = 10, 
          base_width = 20)


# Freq_ratio scatterplots

pair_sample <- sample(unique(pairwise_allele_freqs$pair), 10, replace = F)

scatterplot_data <- left_join(pairwise_allele_freqs, allele_info %>% select(allele, allele_type) %>% unique()) %>%
  filter(t %in% c(10,max(pairwise_allele_freqs$t)), pair %in% pair_sample) %>%
  group_by(t, pair) %>%
  # Spearman correlation ranks (different from allele ranks used for plotting)
  mutate(rank_freq_i = rank(experienced_freq_i, ties.method = 'average'),
         rank_freq_j = rank(experienced_freq_j, ties.method = 'average'),
         rank_freq_ratio_i = rank(freq_ratio_i, ties.method = 'average'),
         rank_freq_ratio_j = rank(freq_ratio_j, ties.method = 'average')) %>%
  ungroup()

scatterplot_data %>%
  ggplot(aes(x = experienced_freq_i, y = experienced_freq_j)) +
  geom_point(size = 3, aes(color = allele_type), alpha = 0.8) +
  geom_smooth(method = 'lm', color = 'black') +
  facet_grid(t~pair, scales = 'free') +
  #scale_color_viridis() +
  theme(legend.position = 'top',
        panel.border = element_rect(color = 'black')) +
  xlab('Experienced frequency in individual i') +
  ylab('Experienced frequency in individual j') 
  #scale_y_log10() +
  #scale_x_log10()

freq_scatterplots <- scatterplot_data %>%
  ggplot(aes(x = experienced_freq_i, y = experienced_freq_j)) +
  geom_point(size = 3, aes(color = allele_type), alpha = 0.8) +
  geom_smooth(method = 'lm', color = 'black') +
  facet_grid(t~pair, scales = 'free') +
  #scale_color_viridis() +
  #scale_y_log10() +
  #scale_x_log10() +
  theme(legend.position = 'top',
        panel.border = element_rect(color = 'black')) +
  xlab('Experienced frequency in individual i') +
  ylab('Experienced frequency in individual j')

rank_freq_scatterplots <- scatterplot_data %>%
  ggplot(aes(x = rank_freq_i, y = rank_freq_j)) +
  geom_point(size = 3, aes(color = allele_type), alpha = 0.8) +
  geom_smooth(method = 'lm', color = 'black') +
  facet_grid(t~pair, scales = 'free') +
  #scale_color_viridis() +
  #scale_y_log10() +
  #scale_x_log10() +
  theme(legend.position = 'top',
        panel.border = element_rect(color = 'black')) +
  xlab('Rank experienced frequency in individual i') +
  ylab('Rank experienced frequency in individual j')



freq_ratio_scatterplots <- scatterplot_data %>%
  ggplot(aes(x = freq_ratio_i, y = freq_ratio_j)) +
  geom_point(size = 3, aes(color = allele_type), alpha = 0.8) +
  geom_smooth(method = 'lm', color = 'black') +
  facet_grid(t~pair, scales = 'free') +
  #scale_color_viridis() +
  theme(legend.position = 'top',
        panel.border = element_rect(color = 'black')) +
  xlab('Experienced-to-naive ratio in individual i') +
  ylab('Experienced-to-naive ratio in individual j') +
  scale_y_log10() +
  scale_x_log10() +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_vline(xintercept = 1, linetype = 2) 

rank_freq_ratio_scatterplots <- scatterplot_data %>%
  ggplot(aes(x = rank_freq_ratio_i, y = rank_freq_ratio_j)) +
  geom_point(size = 3, aes(color = allele_type), alpha = 0.8) +
  geom_smooth(method = 'lm', color = 'black') +
  facet_grid(t~pair, scales = 'free') +
  #scale_color_viridis() +
  theme(legend.position = 'top',
        panel.border = element_rect(color = 'black')) +
  xlab('Rank experienced-to-naive ratio in individual i') +
  ylab('Rank experienced-to-naive ratio in individual j')

save_plot(paste0(results_directory, basename(results_directory), '_freq_ratio_scatterplots.pdf'),
          plot_grid(freq_ratio_scatterplots,
                    rank_freq_ratio_scatterplots, nrow = 2),
          base_height = 16, 
          base_width = 20)

#save_plot(paste0(results_directory, basename(results_directory), '_freq_ratio_scatterplots.pdf'),
#          freq_ratio_scatterplots ,
#          base_height = 8, 
#          base_width = 20)

#save_plot(paste0(results_directory, basename(results_directory), 'rank_freq_ratio_scatterplots.pdf'),
#          rank_freq_ratio_scatterplots ,
#          base_height = 8, 
#          base_width = 20)








# OTHER PLOTS

shared_genes_in_top_10 <- pairwise_allele_freqs %>%
  group_by(pair, t) %>%
  filter(allele_rank_i <= 10 & allele_rank_j <= 10) %>%
  summarise(shared_genes_in_top_10 = n()) %>%
  ungroup()

shared_genes_in_top_10 %>%
  ggplot(aes(x = t, y = shared_genes_in_top_10)) +
  geom_point(size = 2, alpha = 0.05, position = position_jitter(width = 0, height = 0.1)) +
  #geom_smooth(method = 'lm')
  geom_line(data = shared_genes_in_top_10 %>% group_by(t) %>%
              summarise(shared_genes_in_top_10 = mean(shared_genes_in_top_10)), color = 'red')


# Fraction of sequences accounted for by top 10 alleles
repertoire_allele_freqs %>%
  group_by(t, individual,) %>%
  filter(allele_rank <= 10) %>%
  summarise(fraction_seqs_in_top_alleles = sum(experienced_freq)) %>%
  ggplot(aes(x = t, y = fraction_seqs_in_top_alleles, group = individual)) +
  geom_line() +
  scale_x_log10()

# Number of top 10 alleles over time coming from each class of allele
repertoire_allele_freqs %>%
  group_by(t, individual,) %>%
  filter(allele_rank <= 10, experienced_freq > 0) %>% # experienced_freq >0 to only count alleles effectively present
  group_by(t, individual, allele_type) %>%
  count() %>%
  mutate(plotting_group = paste(individual, allele_type, sep = ';')) %>%
  ggplot(aes(x = t, y = n, color = allele_type)) +
  geom_line(aes(group = plotting_group), alpha = 0.3) +
  geom_smooth() +
  xlab('Time') +
  ylab('Number of alleles in top-10')
  

# Combined frequency over time for each type of allele
combined_freq_by_allele_type <- repertoire_allele_freqs %>%
  group_by(individual, t, allele_type) %>%
  summarise(combined_freq = sum(experienced_freq)) %>%
  ungroup()

# Setting naive freqs as time 0 freqs.
combined_freqs_in_naive_rep <- allele_info %>% group_by(allele_type) %>%
  summarise(combined_freq = sum(naive_freq))

combined_freqs_in_naive_rep <- left_join(expand.grid(individual = unique(repertoire_allele_freqs$individual),
                      allele_type = unique(repertoire_allele_freqs$allele_type),
                      t = 0),
          combined_freqs_in_naive_rep)

combined_freq_by_allele_type <- bind_rows(combined_freq_by_allele_type,
                                          combined_freqs_in_naive_rep)


combined_freq_by_allele_type %>%
  mutate(plotting_group = paste(individual, allele_type, sep = ';')) %>%
  ggplot(aes(x = t, y = combined_freq)) +
  geom_line(aes(group = plotting_group, color = allele_type), alpha = 0.3) +
  scale_x_log10() +
  geom_smooth(aes(color = allele_type))

repertoire_allele_freqs %>%
  filter(individual %in% as.character(1:10)) %>%
  ggplot(aes(x = t, y = experienced_freq, group = allele, color = allele_type)) +
  geom_line() +
  scale_x_log10() +
  facet_wrap('individual')







# Old plots (review and delete)

# repertoire_allele_freqs %>%
#   filter(t == max(t)) %>%
#   mutate(change_in_frequency = case_when(
#     freq_ratio ==1 ~ 'none',
#     freq_ratio > 1 ~ 'increase',
#     freq_ratio < 1 ~ 'decrease'
#   )) %>%
#   group_by(allele, allele_type, naive_freq, change_in_frequency) %>%
#   count() %>%
#   group_by(allele, allele_type, naive_freq) %>%
#   mutate(prob = n/sum(n)) %>%
#   pivot_wider(names_from = change_in_frequency, values_from = prob) %>%
#   ggplot(aes(x = naive_freq, y = decrease, color = allele_type)) +
#   geom_point() +
#   geom_smooth() +
#   xlab('Naive frequency') +
#   ylab('Fraction of individuals with frequency decrease') +
#   ylim(0,1)



#example_GC <- example_GCs %>% filter(individual == 1)


# quick_plotting_function(example_GC)  +
#    xlab('Time') +
#    ylab('Number of cells') +
#   theme(legend.position = 'none')


# quick_plotting_function(example_GC) +
#   xlab('Time') +
#   ylab('Number of cells') +
#   theme(legend.position = 'none') +
#   xlim(0,100) +
#   scale_y_log10()

# example_GC_statistics <- example_GC  %>%
#   group_by(t) %>%
#   summarise(n_clones = length(unique(clone_id)),
#             n_alleles = length(unique(allele)),
#             mean_affinity = mean(affinity))

# example_GC_statistics %>% ggplot(aes(x = t, y = mean_affinity)) +
#   geom_line() +
#   xlab('Time') +
#   ylab('Average affinity within GC')

# example_GC_statistics %>% ggplot(aes(x = t, y = n_clones)) + geom_line() +
#   xlab('Time') +
#   ylab('Number of clones in GC')
# 
# GC_statistics %>% ggplot(aes(x = t, y = n_alleles)) + geom_line()  +
#   xlab('Time') +
#   ylab('Number of V alleles in GC')

example_GCs %>% group_by(t, clone_id, individual, GC) %>% count() %>%
  ungroup() %>%
  ggplot(aes(x = t, y = n)) +
  geom_line(aes(group = clone_id)) +
  #xlim(0,50) +
  #scale_y_log10() +
  facet_wrap('individual')


example_GCs %>% group_by(individual, GC, t, clone_id) %>% summarise(mean_affinity = mean(affinity)) %>%
  ungroup() %>%
  ggplot(aes(x = t, y = mean_affinity)) +
  geom_line(aes(group = clone_id)) +
  facet_wrap('individual')


clone_sizes_at_tmax <- example_GCs %>%
  group_by(individual, GC, clone_id) %>%
  filter(t == max(t)) %>%
  count() %>%
  dplyr::rename(clone_final_size = n) %>%
  group_by(individual, GC) %>%
  mutate(clone_final_freq = clone_final_size / sum(clone_final_size)) %>%
  ungroup()

clone_arrival_times <- example_GCs %>%
  group_by(individual, GC, clone_id) %>%
  summarise(clone_arrival = min(t))
  
left_join(clone_sizes_at_tmax, clone_arrival_times) %>%
  ggplot(aes(x = clone_arrival, y = clone_final_freq )) +
  geom_point() +
  facet_wrap('individual')


example_GCs %>% filter(individual == 12) %>%
  filter(t == max(t)) %>%
  group_by(clone_id) %>% count() %>% arrange(desc(n))

example_GCs %>% filter(individual == 12, clone_id %in% c(6,8)) %>%
  group_by(clone_id, t) %>%
  summarise(mean_affinity = mean(affinity)) %>%
  ggplot(aes(x = t, y = mean_affinity, group = clone_id)) +
  geom_line(aes(color = factor(clone_id)))

