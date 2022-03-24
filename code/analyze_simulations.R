library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(cowplot)
library(stringr)
library(viridis)
theme_set(theme_cowplot())
source('simulations_gillespie.R')

# Directory containing simulation results and allele_info.csv file with allele properties
args <- commandArgs(trailingOnly = T) 
results_directory <- args[1] # results_directory <- '../results/simulations/scenario_1/'
facet_wrap_var <- args[2]
color_var <- args[3]



# Read allele info, model parameters and combined simulations file
allele_info <- read_csv(paste0(results_directory, 'allele_info.csv'))
parameter_files <- list.files(paste0(results_directory, 'raw_simulation_files'), recursive = T, pattern = 'model_parameters',
                              full.names = T)
model_parameters <- bind_rows(lapply(as.list(parameter_files), FUN = read_csv))

simulations <- read_csv(paste0(results_directory, 'combined_simulations.csv'))

# Figure out which parameters vary:
n_par_values <- model_parameters %>% summarise(across(everything(), function(x){length(unique(x))})) %>% unlist()
variable_pars <- names(n_par_values)[n_par_values >1]

# Mutate some model parameters into factors
simulations <- simulations %>%
  mutate(mutation_rate = factor(mutation_rate),
         beta = factor(beta))


# Compute allele frequencies over time and metrics of diversity for each GC in each individual

allele_freqs_by_GC <- compute_allele_freqs_per_GC(simulations = simulations, variable_pars = variable_pars)

allele_diversity_by_GC <- compute_allele_diversity_per_GC(allele_freqs_by_GC = allele_freqs_by_GC,
                                                          variable_pars = variable_pars)

clone_diversity_by_GC <- compute_clone_diversity_per_GC(simulations = simulations, variable_pars = variable_pars)

GC_statistics <- left_join(clone_diversity_by_GC, allele_diversity_by_GC, by = c(variable_pars, 'individual','t','GC'))


# Repertoire-wide allele freqs (aggregated across individual GCs)
repertoire_allele_freqs <- compute_repertoire_allele_freqs(allele_freqs_by_GC = allele_freqs_by_GC,
                                                           allele_info = allele_info, variable_pars = variable_pars)

repertoire_allele_diversity <- compute_repertoire_allele_diversity(repertoire_allele_freqs, variable_pars)

# N individuals where each allele is increasing/decreasing relative to naive repertoire over time
n_increasing_decreasing <- count_increases_and_decreases(repertoire_allele_freqs = repertoire_allele_freqs,
                                                         variable_pars = variable_pars)

# Pairwise allele frequency correlations
pairwise_allele_freqs <- get_pairwise_sim_freqs(repertoire_allele_freqs, variable_pars = variable_pars)


pairwise_correlations <- compute_pairwise_correlations(pairwise_allele_freqs = pairwise_allele_freqs, 
                                                       variable_pars = variable_pars)


# Summarise statistics across individual (or pairs of individuals) per time point
summary_repertoire_allele_diversity <- repertoire_allele_diversity %>%
  group_by(across(c(any_of(variable_pars),'t'))) %>%
  summarise_across_individuals(vars_to_summarise = c('repertoire_allele_diversity', 'n_alleles_in_experienced_repertoire')) %>%
  ungroup()

summary_GC_stats <- GC_statistics %>%
  group_by(across(c(any_of(variable_pars),'t'))) %>%
  summarise_across_individuals(vars_to_summarise = c('fraction_biggest_clone', 'fraction_most_common_allele', 'total_GC_pop')) %>%
  ungroup()

summary_pairwise_correlations <- pairwise_correlations %>%
  group_by(across(c(any_of(variable_pars),'t', 'method'))) %>%
  summarise_across_individuals(vars_to_summarise = c('freq_correlation', 'freq_ratio_correlation')) %>%
  ungroup()


# =========== PLOTS ============

base_plotting_function <- function(summary_tibble, y_var, color_var, facet_wrap_var){
  
  summary_tibble %>%
    mutate(across(all_of(facet_wrap_var), function(x){factor(paste0(facet_wrap_var, ' = ',as.character(x)),
                                                             levels = paste0(facet_wrap_var, ' = ',  sort(unique(x))))})) %>%
    ggplot(aes_string(x = 't', y = paste0(y_var, '_median'), color = color_var, group = color_var)) +
    geom_line(size = 1.5) +
    geom_linerange(aes_string(ymin = paste0(y_var, '_lowerq'), ymax = paste0(y_var, '_upperq')),
                   alpha = 0.5, size = 1.2) +
    geom_point(size = 3) +
    facet_wrap(facet_wrap_var, nrow = 1) +
    xlab('Time') +
    theme(legend.position = 'top')
}


color_legend <- get_legend(base_plotting_function(summary_pairwise_correlations %>% filter(method == 'pearson'),
                                                              y_var = 'freq_correlation', color_var = color_var,
                                                              facet_wrap_var = facet_wrap_var))
                           #+ scale_color_discrete(name = '               Mutation rate'))


pairwise_corr_freqs_pl <- base_plotting_function(summary_pairwise_correlations %>% filter(method == 'pearson'),
                                              y_var = 'freq_correlation', color_var = color_var, facet_wrap_var = facet_wrap_var) +
  ylab('Pairwise correlation\nin allele frequencies') +
  geom_hline(yintercept = 0, linetype =2) +
  theme(strip.background = element_blank(), strip.text = element_blank(),
        legend.position = 'none') +
  xlab('')
  
pairwise_corr_freq_ratios_pl <- base_plotting_function(summary_pairwise_correlations %>% filter(method == 'pearson'),
                                                    y_var = 'freq_ratio_correlation', color_var = color_var,
                                                    facet_wrap_var = facet_wrap_var) +
  ylab('Pairwise correlation in\nexperienced-to-naive ratios') +
  geom_hline(yintercept = 0, linetype =2) +
  #ylim(-1,1) +
  theme(strip.background = element_blank(), strip.text = element_blank(),
        legend.position = 'none') +
  xlab('')
  

n_alleles_in_rep_pl <- base_plotting_function(summary_repertoire_allele_diversity,
                                           y_var = 'n_alleles_in_experienced_repertoire', color_var = color_var,
                                           facet_wrap_var = facet_wrap_var) +
  ylab('Number of V alleles\nin the response') +
  theme(strip.background = element_blank(), strip.text = element_blank(),
        legend.position = 'none') +
  xlab('Time')


freq_dominant_clone_pl <- base_plotting_function(summary_GC_stats, y_var = 'fraction_biggest_clone', color_var = color_var,
                                              facet_wrap_var = facet_wrap_var) +
  ylab('Frequency of biggest\nclone within GCs') +
  ylim(0,1) +
  theme(legend.position = 'none') +
  xlab('')
  

main_panel <- plot_grid(color_legend,
                        freq_dominant_clone_pl,
                        pairwise_corr_freqs_pl,
                        pairwise_corr_freq_ratios_pl,
                        n_alleles_in_rep_pl,
                          align = 'v',
                        axis = 'l',
                        nrow = 5,
                        rel_heights = c(0.1,1,1,1,1))
  

save_plot(paste0(results_directory, basename(results_directory), '_main_panel.pdf'),
          main_panel,
          base_height = 12, base_width = 10)


# Other plots
GC_total_pop_pl <- base_plotting_function(summary_GC_stats, y_var = 'total_GC_pop', color_var = color_var,
                                          facet_wrap_var = facet_wrap_var) +
  ylab('Total population per GC') +
  theme(legend.position = 'top') +
  scale_color_discrete(name = 'Mutation rate')

save_plot(paste0(results_directory, basename(results_directory), '_total_GC_population.pdf'),
          GC_total_pop_pl,
          base_height = 4, base_width = 8)


if("high_avg" %in% unique(allele_info$allele_type_affinity)){
  
  # For each high-average allele, plot fraction of individuals where allele has increased / decreased in freq. relative to naive rep. over time.
  fraction_increasing_pl <- n_increasing_decreasing %>%
    filter(allele_type_affinity == 'high_avg') %>%
    ggplot(aes(x = t, y = fraction_increasing)) +
    geom_line(aes(group = allele), alpha= 1) +
    geom_line(data = n_increasing_decreasing %>% 
                filter(allele_type_affinity == 'high_avg') %>%
                group_by(across(c(any_of(variable_pars), 't'))) %>%
                summarise(fraction_increasing = mean(fraction_increasing)),
              color = 'blue', size = 2, alpha = 0.7) +
    #geom_smooth() +
    facet_grid(reformulate(color_var, facet_wrap_var)) +
    geom_hline(yintercept = 0.5, linetype = 2) +
    xlab('Time') +
    ylim(0,1) +
    ylab('Fraction of individuals where allele increased in frequency') +
    theme(panel.border = element_rect(colour = 'black', linetype = 1))
  
  save_plot(paste0(results_directory, basename(results_directory), '_fraction_individuals_high_avg_alleles_increasing.pdf'),
            fraction_increasing_pl,
            base_height = 12, base_width = 10)
  
  combined_freq_of_high_avg_alleles_in_GCs <- left_join(allele_freqs_by_GC, allele_info %>%
              select(individual, allele, allele_type_affinity)) %>%
    filter(s == 0.5) %>%
    group_by(across(c(any_of(variable_pars), 't', 'individual','GC'))) %>%
    summarise(combined_freq_high_avg = sum(allele_freq[allele_type_affinity == 'high_avg']),
              combined_freq_low_avg = sum(allele_freq[allele_type_affinity == 'low_avg'])) %>%
    ungroup() 
  
  total_GCs_across_individuals = length(unique(simulations$individual)) * length(unique(simulations$GC))
  
  combined_freq_of_high_avg_alleles_in_GC_pl <- combined_freq_of_high_avg_alleles_in_GCs %>%
    filter(t %in% c(10, max(t))) %>%
    mutate(s = paste0('allele advantage = ', s),
           t = paste0(t, ' days')) %>%
    ggplot(aes(x = combined_freq_high_avg)) +
    geom_histogram(bins = 20) + 
    facet_grid(reformulate('t', color_var, )) +
    scale_y_continuous(labels = function(x){round(x/total_GCs_across_individuals, 2)}) +
    ylab('Fraction of GCs') +
    xlab('Combined frequency of high-affinity alleles within GC') +
    geom_vline(xintercept = 0.5, alpha = 0.3, linetype = 2)
  
  save_plot(paste0(results_directory, basename(results_directory), '_combined_freq_of_high_avg_alleles_in_GCS_.pdf'),
            combined_freq_of_high_avg_alleles_in_GC_pl,
            base_height = 10, base_width = 10)
  
  # Cool to look at but hard to see patterns
  # left_join(allele_freqs_by_GC, allele_info %>%
  #             select(individual, allele, allele_type_affinity)) %>%
  #   group_by(across(c(any_of(variable_pars), 't', 'individual','GC', 'allele_type_affinity'))) %>%
  #   summarise(combined_freq = sum(allele_freq)) %>%
  #   ungroup() %>%
  #   filter(allele_type_affinity == 'high_avg') %>%
  #   mutate(plotting_group = paste(individual, GC, sep = '_')) %>%
  #   ggplot(aes(x = t, y = combined_freq)) +
  #   geom_line(aes(group = plotting_group)) +
  #   facet_grid(mutation_rate~I_total)
  
}


#simulations %>%
  # filter(mutation_rate == 0.01, I_total == 100) %>%
  # group_by(individual, t,GC) %>%
  # summarise(mean_affinity = sum(clone_freq*mean_affinity)) %>%
  # ungroup() %>%
  # mutate(plotting_group = paste0(individual, GC, sep = ';')) %>%
  # ggplot(aes(x = t, y = mean_affinity)) +
  # geom_line(aes(group = plotting_group), alpha = 0.5) +
  # geom_smooth()





