library(dplyr)
library(tidyr)
library(readr)
library(stringr)
source('simulation_functions.R')

# Directory containing simulation results and allele_info.csv file with allele properties
args <- commandArgs(trailingOnly = T) 
results_directory <- args[1] # e.g. results_directory <- '../results/simulations/neutral_scenario/'

# Read allele info, model parameters and combined simulations file
allele_info <- read_csv(paste0(results_directory, 'allele_info.csv'))
model_parameters <- read_csv(paste0(results_directory, 'combined_model_parameters.csv'))
simulations <- read_csv(paste0(results_directory, 'combined_simulations.csv'))

n_individuals <- length(unique(simulations$individual))
n_GCs_per_individual <- length(unique(simulations$GC))

# Figure out which parameters vary:
variable_pars <- find_variable_parameters(model_parameters)


# Compute allele frequencies over time and metrics of diversity for each GC in each individual (grouped by variable model params)

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
n_inds_increasing_decreasing <- count_increases_and_decreases(repertoire_allele_freqs = repertoire_allele_freqs,
                                                         variable_pars = variable_pars)

# Pairwise allele frequency correlations
pairwise_allele_freqs <- get_pairwise_sim_freqs(repertoire_allele_freqs, variable_pars = variable_pars)


pairwise_correlations <- compute_pairwise_correlations(pairwise_allele_freqs = pairwise_allele_freqs, 
                                                       variable_pars = variable_pars)


# Summarise statistics across individual (or pairs of individuals) per time point (always grouping by variable params.)
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

if("high_avg" %in% unique(allele_info$allele_type_affinity)){
  combined_freq_of_high_avg_alleles_in_GCs <- left_join(allele_freqs_by_GC, allele_info %>%
                                                          select(individual, allele, allele_type_affinity)) %>%
    group_by(across(c(any_of(variable_pars), 't', 'individual','GC'))) %>%
    summarise(combined_freq_high_avg = sum(allele_freq[allele_type_affinity == 'high_avg']),
              combined_freq_low_avg = sum(allele_freq[allele_type_affinity == 'low_avg'])) %>%
    ungroup() 
  
}else{
  combined_freq_of_high_avg_alleles_in_GCs <- NULL
}

if("high_mut" %in% unique(allele_info$allele_type_mutability)){
  combined_freq_of_high_mut_alleles_in_GCs <- left_join(allele_freqs_by_GC, allele_info %>%
                                                          select(individual, allele, allele_type_mutability)) %>%
    group_by(across(c(any_of(variable_pars), 't', 'individual','GC'))) %>%
    summarise(combined_freq_high_mut = sum(allele_freq[allele_type_mutability == 'high_mut']),
              combined_freq_low_mut = sum(allele_freq[allele_type_mutability == 'low_mut'])) %>%
    ungroup() 
    
}else{
  combined_freq_of_high_mut_alleles_in_GCs <- NULL
}


# Export .RData file with a list contaning summaries of current scenario:
summary_object_name = paste0(basename(results_directory), '_summary')
assign(summary_object_name,
       list(n_individuals = n_individuals,
            n_GCs_per_individual = n_GCs_per_individual,
            summary_repertoire_allele_diversity = summary_repertoire_allele_diversity,
            summary_GC_stats = summary_GC_stats,
            summary_pairwise_correlations = summary_pairwise_correlations,
            n_inds_increasing_decreasing = n_inds_increasing_decreasing,
            #allele_diversity_by_GC = allele_diversity_by_GC,
            #allele_freqs_by_GC = allele_freqs_by_GC,
            #clone_diversity_by_GC = clone_diversity_by_GC,
            #GC_statistics = GC_statistics, 
            #repertoire_allele_diversity = repertoire_allele_diversity,
            #repertoire_allele_freqs = repertoire_allele_freqs,
            combined_freq_of_high_avg_alleles_in_GCs = combined_freq_of_high_avg_alleles_in_GCs,
            combined_freq_of_high_mut_alleles_in_GCs = combined_freq_of_high_mut_alleles_in_GCs
       ))


save(list = summary_object_name, 
     file = paste0(results_directory, summary_object_name, '.RData'))

