library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(readr)
theme_set(theme_cowplot())

args <- commandArgs(trailingOnly = T)

source('simulation_functions.R')

# File specifying alleles' affinity distributions and naive frequencies
allele_info_file_path <- args[1] # allele_info_file_path = '../results/simulations/neutral_scenario_1/allele_info.csv'


# File specifying germina center parameters
GC_parameters_file_path <- args[2] # GC_parameters_file_path = '../results/simulations/neutral_scenario_1/GC_parameters.csv'
individual_id <- args[3]

allele_info <- read_csv(allele_info_file_path)
output_directory <- paste0(dirname(allele_info_file_path),'/')

GC_parameters <- read_csv(GC_parameters_file_path)
GC_parameters <- unlist(GC_parameters)
for(i in 1:length(GC_parameters)){
  assign(names(GC_parameters)[i], GC_parameters[i])
}
# nGCs: Number of germinal centers in an individual
# tmax: Number of timesteps observed
# K: carrying capacity of germinal centers
# mu :expected number of newly recruited clones arriving at germinal centers per time step
# lambda_max: 1.5 expected reproductive rate per B cell in an empty germinal center
# mutation_rate: mutation probability per B cell per time step
# mutation_sd: standard deviation for the distribution of mutational effects (mean 0)
# uniform_naive_freqs: TRUE if all alleles have the same naive frequency.

# If randomizing naive frequencies in each individual, export a record of each individual's naive freqs.
#if(randomize_naive_freqs_in_each_individual){
#  allele_info$naive_freq <- sample(allele_info$naive_freq, size = length(allele_info$naive_freq), replace =F)
#  write_csv(allele_info %>% select(allele, naive_freq) %>% mutate(individual = individual_id) %>%
#              select(individual, everything()),
#            file = paste0(output_directory,'randomized_naive_freqs_individual_', individual_id, '.csv'))
#  
#}

individual_simulation <- simulate_repertoire(nGCs = nGCs,
                                             allele_info = allele_info,
                                             lambda_max = lambda_max,
                                             K = K,
                                             mu = mu,
                                             mutation_rate = mutation_rate,
                                             mutation_sd = mutation_sd,
                                             tmax = tmax)
  
write_csv(individual_simulation$allele_counts,
          file = paste0(output_directory,'repertoire_counts_individual_', individual_id, '.csv'))
write_csv(individual_simulation$GC_statistics,
          file = paste0(output_directory,'GC_statistics_individual_', individual_id, '.csv'))

# Export a detailed GC trajectory only for the first 5 individuals
if(as.integer(individual_id) <= 5){
  write_csv(individual_simulation$example_GC_trajectory,
            file = paste0(output_directory,'example_GC_individual_', individual_id, '.csv'))
}

#system.time(
#  simulate_repertoire_allele_counts(nGCs = 10, allele_info, lambda_max, K, mu, mutation_rate, mutation_sd, tmax)
#)  

