library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(readr)
theme_set(theme_cowplot())

args <- commandArgs(trailingOnly = T)

source('simulations_gillespie.R')

# File specifying alleles' affinity distributions and naive frequencies
allele_info_file_path <- args[1] # allele_info_file_path = '../results/simulations/neutral_scenario_1/allele_info.csv'


# File specifying germina center parameters
model_parameters_file_path <- args[2] # model_parameters_file_path = '../results/simulations/neutral_scenario_1/model_parameters.csv'
individual_id <- args[3]

allele_info <- read_csv(allele_info_file_path) %>%
  filter(individual == individual_id) %>% select(-individual)

output_directory <- paste0(dirname(allele_info_file_path),'/')

model_parameters <- read_csv(model_parameters_file_path)
model_parameters <- unlist(model_parameters)
for(i in 1:length(model_parameters)){
  assign(names(model_parameters)[i], model_parameters[i])
}

# See script with simulation functions for parameter definitions

individual_simulation <- run_simulation(K = K,
                                        nGCs = nGCs,
                                        lambda_imm = lambda_imm,
                                        mu_max = mu_max,
                                        delta = delta, 
                                        mutation_rate = mutation_rate,
                                        mutation_sd = mutation_sd,
                                        allele_info = allele_info,
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
