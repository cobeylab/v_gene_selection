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

# File specifying model parameters
model_parameters_file_path <- args[2] # model_parameters_file_path = '../results/simulations/neutral_scenario_1/model_parameters.csv'
# Individual to use naive frequencies and germline gene sets from
individual_id <- args[3]
# 
GC_number <- args[4]

allele_info <- read_csv(allele_info_file_path) %>%
  filter(individual == individual_id) %>% select(-individual)

output_directory <- paste0(dirname(allele_info_file_path),'/')

model_parameters <- read_csv(model_parameters_file_path)
model_parameters <- unlist(model_parameters)
for(i in 1:length(model_parameters)){
  assign(names(model_parameters)[i], model_parameters[i])
}

# See script with simulation functions for parameter definitions

simulation <- master_simulation_function(K = K,
                             I_total = I_total,
                             t_imm = t_imm,
                             mu_max = mu_max,
                             delta = delta, 
                             mutation_rate = mutation_rate,
                             mutation_sd = mutation_sd,
                             allele_info = allele_info,
                             tmax = tmax, 
                             fixed_initial_affinities = fixed_initial_affinities)
  
write_csv(simulation %>%
            mutate(individual = individual_id, GC = GC_number) %>%
            select(individual, GC, everything()),
          file = paste0(output_directory,'simulation_individual_', individual_id, '_GC_', GC_number, '.csv'))

