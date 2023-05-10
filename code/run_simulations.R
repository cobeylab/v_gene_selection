source('simulation_functions.R')

args <- commandArgs(trailingOnly = T)

# File specifying alleles' affinity distributions and naive frequencies
allele_info_file_path <- args[1] # allele_info_file_path = '../results/simulations/neutral_scenario/allele_info.csv'

# Directory for specific parameter values
model_parameters_directory <- args[2]  

# Number of germinal centers to simulate
nGCs <- as.integer(args[3])

# Individual number (specifies a unique individual in simulations. Different sim. individuals can have the same base individual)
individual_number <- args[4]



allele_info <- read_csv(allele_info_file_path)

# Randomly sample a base individual for this simulated individual
base_individual <- sample(unique(allele_info$base_individual), size = 1, replace = F)

allele_info <- allele_info %>%
  filter(base_individual == !!base_individual) %>% select(-base_individual)


model_parameters <- read_csv(paste0(model_parameters_directory,'/model_parameters.csv'))
par_values <- unlist(model_parameters)
for(i in 1:length(par_values)){
  assign(names(par_values)[i], par_values[i])
}

# Assign affinity distributions and relative mutabilities to each allele based on s, sigma r, gamma
allele_info <- assign_allele_properties(allele_info = allele_info, baseline_mean = baseline_mean,
                                        s = s, sigma_r = sigma_r, gamma = gamma)




# Get mutation sd from beta and sigma_r
mutation_sd = sigma_r * beta

# See script with simulation functions for parameter definitions

simulation <- replicate(nGCs,
                        master_simulation_function(K = K,
                             I_total = I_total,
                             t_imm = t_imm,
                             mu_max = mu_max,
                             delta = delta, 
                             mutation_rate = mutation_rate,
                             mutation_sd = mutation_sd,
                             allele_info = allele_info,
                             tmax = tmax),
                        simplify = F
                        )

simulation <- bind_rows(simulation, .id = 'GC') %>%
  mutate(individual = individual_number, base_individual = base_individual)
 
# Add parameter values to simulation tibble
for(i in 1:length(par_values)){
  simulation[,names(par_values)[i]] <- par_values[i]
}

# Store results in directory specific to the combination of parameter values used
write_csv(simulation %>%
            select(individual, base_individual, GC, everything()),
          file = paste0(model_parameters_directory,'simulation_individual_', individual_number, '.csv'))


