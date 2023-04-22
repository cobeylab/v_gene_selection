source('simulation_functions.R')
theme_set(theme_cowplot())

args <- commandArgs(trailingOnly = T)

# File specifying alleles' affinity distributions and naive frequencies
allele_info_file_path <- args[1] # allele_info_file_path = '../results/simulations/scenario_1/allele_info.csv'

# Directory for specific parameter values
model_parameters_directory <- args[2] # 

# Individual to use naive frequencies and germline gene sets from
individual_id <- args[3]

# Arbitrary number assigned to GC being simulated (different GCs simulated by different jobs in the cluster) 
GC_number <- args[4]

allele_info <- read_csv(allele_info_file_path) %>%
  filter(individual == individual_id) %>% select(-individual)


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

simulation <- master_simulation_function(K = K,
                             I_total = I_total,
                             t_imm = t_imm,
                             mu_max = mu_max,
                             delta = delta, 
                             mutation_rate = mutation_rate,
                             mutation_sd = mutation_sd,
                             allele_info = allele_info,
                             tmax = tmax)
 
# Add parameter values to simulation tibble
for(i in 1:length(par_values)){
  simulation[,names(par_values)[i]] <- par_values[i]
}

# Store results in directory specific to the combination of parameter values used
write_csv(simulation %>%
            mutate(individual = individual_id, GC = GC_number) %>%
            select(individual, GC, everything()),
          file = paste0(model_parameters_directory,'simulation_individual_', individual_id, '_GC_', GC_number, '.csv'))


time_end <- Sys.time()

