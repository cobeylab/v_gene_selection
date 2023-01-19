library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(stringr)
theme_set(theme_cowplot())

# ======================== PLOTTING FUNCTIONS  ========================

dir.create('../figures/simulations/', showWarnings = F)

# Base plotting function
base_plotting_function <- function(summary_tibble, y_var, color_var, facet_vars = NULL){
  
  plotting_tibble <- summary_tibble %>%
    mutate(mutation_rate = as.factor(as.character(mutation_rate)))
  
  if('nGCs' %in% names(plotting_tibble)){
    plotting_tibble <- plotting_tibble %>%
      mutate(nGCs = factor(as.character(nGCs),levels = as.character(sort(unique(nGCs)))))
  }

  
  # Stop if parameters not listed as facet or color variables have more than 1 value in summary tibble
  # (Prevents mixing together multiple parameter values in the same panels)
  n_par_values_check <- summary_tibble %>%
    select(any_of(parameter_names)) %>%
    select(-any_of(c(color_var, facet_vars))) %>%
    summarise(n_par_values = across(everything(), function(x){length(unique(x))})) %>%
    unlist()
  
  stopifnot(all(n_par_values_check == 1))

  if(!is.null(facet_vars)){
    stopifnot(length(facet_vars) <= 2)
    plotting_tibble <- plotting_tibble %>%
      mutate(across(all_of(facet_vars[1]), function(x){factor(paste0(facet_vars[1], ' = ', as.character(x)),
                                                              levels = paste0(facet_vars[1], ' = ',  sort(unique(x))))}))
    
    if(length(facet_vars) == 1){
      facet_command <- facet_wrap(facet_vars, ncol = 1)
    }else{
      facet_command <- facet_grid(reformulate(facet_vars[1], facet_vars[2]))
      plotting_tibble <- plotting_tibble %>%
        mutate(across(all_of(facet_vars[2]), function(x){factor(paste0(facet_vars[2], ' = ',as.character(x)),
                                                                levels = paste0(facet_vars[2], ' = ',  sort(unique(x))))}))
      
    }
  }else{
    facet_command <- NULL
  }
  
  pl <- plotting_tibble %>%
    ggplot(aes_string(x = 't', y = paste0(y_var, '_median'), color = color_var, group = color_var)) +
    geom_line(size = 1.5) +
    geom_linerange(aes_string(ymin = paste0(y_var, '_lowerq'), ymax = paste0(y_var, '_upperq')),
                   alpha = 0.5, size = 1.2) +
    geom_point(size = 3) +
    xlab('Time (days)') +
    theme(legend.position = 'top') +
    facet_command
  
  return(pl)
}


# Function for making main fig showing simulated frequency and freq. deviation correlations in different scenarios
generate_freq_cor_panels <- function(summary_list, main_fig_s, main_fig_beta, main_fig_I_total, main_fig_gamma){
  
  subset_tibble <- summary_list$summary_pairwise_correlations %>%
    filter(if_all(any_of('s'), function(x){x == main_fig_s})) %>%
    filter(if_all(any_of('beta'), function(x){x == main_fig_beta})) %>%
    filter(if_all(any_of('I_total'), function(x){x == main_fig_I_total})) %>%
    filter(if_all(any_of('gamma'), function(x){x == main_fig_gamma})) %>%
    filter(method == 'pearson') 
    
  theme_specs <- theme(legend.position = 'none', plot.title = element_text(size = 12,hjust = 0.5)) 
  
  freq_corr <-  base_plotting_function(summary_tibble = subset_tibble,
                                       y_var = 'freq_correlation', color_var = 'nGCs',
                                       facet_vars = 'mutation_rate')  + ylim(-0.2,1) +
    ylab('Pairwise correlation between individuals') +
    theme_specs +
    geom_hline(yintercept = 0, linetype = 2) +
    nGCs_color_scale +
    ggtitle('\n\n\nGermline allele frequencies\nin the response')
  
  freq_ratio_corr <- base_plotting_function(summary_tibble = subset_tibble,
                                            y_var = 'freq_ratio_correlation', color_var = 'nGCs',
                                            facet_vars = 'mutation_rate')  +
    ylab('') + theme(legend.position = 'none') +
    ylim(-0.2,1) +
    theme_specs +
    geom_hline(yintercept = 0, linetype = 2) +
    nGCs_color_scale +
    ggtitle('\n\n\nExperienced-to-naive\nratios')
  
  
  return(plot_grid(freq_corr, freq_ratio_corr, nrow = 1))
}

nGCs_color_scale <- scale_color_manual(name = 'Number of GCs per individual',
                                            values = brewer.pal(n = 9, "BuPu")[c(3,5,7,9)])

mutation_rate_color_scale <- scale_color_manual(name = 'Mutation rate (affinity-changing\nmutations per B cell per division)',
                                                values = brewer.pal(n = 9, "BuPu")[c(3,6,9)])


# ======================== MAKING PLOTS ========================
# Load simulation summaries
scenarios <- c('neutral_scenario', 'high_affinity_scenario', 'high_mutation_scenario')

model_parameters <- bind_rows(
  lapply(as.list(scenarios),
         FUN = function(scen){
           read_csv(paste0('../results/simulations/', scen, '/combined_model_parameters.csv')) %>%
             mutate(scenario = str_remove(scen,'_scenario'))
         })
) %>% select(scenario, everything())

parameter_names <- names(model_parameters)

for(scen in scenarios){
  load(paste0('../results/simulations/', scen, '/', scen, '_summary.RData'))
}


# Main text figure with behavior of frequency correlations
# Across scenarios, where applicable, use simulations with these values
main_fig_s <- 1.5
main_fig_beta <- 4
main_fig_I_total <- 200
main_fig_gamma <- 6


# Get nGC legend from an arbitrary plot
nGC_color_legend <- get_legend(base_plotting_function(summary_tibble = high_affinity_scenario_summary$summary_pairwise_correlations %>%
                                    filter(beta == main_fig_beta,s == main_fig_s,
                                           mutation_rate == max(mutation_rate)),
                                  y_var = 'freq_correlation', color_var = 'nGCs', facet_vars = NULL) +
             nGCs_color_scale)
  

correlations_fig_neutral <- generate_freq_cor_panels(summary_list = neutral_scenario_summary, main_fig_s = main_fig_s, main_fig_beta = main_fig_beta,
                                                     main_fig_I_total = main_fig_I_total, main_fig_gamma = main_fig_gamma)

correlation_fig_high_affinity <- generate_freq_cor_panels(summary_list = high_affinity_scenario_summary, main_fig_s = main_fig_s,
                                                          main_fig_beta = main_fig_beta, main_fig_I_total = main_fig_I_total, 
                                                          main_fig_gamma = main_fig_gamma)


simulated_correlations_main_fig <- plot_grid(correlations_fig_neutral, correlation_fig_high_affinity,
                                             labels = c('A) Germline V alleles are functionally identical',
                                                        'B) Some V alleles encode receptors with higher affinity'),
                                             hjust = -0.05, scale = 0.95,
                                             label_size = 15)

simulated_correlations_main_fig <- plot_grid(
  simulated_correlations_main_fig,
  plot_grid(NULL, nGC_color_legend, nrow =1, rel_widths = c(2,5)),
  nrow = 2,
  rel_heights = c(20,1)
)


save_plot('../figures/simulations/simulated_correlations.pdf',
          simulated_correlations_main_fig,
          base_width = 12, base_height = 9)

# ======= Supplementary figures

# Fraction top lineage in GCs (neutral scenario)
fraction_top_lineage_in_GCs <- base_plotting_function(summary_tibble = neutral_scenario_summary$summary_GC_stats %>%
                                                        filter(beta == main_fig_beta),
                                                      y_var = 'fraction_biggest_clone', color_var = 'mutation_rate', facet_vars = NULL) +
  mutation_rate_color_scale +
  ylab('Frequency of biggest lineage within GCs')

save_plot('../figures/simulations/fraction_top_lineage_in_GCs.pdf',
          fraction_top_lineage_in_GCs,
          base_width = 6, base_height = 5)

# Supplementary correlation panel for high-affinity scenario
# (Assume 15 germinal centers per individual)
high_affinity_scenario_supp_fig <- base_plotting_function(summary_tibble = high_affinity_scenario_summary$summary_pairwise_correlations %>%
                                                            filter(method == 'pearson', nGCs == 15),
                                                          y_var = 'freq_ratio_correlation',
                                                          color_var = 'mutation_rate', facet_vars = c('s', 'beta')) +
  ylim(-0.2,1) + geom_hline(yintercept = 0, linetype = 2) +
  mutation_rate_color_scale +
  ylab('Pairwise correlation in experienced-to-naive frequency ratios')

save_plot('../figures/simulations/high_affinity_scenario_supp_fig.pdf',
          high_affinity_scenario_supp_fig,
          base_width = 8, base_height = 9)


# % GCs occupied by high-affinity alleles in the high-affinity scenario

stopifnot(all(high_affinity_scenario_summary$nGC_values == high_mutation_scenario_summary$nGC_values))

total_GCs_across_individuals <- high_affinity_scenario_summary$n_individuals * max(high_affinity_scenario_summary$nGC_values)


combined_freq_of_high_avg_alleles_in_GC_pl <- high_affinity_scenario_summary$combined_freq_of_high_avg_alleles_in_GCs %>%
  filter(beta == main_fig_beta, s == main_fig_s) %>%
  filter(t %in% c(10, max(t))) %>%
  mutate(s = paste0('allele advantage = ', s),
         t = paste0(t, ' days')) %>%
  ggplot(aes(x = combined_freq_high_avg)) +
  geom_histogram(bins = 20) + 
  facet_grid(reformulate('t', 'mutation_rate')) +
  scale_y_continuous(labels = function(x){round(x/total_GCs_across_individuals, 2)}) +
  ylab('Fraction of GCs') +
  xlab('Combined frequency of high-affinity alleles within GC') +
  geom_vline(xintercept = 0.5, alpha = 0.3, linetype = 2) 

save_plot('../figures/simulations/combined_freq_high_affinity_alleles.pdf',
          combined_freq_of_high_avg_alleles_in_GC_pl,
          base_width = 8, base_height = 10)


# Supplementary correlation panels for high-mutation scenario
high_mutation_scenario_supp_fig <- generate_freq_cor_panels(summary_list = high_mutation_scenario_summary, main_fig_s = main_fig_s,
                                                          main_fig_beta = main_fig_beta, main_fig_I_total = main_fig_I_total, 
                                                          main_fig_gamma = main_fig_gamma)

high_mutation_scenario_supp_fig <- plot_grid(
  high_mutation_scenario_supp_fig,
  plot_grid(NULL, nGC_color_legend, nrow =1, rel_widths = c(2,5)),
  nrow = 2,
  rel_heights = c(20,1)
)

save_plot('../figures/simulations/high_mutation_scenario_supp_fig.pdf',
          high_mutation_scenario_supp_fig,
          base_width = 8, base_height = 9)


# (For the beta x gamma one; assume 15 germinal centers per individual)
# high_mutation_scenario_supp_fig2 <-  base_plotting_function(summary_tibble = high_mutation_scenario_summary$summary_pairwise_correlations %>%
#                                                              filter(method == 'pearson', nGCs == 15), y_var = 'freq_ratio_correlation',
#                                                            color_var = 'mutation_rate', facet_vars = c('gamma', 'beta')) +
#   ylim(-0.2,1) + geom_hline(yintercept = 0, linetype = 2) +
#   mutation_rate_color_scale +
#   ylab('Pairwise correlation in experienced-to-naive frequency ratios')
# 
# save_plot('../figures/simulations/high_mutation_scenario_supp_fig2.pdf',
#           high_mutation_scenario_supp_fig2,
#           base_width = 8, base_height = 9)


# % GCs occupied by high-mutation alleles in the high mutation scenario
combined_freq_of_high_mutation_alleles_in_GC_pl <- high_mutation_scenario_summary$combined_freq_of_high_mut_alleles_in_GCs %>%
  filter(beta == main_fig_beta, gamma == main_fig_gamma) %>%
  filter(t %in% c(10, max(t))) %>%
  mutate(gamma = paste0('Relative mutability = ', gamma),
         t = paste0(t, ' days')) %>%
  ggplot(aes(x = combined_freq_high_mut)) +
  geom_histogram(bins = 20) +
  facet_grid(reformulate('t', 'mutation_rate')) +
  scale_y_continuous(labels = function(x){round(x/total_GCs_across_individuals, 2)}) +
  ylab('Fraction of GCs') +
  xlab('Combined frequency of high-mutability alleles within GC') +
  geom_vline(xintercept = 0.5, alpha = 0.3, linetype = 2)

save_plot('../figures/simulations/combined_freq_high_mut_alleles.pdf',
          combined_freq_of_high_mutation_alleles_in_GC_pl,
          base_width = 8, base_height = 10)



