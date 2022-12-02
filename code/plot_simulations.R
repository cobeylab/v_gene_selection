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
    mutate(across(any_of(color_var), function(x){as.factor(as.character(x))}))
  
  
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
      mutate(across(all_of(facet_vars[1]), function(x){factor(paste0(facet_vars[1], ' = ',as.character(x)),
                                                              levels = paste0(facet_vars[1], ' = ',  sort(unique(x))))}))
    
    if(length(facet_vars) == 1){
      facet_command <- facet_wrap(facet_vars, nrow = 1)
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
  
  fraction_biggest_clone <- base_plotting_function(summary_tibble = summary_list$summary_GC_stats %>%
                                                     filter(across(any_of('s'), function(x){x == main_fig_s})) %>%
                                                     filter(across(any_of('beta'), function(x){x == main_fig_beta})) %>%
                                                     filter(across(any_of('I_total'), function(x){x == main_fig_I_total})) %>%
                                                     filter(across(any_of('gamma'), function(x){x == main_fig_gamma})),
                                                   y_var = 'fraction_biggest_clone', color_var = 'mutation_rate', facet_vars = NULL) +
    xlab('') + theme(legend.position = 'none') + ylim(0,1) +
    mutations_color_scale
  
  freq_corr <-  base_plotting_function(summary_tibble = summary_list$summary_pairwise_correlations %>%
                                         filter(across(any_of('s'), function(x){x == main_fig_s})) %>%
                                         filter(across(any_of('beta'), function(x){x == main_fig_beta})) %>%
                                         filter(across(any_of('I_total'), function(x){x == main_fig_I_total})) %>%
                                         filter(across(any_of('gamma'), function(x){x == main_fig_gamma})) %>%
                                         filter(method == 'pearson'),
                                       y_var = 'freq_correlation', color_var = 'mutation_rate', facet_vars = NULL)  +
    xlab('') + theme(legend.position = 'none') + ylim(-0.2,1) +
    geom_hline(yintercept = 0, linetype = 2) +
    mutations_color_scale
  
  freq_ratio_corr <- base_plotting_function(summary_tibble = summary_list$summary_pairwise_correlations %>%
                                              filter(across(any_of('s'), function(x){x == main_fig_s})) %>%
                                              filter(across(any_of('beta'), function(x){x == main_fig_beta})) %>%
                                              filter(across(any_of('I_total'), function(x){x == main_fig_I_total})) %>%
                                              filter(across(any_of('gamma'), function(x){x == main_fig_gamma})) %>%
                                              filter(method == 'pearson'),
                                            y_var = 'freq_ratio_correlation', color_var = 'mutation_rate', facet_vars = NULL)  +
    xlab('Time (days)') + theme(legend.position = 'none') +
    ylim(-0.2,1) +
    geom_hline(yintercept = 0, linetype = 2) +
    mutations_color_scale
  
  return(list(fraction_biggest_clone = fraction_biggest_clone, freq_corr = freq_corr, freq_ratio_corr = freq_ratio_corr))
}

mutations_color_scale <- scale_color_manual(name = 'Mutation rate (affinity-changing\nmutations per B cell per division)',
                                            values = brewer.pal(n = 9, "BuPu")[c(3,6,9)])



# ======================== MAKING PLOTS ========================
# Load simulation summaries
scenarios <- c('neutral_scenario', 'high_affinity_scenario', 'high_mutation_scenario', 'neutral_uniform_freqs_scenario',
               'test_scenario_neutral_highly_dominant')

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


color_legend <- get_legend(base_plotting_function(summary_tibble = neutral_scenario_summary$summary_GC_stats %>%
                                                    filter(beta == main_fig_beta, I_total == main_fig_I_total),
                                                  y_var = 'fraction_biggest_clone', color_var = 'mutation_rate', facet_vars = NULL) +
                             mutations_color_scale)
  

correlations_fig_neutral <- generate_freq_cor_panels(summary_list = neutral_scenario_summary, main_fig_s = main_fig_s, main_fig_beta = main_fig_beta,
                                                     main_fig_I_total = main_fig_I_total, main_fig_gamma = main_fig_gamma)

correlation_fig_high_affinity <- generate_freq_cor_panels(summary_list = high_affinity_scenario_summary, main_fig_s = main_fig_s,
                                                          main_fig_beta = main_fig_beta, main_fig_I_total = main_fig_I_total, 
                                                          main_fig_gamma = main_fig_gamma)
correlation_fig_high_mutation <- generate_freq_cor_panels(summary_list = high_mutation_scenario_summary, main_fig_s = main_fig_s,
                                                            main_fig_beta = main_fig_beta, main_fig_I_total = main_fig_I_total, 
                                                            main_fig_gamma = main_fig_gamma)

# Adjustments of labels
correlations_fig_neutral$fraction_biggest_clone <- correlations_fig_neutral$fraction_biggest_clone +
  ylab('Frequency of biggest\nlineage within GCs') +
  ggtitle('Affinity distribution and mutation\nrate are the same across alleles') +
  theme(plot.title = element_text(hjust = 0.5, size = 12))
correlations_fig_neutral$freq_corr <- correlations_fig_neutral$freq_corr + ylab('Pairwise correlation\nin allele frequencies')
correlations_fig_neutral$freq_ratio_corr <- correlations_fig_neutral$freq_ratio_corr + ylab('Pairwise correlation in\nexperienced-to-naive ratios') + xlab('')

correlation_fig_high_affinity$fraction_biggest_clone <- correlation_fig_high_affinity$fraction_biggest_clone + ylab('') +
  ggtitle('Some alleles tend to encode\n receptors with higher affinity') +
  theme(plot.title = element_text(hjust = 0.5, size = 12))
correlation_fig_high_affinity$freq_corr <- correlation_fig_high_affinity$freq_corr + ylab('')
correlation_fig_high_affinity$freq_ratio_corr <- correlation_fig_high_affinity$freq_ratio_corr + ylab('')

correlation_fig_high_mutation$fraction_biggest_clone <- correlation_fig_high_mutation$fraction_biggest_clone + ylab('') +
  ggtitle('Some alleles tend to encode\n receptors with higher mutation rate')  +
  theme(plot.title = element_text(hjust = 0.5, size = 12))
correlation_fig_high_mutation$freq_corr <- correlation_fig_high_mutation$freq_corr + ylab('')
correlation_fig_high_mutation$freq_ratio_corr <- correlation_fig_high_mutation$freq_ratio_corr + ylab('') + xlab('')



correlations_fig <- plot_grid(
  plot_grid(plotlist = correlations_fig_neutral, nrow = 3),
  plot_grid(plotlist = correlation_fig_high_affinity, nrow = 3),
  plot_grid(plotlist = correlation_fig_high_mutation, nrow = 3),
  nrow = 1)


correlations_fig <- plot_grid(
  correlations_fig,
  plot_grid(NULL, color_legend, nrow =1, rel_widths = c(2,5)),
  nrow = 2,
  rel_heights = c(15,1)
)

save_plot('../figures/simulations/simulated_correlations.pdf',
          correlations_fig,
          base_width = 11, base_height = 11)

# Main text figure showing % GCs occupied by high-affinity alleles in the high-affinity scenario

total_GCs_across_individuals <- high_affinity_scenario_summary$n_individuals * high_affinity_scenario_summary$n_GCs_per_individual


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

# Supplementary figures
neutral_scenario_supp_fig <- base_plotting_function(summary_tibble = neutral_scenario_summary$summary_pairwise_correlations %>%
                                                      filter(method == 'pearson'), y_var = 'freq_correlation',
                                                    color_var = 'mutation_rate', facet_vars = c('I_total', 'beta')) +
  ylim(-0.2,1) + geom_hline(yintercept = 0, linetype = 2) +
  mutations_color_scale +
  ylab('Pairwise correlation in allele frequencies')

save_plot('../figures/simulations/neutral_scenario_supp_fig.pdf',
          neutral_scenario_supp_fig,
          base_width = 8, base_height = 9)

high_affinity_scenario_supp_fig <- base_plotting_function(summary_tibble = high_affinity_scenario_summary$summary_pairwise_correlations %>%
                                                            filter(method == 'pearson'), y_var = 'freq_ratio_correlation',
                                                          color_var = 'mutation_rate', facet_vars = c('s', 'beta')) +
  ylim(-0.2,1) + geom_hline(yintercept = 0, linetype = 2) +
  mutations_color_scale +
  ylab('Pairwise correlation in experienced-to-naive frequency ratios')

save_plot('../figures/simulations/high_affinity_scenario_supp_fig.pdf',
          high_affinity_scenario_supp_fig,
          base_width = 8, base_height = 9)

high_mutation_scenario_supp_fig <-  base_plotting_function(summary_tibble = high_mutation_scenario_summary$summary_pairwise_correlations %>%
                                                             filter(method == 'pearson'), y_var = 'freq_ratio_correlation',
                                                           color_var = 'mutation_rate', facet_vars = c('gamma', 'beta')) +
  ylim(-0.2,1) + geom_hline(yintercept = 0, linetype = 2) +
  mutations_color_scale +
  ylab('Pairwise correlation in experienced-to-naive frequency ratios')

save_plot('../figures/simulations/high_mutation_scenario_supp_fig.pdf',
          high_mutation_scenario_supp_fig,
          base_width = 8, base_height = 9)


# FOR HIGH MUTATION ALLELES
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



# Spearman correlation in simulations fig.


spearman_freq_ratio_corr_neutral <- base_plotting_function(neutral_scenario_summary$summary_pairwise_correlations %>%
                                                             filter(method == 'spearman', beta == main_fig_beta, I_total == main_fig_I_total),
                                                           y_var = 'freq_ratio_correlation', 
                                                           color_var = 'mutation_rate') +
  ylim(-0.2,1) + geom_hline(yintercept = 0, linetype = 2) +
  mutations_color_scale +
  ylab('Pairwise spearman correlation\nin experienced-to-naive ratios') +
  theme(legend.position = c(0.1,0.8),
        plot.title = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  ggtitle('Alleles have identical affinity distributions\nand mutation rates but different naive frequencies')


spearman_freq_ratio_corr_neutral_uniform_freqs <- base_plotting_function(neutral_uniform_freqs_scenario_summary$summary_pairwise_correlations %>%
                         filter(method == 'spearman', beta == main_fig_beta, I_total == main_fig_I_total),
                       y_var = 'freq_ratio_correlation', 
                       color_var = 'mutation_rate') +
  ylim(-0.2,1) + geom_hline(yintercept = 0, linetype = 2) +
  mutations_color_scale +
  ylab('') +
  theme(legend.position = 'none')  +
  ggtitle('Alleles have identical affinity distributions,\n mutation rates and naive frequencies') +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

freq_ratio_spearman_correlations <- plot_grid(spearman_freq_ratio_corr_neutral,
                                              spearman_freq_ratio_corr_neutral_uniform_freqs,
                                              nrow = 1)

save_plot('../figures/simulations/freq_ratio_spearman_correlations.pdf',
          freq_ratio_spearman_correlations,
          base_width = 12, base_height = 5)

# Spearman correlation in simulations fig., but now with allele freqs
# (Duplicated code; revise)


spearman_freq_corr_neutral <- base_plotting_function(neutral_scenario_summary$summary_pairwise_correlations %>%
                                                             filter(method == 'spearman', beta == main_fig_beta, I_total == main_fig_I_total),
                                                           y_var = 'freq_correlation', 
                                                           color_var = 'mutation_rate') +
  ylim(-0.2,1) + geom_hline(yintercept = 0, linetype = 2) +
  mutations_color_scale +
  ylab('Pairwise spearman correlation\nin allele frequencies') +
  theme(legend.position = c(0.5,0.85),
        plot.title = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  ggtitle('Alleles have identical affinity distributions\nand mutation rates but different naive frequencies')


spearman_freq_corr_neutral_uniform_freqs <- base_plotting_function(neutral_uniform_freqs_scenario_summary$summary_pairwise_correlations %>%
                                                                           filter(method == 'spearman', beta == main_fig_beta, I_total == main_fig_I_total),
                                                                         y_var = 'freq_correlation', 
                                                                         color_var = 'mutation_rate') +
  ylim(-0.2,1) + geom_hline(yintercept = 0, linetype = 2) +
  mutations_color_scale +
  ylab('') +
  theme(legend.position = 'none')  +
  ggtitle('Alleles have identical affinity distributions,\n mutation rates and naive frequencies') +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

freq_spearman_correlations <- plot_grid(spearman_freq_corr_neutral,
                                              spearman_freq_corr_neutral_uniform_freqs,
                                              nrow = 1)

save_plot('../figures/simulations/freq_spearman_correlations.pdf',
          freq_spearman_correlations,
          base_width = 12, base_height = 5)

  



