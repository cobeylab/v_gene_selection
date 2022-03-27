library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(cowplot)
library(viridis)
library(truncnorm)
library(stringr)
theme_set(theme_cowplot())

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
    xlab('Time') +
    theme(legend.position = 'top') +
    facet_command
  
  return(pl)
}


### HAVE TO ADD AN I_TOTAL FILTER FOR NEUTRAL SCENARIO ####


# Load simulatio summaries
scenarios <- c('neutral_scenario', 'high_affinity_scenario')

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
main_fig_s <- 2
main_fig_beta <- 4
sigma_r <- 1 # READ THIS OFF OF IMPORTED COMBINED PARAMETER FILES, ADD CHECK IT IS A SINGLE VALUE





color_legend <- get_legend(base_plotting_function(summary_tibble = neutral_scenario_summary$summary_GC_stats %>%
                                                    filter(beta == main_fig_beta),
                                                  y_var = 'fraction_biggest_clone', color_var = 'mutation_rate', facet_vars = NULL) +
                             scale_color_discrete(name = 'Mutation rate (affinity-changing\nmutations per B cell per division)'))
  



generate_freq_cor_panels <- function(summary_list, main_fig_s, main_fig_beta){
  
  fraction_biggest_clone <- base_plotting_function(summary_tibble = summary_list$summary_GC_stats %>%
                                                           filter(across(any_of('s'), function(x){x == main_fig_s})) %>%
                                                           filter(across(any_of('beta'), function(x){x == main_fig_beta})),
                                                         y_var = 'fraction_biggest_clone', color_var = 'mutation_rate', facet_vars = NULL) +
    xlab('') + theme(legend.position = 'none') + ylim(0,1)
  
  freq_corr <-  base_plotting_function(summary_tibble = summary_list$summary_pairwise_correlations %>%
                                               filter(across(any_of('s'), function(x){x == main_fig_s})) %>%
                                               filter(across(any_of('beta'), function(x){x == main_fig_beta})) %>%
                                               filter(method == 'pearson'),
                                             y_var = 'freq_correlation', color_var = 'mutation_rate', facet_vars = NULL)  +
    xlab('') + theme(legend.position = 'none') + ylim(-0.2,1) +
    geom_hline(yintercept = 0, linetype = 2)
  
  freq_ratio_corr <- base_plotting_function(summary_tibble = summary_list$summary_pairwise_correlations %>%
                                                    filter(across(any_of('s'), function(x){x == main_fig_s})) %>%
                                                    filter(across(any_of('beta'), function(x){x == main_fig_beta})) %>%
                                                    filter(method == 'pearson'),
                                                  y_var = 'freq_ratio_correlation', color_var = 'mutation_rate', facet_vars = NULL)  +
    xlab('Time') + theme(legend.position = 'none') +
    ylim(-0.2,1) +
    geom_hline(yintercept = 0, linetype = 2)
  
  return(list(fraction_biggest_clone = fraction_biggest_clone, freq_corr = freq_corr, freq_ratio_corr = freq_ratio_corr))
}


correlations_fig_neutral <- generate_freq_cor_panels(summary_list = neutral_scenario_summary, main_fig_s = main_fig_s, main_fig_beta = main_fig_beta)
correlation_fig_high_affinity <- generate_freq_cor_panels(summary_list = high_affinity_scenario_summary, main_fig_s = main_fig_s,
                                                          main_fig_beta = main_fig_beta)

# Adjustments of labels
correlations_fig_neutral$fraction_biggest_clone <- correlations_fig_neutral$fraction_biggest_clone +
  ylab('Frequency of biggest\nclone within GCs') +
  ggtitle('All alleles contribute\nequally to fitness') +
  theme(plot.title = element_text(hjust = 0.5, size = 14))
correlations_fig_neutral$freq_corr <- correlations_fig_neutral$freq_corr + ylab('Pairwise correlation\nin allele frequencies')
correlations_fig_neutral$freq_ratio_corr <- correlations_fig_neutral$freq_ratio_corr + ylab('Pairwise correlation in\nexperienced-to-naive ratios')

correlation_fig_high_affinity$fraction_biggest_clone <- correlation_fig_high_affinity$fraction_biggest_clone + ylab('') +
  ggtitle('Some alleles tend to encode\n receptors with higher affinity') +
  theme(plot.title = element_text(hjust = 0.5, size = 14))

correlation_fig_high_affinity$freq_corr <- correlation_fig_high_affinity$freq_corr + ylab('')
correlation_fig_high_affinity$freq_ratio_corr <- correlation_fig_high_affinity$freq_ratio_corr + ylab('')


correlations_fig <- plot_grid(
  plot_grid(plotlist = correlations_fig_neutral, nrow = 3),
  plot_grid(plotlist = correlation_fig_high_affinity, nrow = 3))

correlations_fig <- plot_grid(
  correlations_fig,
  plot_grid(NULL, color_legend, nrow =1, rel_widths = c(1,5)),
  nrow = 2,
  rel_heights = c(15,1)
)

save_plot('../figures/simulated_correlations.pdf',
          correlations_fig,
          base_width = 8, base_height = 10)


freq_correlations_affinity_variation <-
    plot_grid(
      base_plotting_function(summary_tibble = summary_affinity_variation_scenario$summary_GC_stats %>%
                               filter(s == main_fig_s, beta == main_fig_beta),
                             y_var = 'fraction_biggest_clone', color_var = 'mutation_rate', facet_vars = NULL) +
        theme(legend.position = 'none')  + xlab(''),
      base_plotting_function(summary_tibble = summary_affinity_variation_scenario$summary_pairwise_correlations %>%
                               filter(method == 'pearson', s == main_fig_s, beta == main_fig_beta),
                             y_var = 'freq_correlation', color_var = 'mutation_rate', facet_vars = NULL) +
        geom_hline(yintercept = 0, linetype = 2) +
        theme(legend.position = 'none') + xlab(''),
      base_plotting_function(summary_tibble = summary_affinity_variation_scenario$summary_pairwise_correlations %>%
                               filter(method == 'pearson', s == main_fig_s, beta == main_fig_beta),
                             y_var = 'freq_ratio_correlation', color_var = 'mutation_rate', facet_vars = NULL) +
        geom_hline(yintercept = 0, linetype = 2) +
        theme(legend.position = 'none'),
      nrow = 3
    )
  
# Main text figure showing % GCs occupied by high-affinity alleles

total_GCs_across_individuals = length(unique(simulations$individual)) * length(unique(simulations$GC))

combined_freq_of_high_avg_alleles_in_GC_pl <- summary_affinity_variation_scenario$combined_freq_of_high_avg_alleles_in_GCs %>%
  filter(beta == main_fig_beta, s == main_fig_s) %>%
  filter(t %in% c(10, max(t))) %>%
  mutate(s = paste0('allele advantage = ', s),
         t = paste0(t, ' days')) %>%
  ggplot(aes(x = combined_freq_high_avg)) +
  geom_histogram(bins = 20) + 
  facet_grid(reformulate('t', color_var)) +
  scale_y_continuous(labels = function(x){round(x/total_GCs_across_individuals, 2)}) +
  ylab('Fraction of GCs') +
  xlab('Combined frequency of high-affinity alleles within GC') +
  geom_vline(xintercept = 0.5, alpha = 0.3, linetype = 2) 

# Panel showing shape of distributions

bind_rows(tibble(x = seq(0,15,0.01), type = 'low_avg') %>%
            mutate(density = dtruncnorm(x, a = 0, b = Inf, mean = 1, sd = sigma_r)),
          tibble(x = seq(0,15,0.01), type = 'high_avg_same_sd') %>%
            mutate(density = dtruncnorm(x, a = 0, b = Inf, mean = 1 + main_fig_s, sd = sigma_r))) %>%
  ggplot(aes(x = x, y = density, color = type)) +
  geom_line(aes(group = type)) +
  theme(legend.position = 'top')



# =========== OLD PLOT CODE ============



color_legend <- get_legend(base_plotting_function(summary_pairwise_correlations %>% filter(method == 'pearson'),
                                                  y_var = 'freq_correlation', color_var = color_var,
                                                  primary_facet_var = primary_facet_var, secondary_facet_var = secondary_facet_var,
                                                  secondary_facet_value = secondary_facet_value, plot_secondary_facet = F))

pairwise_corr_freqs_pl <- base_plotting_function(summary_pairwise_correlations %>% filter(method == 'pearson'),
                                                 y_var = 'freq_correlation', color_var = color_var, primary_facet_var = primary_facet_var) +
  ylab('Pairwise correlation\nin allele frequencies') +
  geom_hline(yintercept = 0, linetype =2) +
  theme(strip.background = element_blank(), strip.text = element_blank(),
        legend.position = 'none') +
  xlab('')

pairwise_corr_freq_ratios_pl <- base_plotting_function(summary_pairwise_correlations %>% filter(method == 'pearson'),
                                                       y_var = 'freq_ratio_correlation', color_var = color_var,
                                                       primary_facet_var = primary_facet_var) +
  ylab('Pairwise correlation in\nexperienced-to-naive ratios') +
  geom_hline(yintercept = 0, linetype =2) +
  #ylim(-1,1) +
  theme(strip.background = element_blank(), strip.text = element_blank(),
        legend.position = 'none') +
  xlab('')


n_alleles_in_rep_pl <- base_plotting_function(summary_repertoire_allele_diversity,
                                              y_var = 'n_alleles_in_experienced_repertoire', color_var = color_var,
                                              primary_facet_var = primary_facet_var) +
  ylab('Number of V alleles\nin the response') +
  theme(strip.background = element_blank(), strip.text = element_blank(),
        legend.position = 'none') +
  xlab('Time')


freq_dominant_clone_pl <- base_plotting_function(summary_GC_stats, y_var = 'fraction_biggest_clone', color_var = color_var,
                                                 primary_facet_var = primary_facet_var) +
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
                                          primary_facet_var = primary_facet_var) +
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
    facet_grid(reformulate(color_var, primary_facet_var)) +
    geom_hline(yintercept = 0.5, linetype = 2) +
    xlab('Time') +
    ylim(0,1) +
    ylab('Fraction of individuals where allele increased in frequency') +
    theme(panel.border = element_rect(colour = 'black', linetype = 1))
  
  save_plot(paste0(results_directory, basename(results_directory), '_fraction_individuals_high_avg_alleles_increasing.pdf'),
            fraction_increasing_pl,
            base_height = 12, base_width = 10)
  
  
  
  
  
  total_GCs_across_individuals = length(unique(simulations$individual)) * length(unique(simulations$GC))
  
  
  selected_s_value <- 2
  combined_freq_of_high_avg_alleles_in_GC_pl <- combined_freq_of_high_avg_alleles_in_GCs %>%
    filter(s == selected_s_value) %>%
    filter(t %in% c(10, max(t))) %>%
    mutate(s = paste0('allele advantage = ', s),
           t = paste0(t, ' days')) %>%
    ggplot(aes(x = combined_freq_high_avg)) +
    geom_histogram(bins = 20) + 
    facet_grid(reformulate('t', color_var)) +
    scale_y_continuous(labels = function(x){round(x/total_GCs_across_individuals, 2)}) +
    ylab('Fraction of GCs') +
    xlab('Combined frequency of high-affinity alleles within GC') +
    geom_vline(xintercept = 0.5, alpha = 0.3, linetype = 2) +
    ggtitle(paste0('S = ', selected_s_value))
  
  save_plot(paste0(results_directory, basename(results_directory), '_combined_freq_of_high_avg_alleles_in_GCS.pdf'),
            combined_freq_of_high_avg_alleles_in_GC_pl,
            base_height = 10, base_width = 10)
  
  
  # Cool to look at but hard to see patterns
  # combined_freq_of_high_avg_alleles_in_GCs %>%
  #   mutate(plotting_group = paste(individual, GC, sep = '_')) %>%
  #   group_by(mutation_rate, s, plotting_group) %>%
  #   mutate(dominated_by_high_avg = combined_freq_high_avg[t == max(t)] > 0.5) %>%
  #   ungroup() %>%
  #   ggplot(aes(x = t, y = combined_freq_high_avg, color = dominated_by_high_avg)) +
  #   geom_line(aes(group = plotting_group), alpha = 0.5) +
  #   facet_grid(mutation_rate~s)
  
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





