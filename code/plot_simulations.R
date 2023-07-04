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
base_plotting_function <- function(summary_tibble, y_var, color_var, facet_vars = NULL, hline = NULL){
  
  plotting_tibble <- summary_tibble
  
  if('mutation_rate' %in% names(summary_tibble)){
    plotting_tibble <- summary_tibble %>%
      mutate(mutation_rate = as.factor(as.character(mutation_rate)))
  }
  
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
    ggplot(aes_string(x = 't', y = paste0(y_var, '_median'), color = color_var, group = color_var))
  
  if(!is.null(hline)){
    pl <- pl + geom_hline(yintercept = 0, linetype = 2)
  }
  pl <- pl  +
    geom_line(size = 1.5) +
    geom_linerange(aes_string(ymin = paste0(y_var, '_lowerq'), ymax = paste0(y_var, '_upperq')),
                   alpha = 0.5, size = 1.2) +
    geom_point(size = 3) +
    xlab('Time (days)') +
    theme(legend.position = 'top') +
    facet_command
  
  return(pl)
}


# Function for making main fig showing simulated frequency and freq. deviation correlations in the equivalent-alleles scenario
make_neutral_scenario_fig <- function(neutral_scenario_summary, beta, I_total, mutation_rate){
  
  subset_tibble <- neutral_scenario_summary$summary_pairwise_correlations %>%
    filter(beta == !!beta, I_total == !!I_total, mutation_rate == !!mutation_rate) %>%
    filter(method == 'pearson') %>%
    mutate(mutation_rate = paste0('Mutation rate = ', mutation_rate))
    
  theme_specs <- theme(legend.position = 'none', plot.title = element_text(size = 12,hjust = 0.5)) 
  
  freq_corr <-  base_plotting_function(summary_tibble = subset_tibble,
                                       y_var = 'freq_correlation', color_var = 'nGCs',
                                       facet_vars = NULL, hline = 0) + ylim(-0.2,1) +
    ylab('Pairwise correlation between individuals') +
    facet_wrap("mutation_rate") +
    theme_specs +
    nGCs_color_scale +
    ggtitle('\n\n\nGermline allele frequencies\nin the response')
  
  freq_ratio_corr <- base_plotting_function(summary_tibble = subset_tibble,
                                            y_var = 'freq_ratio_correlation', color_var = 'nGCs',
                                            facet_vars = NULL, hline = 0) +
    ylab('') + theme(legend.position = 'none') +
    facet_wrap("mutation_rate") + 
    ylim(-0.2,1) +
    theme_specs +
    nGCs_color_scale +
    ggtitle('\n\n\nExperienced-to-naive\nratios')
  
  
  return(plot_grid(freq_corr, freq_ratio_corr, nrow = 1))
}

nGCs_color_scale <- scale_color_manual(name = 'Number of GCs per individual',
                                            values = brewer.pal(n = 9, "BuPu")[c(3,5,7,9)])

mutation_rate_color_scale <- scale_color_manual(name = 'Mutations per B cell per division',
                                                values = brewer.pal(n = 9, "OrRd")[c(3,6,9)])


# ======================== LOADING SIMULATION SUMMARIES ========================
# Load simulation summaries
scenarios <- c('neutral_scenario', 'high_affinity_scenario', 'high_mutation_scenario',
               'high_affinity_scenario_lower_sigma_r', 'high_affinity_scenario_20inds',
               'high_affinity_scenario_100inds', 'high_affinity_scenario_1000inds',
               'high_affinity_scenario_1allele', 'high_affinity_scenario_10alleles')

# scenarios <- c('high_affinity_scenario_20inds', 'high_affinity_scenario_100inds', 'high_affinity_scenario_1000inds')

model_parameters <- bind_rows(
  lapply(as.list(scenarios),
         FUN = function(scen){
           comb_pars_path <- paste0('../results/simulations/', scen, '/combined_model_parameters.csv')
           if(file.exists(comb_pars_path)){
             output <- read_csv(comb_pars_path) %>% mutate(scenario = str_remove(scen,'_scenario'))
           }else{
             output <- c()
           }
           return(output)
         })
) %>% select(scenario, everything())

parameter_names <- names(model_parameters)

for(scen in scenarios){
  summary_path <- paste0('../results/simulations/', scen, '/', scen, '_summary.RData')
  if(file.exists(summary_path)){
    load(summary_path)
  }
}


# Main text figure with behavior of frequency correlations
# Across scenarios, where applicable, use simulations with these values
default_s <- 1.5
default_gamma <- 6


# ==================================================== Power analyses plots ===========================================
nGC_color_legend <- get_legend(base_plotting_function(summary_tibble = high_affinity_scenario_20inds_summary$summary_pairwise_correlations,
                                                      y_var = 'freq_correlation', color_var = 'nGCs', facet_vars = NULL) +
                                 nGCs_color_scale)

power_analysis_freq_cors <- plot_grid(
  base_plotting_function(summary_tibble = high_affinity_scenario_20inds_summary$summary_pairwise_correlations %>%
                           filter(method == 'pearson'), y_var = 'freq_correlation', color_var = 'nGCs', facet_vars = NULL) +
    nGCs_color_scale + ylab('Pairwise correlation\nin allele frequencies') + xlab('') + theme(legend.position = 'none') ,
  base_plotting_function(summary_tibble = high_affinity_scenario_100inds_summary$summary_pairwise_correlations %>%
                           filter(method == 'pearson'), y_var = 'freq_correlation', color_var = 'nGCs', facet_vars = NULL) +
    nGCs_color_scale + ylab('') + xlab('') + theme(legend.position = 'none'),
  base_plotting_function(summary_tibble = high_affinity_scenario_1000inds_summary$summary_pairwise_correlations %>%
                           filter(method == 'pearson'), y_var = 'freq_correlation', color_var = 'nGCs', facet_vars = NULL) +
    nGCs_color_scale + ylab('') + xlab('') + theme(legend.position = 'none'),
  nrow = 1, 
  labels = c('20 inds', '100 inds', '1000 inds'),
  label_x = 0.5
)

power_analysis_freq_ratio_cors <- plot_grid(
  base_plotting_function(summary_tibble = high_affinity_scenario_20inds_summary$summary_pairwise_correlations %>%
                           filter(method == 'pearson'), y_var = 'freq_ratio_correlation', color_var = 'nGCs', facet_vars = NULL) +
    nGCs_color_scale + ylab('Pairwise correlation\nin exp/naive ratios') + xlab('') + theme(legend.position = 'none') ,
  base_plotting_function(summary_tibble = high_affinity_scenario_100inds_summary$summary_pairwise_correlations %>%
                           filter(method == 'pearson'), y_var = 'freq_ratio_correlation', color_var = 'nGCs', facet_vars = NULL) +
    nGCs_color_scale + ylab('') + theme(legend.position = 'none'),
  base_plotting_function(summary_tibble = high_affinity_scenario_1000inds_summary$summary_pairwise_correlations %>%
                           filter(method == 'pearson'), y_var = 'freq_ratio_correlation', color_var = 'nGCs', facet_vars = NULL) +
    nGCs_color_scale + ylab('') + xlab('') + theme(legend.position = 'none'),
  nrow = 1, 
  label_x = 0.5
)
      
sim_power_analysis <- plot_grid(plot_grid(NULL,nGC_color_legend, nrow = 1, rel_widths = c(1,3)),
                                power_analysis_freq_cors, power_analysis_freq_ratio_cors,
                                rel_heights = c(1,10,10), nrow = 3)

save_plot('../figures/simulations/sim_power_analysis.pdf', sim_power_analysis, base_height = 6, base_width = 12)


# =========================================== Equivalent alleles scenario plots ===========================================
default_beta <- 2
default_I_total <- 200
default_mutation_rate <- 0.01

# Fraction top lineage in GCs (neutral scenario)
fraction_top_lineage_in_GCs <- base_plotting_function(summary_tibble = neutral_scenario_summary$summary_GC_stats %>%
                                                        filter(beta == default_beta),
                                                      y_var = 'fraction_biggest_clone', color_var = 'mutation_rate', facet_vars = NULL) +
  mutation_rate_color_scale +
  ylab('Frequency of biggest lineage within GCs')


mutation_color_legend <- get_legend(fraction_top_lineage_in_GCs)


correlations_fig_neutral <- make_neutral_scenario_fig(neutral_scenario_summary = neutral_scenario_summary, beta = default_beta,
                                                     I_total = default_I_total, mutation_rate = default_mutation_rate)

neutral_scenario_main_fig <- plot_grid(plot_grid(NULL,
                                                fraction_top_lineage_in_GCs + theme(legend.position = 'None'),
                                                mutation_color_legend, 
                                                nrow = 3, rel_heights = c(5,10,1),
                                                labels = c('','A)',''),
                                                hjust = -0.05, scale = 0.95,
                                                vjust = -6,
                                                label_size = 15),
                                      plot_grid(
                                        correlations_fig_neutral,
                                        plot_grid(NULL, nGC_color_legend, nrow = 1, rel_widths = c(2,5)),
                                        nrow = 2,
                                        rel_heights = c(19,1),
                                        labels = c('B)',''),
                                        hjust = -0.05, scale = 0.95,
                                        vjust = 6.5,
                                        label_size = 15
                                       ),
                                      nrow = 1, rel_widths = c(1,2)
                                      )

save_plot('../figures/simulations/neutral_scenario_main_fig.pdf',
          neutral_scenario_main_fig,
          base_width = 14, base_height = 6)

# =========================================== High affinity scenario plots ===========================================

# Main fig correlation panel for high-affinity scenario
# (Assume 15 germinal centers per individual)
high_affinity_scenario_main_fig <- base_plotting_function(summary_tibble = high_affinity_scenario_summary$summary_pairwise_correlations %>%
                                                            # Excluding one of the betas to keep main text fig small
                                                            filter(method == 'pearson', nGCs == 15, beta != 3),
                                                          y_var = 'freq_ratio_correlation',
                                                          color_var = 'mutation_rate', facet_vars = c('s', 'beta')) +
  ylim(-0.2,1) + geom_hline(yintercept = 0, linetype = 2) +
  mutation_rate_color_scale +
  ylab('Pairwise correlation in experienced-to-naive frequency ratios') +
  theme(legend.position = 'bottom')

save_plot('../figures/simulations/high_affinity_scenario_main_fig.pdf',
          high_affinity_scenario_main_fig,
          base_width = 9.5, base_height = 9)

# Supplemental fig. showing correlation in freq ratios in affinity scenario (15 GCs)
high_affinity_supp_freq_corr <- base_plotting_function(summary_tibble = high_affinity_scenario_summary$summary_pairwise_correlations %>%
                                                            filter(method == 'pearson', nGCs == 15),
                                                          y_var = 'freq_correlation',
                                                          color_var = 'mutation_rate', facet_vars = c('s', 'beta')) +
  ylim(-0.2,1) + geom_hline(yintercept = 0, linetype = 2) +
  mutation_rate_color_scale +
  ylab('Pairwise correlation in germline allele frequencies') +
  theme(legend.position = 'bottom')

save_plot('../figures/simulations/high_affinity_supp_freq_corr.pdf',
          high_affinity_supp_freq_corr,
          base_width = 8.5, base_height = 9)


# Supp fig showing % GCs occupied by high-affinity alleles in the high-affinity scenario
stopifnot(all(high_affinity_scenario_summary$nGC_values == high_mutation_scenario_summary$nGC_values))
total_GCs_across_individuals <- high_affinity_scenario_summary$n_individuals * max(high_affinity_scenario_summary$nGC_values)

combined_freq_of_high_avg_alleles_in_GC_pl <- high_affinity_scenario_summary$combined_freq_of_high_avg_alleles_in_GCs %>%
  filter(beta == 4, s == 2) %>%
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

# Supp fig. showing freq ratio correlations in high-affinity scenario with 1 or 30 GCs
high_affinity_freq_ratio_corr_1GC <- base_plotting_function(summary_tibble = high_affinity_scenario_summary$summary_pairwise_correlations %>%
                                                         filter(method == 'pearson', nGCs == 1),
                                                       y_var = 'freq_ratio_correlation',
                                                       color_var = 'mutation_rate', facet_vars = c('s', 'beta')) +
  ylim(-0.2,1) + geom_hline(yintercept = 0, linetype = 2) +
  mutation_rate_color_scale +
  ylab('Pairwise correlation in experienced-to-naive frequency ratios') +
  theme(legend.position = 'None')

high_affinity_freq_ratio_corr_30GCs <- base_plotting_function(summary_tibble = high_affinity_scenario_summary$summary_pairwise_correlations %>%
                                                        filter(method == 'pearson', nGCs == 30),
                                                      y_var = 'freq_ratio_correlation',
                                                      color_var = 'mutation_rate', facet_vars = c('s', 'beta')) +
  ylim(-0.2,1) + geom_hline(yintercept = 0, linetype = 2) +
  mutation_rate_color_scale +
  ylab('') +
  theme(legend.position = 'None')


high_affinity_supp_other_nGCs <- plot_grid(plot_grid(high_affinity_freq_ratio_corr_1GC + ggtitle('1 GC per individual'),
                    high_affinity_freq_ratio_corr_30GCs + ggtitle('30 GCs per individual'),
                    nrow = 1),
          plot_grid(NULL,mutation_color_legend, nrow = 1, rel_widths = c(4,8)),
          nrow = 2, rel_heights = c(19,1))

save_plot('../figures/simulations/high_affinity_supp_other_nGCs.pdf',
                    high_affinity_supp_other_nGCs,
                    base_width = 14, base_height = 11)

# Supp fig showing high affinity scenario freq ratio correlations with 1 or 10 high-affinity alleles

high_affinity_freq_ratio_1allele <-  base_plotting_function(summary_tibble = high_affinity_scenario_1allele_summary$summary_pairwise_correlations %>%
                                                              filter(method == 'pearson', nGCs == 15),
                                                            y_var = 'freq_ratio_correlation',
                                                            color_var = 'mutation_rate', facet_vars = c('s', 'beta')) +
  ylim(-0.2,1) + geom_hline(yintercept = 0, linetype = 2) +
  mutation_rate_color_scale +
  ylab('Pairwise correlation in experienced-to-naive frequency ratios') +
  theme(legend.position = 'None')

high_affinity_freq_ratio_10alleles <-  base_plotting_function(summary_tibble = high_affinity_scenario_10alleles_summary$summary_pairwise_correlations %>%
                                                              filter(method == 'pearson', nGCs == 15),
                                                            y_var = 'freq_ratio_correlation',
                                                            color_var = 'mutation_rate', facet_vars = c('s', 'beta')) +
  ylim(-0.2,1) + geom_hline(yintercept = 0, linetype = 2) +
  mutation_rate_color_scale +
  ylab('Pairwise correlation in experienced-to-naive frequency ratios') +
  theme(legend.position = 'None')

high_affinity_supp_other_nAlleles <- plot_grid(plot_grid(high_affinity_freq_ratio_1allele + ggtitle('1 high-affinity allele'),
                                                         high_affinity_freq_ratio_10alleles + ggtitle('10 high-affinity allele'),
                                                     nrow = 1),
                                           plot_grid(NULL,mutation_color_legend, nrow = 1, rel_widths = c(4,8)),
                                           nrow = 2, rel_heights = c(19,1))

save_plot('../figures/simulations/high_affinity_supp_other_nAlleles.pdf',
          high_affinity_supp_other_nAlleles,
          base_width = 14, base_height = 11)

# Supplementary fig showing high-affinity scenario with lower sigma r
high_affinity_supp_lower_sigma_r <- base_plotting_function(summary_tibble = high_affinity_scenario_lower_sigma_r_summary$summary_pairwise_correlations %>%
                                                             filter(method == 'pearson', nGCs == 15),
                                                           y_var = 'freq_ratio_correlation',
                                                           color_var = 'mutation_rate', facet_vars = c('s', 'beta')) +
  ylim(-0.2,1) + geom_hline(yintercept = 0, linetype = 2) +
  mutation_rate_color_scale +
  ylab('Pairwise correlation in experienced-to-naive frequency ratios') +
  theme(legend.position = 'bottom')

save_plot('../figures/simulations/high_affinity_supp_lower_sigma_r.pdf',
          high_affinity_supp_lower_sigma_r,
          base_width = 8.5, base_height = 9)

# =========================================== High mutation scenario plots ===========================================
high_mutation_scenario_main_fig <- base_plotting_function(summary_tibble = high_mutation_scenario_summary$summary_pairwise_correlations %>%
                                                            filter(method == 'pearson', nGCs == 15),
                                                          y_var = 'freq_ratio_correlation',
                                                          color_var = 'mutation_rate', facet_vars = c('gamma', 'beta')) +
  ylim(-0.2,1) + geom_hline(yintercept = 0, linetype = 2) +
  mutation_rate_color_scale +
  ylab('Pairwise correlation in experienced-to-naive frequency ratios') +
  theme(legend.position = 'bottom')

save_plot('../figures/simulations/high_mutation_scenario_main_fig.pdf',
          high_mutation_scenario_main_fig,
          base_width = 8, base_height = 9)

# % GCs occupied by high-mutation alleles in the high mutation scenario
combined_freq_of_high_mutation_alleles_in_GC_pl <- high_mutation_scenario_summary$combined_freq_of_high_mut_alleles_in_GCs %>%
  filter(beta == 4, gamma == 6) %>%
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



