library(dplyr)
library(fitdistrplus)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
source('gene_frequency_functions.R')
library(viridis)

frequency_type <- 'all_seqs'
use_Greiff2017_naive_freqs <- F

results_directory <- '../results/'
#results_directory <- '~/Desktop/v_gene_selection/results/'

processed_data_directory <- '../processed_data/'
#processed_data_directory <- '~/Desktop/v_gene_selection/processed_data/'

figure_directory <- paste0('../figures/', frequency_type, '_freqs/')
#figure_directory <- paste0('~/Desktop/v_gene_selection/figures/',frequency_type, '_freqs/')

precomputed_freqs_file <- paste0('precomputed_gene_freqs_', frequency_type, '.RData')

if(use_Greiff2017_naive_freqs){
  stopifnot(frequency_type == 'all_seqs')
  figure_directory <- '../figures/all_seqs_freqs_Greiff2017_naive_freqs/'
  precomputed_freqs_file <- 'precomputed_gene_freqs_all_seqs_Greiff2017_naive_freqs.RData'
}

# Load precomputed gene frequencies, neutral realizations, pairwise correlations 
load(paste0(results_directory, precomputed_freqs_file))

# General function for fitting a distribution to naive frequencies
fit_distribution_to_naive_freqs <- function(mouse_id, naive_freqs, distribution){
  values <- naive_freqs %>% filter(mouse_id == !!mouse_id) %>% pull(naive_vgene_seq_freq)
  total_mouse_naive_seqs <- unique(naive_freqs %>% filter(mouse_id == !!mouse_id) %>% pull(total_mouse_naive_seqs))
  stopifnot(length(total_mouse_naive_seqs) == 1)
  # Remove artificial zeros (genes that failed to be sampled in naive rep but were found in experienced cells)
  fitted_dist <- fitdist(values[values!=0], distr = distribution, method = "mle")
  fitted_dist$mouse_id = mouse_id
  fitted_dist$total_mouse_naive_seqs <- total_mouse_naive_seqs 
  fitted_dist$distribution = distribution
      
  return(fitted_dist)
}

# Gamma fits
fitted_gammas_list <- lapply(as.list(unique(naive_freqs$mouse_id)),
                             FUN = fit_distribution_to_naive_freqs,
                             naive_freqs = naive_freqs,
                             distribution = 'gamma')

# Exponential fits
fitted_exponentials_list <- lapply(as.list(unique(naive_freqs$mouse_id)),
                                   FUN = fit_distribution_to_naive_freqs,
                                   naive_freqs = naive_freqs,
                                   distribution = 'exp')

# For sanity, normal fits which are obviously bad:
fitted_normals_list <- lapply(as.list(unique(naive_freqs$mouse_id)),
                             FUN = fit_distribution_to_naive_freqs,
                             naive_freqs = naive_freqs,
                             distribution = 'norm')

# Compiles fits across mice
compile_fit_results <- function(fits_list){

  par_names = names(fits_list[[1]]$estimate)
  compiled_results <- bind_rows(lapply(fits_list,
                   FUN = function(fit){
                     row_values <- unlist(fit[c('mouse_id','total_mouse_naive_seqs',
                                         'distribution', 'estimate', 'loglik','aic')])
                     names(row_values) <- str_remove(names(row_values), 'estimate\\.')
                     return(as_tibble(t(data.frame(row_values))))
                   })) %>%
    mutate(across(any_of(c(par_names,'loglik','aic')), as.numeric)) %>%
    mutate(total_mouse_naive_seqs = as.integer(total_mouse_naive_seqs))
  return(compiled_results)

}

gamma_fits <- compile_fit_results(fitted_gammas_list)
exponential_fits <- compile_fit_results(fitted_exponentials_list)
normal_fits <- compile_fit_results(fitted_normals_list)

# Model selection:
model_selection <- bind_rows(gamma_fits %>% dplyr::select(mouse_id, total_mouse_naive_seqs, distribution, loglik, aic),
                             exponential_fits %>% dplyr::select(mouse_id, total_mouse_naive_seqs, distribution, loglik, aic),
                             normal_fits %>% dplyr::select(mouse_id, total_mouse_naive_seqs, distribution, loglik, aic)) %>%
  arrange(mouse_id) %>%
  group_by(mouse_id, total_mouse_naive_seqs) %>%
  mutate(delta_aic = aic - min(aic)) %>%
  arrange(mouse_id, delta_aic) %>%
  dplyr::summarise(best_model = distribution[1],
            second_best_model = distribution[2],
            delta_aic_for_second_best = delta_aic[2])
  
# In almost all mice, the exponential is either the best model by AIC or has a deltaAIC under 2:
nrow(model_selection %>% 
       filter(best_model == 'exp' |
              (second_best_model == 'exp' & delta_aic_for_second_best <2)))

model_selection %>%
  filter(best_model != 'exp' & delta_aic_for_second_best >= 2)



# Plot showing what the best exponential fits look like for each mouse
exponential_plots <- left_join(tibble(x = seq(0,0.1,0.001)),
                            exponential_fits, by = character()) %>%
  mutate(density = dexp(x = x, rate = rate)) %>%
  ggplot(aes(x = x, y = density, group = mouse_id, color = log10(total_mouse_naive_seqs))) +
  geom_line() +
  xlab('Allele frequency in the naive repertoire') +
  ylab('Density') +
  scale_color_viridis(name = 'Number of\nnaive sequences (log10)') +
  theme(legend.position = c(0.65,0.8)) +
  background_grid()
  
plot(exponential_plots)

# Compute average exponential rate excluding mice with < 100 sequences
fitted_exp_rate_of_naive_freq_distributions <- 
  exponential_fits %>%
  filter(total_mouse_naive_seqs >= 100) %>%
  summarise(M = mean(rate)) %>%
  pull(M)

# export this value for simulation model
save(fitted_exp_rate_of_naive_freq_distributions,
     file = paste0(results_directory,'fitted_exp_rate_of_naive_freq_distributions.RData'))

# Export plot with estimated distributions for each mouse
save_plot(paste0(figure_directory,'fitted_naive_frequency_distributions.pdf'),
          exponential_plots,
          base_width = 7,
          base_height = 6)
# And a detailed example showing the quality of fit for a single mouse
pdf(paste0(figure_directory,'fitted_naive_frequency_example.pdf'), width = 8,
    height = 7)

plot(fitted_exponentials_list[[1]])
dev.off()











