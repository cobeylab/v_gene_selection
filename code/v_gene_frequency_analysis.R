library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

source('gene_frequency_functions.R')
 
# Basic info for each clone (germline genes, CDR lenght, naive CDR seq)
clone_info <- read_csv('../processed_data/clone_info.csv') %>%
        dplyr::rename(clone_id = clone_id_partis, v_gene = v_segment_partis, j_gene = j_segment_partis,
                      d_gene = d_segment_partis)

# Some clone summary statistics (mean n mutations, consensus CDR3 sequence for productive sequences)
clone_summary_statistics <- read_csv('../results/clone_summary_statistics.csv')

clone_info <- left_join(clone_info, clone_summary_statistics)
 

# File in old repo. Update with new one
unique_seq_counts <- read_csv('../processed_data/unique_seq_counts.csv')
 
 
unique_seq_counts <- left_join(unique_seq_counts, clone_info)
unique_seq_counts <- get_info_from_mouse_id(unique_seq_counts)

 
 
# Get naive frequencies excluding the LN, tissue-specific frequencies for other cell types
naive_from_tissue <- c('spleen','BM')
naive_freqs <- (calc_gene_freqs(unique_seq_counts %>% filter(tissue %in% naive_from_tissue), long_format = F, by_tissue = F))$naive_freqs
 
exp_freqs <- (calc_gene_freqs(unique_seq_counts, long_format = F, by_tissue = T))$exp_freqs
 
gene_freqs <- left_join(exp_freqs,naive_freqs) %>%
        mutate(tissue = factor(tissue, levels = c('LN','spleen','BM')),
               group_controls_pooled = factor(group_controls_pooled, levels = group_controls_pooled_factor_levels))


# Example of gene rank frequency plot
LN_PC_freqs_primary8 <- gene_freqs %>% 
        filter(group_controls_pooled == 'primary-8', tissue == 'LN', cell_type == 'PC',
               total_mouse_cell_type_seqs >= 100)

vgene_order <- LN_PC_freqs_primary8 %>%
        group_by(v_gene) %>%
        summarise(median_freq_across_mice = median(vgene_seq_freq),
                  median_naive_freq_across_mice = median(naive_vgene_seq_freq)) %>%
        arrange(desc(median_freq_across_mice)) %>%
        pull(v_gene)


LN_PC_freqs_primary8  %>%
        mutate(v_gene = factor(v_gene, levels = vgene_order)) %>%
        ggplot(aes(x = v_gene)) +
        geom_point(aes(y = vgene_seq_freq)) +
        facet_wrap('mouse_id', nrow = 3) +
        scale_y_log10() +
        background_grid(minor = 'none', major = 'y') +
        theme(axis.text.x = element_blank()) +
        xlab('V genes\n(ordered by median frequency across mice)') +
        ylab('Frequency in lymph node plasma cells')

# Select scatterplots showing exp. vs. naive freq correlations

plot_naive_freq_corr <- function(group_controls_pooled, tissue, cell_type){
        stopifnot(tissue == 'LN') # If looking at other tissues have to change y axis
        
        cell_type_label <- switch(cell_type,
                                  'GC' = 'germinal center cells',
                                  'PC' = 'plasma cells',
                                  'mem' = 'memory cells')
        
        gene_freqs %>% 
                filter(group_controls_pooled == !!group_controls_pooled, tissue == !!tissue, cell_type == !!cell_type) %>%
                #filter(total_mouse_cell_type_seqs >= 100, total_mouse_naive_seqs >= 100) %>%
                ggplot(aes(x = naive_vgene_seq_freq, vgene_seq_freq)) +
                geom_point() +
                facet_wrap('mouse_id', nrow = 2) +
                scale_y_log10() +
                scale_x_log10() + 
                geom_abline(slope = 1, intercept = 0, linetype = 2) +
                background_grid() +
                xlab('Frequency in naive repertoire') +
                ylab(paste0('Frequency in lymph node ', cell_type_label))
}
plot_naive_freq_corr('primary-8', 'LN','PC')
plot_naive_freq_corr('primary-16', 'LN','GC')
plot_naive_freq_corr('secondary-40', 'LN','GC')


# Plot of naive-exp correlation coeffs. over timex.
# Exclude mouse 40-7
naive_exp_correlations <- gene_freqs %>% 
        filter(mouse_id != '40-7') %>%
        filter(total_mouse_cell_type_seqs > 100, total_mouse_naive_seqs > 100) %>%
        group_by(mouse_id, day, infection_status, group_controls_pooled, cell_type,
                 tissue, total_mouse_cell_type_seqs) %>%
        summarise(naive_exp_corr = cor.test(vgene_seq_freq,naive_vgene_seq_freq,
                                            method = 'spearman')$estimate) %>%
        ungroup()

naive_exp_corr_weighted_means <- naive_exp_correlations %>%
        group_by(group_controls_pooled, infection_status, cell_type, tissue) %>%
        mutate(weights = total_mouse_cell_type_seqs/sum(total_mouse_cell_type_seqs)) %>%
        summarise(naive_exp_corr_weighted_mean = sum(naive_exp_corr*weights)) %>%
        ungroup()

naive_exp_correlations %>%
        filter(tissue == 'LN', group_controls_pooled != 'control') %>%
        ggplot(aes(x = group_controls_pooled, y = naive_exp_corr, color = infection_status)) +
        geom_boxplot(outlier.alpha = 0) +
        geom_point(aes(size = total_mouse_cell_type_seqs), shape = 1,
                   position = position_jitter(width = 0.1)) +
        geom_point(data = naive_exp_corr_weighted_means %>%
                           filter(tissue == 'LN', group_controls_pooled != 'control'),
                   aes(y = naive_exp_corr_weighted_mean),
                   shape = 4, size = 4, stroke = 2, show.legend = F) +
        facet_grid(.~cell_type, scales = 'free') +
        theme(axis.text.x = element_text(angle = 20, vjust = 0.5), legend.position = 'top') +
        xlab('Group') +
        ylab('Correlation between naive and experienced gene frequencies') +
        geom_hline(yintercept = 0, linetype = 2) +
        scale_color_manual(values = c('green3','dodgerblue2')) +
        scale_size_continuous(name = 'Number of unique sequences',
                              breaks = c(1000,10000,20000)) +
        guides(color = 'none') +
        background_grid()

naive_exp_correlations %>%
        filter(tissue == 'spleen') %>%
        ggplot(aes(x = group_controls_pooled, y = naive_exp_corr, color = infection_status)) +
        geom_boxplot(outlier.alpha = 0) +
        geom_point(aes(size = total_mouse_cell_type_seqs), shape = 1,
                   position = position_jitter(width = 0.1)) +
        geom_point(data = naive_exp_corr_weighted_means %>%
                           filter(tissue == 'spleen'),
                   aes(y = naive_exp_corr_weighted_mean),
                   shape = 4, size = 4, stroke = 2, show.legend = F) +
        facet_grid(.~cell_type, scales = 'free') +
        theme(axis.text.x = element_text(angle = 20, vjust = 0.5), legend.position = 'top') +
        xlab('Group') +
        ylab('Correlation between naive and experienced gene frequencies') +
        geom_hline(yintercept = 0, linetype = 2) +
        scale_size_continuous(name = 'Number of unique sequences') +
        guides(color = 'none') +
        background_grid()

# Test of overrepresentation relative to naive using bootstrap
# For these analyses, adjust observed zeros in naive frequencies for genes known to be present in the mouse
adj_gene_freqs <- left_join(exp_freqs, adjust_zero_naive_freqs(naive_freqs)) %>%
        mutate(tissue = factor(tissue, levels = c('LN','spleen','BM')),
               group_controls_pooled = factor(group_controls_pooled,
                                              levels = group_controls_pooled_factor_levels))

neutral_realizations <- simulate_selection_freq_changes(unique_seq_counts,
                                                        synth_data_input_tibble = 'neutral',
                                                        min_seqs = 0,
                                                        extra_mice = 0,
                                                        by_tissue = T, 
                                                        naive_from_tissue = naive_from_tissue,
                                                        n_reps = 100)

# Rho: ratio of experienced to naive frequency
obs_rhos <- adj_gene_freqs %>%
        mutate(log_rho = log(vgene_seq_freq) - log(naive_vgene_seq_freq)) %>%
        mutate(obs_rho = exp(log_rho)) %>%
        rename(obs_n_vgene_seqs = n_vgene_seqs) %>%
        select(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type, v_gene,
               obs_n_vgene_seqs,obs_rho)

neutral_rhos <- neutral_realizations %>%
        mutate(log_rho = log(vgene_seq_freq) - log(naive_vgene_seq_freq)) %>%
        mutate(rho = exp(log_rho)) %>%
        group_by(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type, v_gene) %>%
        summarise(mean_sim_rho = mean(rho),
                  lbound_sim_rho = quantile(rho, 0.025),
                  ubound_sim_rho = quantile(rho, 0.975)) %>%
        ungroup()

deviation_from_naive <- left_join(obs_rhos, neutral_rhos) %>%
        mutate(deviation_from_naive = case_when(
                obs_rho > ubound_sim_rho ~ 'positive',
                obs_rho < lbound_sim_rho ~ 'negative',
                (obs_rho >= lbound_sim_rho) & (obs_rho <= ubound_sim_rho) ~ 'neutral'
        )) %>%
        mutate(group_controls_pooled = factor(group_controls_pooled, levels = group_controls_pooled_factor_levels)) %>%
        #filter(!is.na(deviation_from_naive))

n_genes_by_deviation <- deviation_from_naive %>%
        group_by(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type) %>%
        mutate(total_genes = length(unique(v_gene))) %>%
        ungroup() %>%
        group_by(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type, total_genes, deviation_from_naive) %>%
        summarise(n_genes = length(unique(v_gene))) %>%
        ungroup() %>%
        mutate(fraction_genes = n_genes / total_genes)

n_genes_by_deviation %>% filter(group_controls_pooled != 'control', tissue == 'LN') %>%
        #filter(!is.na(deviation_from_naive)) %>%
        ggplot(aes(x = group_controls_pooled, y = fraction_genes, color = infection_status)) +
        geom_boxplot(outlier.alpha = 0) +
        geom_point() +
        facet_grid(deviation_from_naive~cell_type) +
        background_grid() + 
        scale_color_manual(values = c('green3','dodgerblue2')) +
        theme(axis.text.x = element_text(angle = 20, vjust = 0.5), legend.position = 'top') +
        xlab('Group') + ylab('Fraction of genes')


# Genes deviating somewhat consistently across mice 
genes_by_n_mice_called <- 
        left_join(deviation_from_naive %>% filter(!is.na(deviation_from_naive)), 
                  gene_freqs %>% select(mouse_id, cell_type, tissue, total_mouse_cell_type_seqs, total_mouse_naive_seqs) %>%
                          unique()) %>%
        filter(total_mouse_cell_type_seqs >= 100, total_mouse_naive_seqs >= 100, mouse_id != '40-7') %>%
        group_by(group_controls_pooled, v_gene) %>% 
        mutate(n_mice_gene_occurs_in_this_group = length(unique(mouse_id))) %>%
        group_by(group_controls_pooled, v_gene, n_mice_gene_occurs_in_this_group, cell_type, tissue, deviation_from_naive) %>%
        summarise(n_mice_with_call = length(unique(mouse_id))) %>% ungroup()

get_consistent_genes_table <- function(candidate_genes_tissue, candidate_genes_cell_type, candidate_genes_group){
        candidate_genes <- genes_by_n_mice_called %>% 
                filter(tissue == candidate_genes_tissue, group_controls_pooled == candidate_genes_group, cell_type == candidate_genes_cell_type,
                       deviation_from_naive == 'positive') %>%
                arrange(desc(n_mice_with_call)) %>%
                filter(n_mice_with_call >= 2) %>% pull(v_gene)
        genes_by_n_mice_called %>% 
                filter(tissue == candidate_genes_tissue , group_controls_pooled == candidate_genes_group, cell_type == candidate_genes_cell_type) %>%
                filter(v_gene %in% candidate_genes) %>% pivot_wider(names_from = deviation_from_naive, values_from = n_mice_with_call) %>%
                arrange(desc(positive)) %>% select(v_gene, positive, neutral, negative) %>%
                mutate(positive = ifelse(is.na(positive), 0, positive),
                       negative = ifelse(is.na(negative), 0, negative),
                       neutral = ifelse(is.na(neutral), 0, neutral)
                ) %>%
                rename(nonsignificant = neutral)
}
 
get_consistent_genes_table('LN','PC', 'primary-8')


# Plots showing deviation from naive freq for major genes
gene_freqs <- left_join(gene_freqs, 
                        deviation_from_naive %>%
                                select(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type, v_gene, matches('rho'), deviation_from_naive))

gene_freqs <- gene_freqs %>%
        group_by(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type) %>%
        mutate(v_gene_rank = rank(-vgene_seq_freq, ties.method = 'first')) %>%
        ungroup()

plot_most_common_genes <- function(plot_cell_type, plot_tissue, plot_group){
        gene_freqs %>% 
                filter(cell_type == plot_cell_type, tissue == plot_tissue, group_controls_pooled == plot_group, v_gene_rank <= 20) %>%
                filter(total_mouse_cell_type_seqs >= 100, total_mouse_naive_seqs >= 100,
                       mouse_id != '40-7') %>%
                rowwise() %>%
                mutate(label_position = ifelse(vgene_seq_freq > naive_vgene_seq_freq, 1.05*vgene_seq_freq, 1.05*naive_vgene_seq_freq)) %>%
                ggplot() +
                geom_point(aes(x = v_gene_rank, y = naive_vgene_seq_freq, color = deviation_from_naive)) +
                geom_segment(aes(x = v_gene_rank, xend = v_gene_rank, y = naive_vgene_seq_freq, yend = vgene_seq_freq,
                                 color = deviation_from_naive),
                             arrow = arrow(ends = 'last', length = unit(10, "pt"), type = 'closed')) +
                geom_text(aes(label = str_remove(v_gene, 'IGHV'),
                              x = v_gene_rank, y = label_position), angle = 20, size = 4, alpha = 0.8) +
                facet_wrap('mouse_id', scales = 'free') +
                xlab(paste0('V gene rank in ', plot_group, ' ', plot_tissue, ' ',
                            switch(plot_cell_type, 'PC' = 'plasma cells', 'GC' = 'germinal center cells', 'mem' = 'memory cells'), ' (top 20 genes only)')) +
                ylab('V gene frequency') +
                theme(legend.position = c(0.4,0.25)) +
                scale_color_discrete(name = 'Deviation from naive frequency\n(bootstrap)', labels = c('negative','non-significant','positive'))
        
}

plot_most_common_genes('PC','LN','primary-8')


# ============== CLONE SIZE DISTRIBUTIONS

# Clone size distributions over time
clone_size_dist_whole_mouse <- unique_seq_counts %>%
        group_by(mouse_id, day, infection_status, group_controls_pooled, clone_id, v_gene) %>%
        summarise(clone_size = sum(uniq_prod_seqs)) %>%
        group_by(mouse_id, day, infection_status, group_controls_pooled) %>%
        mutate(total_seqs = sum(clone_size),
               clone_freq = clone_size / total_seqs) %>%
        mutate(clone_rank = rank(-clone_freq, ties.method = 'first')) %>%
        ungroup() %>%
        mutate(group_controls_pooled = factor(group_controls_pooled, levels = group_controls_pooled_factor_levels))

clone_size_dist_by_tissue_and_cell_type <- unique_seq_counts %>%
        group_by(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type, clone_id, v_gene) %>%
        summarise(clone_size = sum(uniq_prod_seqs)) %>%
        group_by(mouse_id, day, infection_status, group_controls_pooled, tissue, cell_type) %>%
        mutate(total_seqs = sum(clone_size),
               clone_freq = clone_size / total_seqs) %>%
        mutate(clone_rank = rank(-clone_freq, ties.method = 'first')) %>%
        ungroup() %>%
        mutate(group_controls_pooled = factor(group_controls_pooled, levels = group_controls_pooled_factor_levels)) %>%
        mutate(tissue = factor(tissue, levels = c('LN','spleen','BM')),
               cell_type = factor(cell_type, levels = c('naive','GC','PC','mem')))

clone_size_dist_whole_mouse <- left_join(clone_size_dist_whole_mouse, clone_info)





