# Exploration of germline V alleles identified by partis
source('gene_frequency_functions.R')
source('plot_options.R')
library(readr)
library(vegan)
theme_set(theme_cowplot())

clone_info <- read_csv('../processed_data/clone_info.csv')

# Load pre-computed gene frequencies
load('../results/precomputed_gene_freqs_all_seqs.RData')

min_compartment_size = 100 # For certain plots, exclude mice with fewer than 100 sequences.

# Function for computing Chao 1 estimate of number of alleles
compute_chao1 <- function(gene_freqs){
  
  gene_counts <- gene_freqs %>%
    select(mouse_id, day, infection_status, v_gene, cell_type, matches('tissue'), n_vgene_seqs) %>%
    pivot_wider(names_from = v_gene, values_from =  n_vgene_seqs, values_fill = 0) %>%
    arrange(mouse_id)
  
  cell_types <- unique(gene_counts$cell_type)
  tissues <- unique(gene_counts$tissue)
  
  base_function <- function(cell_type, tissue, gene_counts){
    gene_counts_subset <- gene_counts %>% 
      filter(cell_type == !!cell_type, tissue == !!tissue)
    if(nrow(gene_counts_subset) > 0){
      stopifnot(length(unique(gene_counts_subset$cell_type)) == 1)
      stopifnot(length(unique(gene_counts_subset$tissue)) == 1)
      
      chao1_est_subset <- gene_counts_subset %>%
        select(mouse_id, day, infection_status) %>%
        mutate(obs = estimateR(gene_counts_subset %>% select(matches('IGHV')))['S.obs',],
               chao1 = estimateR(gene_counts_subset %>% select(matches('IGHV')))['S.chao1',]) %>%
        mutate(cell_type = cell_type, tissue = tissue)
      return(chao1_est_subset)
    }else{
      return(c())
    }
  }
  
  ctype_tissue_combinations <- gene_counts %>% select(cell_type, tissue) %>% arrange(tissue)%>%
    unique()
  
  
  chao1_ests <- mapply(FUN = base_function, cell_type = ctype_tissue_combinations$cell_type,
                       tissue = ctype_tissue_combinations$tissue,
                       MoreArgs = list(gene_counts = gene_counts), SIMPLIFY = F)
  chao1_ests <- bind_rows(chao1_ests)
  
  return(chao1_ests)
}

# ==========================================================================================
# How many alleles do mice use?
# ==========================================================================================
# Put gene freqs back in long format
gene_freqs <- bind_rows(exp_freqs,
                        naive_freqs %>%
                          dplyr::rename(n_vgene_seqs = n_naive_vgene_seqs,
                                        vgene_seq_freq = naive_vgene_seq_freq,
                                        total_compartment_seqs = total_mouse_naive_seqs) %>%
                          mutate(tissue = 'all', cell_type = 'naive'))
  
# Compute allele freqs across all experienced cells across all tissues
exp_freqs_across_all_tissues <- exp_freqs %>%
  group_by(mouse_id, day, infection_status, group, group_controls_pooled, cell_type, v_gene) %>%
  dplyr::summarise(across(c('n_vgene_seqs','total_compartment_seqs'), sum)) %>%
  mutate(vgene_seq_freq = n_vgene_seqs / total_compartment_seqs) %>%
  mutate(tissue = 'all') %>%
  ungroup()

gene_freqs <- bind_rows(gene_freqs,
                        exp_freqs_across_all_tissues)
 
# For each tissue / cell type combination in each mouse, count number of alleles with freq. > 0
obs_n_genes <- gene_freqs %>% filter(n_vgene_seqs > 0) %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, cell_type, tissue,
           total_compartment_seqs) %>%
  summarise(n_genes = length(unique(v_gene))) %>% ungroup() 

# Add Chao 1 estimates
chao1_estimates <- compute_chao1(gene_freqs)

obs_n_genes <- left_join(obs_n_genes,
                         chao1_estimates %>% select(-obs) %>% dplyr::rename(n_genes_chao1 = chao1)) %>%
  mutate(cell_type = factor(cell_type, levels = c('naive','experienced','nonnaive_IgD+B220+','GC','PC','mem')),
         tissue = factor(tissue, levels = c('all','LN','spleen','BM')),
         group_controls_pooled = factor(group_controls_pooled, levels = group_controls_pooled_factor_levels))

# Export csv files
write_csv(obs_n_genes, '../results/n_v_genes_by_mouse.csv')


# PLOTS:

# Number of genes shared by pairs of mice
n_shared_genes <- pairwise_gene_freqs %>%
  filter(cell_type == 'naive') %>%
  filter(vgene_seq_freq_i != 0, vgene_seq_freq_j != 0,
         total_mouse_naive_seqs_i >= min_compartment_size, total_mouse_naive_seqs_j >= min_compartment_size) %>%
  group_by(mouse_pair, pair_type) %>%
  dplyr::summarise(n_genes_shared = n()) %>%
  ungroup() %>%
  ggplot(aes(x = pair_type, y = n_genes_shared, color = pair_type)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.5) +
  xlab('Type of pair') +
  ylab('Number of genes shared by mouse pair') +
  pair_types_color_scale(name = '') +
  background_grid() +
  theme(legend.position = 'none')

total_genes_and_genes_in_LN_pops <- 
  bind_rows(obs_n_genes %>%
              filter(tissue == 'all', cell_type %in% c('naive', 'experienced')) %>%
              mutate(cell_type = case_when(
                cell_type == 'naive' ~ 'Naive cells (all tissues)',
                cell_type == 'experienced' ~ 'Experienced cells (all tissues)'
              )) %>%
              select(-tissue),
            obs_n_genes %>%
              filter(tissue == 'LN', cell_type %in% c('GC','PC','mem')) %>%
              mutate(cell_type = case_when(
                cell_type == 'GC' ~ 'Lymph node GC cells',
                cell_type == 'PC' ~ 'Lymph node plasma cells',
                cell_type == 'mem' ~ 'Lymph node memory cells'
              )) %>%
              select(-tissue)) %>%
  mutate(cell_type = factor(cell_type,
         levels = c('Naive cells (all tissues)',
                    'Experienced cells (all tissues)',
                    'Lymph node GC cells',
                    'Lymph node plasma cells',
                    'Lymph node memory cells'))) %>%
  ggplot(aes(x = group_controls_pooled, y = n_genes_chao1, color = infection_status)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(aes(size = total_compartment_seqs),
             position = position_jitter(height = 0, width = 0.1),
             alpha = point_alpha) +
  facet_grid(.~cell_type) +
  theme(legend.position = 'top',
        axis.text.x = axis_text_x) +
  groups_color_scale(name = 'Infection') +
  scale_size(name = 'Number of sequences') +
  background_grid() +
  xlab('Group') +
  ylab('Number of V genes (Chao1 estimate)')

# Put these plots in the folder with figures based on analysis of all sequences (main analysis; as opposed to unique seqs only)

figure_output_dir = '../figures/all_seqs_freqs/exported_ggplot_objects/'
dir.create(figure_output_dir, recursive = T, showWarnings = F)

save(total_genes_and_genes_in_LN_pops, n_shared_genes,
     file = paste0(figure_output_dir,'n_genes.RData'))
