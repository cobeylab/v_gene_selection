# Exploration of germline V alleles identified by partis
source('gene_frequency_functions.R')
source('plot_options.R')
library(readr)
library(vegan)
theme_set(theme_cowplot())

args <- commandArgs(trailingOnly = T)
precomputed_freqs_file <- as.character(args[1])

# Load pre-computed gene frequencies, define and create fig directory (if non-existent)

precomputed_file_name <- basename(precomputed_freqs_file)

output_label <- precomputed_files_labeller(precomputed_file_name)

figure_directory <- paste0('../figures/', output_label, '/exported_ggplot_objects/')
dir.create(figure_directory, recursive = T, showWarnings = F)

load(precomputed_freqs_file)

n_v_genes_by_mouse_path <- paste0('../results/n_vgenes_by_mouse_', output_label, '.csv')

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

build_rarefaction_curves <- function(all_tissue_all_cell_freqs){

  sample_sizes <- seq(1000, max(all_tissue_all_cell_freqs$total_compartment_seqs), 1000)
  
  sample_sizes <- unique(c(sample_sizes, unique(all_tissue_all_cell_freqs$total_compartment_seqs)))
  
  wide_format_tibble <- all_tissue_all_cell_freqs %>%
    select(mouse_id, v_gene, n_vgene_seqs) %>%
    pivot_wider(names_from = 'v_gene', values_from = 'n_vgene_seqs',
                values_fill = 0)
  
  count_matrix <- as.matrix(wide_format_tibble %>% select(-mouse_id))
  
  rarefaction_curves <- rarefy(count_matrix, sample_sizes)
  
  colnames(rarefaction_curves) <- str_remove(colnames(rarefaction_curves), 'N')
  
  rarefaction_curves <- as_tibble(rarefaction_curves) %>%
    mutate(mouse_id = wide_format_tibble$mouse_id) %>%
    select(mouse_id, everything()) %>%
    pivot_longer(cols = as.character(sample_sizes),
                 names_to = 'sample_size', 
                 values_to = 'n_v_genes')  %>%
    mutate(sample_size = as.numeric(as.character(sample_size)))
  
  rarefaction_curves <- get_info_from_mouse_id(rarefaction_curves)
  
  # Exclude values for sample size larger than n. seqs in mouse
  rarefaction_curves <- left_join(rarefaction_curves,
                                  all_tissue_all_cell_freqs %>% 
                                    select(mouse_id, total_compartment_seqs) %>%
                                    unique()) %>%
    filter(sample_size <= total_compartment_seqs)
  
  return(rarefaction_curves)

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

gene_freqs <- bind_rows(gene_freqs, exp_freqs_across_all_tissues)

# Allele frequencies across all cell types and tissues for each mouse (for rarefaction curves)
all_tissue_all_cell_freqs <- gene_freqs %>%
  filter(cell_type %in% c('naive','experienced'), tissue == 'all') %>%
  group_by(mouse_id, day, infection_status, group, group_controls_pooled, v_gene) %>%
  dplyr::summarise(across(c('n_vgene_seqs','total_compartment_seqs'), sum)) %>%
  group_by(mouse_id) %>%
  ungroup() %>%
  mutate(vgene_seq_freq = n_vgene_seqs / total_compartment_seqs) %>%
  mutate(cell_type = 'all', tissue = 'all') %>%
  select(mouse_id, day, infection_status, group, group_controls_pooled, tissue, cell_type, v_gene, n_vgene_seqs,
         total_compartment_seqs, vgene_seq_freq)
 
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
write_csv(obs_n_genes, n_v_genes_by_mouse_path)


# PLOTS:

# Rarefaction_curves
mice_with_enough_LN_seqs <- exp_freqs %>% filter(tissue == 'LN') %>%
  filter(total_compartment_seqs >= min_compartment_size) %>%
  select(mouse_id) %>%
  unique()


rarefaction_curves <- build_rarefaction_curves(all_tissue_all_cell_freqs)

rarefaction_curves_pl <- rarefaction_curves %>%
  ggplot(aes(x = sample_size, y = n_v_genes, group = mouse_id, color = infection_status)) +
  geom_line() +
  scale_x_log10() +
  geom_point(data = rarefaction_curves %>%
               group_by(mouse_id) %>% filter(sample_size == max(sample_size)) %>%
               ungroup()) +
  theme(legend.position = 'top') +
  groups_color_scale(name = 'Infection') +
  xlab('Number of sequences') +
  ylab('Number of V genes')



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


save(total_genes_and_genes_in_LN_pops, n_shared_genes,
     rarefaction_curves_pl, file = paste0(figure_directory,'n_genes.RData'))
