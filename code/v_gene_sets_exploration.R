# Exploration of germline V alleles identified by partis
source('gene_frequency_functions.R')
source('plot_options.R')
library(readr)
library(vegan)
theme_set(theme_cowplot())

# ==== Parsing arguments, importing data ====

args <- commandArgs(trailingOnly = T)

assignment <- as.character(args[1])
collapse_novel_alleles <- as.logical(as.character(args[2]))

clone_info <- read_csv(paste0('../processed_data/clone_info_', assignment, '.csv'))

if(collapse_novel_alleles){
  stopifnot(assignment == 'partis') # Other assignments do not look for novel alleles
  clone_info <- clone_info %>% mutate(v_gene = str_remove(v_gene, '\\+.*'))
}

# For the analyses of germline V sets detected in the mice, include unproductive sequences.
seq_counts_incl_unprod <- read_csv(paste0('../processed_data/seq_counts_', assignment, '_incl_unprod.csv')) 

seq_counts_incl_unprod <- left_join(seq_counts_incl_unprod, clone_info %>% select(mouse_id, clone_id, v_gene))


figure_directory <- paste0('../figures/all_seqs_', assignment, 
                           ifelse(collapse_novel_alleles, '_collapsed_novel_alleles/','/'),
                           'exported_ggplot_objects/')

n_vgenes_by_mouse_path <- paste0('../results/n_vgenes_by_mouse_', assignment,
                                 ifelse(collapse_novel_alleles, '_collapsed_novel_alleles.csv','.csv'))

# ==== Defining additional functions ====

# Function for computing Chao 1 estimate of number of alleles
compute_chao1 <- function(gene_counts){
  
  gene_counts <- gene_counts %>%
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

# Function for building rarefaction curve
build_rarefaction_curves <- function(whole_mice_gene_counts){
  
  
  
  sample_sizes <- seq(1000, max(whole_mice_gene_counts$total_compartment_seqs), 1000)
  
  sample_sizes <- unique(c(sample_sizes, unique(whole_mice_gene_counts$total_compartment_seqs)))
  
  wide_format_tibble <- whole_mice_gene_counts %>%
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
                                  whole_mice_gene_counts %>% 
                                    select(mouse_id, total_compartment_seqs) %>%
                                    unique()) %>%
    filter(sample_size <= total_compartment_seqs)
  
  return(rarefaction_curves)
  
}

count_shared_genes <- function(whole_mice_gene_counts){
  
  unique_pairs <- get_unique_pairs(whole_mice_gene_counts, within_groups_only = F)
  
  base_function <- function(mouse_pair, whole_mice_gene_counts){
    mouse_ids <- str_split(mouse_pair, ';')[[1]] 
    
    mouse_pair_counts <- whole_mice_gene_counts %>% 
      filter(mouse_id %in% mouse_ids)
    
    pair_type <- get_pair_type(unique(mouse_pair_counts$infection_status))
    
    mouse_i_genes <- mouse_pair_counts %>% filter(mouse_id == mouse_ids[1]) %>%
      pull(v_gene)
    
    mouse_j_genes <- mouse_pair_counts %>% filter(mouse_id == mouse_ids[2]) %>%
      pull(v_gene)
    
    n_shared_genes <- length(intersect(mouse_i_genes, mouse_j_genes))
    
    return(tibble(
      mouse_pair = mouse_pair,
      pair_type = pair_type,
      n_shared_genes
    ))
  }
  
  shared_genes <- bind_rows(lapply(as.list(unique_pairs),
                                   FUN = base_function,
                                   whole_mice_gene_counts = whole_mice_gene_counts))
  return(shared_genes)
  
}

# ==== Running the analyses ====

# Sequence counts per V gene (by tissue/cell type & for the whole mouse)

vgene_counts_by_tissue_cell_type <- seq_counts_incl_unprod %>%
  group_by(mouse_id, tissue, cell_type, v_gene) %>%
  summarise(n_vgene_seqs = sum(seqs)) %>%
  ungroup() 

vgene_counts_whole_mice <- seq_counts_incl_unprod %>%
  group_by(mouse_id, v_gene) %>%
  summarise(n_vgene_seqs = sum(seqs)) %>%
  ungroup() %>%
  mutate(tissue = 'all', cell_type = 'all') 

vgene_counts <- bind_rows(vgene_counts_by_tissue_cell_type, vgene_counts_whole_mice) %>%
  filter(!is.na(v_gene)) %>%
  group_by(mouse_id, tissue, cell_type) %>%
  mutate(total_compartment_seqs = sum(n_vgene_seqs)) %>%
  ungroup() %>%
  get_info_from_mouse_id()

rm(vgene_counts_by_tissue_cell_type)
rm(vgene_counts_whole_mice)


# Number of V genes (whole mouse and by tissue/cell type)
obs_n_genes <- vgene_counts %>% 
  group_by(mouse_id, day, infection_status, group_controls_pooled, cell_type, tissue,
           total_compartment_seqs) %>%
  summarise(n_genes = length(unique(v_gene))) %>% ungroup() 

# Add in Chao 1 estimates
chao1_estimates <- compute_chao1(gene_counts = vgene_counts)

obs_n_genes <- left_join(obs_n_genes,
                         chao1_estimates %>% select(-obs) %>% dplyr::rename(n_genes_chao1 = chao1)) %>%
  mutate(cell_type = factor(cell_type, levels = c('all','GC','PC','mem','naive','nonnaive_IgD+B220+','unassigned_IgD+B220+')),
         tissue = factor(tissue, levels = c('all','LN','spleen','BM')),
         group_controls_pooled = factor(group_controls_pooled, levels = group_controls_pooled_factor_levels))


# Number of shared genes (whole mice only)
shared_genes <- count_shared_genes(whole_mice_gene_counts = vgene_counts %>% filter(cell_type == 'all', tissue == 'all'))

# Rarefaction curves (whole mice only)
rarefaction_curves <- build_rarefaction_curves(vgene_counts %>% filter(tissue == 'all', cell_type == 'all'))

# Export csv file with number of alleles
write_csv(obs_n_genes, n_vgenes_by_mouse_path)


# PLOTS:

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


total_genes_and_genes_in_LN_pops <- 
  bind_rows(obs_n_genes %>%
              filter(tissue == 'all', cell_type == 'all') %>%
              mutate(cell_type = 'All cell types and tissues') %>%
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
         levels = c('All cell types and tissues',
                    'Lymph node GC cells',
                    'Lymph node plasma cells',
                    'Lymph node memory cells'))) %>%
  ggplot(aes(x = group_controls_pooled, y = n_genes, color = infection_status)) +
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
  ylab('Number of V genes')

# Number of genes shared by pairs of mice
n_shared_genes_pl <- shared_genes %>%
  mutate(pair_type = factor(pair_type, levels = c('control','control/primary','primary',
                                                  'control/secondary','secondary','primary/secondary'))) %>%
  ggplot(aes(x = pair_type, y = n_shared_genes, color = pair_type)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.5) +
  xlab('Type of pair') +
  ylab('Number of genes shared by mouse pair') +
  pair_types_color_scale(name = '') +
  background_grid() +
  ylim(c(0,NA)) +
  theme(legend.position = 'none')

save(total_genes_and_genes_in_LN_pops, n_shared_genes_pl,
     rarefaction_curves_pl, file = paste0(figure_directory,'n_genes.RData'))
