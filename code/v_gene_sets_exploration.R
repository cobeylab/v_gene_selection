# Exploration of V genes identified by partis
# - How many V genes each mouse uses
# - How many undescribed alleles were found in each mouse
library(dplyr)
library(readr)
library(ggplot2)
library(cowplot)
library(vegan)
theme_set(theme_cowplot())
source('gene_frequency_functions.R')


clone_info <- read_csv('../processed_data/clone_info.csv') %>%
  dplyr::rename(clone_id = clone_id_partis, v_gene = v_segment_partis, j_gene = j_segment_partis,
                d_gene = d_segment_partis)

unique_seq_counts <- read_csv('../processed_data/unique_seq_counts.csv')

unique_seq_counts <- left_join(unique_seq_counts, clone_info)

unique_seq_counts <- get_info_from_mouse_id(unique_seq_counts)

unique_seq_counts$mouse_id = factor(unique_seq_counts$mouse_id, levels = mouse_id_factor_levels)
unique_seq_counts$group <- factor(unique_seq_counts$group, 
                           levels = group_factor_levels)
unique_seq_counts$group_controls_pooled <- factor(unique_seq_counts$group_controls_pooled,
                                           levels = group_controls_pooled_factor_levels)

# ==========================================================================================
# How many genes do mice use?
# ==========================================================================================
gene_freqs <- calc_gene_freqs(unique_seq_counts, long_format = T)
gene_freqs_by_tissue <- calc_gene_freqs(unique_seq_counts, long_format = T, by_tissue = T)
gene_freqs <- bind_rows(gene_freqs %>% mutate(tissue = 'all'),
                        gene_freqs_by_tissue)


obs_n_genes <- gene_freqs %>% filter(n_vgene_seqs > 0) %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, cell_type, tissue,
           total_mouse_cell_type_seqs) %>%
  summarise(n_genes = length(unique(v_gene))) %>% ungroup() 

#obs_n_genes <- left_join(obs_n_genes, uniq_seqs_by_tissue)

# Chao1 estimates
compute_chao1 <- function(gene_freqs){
  
  gene_counts <- gene_freqs %>%
    select(mouse_id, day, infection_status, v_gene, cell_type, matches('tissue'), n_vgene_seqs) %>%
    pivot_wider(names_from = v_gene, values_from =  n_vgene_seqs, values_fill = 0) %>%
    arrange(mouse_id)
  
  if(('tissue' %in% names(gene_counts) == F)){
    gene_counts <- gene_counts %>% mutate(tissue = 'all')
  }
  
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
chao1_estimates <- compute_chao1(gene_freqs)

obs_n_genes <- left_join(obs_n_genes,
                         chao1_estimates %>% select(-obs) %>% dplyr::rename(n_genes_chao1 = chao1)) %>%
  mutate(cell_type = factor(cell_type, levels = c('naive','experienced','GC','PC','mem')),
         tissue = factor(tissue, levels = c('all','LN','spleen','BM')))
write_csv(obs_n_genes, '../results/n_v_genes_by_mouse.csv')


n_vgenes_by_group <- obs_n_genes %>%
  filter(tissue == 'all') %>%
  pivot_longer(cols = c('n_genes','n_genes_chao1')) %>%
  dplyr::rename(value_type = name) %>%
  rowwise() %>%
  mutate(point_label = ifelse(value < 50, as.character(mouse_id), '')) %>%
  ungroup() %>%
  ggplot(aes(x = group_controls_pooled, y = value, color = infection_status)) +
  geom_boxplot(outlier.alpha = 0, show.legend = F) +
  geom_point(aes(size = total_mouse_cell_type_seqs),
             position = position_jitter(height = 0, width = 0.2),
             shape = 1) +
  geom_text(aes(label = point_label), show.legend = F, size = 4,
            position = position_nudge(y = 0, x = 0.5)) +
  facet_grid(value_type~cell_type) +
  theme(legend.position = 'top',
        axis.text.x = element_text(size = 10, angle = 40, vjust = 0.5)) +
  scale_color_discrete(guide = 'none') +
  scale_size(name = 'Number of unique sequences') +
  background_grid() +
  xlab('Group') +
  ylab('Number of V genes')

save_plot('../figures/v_gene_sets/n_vgenes_by_group.pdf',
          n_vgenes_by_group, base_width = 18, base_height = 9)

# Number of V genes by group by tissue
chao1_genes_by_group_by_tissue <- obs_n_genes %>%
  ggplot(aes(x = group_controls_pooled, y = n_genes_chao1, color = infection_status)) +
  geom_boxplot(outlier.alpha = 0, show.legend = F) +
  geom_point(aes(size = total_mouse_cell_type_seqs),
             position = position_jitter(height = 0, width = 0.2),
             shape = 1) +
  facet_grid(tissue~cell_type) +
  theme(legend.position = 'top',
        axis.text.x = element_text(size = 10, angle = 40, vjust = 0.5)) +
  scale_color_discrete(guide = 'none') +
  scale_size(name = 'Number of unique sequences') +
  background_grid() +
  xlab('Group') +
  ylab('Number of V genes (chao1 estimate)')

save_plot('../figures/v_gene_sets/chao1_genes_by_group_by_tissue.pdf',
          chao1_genes_by_group_by_tissue, base_height = 8, base_width = 14)

# Number of V genes per pouse
n_vgenes_by_mouse <- obs_n_genes %>%
  filter(tissue == 'all') %>%
  pivot_longer(cols = c('n_genes','n_genes_chao1')) %>%
  dplyr::rename(value_type = name) %>%
  mutate(total_mouse_cell_type_seqs = ifelse(value_type == 'n_genes',total_mouse_cell_type_seqs,'')) %>%
  ggplot(aes(x = mouse_id, y = value, color = infection_status)) +
  geom_point(aes(shape = value_type), size = 4) +
  geom_text(aes(label = total_mouse_cell_type_seqs),
            position = position_nudge(x = 0, y = -5), size = 2) +
  facet_wrap('cell_type', nrow = 3) +
  theme(legend.position = c(0.55,0.12),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5)) +
  scale_color_discrete(name = 'Infection status') +
  background_grid() +
  xlab('Mouse') +
  ylab('Number of V genes') +
  scale_shape_manual(values = c(1,4), 
                     name = '',
                     labels = c('Observed','Chao1'))

save_plot(paste0('../figures/', clonal_method, '_assignment/v_gene_sets_exploration/n_vgenes_by_mouse.pdf'),
          n_vgenes_by_mouse, base_width = 18, base_height = 9)

# Plot of Chao1 estimates
chao1_estimates_pl <- chao1_estimates %>%
  ggplot(aes(x = obs, y = chao1)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) + 
  geom_point(aes(color = infection_status), size = 2, alpha = 0.5) +
  facet_grid(cell_type~tissue) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_x_continuous(limits = c(0, NA)) +
  xlab('Observed number of V genes') +
  ylab('Chao1 estimate') +
  scale_color_discrete(name = 'Infection group') +
  theme(legend.position = 'top',
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10)) +
  background_grid()

save_plot(paste0('../figures/', clonal_method, '_assignment/v_gene_sets_exploration/chao1_estimates_pl.pdf'),
          chao1_estimates_pl, base_width = 13, base_height = 7)

# Number of genes as a function of the number of sequences in each compartment.
obs_n_genes %>%
  filter(tissue == 'all') %>%
  pivot_longer(cols = c('n_genes','n_genes_chao1')) %>%
  ggplot(aes(x = total_mouse_cell_type_seqs, y = value)) +
  geom_text(aes(label = mouse_id, color = infection_status)) +
  facet_grid(cell_type~name) +
  ylab('Number of V genes') +
  xlab('Number of unique sequences in compartment') +
  scale_x_log10() +
  #geom_smooth(method = 'loess', color = 'black') +
  theme(legend.position = 'top') +
  background_grid()



###### Rarefaction analysis
# (using mouse specific gene frequencies in the naive or experienced repertoire for the randomizations.)

rarefact_v_genes <- function(clone_info, min_seqs = 0, n_resamplings = 50,
                             n_sample_size_points = 20){
  
  gene_freqs <- calc_gene_freqs(clone_info, long_format = T)
  
  

  cell_types <- unique(gene_freqs$cell_type)
  sample_sizes <- lapply(as.list(cell_types),
                             FUN = function(ctype, gene_freqs){
                               ctype_freqs <- gene_freqs %>% filter(cell_type == ctype) 
                               min_sample_size <- min(ctype_freqs$total_mouse_cell_type_seqs)
                               max_sample_size <- max(ctype_freqs$total_mouse_cell_type_seqs)
                               stopifnot(max_sample_size > 1000)
                               
                               sample_sizes <- c(
                                 seq(100, max(min_sample_size, 1000), length.out = round(n_sample_size_points/2)),
                                 seq(max(min_sample_size, 1000), max_sample_size, length.out = round(n_sample_size_points/2)))
                               return(unique(sample_sizes))
                             }, gene_freqs = gene_freqs)
  names(sample_sizes) <-  cell_types
  
  
  # Function for resampling V genes
  resample_v_genes <- function(gene_freqs, sample_sizes){
    resampling_tibble <- c()
    
    for(ssize in sample_sizes){
      resampling_tibble <- bind_rows(resampling_tibble,
                gene_freqs %>% group_by(mouse_id, day, infection_status, group_controls_pooled) %>%
                  filter(total_mouse_cell_type_seqs >= ssize) %>%
                  summarise(obs_v_genes = length(unique(v_gene)),
                            total_mouse_cell_type_seqs = unique(total_mouse_cell_type_seqs),
                            resampled_v_genes = length(unique(sample(unique(v_gene), size = ssize,
                                                                     prob = vgene_seq_freq,
                                                                     replace = T)))) %>%
                  mutate(sample_size = ssize)
                
                )
    }
    return(resampling_tibble)
  }
  
  # Perform n_resamplings for the naive and experienced repertoires of each mouse
  resampled_freqs <- c()
  for(ctype in cell_types){
    
    ctype_replicates <- bind_rows(replicate(n_resamplings,
                                            resample_v_genes(gene_freqs %>%
                                                               filter(cell_type == ctype),
                                                                      sample_sizes[[ctype]]),
                                            simplify = F)) %>%
      mutate(cell_type = ctype)
    
    resampled_freqs <- bind_rows(resampled_freqs,
                                             ctype_replicates)
    
  }
  
  
  # Test: the resampled number of genes should never be greater than the observed number of genes
  stopifnot(nrow(resampled_freqs %>% filter(resampled_v_genes > obs_v_genes)) == 0)
 
  
  # Summarise
  v_gene_rarefaction <- resampled_freqs %>%
    group_by(mouse_id, day, infection_status, group_controls_pooled, cell_type, sample_size) %>%
    summarise(v_genes_mean = mean(resampled_v_genes),
              v_genes_lbound = quantile(resampled_v_genes, 0.05),
              v_genes_ubound = quantile(resampled_v_genes, 0.95)) %>%
    ungroup() %>%
    mutate(cell_type = factor(cell_type,
                               levels = c('naive','experienced','GC','mem','PC')))

  
  return(v_gene_rarefaction)
}

v_gene_rarefaction <- rarefact_v_genes(clone_info) 
total_n_mice <- length(unique(v_gene_rarefaction$mouse_id))

v_gene_rarefaction <- bind_rows(v_gene_rarefaction %>% mutate(is_rarefaction = T),
          obs_n_genes %>% filter(tissue == 'all') %>%
            select(mouse_id, day, infection_status, group_controls_pooled, cell_type,
                                 total_mouse_cell_type_seqs, n_genes) %>% 
            dplyr::rename(v_genes_mean = n_genes, sample_size = total_mouse_cell_type_seqs) %>% 
            mutate(v_genes_lbound = NA, v_genes_ubound = NA, is_rarefaction = F))

rarefaction_curves <- v_gene_rarefaction %>%
  ggplot(aes(x = sample_size, y = v_genes_mean, group = mouse_id, color = is_rarefaction)) +
  geom_line()  +
  geom_linerange(data = v_gene_rarefaction %>% filter(is_rarefaction),
                   aes(ymin = v_genes_lbound,
                                    ymax = v_genes_ubound)) +
  geom_point(data = v_gene_rarefaction %>%
                 filter(!is_rarefaction)) +
  facet_wrap('cell_type', nrow = 2, scales = 'free') + 
  xlab('Number of unique sequences') +
  ylab('Number of V genes') +
  scale_color_discrete(name = '', labels = c('Observed number', 'Rarefaction')) +
  theme(legend.position = c(0.75,0.15)) +
  scale_x_log10() +
  background_grid()

save_plot(paste0('../figures/', clonal_method,
                 '_assignment/v_gene_sets_exploration/rarefaction_curves.pdf'),
          rarefaction_curves, base_height = 10, base_width = 16)


# ==========================================================================================
# Matrix of V gene usage
# ==========================================================================================
create_gene_by_mouse_matrix <- function(clone_info, rep_subset = NULL){
  # rep_subset: 'naive' repertoire only or 'experienced' repertoire only. Whole mouse if left as NULL
  gene_by_mouse_matrix <- as_tibble(expand.grid(mouse_id = unique(clone_info$mouse_id),
                                                v_gene = unique(clone_info$v_gene)))
  
  non_zero_entries <- clone_info %>%
    filter(uniq_prod_seqs > 0)
  
  if(!is.null(rep_subset)){
    
    stopifnot(rep_subset %in% c('naive','experienced','mem','PC','GC'))
    if(rep_subset == 'experienced'){
      non_zero_entries <- non_zero_entries %>%
        filter(cell_type != 'naive')
    }else{
      non_zero_entries <- non_zero_entries %>%
        filter(cell_type == rep_subset)
    }
  }
  
  non_zero_entries <- non_zero_entries %>%
    select(mouse_id, v_gene) %>%
    unique() %>%
    group_by(mouse_id, v_gene) %>%
    dplyr::count() %>%
    dplyr::rename(has_v_gene = n)
  
  gene_by_mouse_matrix <- left_join(gene_by_mouse_matrix, non_zero_entries)
  gene_by_mouse_matrix <- left_join(gene_by_mouse_matrix, 
                                    clone_info %>% select(mouse_id, day, infection_status, group_controls_pooled) %>%
                                      unique(), by = 'mouse_id') %>%
    mutate(cell_type = ifelse(is.null(rep_subset), 'whole_mouse', rep_subset)) %>%
    select(cell_type, everything())
  

  # Add IGMT locus and IMGT allele columns (so we know which IMGT alleles new alleles come from)
  gene_by_mouse_matrix <- gene_by_mouse_matrix %>%
    mutate(imgt_allele = str_extract(v_gene, '[^\\+]*'), # Gene name up to '+' sign partis uses to indicate new alleles
           imgt_locus = str_extract(v_gene, '[^\\*]*')) # Gene name up to * IGMT uses to indicate allele
  
  return(gene_by_mouse_matrix)
  
}


gene_by_mouse_matrix <- lapply(list(NULL, 'naive','experienced','GC','PC','mem'), FUN = create_gene_by_mouse_matrix,
                               clone_info = clone_info)
gene_by_mouse_matrix <- bind_rows(gene_by_mouse_matrix)


# Genes missing from naive rep. but present in experienced rep:
gene_by_mouse_matrix %>%
  select(mouse_id, v_gene, cell_type, has_v_gene) %>%
  pivot_wider(names_from = cell_type, values_from = has_v_gene) %>%
  #filter(grepl('\\+', v_gene)) %>%
  filter(is.na(naive), !is.na(experienced))

gene_freqs %>%
  select(mouse_id, v_gene, cell_type, n_vgene_seqs) %>%
  pivot_wider(names_from = cell_type, values_from = n_vgene_seqs) %>%
  filter(naive == 0, experienced > 0) %>%
  group_by(mouse_id) %>% 
  dplyr::count()


# Export table with genes and how often they occur in different groups
n_mice_by_v_gene_table <- gene_by_mouse_matrix %>%
  group_by(infection_status) %>%
  mutate(total_mice_in_infection_status = length(unique(mouse_id))) %>%
  ungroup() %>%
  group_by(v_gene, infection_status, total_mice_in_infection_status, cell_type, .drop = F) %>%
  filter(!is.na(has_v_gene)) %>%
  summarise(n_mice_gene_occurs = length(unique(mouse_id))) %>%
  select(cell_type, v_gene, infection_status, n_mice_gene_occurs, total_mice_in_infection_status) %>%
  ungroup()

stopifnot(all(n_mice_by_v_gene_table$n_mice_gene_occurs <= n_mice_by_v_gene_table$total_mice_in_infection_status))
write_csv(n_mice_by_v_gene_table, '../results/germline_sets/n_mice_by_v_gene_table.csv')


# CSV file with genes present in all controls
genes_in_all_controls <- n_mice_by_v_gene_table %>% 
  filter(cell_type == 'whole_mouse') %>%
  filter(n_mice_gene_occurs == total_mice_in_infection_status, infection_status == 'control') %>%
  pull(v_gene)

n_mice_by_v_gene_table %>%
  filter(v_gene %in% genes_in_all_controls, cell_type == 'whole_mouse') %>%
  dplyr::rename
  pivot_wider(id_cols = 'v_gene', names_from = 'infection_status', values_from = c('n_mice_gene_occurs','total_mice')) %>%
  write_csv('../results/germline_sets/genes_present_in_all_controls.csv')

# CSV file with genes present in all mice
gene_by_mouse_matrix %>%
  filter(cell_type == 'whole_mouse') %>%
  mutate(total_n_mice = length(unique(mouse_id))) %>%
  filter(has_v_gene == 1) %>%
  group_by(v_gene, total_n_mice) %>%
  summarise(n_mice_gene_is_present = length(unique(mouse_id))) %>%
  filter(n_mice_gene_is_present == total_n_mice) %>%
  select(v_gene) %>%
  write_csv('../results/germline_sets/genes_present_in_all_mice.csv')


gene_by_mouse_matrix %>%
  filter(cell_type == 'whole_mouse') %>%
  mutate(is_new_allele = grepl('\\+', v_gene)) %>%
  group_by(is_new_allele) %>%
  summarise(n_genes = length(unique(v_gene))) %>%
  write_csv(paste0('../results/unique_v_genes_', clonal_method,'.csv'))

n_mice_by_v_gene <- gene_by_mouse_matrix %>% 
  group_by(v_gene, cell_type) %>%
  summarise(n_mice = sum(has_v_gene, na.rm = T)) %>%
  arrange(n_mice)

gene_by_mouse_matrix <- left_join(gene_by_mouse_matrix, n_mice_by_v_gene)

# Order genes by number of mice that have that gene at all
v_gene_order <- n_mice_by_v_gene %>% 
  filter(cell_type == 'whole_mouse') %>%
  pull(v_gene)

n_mice_by_v_gene <- n_mice_by_v_gene %>%
  group_by(n_mice, cell_type) %>%
  mutate(label = ifelse(row_number() == 1,
                        as.character(n_mice),''))
  
gene_by_mouse_matrix$v_gene <- factor(gene_by_mouse_matrix$v_gene,
                                      levels = v_gene_order)

gene_by_mouse_matrix_pl <-  gene_by_mouse_matrix %>% 
  filter(cell_type == 'whole_mouse') %>%
  ggplot(aes(x = v_gene, y = mouse_id)) +
  geom_tile(aes(fill = factor(has_v_gene), color = factor(has_v_gene))) +
  theme(legend.position = 'none',
        axis.text.x = element_text(size = 8)) +
  #theme(axis.text.y = element_text(size = 2)) +
  scale_x_discrete(labels = n_mice_by_v_gene %>% filter(cell_type == 'whole_mouse') %>% pull(label)) +
  ylab('Mouse') +
  xlab('V genes grouped and ordered by the number (n) of mice in which they occur') +
  scale_fill_brewer(type = 'qual') +
  scale_color_manual(values = c('black'))

save_plot(paste0('../figures/', clonal_method, '_assignment/v_gene_sets_exploration/gene_by_mouse_matrix_pl.pdf'),
          gene_by_mouse_matrix_pl, base_width = 16, base_height = 8)

# ==========================================================================================
# New alleles
# ==========================================================================================
new_alleles <- gene_by_mouse_matrix %>% 
  filter(cell_type == 'whole_mouse') %>%
  filter(grepl('\\+',v_gene)) %>%
  group_by(v_gene, imgt_allele, imgt_locus) %>%
  summarise(n_mice = sum(has_v_gene, na.rm = T)) %>%
  arrange(desc(n_mice)) %>%
  mutate(fraction_mice = n_mice / !!total_n_mice) # Divide n mice by (total) n

# Distribution of allele frequencies
allele_freq_distribution <- new_alleles %>%
  group_by(fraction_mice) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  mutate(freq = n/sum(n)) %>%
  ggplot(aes(x = fraction_mice, y = freq)) +
  geom_col() +
  scale_x_continuous(breaks = seq(0,0.5,0.1)) +
  xlab('Alelle frequency') +
  ylab('Fraction of alleles with that frequency')
save_plot('../figures/partis_assignment/v_gene_sets_exploration/allele_freq_distribution.pdf',
          allele_freq_distribution, base_width = 5, base_height = 4)

# Each individual new allele and how many mice it occurs in.
new_alleles_pl <- new_alleles %>%
  ggplot(aes(y = v_gene, x = n_mice)) +
  geom_point() +
  theme(axis.text.y = element_text(size = 7)) +
  xlab('Number of mice') +
  ylab('Allele')

save_plot('../figures/partis_assignment/v_gene_sets_exploration/new_alleles.pdf',
         new_alleles_pl, base_height = 10, base_width = 8)

# For each IMGT locus, how many new alleles were found
# And the range of the number of mice in which each new allele occurs
new_alleles_by_IMGT_locus <- new_alleles %>% group_by(imgt_locus) %>%
  summarise(n_new_alleles = n(),
            min_n_mice = min(n_mice),
            max_n_mice = max(n_mice)) %>%
  mutate(min_max_equal = min_n_mice == max_n_mice) %>%
  ungroup() %>%
  mutate(n_mice_range = ifelse(min_max_equal, 
                               paste0("(",min_n_mice,")"),
                               paste0("(", min_n_mice, "-",max_n_mice,")"))) %>%
  arrange(n_new_alleles) %>%
  mutate(imgt_locus = factor(imgt_locus, levels = unique(imgt_locus))) %>%
  ggplot(aes(x = n_new_alleles, y = imgt_locus)) +
  geom_point() +
  scale_x_continuous(breaks = 1:10) +
  geom_text(aes(label = n_mice_range),
            position = position_nudge(x = 0.6),
            size = 3) +
  background_grid(major = 'x') +
  xlab('Number of new alleles\n(each occuring in x-y mice)') +
  ylab('IMGT locus')

save_plot('../figures/partis_assignment/v_gene_sets_exploration/new_alleles_by_IMGT_locus.pdf',
          new_alleles_by_IMGT_locus, base_height =7, base_width = 6)

# For each IMGT allele, number of mice that have at least one new allele of that IMGT allele
# (and the range x-y of the number of alleles of that IMGT allele in a single mouse)

n_mice_with_new_alleles_by_IMGT_locus <- gene_by_mouse_matrix %>% 
  filter(cell_type == 'whole_mouse') %>%
  filter(grepl('\\+',v_gene)) %>%
  filter(!is.na(has_v_gene)) %>%
  group_by(imgt_locus, mouse_id) %>%
  summarise(n_new_alleles = n()) %>%
  group_by(imgt_locus) %>%
  summarise(n_mice = n(),
            min_alleles_in_a_mouse = min(n_new_alleles),
            max_alleles_in_a_mouse = max(n_new_alleles)) %>%
  mutate(min_max_equal = min_alleles_in_a_mouse == max_alleles_in_a_mouse) %>%
  ungroup() %>%
  mutate(n_alleles_range = ifelse(min_max_equal, 
                               paste0("(",min_alleles_in_a_mouse,")"),
                               paste0("(", min_alleles_in_a_mouse, "-",max_alleles_in_a_mouse,")"))) %>%
  arrange(desc(n_mice)) %>%
  mutate(imgt_locus = factor(imgt_locus, levels = unique(imgt_locus))) %>%
  ggplot(aes(x = n_mice, y = imgt_locus)) +
  geom_point() +
  geom_text(aes(label = n_alleles_range),
            position = position_nudge(x = 1.5),
            size = 3) +
  background_grid(major = 'x') +
  xlab('Number of mice with at least one new alelle\n(range of new alleles per mouse)') +
  ylab('IMGT locus')
save_plot('../figures/partis_assignment/v_gene_sets_exploration/n_mice_with_new_alleles_by_IMGT_locus.pdf',
          n_mice_with_new_alleles_by_IMGT_locus, base_height =7, base_width = 6)


# ==========================================================================================
# How many loci does each mouse have?
# How many alleles in each locus (more than 2 in a locus could mean a duplication)
# ==========================================================================================
n_imgt_loci_distribution <- gene_by_mouse_matrix %>% 
  filter(cell_type == 'whole_mouse') %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, imgt_locus) %>%
  summarise(n_alleles = sum(has_v_gene, na.rm = T)) %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled) %>%
  summarise(n_loci = sum(n_alleles>0)) %>%
  ggplot(aes(x = n_loci)) + 
  geom_histogram(binwidth = 1, color = 'white') +
  xlab('Number of IMGT loci') +
  ylab('Number of mice')

save_plot('../figures/partis_assignment/v_gene_sets_exploration/n_imgt_loci_distribution.pdf',
          n_imgt_loci_distribution, base_height =5, base_width = 6)

# Distribution of the number of heterozygous loci
n_heterozygous_loci_distribution <- gene_by_mouse_matrix %>%
  filter(cell_type == 'whole_mouse') %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, imgt_locus) %>%
  summarise(n_alleles = sum(has_v_gene, na.rm = T)) %>%
  filter(n_alleles > 1) %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, n_alleles) %>%
  summarise(n_loci_with_more_than_two_alleles = n()) %>%
  mutate(n_alleles = paste0(n_alleles, ' alleles')) %>%
  ggplot(aes(x = n_loci_with_more_than_two_alleles)) +
  facet_grid(.~n_alleles) +
  geom_histogram(binwidth = 1, color = 'white') +
  xlab('Number of IMGT loci') +
  ylab('Number of mice') +
  scale_y_continuous(breaks = 0:12) +
  scale_x_continuous(breaks = 1:12)
  
save_plot('../figures/partis_assignment/v_gene_sets_exploration/n_heterozygous_loci_distribution.pdf',
          n_heterozygous_loci_distribution, base_height =5, base_width = 10)

# Similar analyses of heterozygosity, but broken down by naive and experienced repertoires
gene_by_mouse_matrix_naive <- gene_by_mouse_matrix %>% filter(cell_type == 'naive')
gene_by_mouse_matrix_experienced <- gene_by_mouse_matrix %>% filter(cell_type == 'experienced')

n_heterozygous_loci_distribution_by_subset <- bind_rows(gene_by_mouse_matrix_experienced %>% 
                                                  mutate(repertoire = 'Experienced repertoire'),
                                                gene_by_mouse_matrix_naive %>%
                                                  mutate(repertoire = 'Naive repertoire')) %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, imgt_locus, repertoire) %>%
  summarise(n_alleles = sum(has_v_gene, na.rm = T)) %>%
  filter(n_alleles > 1) %>%
  group_by(mouse_id, day, infection_status, group_controls_pooled, repertoire, n_alleles) %>%
  summarise(n_loci_with_more_than_two_alleles = n()) %>%
  mutate(n_alleles = paste0(n_alleles, ' alleles')) %>%
  ggplot(aes(x = n_loci_with_more_than_two_alleles)) +
  facet_grid(repertoire~n_alleles) +
  geom_histogram(binwidth = 1, color = 'white') +
  xlab('Number of IMGT loci') +
  ylab('Number of mice') +
  scale_y_continuous(breaks = 0:13) +
  scale_x_continuous(breaks = 1:12) +
  background_grid()
  
save_plot('../figures/partis_assignment/v_gene_sets_exploration/n_heterozygous_loci_distribution_by_subset.pdf',
          n_heterozygous_loci_distribution_by_subset, base_height =6, base_width = 10)

          

# ==========================================================================================
# OLD CODE BELOW - REVIEW OR DELETE
# ==========================================================================================
# ==========================================================================================
# How many V genes do mice use in the experienced and naive repertoires?
# ==========================================================================================

n_v_genes_naive <- gene_freq_changes %>%
  filter(cell_type == 'experienced') %>%
  group_by(mouse_id,infection_status) %>%
  filter(naive_vgene_seq_freq > 0) %>%
  summarise(n_v_genes = length(unique(v_gene))) %>%
  ungroup()

n_v_genes_naive_test <- glm(n_v_genes~infection_status, data = n_v_genes_naive,
                            family = poisson())
plot(n_v_genes_naive_test)
summary(n_v_genes_naive_test)
n_v_genes_naive %>% 
  summarise(median_n_v_genes = median(n_v_genes),
            min_n_v_genes = min(n_v_genes),
            max_n_v_genes = max(n_v_genes))
n_v_genes_naive %>% group_by(infection_status) %>%
  summarise(median_n_v_genes = median(n_v_genes),
            min_n_v_genes = min(n_v_genes),
            max_n_v_genes = max(n_v_genes))


n_v_genes_exp <- gene_freq_changes %>%
  filter(cell_type == 'experienced') %>%
  group_by(mouse_id,infection_status) %>%
  filter(vgene_seq_freq > 0) %>%
  summarise(n_v_genes = length(unique(v_gene))) %>%
  ungroup()
n_v_genes_exp_test <- glm(n_v_genes~infection_status, data = n_v_genes_exp,
                          family = poisson())
plot(n_v_genes_exp_test)
summary(n_v_genes_exp_test)
n_v_genes_exp %>% 
  summarise(median_n_v_genes = median(n_v_genes),
            min_n_v_genes = min(n_v_genes),
            max_n_v_genes = max(n_v_genes))
n_v_genes_exp %>% 
  group_by(infection_status) %>%
  summarise(median_n_v_genes = median(n_v_genes),
            min_n_v_genes = min(n_v_genes),
            max_n_v_genes = max(n_v_genes))

shared_v_genes_exp <- paired_gene_seq_freq_changes %>% 
  filter(vgene_seq_freq_i > 0 & vgene_seq_freq_j > 0) %>%
  group_by(mouse_pair, pair_type, cell_type) %>%
  summarise(n_genes_shared = n()) %>%
  separate(mouse_pair, c('mouse_i','mouse_j'), sep = ';') %>%
  filter(mouse_i != mouse_j) %>%
  ungroup()

shared_v_genes_exp %>% 
  ggplot(aes(x = pair_type, y = n_genes_shared)) +
  geom_boxplot() +
  geom_point()



# Number of unique V genes (V genes only a mouse has) for each mouse
count_unique_mouse_genes(gene_freq_changes) %>%
  summarise(median_unique_genes = median(mouse_unique_genes),
            first_quartile = quantile(mouse_unique_genes, 0.25),
            third_quartile = quantile(mouse_unique_genes, 0.75),
            min = min(mouse_unique_genes),
            max = max(mouse_unique_genes))










# Gets distribution of the number of v genes in each mouse when number of unique sequences is downsampled to sample_size
sample_size <- 1000


downsample_mouse <- function(mouse_data, sample_size, tolerance = 0.05){
    n_seqs <- 0
    v_genes <- c()
    while(n_seqs < sample_size){
      clone <- sample(unique(mouse_data$clone_id),1)
      new_seqs <- mouse_data %>% filter(clone_id == clone) %>% pull(clone_size)
      v_genes <- c(v_genes,
                   unique(as.character(mouse_data %>% filter(clone_id == clone) %>% pull(v_gene)))
                   )
      
      if((n_seqs + new_seqs) < (sample_size * (1 + tolerance))){
        n_seqs <- n_seqs + new_seqs
        #print(n_seqs)
      }
    }
    return(length(unique(v_genes)))
}


mouse_data <- mouse_data %>% group_by(mouse_id, v_gene, clone_id) %>%
  summarise(clone_size = sum(uniq_prod_seqs)) %>%
  ungroup()
