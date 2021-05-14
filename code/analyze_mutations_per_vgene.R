library(tidyverse)
library(cowplot)
library(shazam)
library(stringr)
theme_set(theme_cowplot())
source('gene_frequency_functions.R')

mutations_per_vgene_base <- read_csv('../results/mutations_per_vgene_base.csv')
germline_v_genes <- read_csv('../results/germline_genes.csv') %>%
  mutate(v_gene_length = nchar(v_gene_seq))


# Complete mutations per base tibble with unmutated bases represented explicitly as zeros.
complete_mutations_tibble <- left_join(mutations_per_vgene_base, germline_v_genes %>% select(v_gene, v_gene_length)) %>%
  group_by(mouse_id, cell_type, tissue, v_gene) %>%
  complete(position = 1:v_gene_length[1],
           fill = list(n_vgene_seqs_in_compartment = NA,
                       n_mutations = 0, mutation_freq = 0, v_gene_length = NA)) %>%
  mutate(n_vgene_seqs_in_compartment = unique(n_vgene_seqs_in_compartment[!is.na(n_vgene_seqs_in_compartment)]),
         v_gene_length = unique(v_gene_length[!is.na(v_gene_length)])) %>%
  ungroup()
  

# Calculate the mutability of each germline V gene
get_position_mutability <- function(germline_v_genes){

  # For each gene, get 5-mer centered at
  get_position_fivemers <- function(selected_v_gene, germline_v_genes){
    v_gene_seq = germline_v_genes %>% filter(v_gene == selected_v_gene) %>% pull(v_gene_seq)
    defined_positions <- seq(3, nchar(v_gene_seq) - 2)
    fivemers <- lapply(as.list(defined_positions),
           FUN = function(position){
             fivemer = str_sub(v_gene_seq, position - 2, position + 2)
           })
    
    return(tibble(v_gene = selected_v_gene,
                  position = defined_positions,
                  fivemer = as.character(fivemers)))
  }
  
  position_fivemers <- lapply(as.list(unique(germline_v_genes$v_gene)),
         FUN = get_position_fivemers, germline_v_genes = germline_v_genes)
  position_fivemers <- bind_rows(position_fivemers)
  
  # Get predicted mutability values for each fivemer using the 'RS5NF' model from Cui et al. (2016)
  predicted_mutability <- tibble(fivemer = names(MK_RS5NF@mutability),
                                 mutability = MK_RS5NF@mutability)
  
  return(left_join(position_fivemers, predicted_mutability, by = 'fivemer'))
}

per_position_mutability <- get_position_mutability(germline_v_genes)

complete_mutations_tibble <- left_join(complete_mutations_tibble,
                                      per_position_mutability, by = c('v_gene','position'))

complete_mutations_tibble <- get_info_from_mouse_id(complete_mutations_tibble) %>%
  mutate(group_controls_pooled = factor(group_controls_pooled, levels = group_controls_pooled_factor_levels))

naive_cell_correlation_with_mutability <- complete_mutations_tibble %>% 
  filter(cell_type == 'naive', n_vgene_seqs_in_compartment >= 100) %>%
  group_by(mouse_id, day, infection_status, group, group_controls_pooled, tissue, v_gene) %>%
  summarise(cor_coef = cor.test(mutability, mutation_freq, method = 'spearman')$estimate) %>%
  ungroup()

naive_cell_correlation_with_mutability %>%
  ggplot(aes(x = group_controls_pooled, y = cor_coef)) +
  geom_boxplot(outlier.alpha = 0, color = 'red') +
  #geom_violin() +
  geom_point(position = position_jitter(width = 0.05),
             alpha = 0.5) +
  facet_wrap('tissue') +
  ylab('Correlation between observed site mutation frequency and \n predicted site mutability (naive cells only)') +
  xlab('Group') +
  scale_x_discrete(labels = function(x){str_replace(x,'-','\n')}) +
  ggtitle('Each point is a V gene in a mouse. V genes with < naive 100 seqs. in a mouse we excluded')
  
complete_mutations_tibble %>%
  filter(n_vgene_seqs_in_compartment > 100) %>%
  filter(mouse_id == '16-1', tissue == 'spleen', cell_type == 'naive') %>%
  ggplot(aes(x = mutability, y = mutation_freq)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = 'lm') +
  facet_wrap('v_gene') +
  xlab('Predicted (RS5NF) site mutability') +
  ylab('Observed site mutation frequency') +
  ggtitle('Naive B cells from the spleen of mouse 16-1')

# Considering only third base in each codon
naive_cell_correlation_with_mutability_third_base <- complete_mutations_tibble %>% 
  filter(cell_type == 'naive', n_vgene_seqs_in_compartment >= 100) %>%
  filter(position%%3 == 0) %>%
  group_by(mouse_id, day, infection_status, group, group_controls_pooled, tissue, v_gene) %>%
  summarise(cor_coef = cor.test(mutability, mutation_freq, method = 'spearman')$estimate) %>%
  ungroup()

naive_cell_correlation_with_mutability_third_base %>%
  ggplot(aes(x = group_controls_pooled, y = cor_coef)) +
  geom_boxplot(outlier.alpha = 0, color = 'red') +
  #geom_violin() +
  geom_point(position = position_jitter(width = 0.05),
             alpha = 0.5) +
  facet_wrap('tissue') +
  ylab('Correlation between observed site mutation frequency and \n predicted site mutability (naive cells only, third codon bases only)') +
  xlab('Group') +
  scale_x_discrete(labels = function(x){str_replace(x,'-','\n')}) +
  ggtitle('Each point is a V gene in a mouse. V genes with < naive 100 seqs. in a mouse we excluded')

complete_mutations_tibble %>%
  filter(n_vgene_seqs_in_compartment > 100) %>%
  filter(position%%3 == 0) %>%
  filter(mouse_id == '16-1', tissue == 'spleen', cell_type == 'naive') %>%
  ggplot(aes(x = mutability, y = mutation_freq)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = 'lm') +
  facet_wrap('v_gene') +
  xlab('Predicted (RS5NF) site mutability') +
  ylab('Observed site mutation frequency') +
  ggtitle('Naive B cells from the spleen of mouse 16-1')



