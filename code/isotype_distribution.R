# Looks at the distribution of isotypes in each mouse / tissue / cell type

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

source('gene_frequency_functions.R')

annotated_seqs <- read_csv('../processed_data/annotated_seqs.csv') 
# annotated_seqs <- read_csv('~/Desktop/v_gene_selection_files/annotated_seqs.csv')

annotated_seqs <- annotated_seqs %>%
  mutate(across(c('clone_id_partis','partis_uniq_ref_seq','seq_id'), as.character))

annotated_seqs$specimen_cell_subset[annotated_seqs$specimen_cell_subset == 'naïve'] <- 'naive'

annotated_seqs <- annotated_seqs %>% filter(!is.na(n_mutations_partis_nt), !is.na(vgene_mutations_partis_nt))
annotated_seqs <- get_info_from_mouse_id(annotated_seqs)
annotated_seqs <- annotated_seqs %>% dplyr::rename(tissue = specimen_tissue, cell_type = specimen_cell_subset)

# Ignore sequences inferred to be unproductive by partis
annotated_seqs <- annotated_seqs %>% filter(productive_partis)


isotype_distribution <- annotated_seqs %>%
  group_by(mouse_id, day, infection_status, group, group_controls_pooled, cell_type, tissue, isotype) %>%
  summarise(n_seqs = n()) %>%
  group_by(mouse_id, day, infection_status, group, group_controls_pooled, cell_type, tissue) %>%
  mutate(compartment_seqs = sum(n_seqs),
         fraction_seqs = n_seqs/compartment_seqs) %>%
  ungroup()

isotype_distribution_pooling_mice_within_groups <- annotated_seqs %>%
  group_by(group_controls_pooled, cell_type, tissue, isotype) %>%
  summarise(n_seqs = n()) %>%
  group_by(group_controls_pooled, cell_type, tissue) %>%
  mutate(compartment_seqs = sum(n_seqs),
         fraction_seqs = n_seqs/compartment_seqs) %>%
  ungroup() %>%
  mutate(cell_type = factor(cell_type, levels = c('naive','GC','PC','mem')),
         tissue = factor(tissue, levels = c('LN','spleen','BM')),
         group_controls_pooled = factor(group_controls_pooled,
                                        levels = group_controls_pooled_factor_levels))

isotype_distribution_pooling_mice_within_groups  %>%
  ggplot(aes(x = group_controls_pooled, y = fraction_seqs, fill = isotype)) +
  facet_grid(tissue~cell_type) +
  geom_col(color = 'black', size = 0.1) +
  scale_fill_brewer(type = 'qual', direction = 1, palette = 1) +
  scale_x_discrete(labels = function(x){str_replace(x,'-','\n')}) +
  ylab("Fraction of sequences (pooled across mice in group)") +
  xlab("Group") +
  theme(legend.position = c(0.85, 0.17),
        axis.text.x = element_text(size = 11))





