library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(scales)
library(stringr)
source('gene_frequency_functions.R')
theme_set(theme_cowplot())

# ======= Number of sorted cells =======
# Read data
n_sorted_cells <- read_csv('../data/n_sorted_cells.csv') %>%
  mutate(cell_type = as.character(cell_type))
n_sorted_cells$cell_type[n_sorted_cells$cell_type == 'naive'] <- 'IgD+B220+'

# Label tissues and cell types, set plotting order using factor()
n_sorted_cells <- get_info_from_mouse_id(n_sorted_cells) %>%
  mutate(mouse_id = factor(mouse_id, levels = mouse_id_factor_levels),
         group = factor(group, levels = group_factor_levels),
         group_controls_pooled = factor(group_controls_pooled,
                                        levels = group_controls_pooled_factor_levels)) %>%
  mutate(cell_type = case_when(
    cell_type == 'IgD+B220+' ~ 'IgD+B220+ (putatively naive) cells',
    cell_type == 'GC' ~ 'Germinal center cells',
    cell_type == 'PC' ~ 'Plasma cells',
    cell_type == 'mem' ~ 'Memory cells'),
    tissue = case_when(
      tissue == 'LN' ~ 'Mediastinal LN',
      tissue == 'spleen' ~ 'Spleen',
      tissue == 'BM' ~ 'Bone marrow'
    )
  ) %>%
  mutate(cell_type = factor(cell_type, levels = c('IgD+B220+ (putatively naive) cells',
                                                  'Germinal center cells','Plasma cells',
                                                  'Memory cells')),
         tissue = factor(tissue, levels = c('Mediastinal LN',
                                            'Spleen','Bone marrow')))

# Make plot
n_sorted_cells_plot <- n_sorted_cells %>% 
  ggplot(aes(x = group_controls_pooled, y = sorted_cells, color = infection_status)) +
  facet_grid(tissue~cell_type, scales = 'free') +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position = position_jitter(width = 0.1)) +
  scale_y_log10(breaks = c(10,100,1000,10000,100000,1000000),
                labels = c(10,100,1000,10000,100000,1000000)) +
  background_grid() +
  xlab('Group') +
  ylab('Number of cells sorted') +
  theme(legend.position = 'top',
        axis.text.x = element_text(size = 8, angle = 40,
                                   vjust = 0.5)) +
  scale_color_discrete(name = 'Infection')

# Export 
save(n_sorted_cells_plot,
     file = '../figures/all_seqs_freqs/exported_ggplot_objects/n_cells_sorted.RData')

# ======= ELISA TITERS =======

titers <- read_csv('../data/ELISA_titers.csv')
titers <- get_info_from_mouse_id(titers) %>%
  mutate(group_controls_pooled = factor(group_controls_pooled, levels = group_controls_pooled_factor_levels))


# For plotting, set titers below measurement threshold as 0
titers <- titers %>%
  mutate(below_threshold = str_detect(titer, '<')) %>%
  mutate(titer = ifelse(below_threshold, 0, titer),
         titer = as.integer(titer))

# Make plot
titers_against_NL09 <- titers %>% 
  filter(antigen == "NL09_whole_virus") %>%
  ggplot(aes(x = group_controls_pooled, y = titer, color = infection_status)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position = position_jitter(width = 0.2, height = 0),
             size = 3, alpha =0.8) +
  scale_y_log10(oob = scales::squish_infinite,
                labels=trans_format('log10',math_format(10^.x))) + 
  background_grid() +
  geom_hline(yintercept = -100, linetype =2) +
  ylab('IgG ELISA titer against the infecting strain') +
  xlab('Group')

# Export
save(titers_against_NL09,
     file = '../figures/all_seqs_freqs/exported_ggplot_objects/titers.RData')


