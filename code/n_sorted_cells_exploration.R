library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(stringr)
source('gene_frequency_functions.R')
theme_set(theme_cowplot())

n_sorted_cells <- read_csv('../data/n_sorted_cells.csv')

n_sorted_cells <- get_info_from_mouse_id(n_sorted_cells) %>%
  mutate(mouse_id = factor(mouse_id, levels = mouse_id_factor_levels),
         group = factor(group, levels = group_factor_levels),
         group_controls_pooled = factor(group_controls_pooled,
                                        levels = group_controls_pooled_factor_levels)) %>%
  mutate(cell_type = factor(cell_type, levels = c('naive','GC','PC','mem')),
         tissue = factor(tissue, levels = c('LN','spleen','BM')))

n_sorted_cells_pl <- n_sorted_cells %>% 
  ggplot(aes(x = group_controls_pooled, y = sorted_cells, color = infection_status)) +
  facet_grid(tissue~cell_type, scales = 'free') +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position = position_jitter(width = 0.1)) +
  #geom_text(aes(label = mouse_id), position = position_jitter(width = 0.1)) +
  scale_y_log10(breaks = c(10,100,1000,10000,100000,1000000),
                labels = c(10,100,1000,10000,100000,1000000)) +
  background_grid() +
  scale_x_discrete(labels = function(x){str_replace(x,'-','\n')}) +
  xlab('Group') +
  ylab('Number of cells sorted') +
  theme(legend.position = 'top',
        axis.text.x = element_text(size = 10, angle = 0,
                                   vjust = 0.5)) +
  scale_color_discrete(name = 'Infection status')

#save_plot('../figures/partis_assignment/n_sorted_cells.pdf',
#          n_sorted_cells_pl,
#          base_height = 8, base_width = 14)

  
  