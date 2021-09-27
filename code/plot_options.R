library(ggplot2)
library(cowplot)
library(dplyr)
library(RColorBrewer)

point_jitter_width <- 0.2
point_alpha <- 0.5
axis_text_x <- element_text(size = 10, angle = 40, vjust = 0.5)


group_order <- c('primary-8','primary-16','primary-24',
                 'secondary-40', 'secondary-56')

group_controls_pooled_palette <- 
  tibble(group_controls_pooled = factor(group_order, levels = group_order),
         group_color = c(#brewer.pal(3, name = 'Reds')[3],
                   brewer.pal(3, name = 'Greens')[c(3,2,1)],
                   brewer.pal(3, name = 'Blues')[c(3,2)]),
         y = 1)

group_controls_pooled_palette_plot <- group_controls_pooled_palette %>%
  ggplot(aes(x = group_controls_pooled, y = 1)) + geom_col(aes(fill = group_controls_pooled)) +
  scale_fill_manual(name = 'Days after\nprimary\ninfection', values = group_controls_pooled_palette$group_color[],
                    labels = function(x){str_extract(x,'[0-9]+')}) +
  theme(legend.position = 'right',
        legend.text = element_text(margin = margin(b = 10, t = 10))) +
  guides(fill = guide_legend(ncol = 1))

group_controls_pooled_legend <- get_legend(group_controls_pooled_palette_plot)
