library(ggplot2)
library(cowplot)
library(magick)
library(stringr)
theme_set(theme_cowplot())

source('plot_options.R')

exported_objecs_dir <- '../figures/exported_ggplot_objects/'
#exported_objecs_dir <- '~/Desktop/v_gene_selection_files/figures/exported_ggplot_objects/'

final_figures_dir <- paste0(dirname(exported_objecs_dir),'/')


ggplot_object_files <- list.files(exported_objecs_dir, pattern = 'RData',
                                  full.names = T)
for(f in ggplot_object_files){load(f)}

# Gene sets and naive freqs

top_2_rows <- plot_grid(total_genes_and_genes_in_LN_pops +
                          ylab('Number of V genes\n(Chao1 estimate)'), 
                        plot_grid(n_shared_genes + 
                                    ylab ('Number of V genes\nshared by mouse pair') +
                                    theme(axis.text.x = element_text(size = 10, angle = 15, vjust = 0.5),
                                          axis.title.y = element_text(size = 12)),
                                  pairwise_naive_correlations_plot +
                                    ylab('Correlation in naive V gene\nfrequencies between mice') +
                                    theme(axis.text.x = element_text(size = 9.5, angle = 15, vjust = 0.5),
                                          axis.title.y = element_text(size = 12)),
                                  align = 'h',
                                  nrow = 1),
                        ncol = 1, labels = c('A','B'), label_size = 16) 

bottom_row <- plot_grid(naive_exp_correlations_plot +
                          ylab('Correlation with naive\nfrequencies within each mouse'),
                        labels = c('C'), label_size = 16)

gene_sets_and_naive_freqs <- plot_grid(top_2_rows, bottom_row, nrow = 2, rel_heights = c(1.7,1)) 
save_plot(paste0(final_figures_dir, 'gene_sets_and_naive_freqs.pdf'),
          gene_sets_and_naive_freqs, 
          base_height = 13, base_width = 15)


# Correlations in V gene freqs and freq deviations
top_row <- plot_grid(pairwise_freq_correlations_plot +
                       ylab('Mouse-pair correlations in V gene\nfrequencies (mice with >= 100 seqs.)') +
                       theme(axis.title.y = element_text(size = 12),
                             plot.margin = margin(l = 20, r = 10,t = 5, b = 30)),
                     pairwise_freq_deviations_plot +
                       ylab('Mouse-pair correlations in frequency\ndeviations (mice with >= 100 seqs.)') +
                       theme(axis.title.y = element_text(size = 12),
                             plot.margin = margin(r = 20, l = 10, t = 5, b = 30)),
                     align = 'h')



top_row_legend <- get_legend(pairwise_freq_correlations_plot + 
                               guides(color = guide_legend('Infection')) +
                               theme(legend.position = 'top',
                                     legend.box.margin = margin(0, 0, 0, 500)))
top_row <- plot_grid(top_row_legend, top_row, nrow = 2, rel_heights = c(0.1,1))


freqs_heatmap_LN_PCs <- image_read_pdf(paste0(exported_objecs_dir,'freqs_heatmap_LN_PCs.pdf'),
                                       density = 600) %>%
  image_resize("1688x1125")

deviations_heatmap_LN_PCs <- image_read_pdf(paste0(exported_objecs_dir,'deviations_heatmap_LN_PCs.pdf'),
                                            density = 600) %>%
  image_resize("1688x1125")

middle_row <- plot_grid(ggdraw() + draw_image(freqs_heatmap_LN_PCs),
                        group_controls_pooled_legend,
                        ggdraw() + draw_image(deviations_heatmap_LN_PCs),
                        nrow = 1, rel_widths = c(1,0.02,1))

bottom_row <- top_genes_LN_PC_day8_plot +
  theme(plot.margin = margin(t = 15, b = 5, l = 20, r = 20),
        legend.position = c(0.200,0.9),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 11, margin = margin(l = 4, r = 4, t = 2)),
        legend.background = element_rect(fill = 'white', color = 'black')) +
  #scale_y_continuous(expand = c(0.001,0.1),
  #                   limits = c(-0.05,NA)) +
  xlab('Top 20 genes in lymph node plasma cells from each infected mouse on day 8')




freq_and_deviation_correlations <- plot_grid(top_row, middle_row, bottom_row, nrow = 3,
                                             rel_heights = c(1,1.5,2),
                                             labels = c('A','B','C'),
                                             label_size = 16)


save_plot(paste0(final_figures_dir, 'freq_and_deviation_correlations.pdf'),
          freq_and_deviation_correlations, 
          base_height = 16, base_width = 16)


# Increasing titers, clonality, mutations over time

titers_figure_path <- '../figures/titers.png'
titers_figure <- image_read(titers_figure_path)

left_panel <- plot_grid(NULL,
                        ggdraw() +
                        draw_image(titers_figure),
                        NULL,
                        nrow = 3,
                        rel_heights = c(1,2,1),
                        labels = c('','A [placeholder]',''), label_size = 16,
                        hjust =  -1)

right_panels_top_row <- fraction_in_top_10_clones_plot +
  theme(legend.box.spacing = unit(1,'pt'),
        legend.spacing = unit(50,'pt'),
        plot.margin = margin(l = 20, r = 10, b = 20)) +
  guides(color = guide_legend(keywidth = 1),
         size = guide_legend(keywidth = 1)) +
  ylab('Fraction of sequences\nin the 10 largest clones')

right_panels_bottom_row <- plot_grid(fraction_clones_with_high_freq_muts_LN_plot +
                                  xlab('') +
                                  theme(legend.position = 'none',
                                        plot.margin = margin(l = 20, r = 10)) +
                                  ylab('Fraction of clones (10+ seqs.) with at\nleast one high-frequency mutation'),
                                mean_n_mutations_LN_plot +
                                  theme(strip.background = element_rect(fill = 'white'),
                                        strip.text = element_text(color = 'white'),
                                        legend.position = 'none',
                                        plot.margin = margin(l = 20, r = 10)) +
                                  ylab('Average number of high-frequency\nmutations (clones with 10+ seqs)'),

          nrow = 2, 
          align = 'v')
bottom_row_legend <- get_legend(fraction_clones_with_high_freq_muts_LN_plot +
             theme(legend.box.margin = margin(l = 70)))
right_panels_bottom_row <- plot_grid(bottom_row_legend, right_panels_bottom_row, nrow = 2,
                                rel_heights = c(0.05,1))

right_panels <- plot_grid(right_panels_top_row, right_panels_bottom_row, nrow = 2, rel_heights = c(1,2),
          labels = c('B','C'), label_size = 16)

evidence_of_clonal_evolution <- plot_grid(left_panel, right_panels, nrow =1, rel_widths = c(1,2))

save_plot('../figures/evidence_of_clonal_evolution.pdf',
          evidence_of_clonal_evolution,
          base_height = 12, base_width = 16)

# Supplementary figs with GC heat maps, GC top genes on days 8 and 16, 

freqs_heatmap_LN_GCs <- image_read_pdf(paste0(exported_objecs_dir,'freqs_heatmap_LN_GCs.pdf'),
                                       density = 600) %>%
  image_resize("1688x1125")

deviations_heatmap_LN_GCs <- image_read_pdf(paste0(exported_objecs_dir,'deviations_heatmap_LN_GCs.pdf'),
                                            density = 600) %>%
  image_resize("1688x1125")

LN_GC_day16_heatmaps <- plot_grid(ggdraw() + draw_image(freqs_heatmap_LN_GCs),
                                  NULL,
                                  group_controls_pooled_legend,
                                  ggdraw() + draw_image(deviations_heatmap_LN_GCs),
                                  nrow = 1, rel_widths = c(1,0.1,0.2,1))

save_plot('../figures/LN_GC_day16_heatmaps.pdf',
          LN_GC_day16_heatmaps,
          base_width = 15, base_height = 8)



top_genes_LN_GC_day8_plot$layers[[3]]$aes_params$size <- 2.5
top_genes_LN_GC_day16_plot$layers[[3]]$aes_params$size <- 3
top_genes_LN_GC_day16_plot$layers[[3]]$aes_params$angle <- -25

LN_GC_top_genes <- plot_grid(top_genes_LN_GC_day8_plot + background_grid() +
                               xlab('Top 20 genes in lymph node GC cells from each infected mouse on day 8') +
                               theme(plot.margin = margin(l = 20, t = 10, b = 10, r = 10)),
                             top_genes_LN_GC_day16_plot + background_grid() +
                               xlab('Top 20 genes in lymph node GC cells from each infected mouse on day 16') + 
                               theme(plot.margin = margin(l = 20, t = 10, b = 10, r = 10)),
                             nrow = 2, rel_heights = c(2,3), labels = c('A','B'), label_size = 16)


save_plot('../figures/LN_GC_top_genes.pdf',
          LN_GC_top_genes,
          base_height = 15, base_width = 16)

# Supplementary fig with top genes on day 16 LN plasma cells

top_genes_LN_PC_day16_plot <- top_genes_LN_PC_day16_plot + background_grid() +
  xlab('Top 20 genes in lymph node plasma cells from each infected mouse on day 16')
save_plot('../figures/top_genes_LN_PC_day16_plot.pdf',
          top_genes_LN_PC_day16_plot,
          base_height = 8, base_width = 14)






