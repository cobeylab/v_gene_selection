library(ggplot2)
library(cowplot)
library(magick)
library(stringr)
theme_set(theme_cowplot())


source('plot_options.R')

exported_objecs_dir <- '../figures/all_seqs_freqs/exported_ggplot_objects/'

final_figures_dir <- paste0(dirname(exported_objecs_dir),'/')


ggplot_object_files <- list.files(exported_objecs_dir, pattern = 'RData',
                                  full.names = T)
for(f in ggplot_object_files){load(f)}

# Gene sets and naive freqs

top_2_rows <- plot_grid(total_genes_and_genes_in_LN_pops +
                          ylab('Number of V genes\n(Chao1 estimate)'), 
                        plot_grid(n_shared_genes + 
                                    ylab ('Number of V alleles\nshared by mouse pair') +
                                    theme(axis.text.x = element_text(size = 10, angle = 15, vjust = 0.5),
                                          axis.title.y = element_text(size = 12)),
                                  pairwise_naive_correlations_plot +
                                    ylab('Correlation in naive V allele\nfrequencies between mice') +
                                    theme(axis.text.x = element_text(size = 9.5, angle = 15, vjust = 0.5),
                                          axis.title.y = element_text(size = 12)),
                                  nrow = 1),
                        ncol = 1, labels = c('A','B'), label_size = 16) 

bottom_row <- plot_grid(naive_exp_pearson_corr_plot  +
                          ylab('Correlation with naive\nfrequencies within each mouse'),
                        labels = c('C'), label_size = 16)

gene_sets_and_naive_freqs <- plot_grid(top_2_rows, bottom_row, nrow = 2, rel_heights = c(1.7,1)) 
save_plot(paste0(final_figures_dir, 'gene_sets_and_naive_freqs.pdf'),
          gene_sets_and_naive_freqs, 
          base_height = 13, base_width = 15)


# Correlations in V gene freqs and freq deviations
top_row <- plot_grid(pairwise_freq_correlations_plot +
                       ylab('Correlation between pairs of mice') +
                       ylim(-0.2,0.9) +
                       theme(axis.title = element_text(size = 16),
                             plot.margin = margin(l = 20, r = 10,t = 5, b = 30),
                             plot.title = element_text(hjust = 0.5)) +
                       ggtitle('Allele frequencies'),
                     pairwise_freq_deviations_plot +
                       ylab('Correlation between pairs of mice') +
                       ylim(-0.2,0.9) +
                       theme(axis.title = element_text(size = 16),
                             plot.margin = margin(r = 20, l = 10, t = 5, b = 30),
                             plot.title = element_text(hjust = 0.5)) +
                       ggtitle('Frequency deviations from the naive repertoire'),
                     align = 'h')



top_row_legend <- get_legend(pairwise_freq_correlations_plot + 
                               guides(color = guide_legend('Infection')) +
                               theme(legend.position = 'top',
                                     legend.box.margin = margin(0, 0, 0, 500)))
top_row <- plot_grid(top_row_legend, top_row, nrow = 2, rel_heights = c(0.1,1))



bottom_row <- top_genes_LN_PC_day8_plot +
  theme(plot.margin = margin(t = 15, b = 5, l = 20, r = 20),
        legend.position = c(0.200,0.9),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 11, margin = margin(l = 4, r = 4, t = 2)),
        legend.background = element_rect(fill = 'white', color = 'black'),
        axis.title = element_text(size = 16)) +
  #scale_y_continuous(expand = c(0.001,0.1),
  #                   limits = c(-0.05,NA)) +
  xlab('Top 20 alleles in lymph node plasma cells from each infected mouse on day 8') +
  ylab('V allele frequency')




freq_and_deviation_correlations <- plot_grid(top_row, bottom_row, nrow = 2,
                                             rel_heights = c(2,2),
                                             labels = c('A','B'),
                                             label_size = 16)


save_plot(paste0(final_figures_dir, 'freq_and_deviation_correlations.pdf'),
          freq_and_deviation_correlations, 
          base_height = 16, base_width = 16)


# Heatmaps go to the supplement

freqs_heatmap_LN_PCs_pearson <- image_read_pdf(paste0(exported_objecs_dir,'freqs_heatmap_LN_PCs_pearson.pdf'),
                                       density = 600) %>% image_resize("1688x1125")

freqs_heatmap_LN_PCs_spearman <- image_read_pdf(paste0(exported_objecs_dir,'freqs_heatmap_LN_PCs_spearman.pdf'),
                                               density = 600) %>% image_resize("1688x1125")

deviations_heatmap_LN_PCs_pearson <- image_read_pdf(paste0(exported_objecs_dir,'deviations_heatmap_LN_PCs_pearson.pdf'),
                                            density = 600) %>% image_resize("1688x1125")

deviations_heatmap_LN_PCs_spearman <- image_read_pdf(paste0(exported_objecs_dir,'deviations_heatmap_LN_PCs_spearman.pdf'),
                                                    density = 600) %>% image_resize("1688x1125")

top_row <- plot_grid(ggdraw() + draw_image(freqs_heatmap_LN_PCs_pearson),
                        group_controls_pooled_legend,
                        ggdraw() + draw_image(deviations_heatmap_LN_PCs_pearson),
                        nrow = 1, rel_widths = c(1,0.02,1))
bottom_row <- plot_grid(ggdraw() + draw_image(freqs_heatmap_LN_PCs_spearman),
                        NULL,
                        ggdraw() + draw_image(deviations_heatmap_LN_PCs_spearman),
                        nrow = 1, rel_widths = c(1,0.02,1))

heatmaps_figure <- plot_grid(top_row, bottom_row, nrow = 2, label_x = c(0.5,0.5))

save_plot(paste0(final_figures_dir, 'heatmaps_figure.pdf'),
          heatmaps_figure, 
          base_height = 16, base_width = 20)




# Increasing titers, clonality, mutations over time

left_panel <- plot_grid(NULL,
                        titers_against_NL09 +
                          theme(legend.position = 'top') +
                          scale_color_discrete(name = 'Infection') +
                          scale_x_discrete(labels = function(x){str_replace(x,'-','\n')}),
                        NULL,
                        nrow = 3,
                        rel_heights = c(1,2,1),
                        labels = c('','A',''), label_size = 16,
                        hjust =  -1)

right_panels_top_row <- fraction_in_top_10_clones_plot +
  theme(legend.box.spacing = unit(1,'pt'),
        legend.spacing = unit(50,'pt'),
        plot.margin = margin(l = 20, r = 10, b = 20)) +
  guides(color = guide_legend(keywidth = 1),
         size = guide_legend(keywidth = 1)) +
  ylab('Fraction of sequences\nin the 10 largest lineages')

right_panels_bottom_row <- plot_grid(fraction_clones_with_high_freq_muts_LN_plot +
                                  xlab('') +
                                  theme(legend.position = 'none',
                                        plot.margin = margin(l = 20, r = 10)) +
                                  ylab('Fraction of lineages (10+ seqs.) with at\nleast one high-frequency mutation'),
                                mean_n_mutations_LN_plot +
                                  theme(strip.background = element_rect(fill = 'white'),
                                        strip.text = element_text(color = 'white'),
                                        legend.position = 'none',
                                        plot.margin = margin(l = 20, r = 10)) +
                                  ylab('Average number of high-frequency\nmutations (lineages with 10+ seqs)'),

          nrow = 2, 
          align = 'v')
bottom_row_legend <- get_legend(fraction_clones_with_high_freq_muts_LN_plot +
             theme(legend.box.margin = margin(l = 70)))
right_panels_bottom_row <- plot_grid(bottom_row_legend, right_panels_bottom_row, nrow = 2,
                                rel_heights = c(0.05,1))

right_panels <- plot_grid(right_panels_top_row, right_panels_bottom_row, nrow = 2, rel_heights = c(1,2),
          labels = c('B','C'), label_size = 16)

evidence_of_clonal_evolution <- plot_grid(left_panel, right_panels, nrow =1, rel_widths = c(1.2,2))

save_plot(paste0(final_figures_dir, 'evidence_of_clonal_evolution.pdf'),
          evidence_of_clonal_evolution,
          base_height = 12, base_width = 16)

# Supp fig. with patterns for all tissues:
  
joint_legend <- get_legend(fraction_clones_with_high_freq_muts_all_tissues_plot +
                             theme(legend.box.margin = margin(l = 50)))

evidence_of_clonal_evolution_all_tissues <-  plot_grid(joint_legend,
           plot_grid(fraction_clones_with_high_freq_muts_all_tissues_plot +
                       theme(legend.position = 'none',
                             axis.text.x = element_blank()) +
                       xlab('') +
                       ylab('Fraction of lineages (10+ seqs.) with at\nleast one high-frequency mutation'),
                     mean_n_mutations_all_tissues_plot +
                       theme(legend.position = 'none',
                             strip.background.x = element_blank(),
                             strip.text.x = element_blank(),
                             axis.text.x = element_text(size = 10, angle = 20)) +
                       ylab('Average number of high-frequency\nmutations (lineages with 10+ seqs)'),
                     nrow = 2, align = 'v', rel_heights = c(1,1.1)),
           nrow = 2, rel_heights = c(0.1,2))

save_plot(paste0(final_figures_dir, 'evidence_of_clonal_evolution_all_tissues.pdf'),
          evidence_of_clonal_evolution_all_tissues,
          base_height = 14, base_width = 13)           
           

# Supplementary fig with numbers of sorted cells
save_plot(paste0(final_figures_dir,'n_sorted_cells.pdf'),
          n_sorted_cells_plot,
          base_width = 12, base_height = 8 )


# Supplementary fig with top genes on day 16 LN plasma cells

top_genes_LN_PC_day16_plot <- top_genes_LN_PC_day16_plot + background_grid() +
  xlab('Top 20 genes in lymph node plasma cells from each infected mouse on day 16')

save_plot(paste0(final_figures_dir,'top_genes_LN_PC_day16_plot.pdf'),
          top_genes_LN_PC_day16_plot,
          base_height = 8, base_width = 14)

# Supplementary figs with GC top genes on days 8 and 16, 


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


save_plot(paste0(final_figures_dir,'LN_GC_top_genes.pdf'),
          LN_GC_top_genes,
          base_height = 15, base_width = 16)


# Correlation between germline mutability and freq-deviations from naive repertoire
save_plot(paste0(final_figures_dir,'germline_mutability_fig.pdf'),
          plot_grid(germline_mutability_by_region_pl +
                      theme(axis.text.x = element_text(angle = -45)) +
                      ylab('\nNumber of alleles'),
                    freq_ratio_mutability_correlations_pl ,
                    nrow = 2,
                    labels = c('A','B'),
                    align = 'v'),
          base_height = 13, base_width = 12)
          


# Clone rank vs. number of high frequency mutations.
clone_rank_vs_high_freq_muts <- plot_grid(clone_rank_vs_high_freq_muts_LN_GCs_plot +
                                            background_grid() +
                                            theme(title = element_text(size = 10)),
                                          clone_rank_vs_high_freq_muts_LN_PCs_plot +
                                            background_grid() +
                                            theme(title = element_text(size = 10)) ,
                                          nrow=2)

save_plot(paste0(final_figures_dir,'clone_rank_vs_high_freq_muts.pdf'),
          clone_rank_vs_high_freq_muts,
          base_width = 10, base_height = 12)
  

# Shared mutations:
save_plot(paste0(final_figures_dir,'shared_mutations_in_LN_clones.pdf'),
          shared_mutations_in_LN_clones_pl + background_grid() +
          ylab('Probability that two lineages (10+ seqs.) sharing the same\nV allele have high-frequency mutations in common'),
          base_width = 15, base_height = 7)



# Fraction of clones dominated by a single tissue or a single cell type
joint_legend <- get_legend(fraction_clones_dominated_by_single_tissue_plot +
                             theme(legend.box.margin = margin(l = 250)))

clonal_composition <- plot_grid(joint_legend,
          plot_grid(fraction_clones_dominated_by_single_tissue_plot +
                      theme(legend.position = 'none'),
                    fraction_clones_dominated_by_single_cell_type_plot +
                      theme(legend.position = 'none')),
          nrow = 2, rel_heights = c(0.1,1))

save_plot(paste0(final_figures_dir,'clonal_compostion.pdf'),
          clonal_composition,
          base_width = 16, base_height = 7)

# Top genes in day 24 memory cells
save_plot(paste0(final_figures_dir,'day24_LN_mem_top_genes.pdf'),
          top_genes_LN_mem_day24_plot +
            xlab('Top 20 genes in lymph node memory cells from each infected mouse on day 24'),
          base_height = 5,
          base_width = 8)



# CDR3 analyses

cdr3_similarity <- plot_grid(length_matched_CDR3_similarity_plot + ylab('Similarity of sequence pairs\nsampled from different mice') +
            ggtitle('\nCDR3 sequences matched for length') +
            theme(plot.title = element_text(hjust = 0.5, size = 14)) + background_grid(),
          length_and_allele_matched_CDR3_similarity_plot +
            ggtitle('\nCDR3 sequences matched for length and V allele') +
            ylab('Similarity of sequence pairs\nsampled from different mice') +
            theme(plot.title = element_text(hjust = 0.5, size = 14)) + background_grid(),
          nrow = 2
)
save_plot(paste0(final_figures_dir,'cdr3_similarity.pdf'),
          cdr3_similarity,
          base_width = 15, base_height = 13)




convergent_cdr3_seqs <- plot_grid(high_similarity_length_and_allele_matched_seqs_day56_LN_PCs,
                                  plot_grid(allele_usage_day56_LN_PC_convergent_CDRs +
                                              ylab('Fraction of sequences') + ylim(0,1),
                                            combined_freq_of_day56_LN_PC_convergent_CDRs +
                                              theme(legend.position = 'top') , nrow = 2, 
                                            labels = c("B", "C")),
                                  nrow = 1, rel_widths = c(5,10),
                                  labels = c("A",""))

save_plot(paste0(final_figures_dir,'convergent_cdr3_seqs.pdf'),
          convergent_cdr3_seqs,
          base_width = 16, base_height = 11)
                                  


