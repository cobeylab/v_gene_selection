library(stringr)
source('plot_options.R')
theme_set(theme_cowplot())

args <- commandArgs(trailingOnly = T)
figures_dir <- args[1] # Path to folder within figures directory

exported_objecs_dir <- paste0(figures_dir, '/exported_ggplot_objects/')

ggplot_object_files <- list.files(exported_objecs_dir, pattern = 'RData',
                                  full.names = T)
for(f in ggplot_object_files){load(f)}

# Some of the analysis and plots from the main analysis are not present in the sensitivity analysis, 
# (hence the if statements in some of the plots)

# Gene sets and naive freqs

if('total_genes_and_genes_in_LN_pops' %in% ls()){
  
  top_2_rows <- plot_grid(total_genes_and_genes_in_LN_pops, 
                          plot_grid(n_shared_genes_pl + 
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
  save_plot(paste0(figures_dir, 'gene_sets_and_naive_freqs.pdf'),
            gene_sets_and_naive_freqs, 
            base_height = 13, base_width = 15)
  
  # Save rarefaction curve plot
  save_plot(paste0(figures_dir, 'rarefaction_curves.pdf'),
            rarefaction_curves_pl,
            base_width = 7, base_height = 5)
  
}



# Correlations in V gene freqs and freq deviations
top_row <- plot_grid(pairwise_freq_correlations_plot +
                       ylab('Pairwise correlation\nbetween individuals') +
                       ylim(-0.2,0.9) +
                       theme(axis.title = element_text(size = 16),
                             plot.margin = margin(l = 20, r = 10,t = 5, b = 30),
                             plot.title = element_text(hjust = 0.5),
                             legend.position = 'none',
                             axis.text.x = element_text(angle = 45, vjust = 0.7)) +
                       ggtitle('Germline V allele frequencies in the response'),
                     pairwise_freq_deviations_plot +
                       ylab('') +
                       ylim(-0.2,0.9) +
                       theme(axis.title = element_text(size = 16),
                             plot.margin = margin(r = 20, l = 10, t = 5, b = 30),
                             plot.title = element_text(hjust = 0.5),
                             legend.position = 'none',
                             axis.text.x = element_text(angle = 45, vjust = 0.7)) +
                       ggtitle('Experienced-to-naive ratios'),
                     align = 'h')



top_row_legend <- get_legend(pairwise_freq_correlations_plot + 
                               guides(color = guide_legend('Infection')) +
                               theme(legend.position = 'top',
                                     legend.box.margin = margin(0, 0, 0, 500)))
top_row <- plot_grid(top_row_legend, top_row, nrow = 2, rel_heights = c(0.1,1))



bottom_row <- top_20_genes_day8_LN_PC +
  facet_wrap('mouse_id', scales = 'free') +
  theme(plot.margin = margin(t = 15, b = 5, l = 20, r = 20),
        legend.position = c(0.200,0.9),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 11, margin = margin(l = 4, r = 4, t = 2)),
        legend.background = element_rect(fill = 'white', color = 'black'),
        axis.title = element_text(size = 16)) +
  #scale_y_continuous(expand = c(0.001,0.1),
  #                   limits = c(-0.05,NA)) +
  xlab('\nTop 20 germline V alleles in lymph node plasma cells from each infected mouse on day 8') +
  ylab('V allele frequency')

freq_and_deviation_correlations <- plot_grid(top_row, bottom_row, nrow = 2,
                                             rel_heights = c(1.3,2),
                                             labels = c('A','B'),
                                             label_size = 16)


save_plot(paste0(figures_dir, 'freq_and_deviation_correlations.pdf'),
          freq_and_deviation_correlations, 
          base_height = 13, base_width = 19)


# Increasing titers, clonality, mutations over time

if('titers_against_NL09' %in% ls()){
  left_panel <- plot_grid(NULL,
                          titers_against_NL09 +
                            theme(legend.position = 'top') +
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
  
  save_plot(paste0(figures_dir, 'evidence_of_clonal_evolution.pdf'),
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
  
  save_plot(paste0(figures_dir, 'evidence_of_clonal_evolution_all_tissues.pdf'),
            evidence_of_clonal_evolution_all_tissues,
            base_height = 14, base_width = 13)           
  
  
  # Supplementary fig with numbers of sorted cells
  save_plot(paste0(figures_dir,'n_sorted_cells.pdf'),
            n_sorted_cells_plot,
            base_width = 12, base_height = 8 )
  
  
}




# Correlation between germline mutability and freq-deviations from naive repertoire
if(!is.na(germline_mutability_by_region_pl)){
  save_plot(paste0(figures_dir,'germline_mutability_fig.pdf'),
            plot_grid(germline_mutability_by_region_pl +
                        theme(axis.text.x = element_text(angle = -45)) +
                        ylab('\nNumber of alleles'),
                      freq_ratio_mutability_correlations_pl ,
                      nrow = 2,
                      labels = c('A','B'),
                      align = 'v'),
            base_height = 13, base_width = 12)
}
          


# Clone rank vs. number of high frequency mutations.
if('clone_rank_vs_high_freq_muts_LN_GCs_plot' %in% ls()){
  clone_rank_vs_high_freq_muts <- plot_grid(clone_rank_vs_high_freq_muts_LN_GCs_plot +
                                              background_grid() +
                                              theme(title = element_text(size = 10)),
                                            clone_rank_vs_high_freq_muts_LN_PCs_plot +
                                              background_grid() +
                                              theme(title = element_text(size = 10)) ,
                                            nrow=2)
  
  save_plot(paste0(figures_dir,'clone_rank_vs_high_freq_muts.pdf'),
            clone_rank_vs_high_freq_muts,
            base_width = 10, base_height = 12)
  
}
  

# Shared mutations:
if('shared_mutations_in_LN_clones_pl' %in% ls()){
  save_plot(paste0(figures_dir,'shared_mutations_in_LN_clones.pdf'),
            shared_mutations_in_LN_clones_pl + background_grid() +
              ylab('Probability that two lineages (10+ seqs.) sharing the same\nV allele have high-frequency mutations in common'),
            base_width = 15, base_height = 7)
}



# Fraction of clones dominated by a single tissue or a single cell type
if('fraction_clones_dominated_by_single_tissue_plot' %in% ls()){
  joint_legend <- get_legend(fraction_clones_dominated_by_single_tissue_plot +
                               theme(legend.box.margin = margin(l = 250)))
  
  clonal_composition <- plot_grid(joint_legend,
                                  plot_grid(fraction_clones_dominated_by_single_tissue_plot +
                                              theme(legend.position = 'none'),
                                            fraction_clones_dominated_by_single_cell_type_plot +
                                              theme(legend.position = 'none')),
                                  nrow = 2, rel_heights = c(0.1,1))
  
  save_plot(paste0(figures_dir,'clonal_compostion.pdf'),
            clonal_composition,
            base_width = 16, base_height = 7)
}


# CDR3 analyses

if('length_matched_CDR3_similarity_plot' %in% ls()){
  cdr3_similarity <- plot_grid(length_matched_CDR3_similarity_plot + 
                                 ylab('Similarity of sequence pairs\nsampled from different mice') +
                                 ggtitle('\nCDR3 sequences matched for length') +
                                 theme(plot.title = element_text(hjust = 0.5, size = 14)) + background_grid(),
                               length_and_allele_matched_CDR3_similarity_plot +
                                 ggtitle('\nCDR3 sequences matched for length and V allele') +
                                 ylab('Similarity of sequence pairs\nsampled from different mice') +
                                 theme(plot.title = element_text(hjust = 0.5, size = 14)) + background_grid(),
                               nrow = 2
  )
  
  save_plot(paste0(figures_dir,'cdr3_similarity.pdf'),
            cdr3_similarity,
            base_width = 15, base_height = 13)
  
  # In CDR3 diversity per V gene plot, color
  
  
  
  cdr3_diversity_per_vgene <- CDR3_similarity_NAIVE_same_mouse_pl +
    scale_x_discrete(labels = function(x){str_remove(x, 'IGHV')}) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 5)) +
    xlab('V allele') +
    ylab('Dissimilarity between length-matched CDR3 sequences') 
  
  
}




# Null model with fixed lineage abundances but randomized lineage V alleles
save_plot(
  paste0(figures_dir,'pw_cors_randomized_lineage_V_alleles.pdf'),
  plot_grid(top_row_legend,
            pw_freq_cors_randomized_lineage_V_alleles,
            pw_freq_deviations_randomized_lineage_V_alleles,
            nrow = 3, rel_heights = c(0.1,1,1)),
  base_width = 16,
  base_height = 12
)

# Plot with focal alleles (those consistently overrepresented early in the plasma cell response)
if(!is.na(focal_alleles_plot)){
  save_plot(
    paste0(figures_dir,'focal_alleles_plot.pdf'),
    focal_alleles_plot$LN$freq_ratios,
    base_width = 16,
    base_height = 10
  )
}

# Detailed arrow plots

make_arrow_plots_panel <- function(plot_list){
  
  # Get legend from first plot in list
  stopifnot(names(plot_list) %in% c('primary-8','primary-16','primary-24','secondary-40','secondary-56'))
  pl_legend <- get_legend(plot_list[[1]] + theme(legend.position = 'bottom'))
  
  cell_type <- unique(plot_list[[1]]$data$cell_type)
  tissue <- unique(plot_list[[1]]$data$tissue)
  
  stopifnot(length(cell_type) == 1)
  stopifnot(length(cell_type) == 1 & tissue == 'LN')
  
  cell_type_label <- case_when(
    cell_type == 'GC' ~ 'germinal center cells',
    cell_type == 'PC' ~ 'plasma cells',
    cell_type == 'mem' ~ 'memory cells'
  )
  
  
  adjusted_plots <- plot_list
  for(i in 1:length(adjusted_plots)){
    day <- str_extract(names(plot_list)[i], '[0-9]+')
    
    ylabel <- case_when(
      day == 8 ~ paste0('V allele frequency in lymph node ', cell_type_label),
      T ~ ''
    )
    xlabel <- case_when(
      day == 24 ~ 'Top 10 germline alleles',
      T ~ ''
    )
    x_axis_expansion <- ifelse(day == 8 & cell_type == 'GC', 0.3, 0.2)
    
    
    adjusted_plots[[i]] <- adjusted_plots[[i]] + 
      theme(legend.position = 'none',
            strip.background = element_blank(),
            strip.text = element_blank()) +
      xlab(xlabel) +
      ylab(ylabel) +
      ggtitle(paste0('Day ', day)) +
      scale_y_continuous(expand = expansion(mult = 0.2)) +
      scale_x_continuous(expand = expansion(mult = x_axis_expansion), breaks = 1:10) 
  
  }
  
  return(
    plot_grid(plot_grid(plotlist = adjusted_plots, nrow = 1),
            pl_legend,
            rel_heights = c(20,1), nrow = 2)
  )
 
}


save_plot(
  paste0(figures_dir,'top_genes_LN_GC.pdf'),
  make_arrow_plots_panel(top_genes_LN_GC),
  base_width = 20,
  base_height = 12
)

save_plot(
  paste0(figures_dir,'top_genes_LN_PC.pdf'),
  make_arrow_plots_panel(top_genes_LN_PC),
  base_width = 20,
  base_height = 12
)

save_plot(
  paste0(figures_dir,'top_genes_LN_mem.pdf'),
  make_arrow_plots_panel(top_genes_LN_mem),
  base_width = 20,
  base_height = 12
)



# Heatmaps go to the supplement

#freqs_heatmap_LN_PCs_pearson <- image_read_pdf(paste0(exported_objecs_dir,'freqs_heatmap_LN_PCs_pearson.pdf'),
#                                       density = 600) %>% image_resize("1688x1125")

#freqs_heatmap_LN_PCs_spearman <- image_read_pdf(paste0(exported_objecs_dir,'freqs_heatmap_LN_PCs_spearman.pdf'),
#                                               density = 600) %>% image_resize("1688x1125")

#deviations_heatmap_LN_PCs_pearson <- image_read_pdf(paste0(exported_objecs_dir,'deviations_heatmap_LN_PCs_pearson.pdf'),
#                                            density = 600) %>% image_resize("1688x1125")

#deviations_heatmap_LN_PCs_spearman <- image_read_pdf(paste0(exported_objecs_dir,'deviations_heatmap_LN_PCs_spearman.pdf'),
#                                                    density = 600) %>% image_resize("1688x1125")

#top_row <- plot_grid(ggdraw() + draw_image(freqs_heatmap_LN_PCs_pearson),
#                        group_controls_pooled_legend,
#                        ggdraw() + draw_image(deviations_heatmap_LN_PCs_pearson),
#                        nrow = 1, rel_widths = c(1,0.02,1))
#bottom_row <- plot_grid(ggdraw() + draw_image(freqs_heatmap_LN_PCs_spearman),
#                        NULL,
#                        ggdraw() + draw_image(deviations_heatmap_LN_PCs_spearman),
#                        nrow = 1, rel_widths = c(1,0.02,1))

#heatmaps_figure <- plot_grid(top_row, bottom_row, nrow = 2, label_x = c(0.5,0.5))

#save_plot(paste0(figures_dir, 'heatmaps_figure.pdf'),
#          heatmaps_figure, 
#          base_height = 16, base_width = 20)



                                  


