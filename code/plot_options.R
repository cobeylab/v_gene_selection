library(ggplot2)
library(cowplot)
library(dplyr)
library(RColorBrewer)

point_jitter_width <- 0.2
point_alpha <- 0.5
axis_text_x <- element_text(size = 10, angle = 40, vjust = 0.5) # Not currently used


group_order <- c('primary-8','primary-16','primary-24',
                 'secondary-40', 'secondary-56') # (Not currently used)

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

# Labels cell type panels (for plots made with single-tissue tibbles)
# (Could be made into a pipe function.)
cell_type_facet_labeller <- function(data_tibble){

  if('compartment_cell_type' %in% names(data_tibble)){
    stopifnot('compartment_tissue' %in% names(data_tibble))
    rename_cell_type_var <- T # Handles tibbles where cell type and tissue vars are 'compartment'
    data_tibble <- data_tibble %>% dplyr::rename(cell_type = compartment_cell_type,
                                                 tissue = compartment_tissue)
    # (Changes names back after labeling operation)
  }else{
    rename_cell_type_var <- F
  }
  
  selected_tissue <- unique(data_tibble$tissue)
  stopifnot(length(selected_tissue) == 1)
  
  tissue_label <- case_when(selected_tissue == 'LN' ~ 'Lymph node',
                            selected_tissue == 'BM' ~ 'Bone marrow',
                            selected_tissue == 'spleen' ~ 'Spleen')
    
    
  data_tibble <- data_tibble  %>% 
    mutate(day = as.integer(as.character(day))) %>%
    mutate(cell_type = case_when(
      cell_type == 'GC' ~ paste0(tissue_label, ' GC cells'),
      cell_type == 'PC' ~ paste0(tissue_label, ' plasma cells'),
      cell_type == 'mem' ~ paste0(tissue_label, ' memory cells'),
      cell_type == 'nonnaive_IgD+B220+' ~ paste0(tissue_label, ' non-naive IgD+B220+ cells')
    )) %>%
    mutate(cell_type = factor(cell_type,
                              levels = c(paste0(tissue_label, ' GC cells'),
                                         paste0(tissue_label, ' plasma cells'),
                                         paste0(tissue_label, ' memory cells'),
                                         paste0(tissue_label, ' non-naive IgD+B220+ cells')))) 
  if(rename_cell_type_var){
    data_tibble <- data_tibble %>% dplyr::rename(compartment_cell_type = cell_type,
                                                 compartment_tissue = tissue)
  }
  return(data_tibble)
}

set_controls_as_day_0 <- function(data_tibble){
  data_tibble%>% 
    mutate(day = as.integer(as.character(day))) %>%
    mutate(day = ifelse(group_controls_pooled == 'control', 0, day))
  
}

label_controls_as_day_0 <- scale_x_continuous(breaks = c(0,8,16,24,40,56), labels = c('control', '8','16','24','40','56'))
