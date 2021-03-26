 library(readr)
 library(dplyr)
 library(tidyr)
 library(ggplot2)
 library(cowplot)
 theme_set(theme_cowplot())

 source('gene_frequency_functions.R')
 
 clone_info <- read_csv('../processed_data/clone_info.csv')

 # File in old repo. Update with new one
 #uniq_seq_counts <- read_csv('../../mouse_experiments/results/partis_clone_info.csv')
 
