library(tidyverse)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(readr)
source('gene_frequency_functions.R')

processed_data_dir <- '../processed_data/'
# processed_data_dir <- '~/Desktop/'
annotated_seqs <- read_csv(paste0(processed_data_dir, 'annotated_seqs.csv'))

annotated_seqs
