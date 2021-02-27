library(dplyr)
library(tidyr)
library(readr)

theme_set(theme_cowplot())

seq_level_files <- list.files('../processed_data/seq_level_files/', pattern = 'csv', full.names = T)

seq_level_data <- lapply(seq_level_files, read_csv)
seq_level_data <- bind_rows(seq_level_data)
write_csv(seq_level_data, '../processed_data/seq_level_files/seq_level_data.csv')
