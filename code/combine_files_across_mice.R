library(dplyr)
library(tidyr)
library(readr)

seq_level_files <- list.files('../processed_data/seq_level_files/', pattern = 'csv', full.names = T)
clone_info_files <- list.files('../processed_data/clone_info_files/', pattern = 'csv', full.names = T)

seq_level_data <- lapply(seq_level_files, read_csv)
seq_level_data <- bind_rows(seq_level_data)
write_csv(seq_level_data, '../processed_data/seq_level_data.csv')

clone_info <- lapply(clone_info_files, read_csv)
clone_info <- bind_rows(clone_info)
write_csv(clone_info, '../processed_data/clone_info.csv')

