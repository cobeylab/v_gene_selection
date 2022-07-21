library(readr)
source('partis_output_functions.R')
source('gene_frequency_functions.R')

args <- commandArgs(trailingOnly = T)
yaml_file_path <- args[1]

yaml_object <- read_yaml(yaml_file_path)

annotated_seqs <- format_partis_info(yaml_object)

base_path <- str_remove(yaml_file_path,'\\.yaml')

write_csv(annotated_seqs, paste0(base_path, '_annotated_seqs.csv'))

write.fasta(yaml_object$`germline-info`$seqs$v, names = names(yaml_object$`germline-info`$seqs$v),
            file.out = paste0(base_path,'_v_genes.fasta'))
write.fasta(yaml_object$`germline-info`$seqs$d, names = names(yaml_object$`germline-info`$seqs$d),
            file.out = paste0(base_path,'_d_genes.fasta'))
write.fasta(yaml_object$`germline-info`$seqs$j, names = names(yaml_object$`germline-info`$seqs$j),
            file.out = paste0(base_path,'_j_genes.fasta'))
