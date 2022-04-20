library(readr)
library(dplyr)
library(tidyr)
library(shazam)


clone_info <- read_csv('../processed_data/clone_info.csv') %>% mutate(clone_id = as.character(clone_id))
# clone_info <- read_csv('~/Desktop/v_gene_selection/processed_data/clone_info.csv')  %>% mutate(clone_id = as.character(clone_id))

annotated_seqs <- read_csv('../processed_data/annotated_seqs.csv') %>%
  mutate(seq_id = as.character(seq_id),
         clone_id = as.character(clone_id))

# annotated_seqs <- read_csv('~/Desktop/v_gene_selection/processed_data/annotated_seqs.csv') %>% mutate(seq_id = as.character(seq_id), clone_id = as.character(clone_id))


naive_ancestor_annotations <- read_csv('../results/clone_naive_seqs_igblast.tsv')
# naive_ancestor_annotations <- read_tsv('../results/clone_naive_seqs_igblast_SAMPLE.tsv')

naive_ancestor_annotations <- naive_ancestor_annotations %>%
  separate(sequence_id, into = c('mouse_id', 'clone_id'), sep = '_')

focal_seqs <- annotated_seqs %>% filter(tissue == 'LN', cell_type == 'PC', productive_partis == T) %>%
            select(mouse_id, clone_id, seq_id, tissue, cell_type, vgene_region_seq_partis)

focal_seqs <- left_join(focal_seqs, 
                        clone_info %>% select(mouse_id, clone_id, clone_naive_seq_nt_partis)) %>%
  filter(nchar(vgene_region_seq_partis) != nchar(clone_naive_seq_nt_partis))


clone_effective_seqs <- collapseClones(focal_seqs, sequenceColumn = "vgene_region_seq_partis", cloneColumn = "clone_id",
               germlineColumn = 'clone_naive_seq_nt_partis',
               method="thresholdedFreq", minimumFrequency=0.6,includeAmbiguous=FALSE,
               breakTiesStochastic=FALSE)
