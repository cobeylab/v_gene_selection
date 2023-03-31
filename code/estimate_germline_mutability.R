library(shazam)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)

# Read germline sequences annotated with FRs and CDRS
# Uses partis inferred genes but runs them through IgBlast (since partis annotates only CDR3, not the other CDRs and FRs)
annotated_germline_seqs_file <- '../results/germline_genes_partis_igblast.tsv'

annotated_germline_seqs <- read_tsv(annotated_germline_seqs_file)

germline_mutability_by_region  <- annotated_germline_seqs %>% select(sequence_id, sequence, fwr1, fwr2, fwr3, cdr1, cdr2) %>%
  dplyr::rename(whole_sequence = sequence) %>%
  pivot_longer(cols = c('whole_sequence','fwr1','fwr2','fwr3','cdr1','cdr2'),
               names_to = 'region', values_to = 'region_seq') %>%
  mutate(region_type = str_remove(region,'[0-9]+')) %>%
  select(sequence_id, region_type, region, region_seq) %>%
  mutate(region_length = nchar(region_seq)) %>%
  mutate(average_S5F_mutability = calculateMutability(sequences = region_seq, model = HH_S5F)/region_length,
         average_RS5NF_mutability = calculateMutability(sequences = region_seq, model = MK_RS5NF)/region_length) %>%
  mutate(region = factor(region, levels = c('fwr1','cdr1','fwr2','cdr2','fwr3','whole_sequence'))) %>%
  dplyr::rename(v_gene = sequence_id)
  

write_csv(germline_mutability_by_region, file = '../results/germline_mutability_by_region.csv')

# Compute average mutability for type of region (all FRs, all CDRs, instead of individual ones)
# doing a weighted average with the region lengths of the same type as weights
germline_mutability_by_region_type <- germline_mutability_by_region %>%
  group_by(v_gene, region_type) %>%
  summarise(average_S5F_mutability = sum(region_length*average_S5F_mutability)/sum(region_length),
            average_RS5NF_mutability = sum(region_length*average_RS5NF_mutability)/sum(region_length)) %>%
  ungroup() %>%
  pivot_wider(names_from = 'region_type', values_from = matches('_mutability')) %>%
  mutate(average_S5F_mutability_cdr_minus_fwr = average_S5F_mutability_cdr - average_S5F_mutability_fwr,
         average_RS5NF_mutability_cdr_minus_fwr = average_RS5NF_mutability_cdr - average_RS5NF_mutability_fwr)

write_csv(germline_mutability_by_region_type, file = '../results/germline_mutability_by_region_type.csv')

# Some plots for visual inspection:
germline_mutability_by_region %>%
  filter(region != 'whole_sequence') %>%
  ggplot(aes(x = region, y = average_RS5NF_mutability)) +
  geom_point() +
  geom_boxplot()

germline_mutability_by_region %>%
  filter(region != 'whole_sequence') %>%
  ggplot(aes(x = average_S5F_mutability, y = average_RS5NF_mutability)) +
  geom_point() +
  facet_wrap('region')

  

