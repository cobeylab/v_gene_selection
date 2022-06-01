# v_gene_selection
V gene usage in mice infected with flu. 

We provide instructions for reproducing the analyses in two steps:

1. Pre-processing and annotation of sequence data.
2. Subsequent analyses of processed data.

Because much of step 1 is computationally expensive and assumes access to a computing cluster, we provide the processed data in this Dryad repository so users can choose to reproduce step 2 without having to reproduce step 1.

## 1. Pre-processing and annotation of sequence data ##
The full sequence data are stored in the MAIN_SEQ_FILE. Executing `split_full_data.sh` breaks this file into a separate csv file for each mouse. Executing `run_partis.sh` runs [*partis*](https://github.com/psathyrella/partis) v0.15.0 for all mice from the specified time point (passed as an argument: 8, 16, 24, 40 or 56). It does so by generating a sbatch file for each mouse and using it to submit a job to a SLURM-based cluster (precise sbatch configurations need to be modified by the user).

`process_partis_output.sh` then processes the yaml files produced by partis (running one job per mouse). Alternatively, the user can run the associated R script, `process_partis_output.R` with the paths to the yaml and csv files for a single mouse as arguments. For each mouse, this processing produces the following files

- `[mouse_id]_annotated_seqs.csv`: sequence-level annotations (one sequence per row).
- `[mouse_id]_clone_info.csv`: clone-level annotation (one B cell clone per row).
- `[mouse_id]_mutations_per_vgene_base.csv`: frequency of mutations for each position in each germline V allele (one row per position).
- `[v/d/j]_genes_[mouse_id].fasta`: fasta file with the sequences of germline alleles detected in each mouse.


Python 2.7.15 with packages sys, csv and os is assumed. LIST R DEPENDENCIES


