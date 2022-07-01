# v_gene_selection
V gene usage in mice infected with flu. 

1. Pre-processing and annotation of sequence data.
2. Pre-calculation of germline allele frequencies, lineage sizes, mutation frequencies, and randomization-based null distributions.
3. Other steps [TODO]

Because parts of steps 1 and 2 are computationally expensive and assume access to a computing cluster, we provide the output of those steps data in this Dryad repository so users can choose to skip them and start from subsequent steps.

## 1. Pre-processing and annotation of sequence data ##

[After setting up directory structure...]

1. Run `split_full_data.sh` to break the file containing the main BCR sequence dataset ([NAME_OF_FILE]) into a separate csv file for each mouse.
2. Run `run_partis.sh` to run [partis](https://github.com/psathyrella/partis) v0.15.0 for all mice from the specified time point (passed as an argument: 8, 16, 24, 40 or 56). This script generates a sbatch file for each mouse and uses it to submit a job to a SLURM-based cluster (precise sbatch configurations need to be modified by the user).
3. Run `process_partis_output.sh` to process the yaml files produced by partis by running one SLURM job per mouse. Alternatively, the user can run the associated R script, `process_partis_output.R`, with the paths to the yaml and csv files for a single mouse as arguments (in this order). For each mouse, this processing produces the following files

- `[mouse_id]_annotated_seqs.csv`: sequence-level annotations (one sequence per row).
- `[mouse_id]_clone_info.csv`: clone-level annotation (one B cell clone per row).
- `[mouse_id]_mutations_per_vgene_base.csv`: frequency of mutations for each position in each germline V allele (one row per position).
- `[v/d/j]_genes_[mouse_id].fasta`: fasta file with the sequences of germline alleles detected in each mouse.

4. Run `combine_files_across_mice.R` to combine those files across mice to produce a single file of each type, while also exporting counts of sequences per mouse/cell type/tissue/clone.
5. Run `process_Greiff2017_reads.sh` to process paired-end reads from an independent naive B cell data set from [Greiff et al.(2017)](https://www.sciencedirect.com/science/article/pii/S221112471730565X) using [pRESTO](https://presto.readthedocs.io/en/stable/) v0.6.2.
6. Run  `run_partis_seq_data_Greiff2017.sh` to run partis on the processed reads from this second dataset and process the resulting yaml files.

Python 2.7.15 with packages sys, csv and os is assumed. LIST R DEPENDENCIES

## 2. Pre-calculation of germline allele frequencies, lineage sizes, mutation frequencies, and randomization-based null distributions. ##



