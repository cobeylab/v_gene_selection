# v_gene_selection
V gene usage in mice infected with flu. 

1. Pre-processing and annotation of sequence data.
2. Pre-calculation of germline allele frequencies, lineage sizes, mutation frequencies, and randomization-based null distributions.
3. Empirical analyses
4. Running simulations
5. Analyzing simulation results

[WILL JUST SAY WE PROVIDE ALL INTERMEDIATE FILES]

Because parts of steps 1 and 2 are computationally expensive and assume access to a computing cluster, we provide the output of those steps data in this Dryad repository so users can choose to skip them and start from subsequent steps. Similarly, we provide the output of simulations so users can replicate step 5 without having to run the computationally intensive step 4 (MUST ALSO PROVIDE GERMLINE ALLELE IGBLAST TSV WITH FR/CDR3 partitions).

## 1. Pre-processing and annotation of sequence data ##

[After setting up directory structure...]

1.1. Run `split_full_data.sh` to break the file containing the main BCR sequence dataset ([NAME_OF_FILE]) into a separate csv file for each mouse.

1.2. Run `run_partis.sh` to run [partis](https://github.com/psathyrella/partis) v0.15.0 for all mice from the specified time point (passed as an argument: 8, 16, 24, 40 or 56). This script generates a sbatch file for each mouse and uses it to submit a job to a SLURM-based cluster (precise sbatch configurations need to be modified by the user).

1.3. Run `process_partis_output.sh` to process the yaml files produced by partis by running one SLURM job per mouse. Alternatively, the user can run the associated R script, `process_partis_output.R`, with the paths to the yaml and csv files for a single mouse as arguments (in this order). For each mouse, this processing produces the following files

- `[mouse_id]_annotated_seqs.csv`: sequence-level annotations (one sequence per row).
- `[mouse_id]_clone_info.csv`: clone-level annotation (one B cell clone per row).
- `[mouse_id]_mutations_per_vgene_base.csv`: frequency of mutations for each position in each germline V allele (one row per position).
- `[v/d/j]_genes_[mouse_id].fasta`: fasta file with the sequences of germline alleles detected in each mouse.

1.4. Run `combine_files_across_mice.R` to combine those files across mice to produce a single file of each type, while also exporting counts of sequences per mouse/cell type/tissue/clone.

1.5. Run `process_Greiff2017_reads.sh` to process paired-end reads from an independent naive B cell data set from [Greiff et al.(2017)](https://www.sciencedirect.com/science/article/pii/S221112471730565X) using [pRESTO](https://presto.readthedocs.io/en/stable/) v0.6.2.

1.6. Run  `run_partis_seq_data_Greiff2017.sh` to run partis on the processed reads from this second dataset and process the resulting yaml files.

1.7 Run `estimate_error_rate.sh` to estimate the sequencing/amplification error rate based on mutated bases in the constant region, using a SLURM-based cluster (with user-specific configurations) to run `estimate_error_rate.py` for multiple mouse-specific `.csv` files.

Python 2.7.15 with packages sys, csv and os is assumed. LIST R DEPENDENCIES

## 2. Pre-calculation of germline allele frequencies, lineage sizes, mutation frequencies, and randomization-based null distributions. ##

`precompute_gene_and_mutation_frequencies.R` runs the bulk of the computations used by subsequent empirical analyses. It has 3 positional arguments. The first is either `all_seqs` (to compute frequencies based on all productive sequences) or `unique_seqs` (to compute frequencies using only unique productive sequences). The second (`TRUE` or `FALSE`) determines whether naive frequencies are to be estimated from the alternative dataset by Greiff et al. 2017, and the third (`TRUE` or `FALSE`) determines whether novel alleles identified by partis are to be counted together with their inferred parent alleles. To perform the analyses, run:

2.1. `Rscript precompute_gene_and_mutation_frequencies.R all_seqs FALSE FALSE`: pre-computations for core analyses presented in the main text.

2.2. `Rscript precompute_gene_and_mutation_frequencies.R unique_seqs FALSE FALSE`: sensitivity analysis for using unique sequences only.

2.3. `Rscript compute_naive_freqs_seq_data_Greiff2017.R`: computes naive frequencies based on the alternative naive B cell dataset (requires the output of 2.1).

2.4. `Rscript precompute_gene_and_mutation_frequencies.R all_seqs TRUE FALSE`: sensitivity analysis for using naive alleles frequencies from the independent naive B cell dataset. (requires the output of 2.3).

2.5. `Rscript precompute_gene_and_mutation_frequencies.R unique_seqs FALSE TRUE`: sensitivity analysis for collapsing novel alleles.

Because of bootstrapping and replicated randomizations, `precompute_gene_and_mutation_frequencies.R` takes several hours to run. It could be modified to run the randomizations in parallel, but we did not find it necessary because it only needs to be run once for each case (main analysis or sensitivity analysis). Each run of `precompute_gene_and_mutation_frequencies.R` produces an `.RData` object that to be used by downstream scripts.

## 3. Empirical analyses
Different scripts execute different parts of the analysis, exporting plots as `.RData` objects to be subsequently combined by `make_MS_figures.R`. By default, figures are exported to `figures/all_seqs_freqs/`. Figures specificic to sensitivity analyses are exported to the other directories in `figures` (see below).

3.1. *Sorted cells and ELISA titers*

Run `sorted_cells_and_ELISA_titers.R` 

3.2. *Number of V alleles in each mouse, and number of alleles shared between mice*

Run `v_gene_sets_exploration.R`

3.3. *Size and composition of B cell lineages (clones)*

Run `clone_size_and_composition.R`

3.3.*Mutability of germline alleles*
 
 `annotate_germline_FRs_CDRs.sbatch` annotates germline allele sequences with FR and CDR positions using the [Immcantation wrapper for IgBlast](https://changeo.readthedocs.io/en/stable/examples/igblast.html). Because this script is specific to our cluster configuration, we provide the output file (`germline_genes_igblast.tsv`) in the results directory via the Dryad repository. Run `estimate_germline_mutability.R` to estimate the mutability of germline V alleles.
 
3.4. *Analysis of high-frequency mutations*

Run `SHM_analysis.R`

3.5 *CDR3 analysis*

Run `CDR3_analysis.R`.

3.6 *Analysis of germline allele frequencies*

 `allele_frequency_analysis.R` follows the same argument structure outlined in step 2. For instance, `Rscript allele_frequency_analysis.R all_seqs FALSE FALSE` runs the frequency analyses presented in the main text, based on `.RData` object from 2.1.

3.7 *Make figures*

From the `code` directory, run `make_MS_figures.R` with one of the following paths as an argument:

`../figures/all_seqs_freqs`: main analysis presented in the main text.

`../figures/unique_seqs_freqs`: sensitivity analysis for using unique sequences only.

`../figures/all_seqs_freqs_Greiff2017_naive_freqs`: analysis based on alternative naive sequence data.

`../figures/all_seqs_freqs_collapsed_novel_alleles`: sensitivity analysis for collapsing novel alleles.
 
## 4. Running simulations

4.1 *Create input files specifying different scenarios*

Run `simulation_scenarios.R` to create input files for the different scenarios. Input files and results are stored in the directory corresponding to each scenario in `results/simulations/`. For each scenario, there are two types of input files: First, `allele_info.csv` specificies alleles' affinity distributions, mutabilities and naive frequencies. Second `model_parameters.csv` specificies model parameters (different parameter combinations are kept in different directories within `raw_simulation_files/`).

4.2 *Run simulations*

`run_simulations.sh [scenario_directory] [n. individuals] [n. GCs]` runs simulations for multiple germinal centers, individuals and parameter combinations, assuming a (user-specific) SLURM cluster configuration. A single germinal center in a single individual is simulated for a single parameter combination by running `run_simulations.R` with positional arguments:

`Rscript run_simulations.R [path to allele_info.csv] [path to model_paramters.csv] [individual id] [GC_number]`

where `individual_id` is an integer specifying a single individual represented in `allele_info.csv` and GC number is an arbitrary integer id for the germinal center being simulated.

4.3 *Combine simulated GCs*
`combine_simulated_GCs.sh [scenario directory]` combines output files in the `raw_simulation_files` subdirectory for a given scenario.

## 5. Analyzing simulations

5.1 *Summarize simulation results* 

Run `Rscript summarize_simulations [scenario directory]` to produce an `.RData` object with summary statistics for the chosen scenario.

5.2 *Plot simulation results*

Run `plot_simulations.R` to make figures for all scenarios combined (exported to `figures/simulations/`).

