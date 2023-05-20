#!/bin/bash
#SBATCH --job-name=precompute_gene_frequencies
#SBATCH --output=out_err_files/precompute_gene_frequencies.out
#SBATCH --error=out_err_files/precompute_gene_frequencies.err
#SBATCH --partition=cobey
#SBATCH --account=pi-cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20:00:00
#SBATCH --mem-per-cpu=128G

module load R/3.6.1

Rscript precompute_gene_and_mutation_frequencies.R all_seqs FALSE FALSE igblast

