#!/bin/bash
#SBATCH --job-name=get_cluster_level_stats
#SBATCH --output=out_err_files/get_cluster_level_stats.out
#SBATCH --error=out_err_files/get_cluster_level_stats.err
#SBATCH --time=100:00:00
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=15000

module load R
Rscript get_seq_cluster_level_statistics.R