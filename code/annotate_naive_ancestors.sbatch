#!/bin/bash
#SBATCH --job-name=annotate_naive_ancestors
#SBATCH --output=out_err_files/annotate_naive_ancestors.out
#SBATCH --error=out_err_files/annotate_naive_ancestors.err
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --time=100:00:00
#SBATCH --mem-per-cpu=1000


module load python/anaconda-2020.02

# Uses immcantation wrapper to run Igblast on clone's naive ancestors
# (to get their FRs and CDR partitions)
AssignGenes.py igblast -s ../results/clone_naive_seqs.fasta -b ~/share/igblast --organism mouse --loci ig --format airr --nproc 27 --outdir ../results

