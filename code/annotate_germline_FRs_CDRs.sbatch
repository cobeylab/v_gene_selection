#!/bin/bash
#SBATCH --job-name=annotate_germline_FRs_CDRs
#SBATCH --output=out_err_files/annotate_germline_FRs_CDRs.out
#SBATCH --error=out_err_files/annotate_germline_FRs_CDRs.err
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=100:00:00
#SBATCH --mem-per-cpu=16000


# Create fasta file with germline seqs (from csv) to pass to IgBlast 
tail -n+2 ../results/germline_genes_partis.csv | sed 's/^/>/' | tr ',' '\n' > ../results/germline_genes_partis.fasta

tail -n+2 ../results/germline_genes_partis_ogrdb.csv | sed 's/^/>/' | tr ',' '\n' > ../results/germline_genes_partis_ogrdb.fasta


module load python/anaconda-2020.02

# Uses immcantation wrapper to run Igblast on germline alleles identified by partis
# (just to get their FRs and CDR partitions)
/home/mvieira/.local/bin/AssignGenes.py igblast -s ../results/germline_genes_partis.fasta -b ~/share/igblast --organism mouse --loci ig --format airr --outdir ../results

/home/mvieira/.local/bin/AssignGenes.py igblast -s ../results/germline_genes_partis_ogrdb.fasta -b ~/share/igblast --organism mouse --loci ig --format airr --outdir ../results
