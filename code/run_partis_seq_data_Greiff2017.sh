#!/bin/bash

# Create output directory if it doesn't exist
mkdir -p ../results/partis/seq_data_Greiff2017/

for mouse_dir in ../data/seq_data_Greiff2017/ERR*
do

    dataset_name=$(basename $mouse_dir)
    sbatch_file=run_partis_${dataset_name}.sbatch
    
    echo "#!/bin/bash
#SBATCH --job-name=run_partis_${dataset_name}_20k
#SBATCH --output=out_err_files/run_partis_${dataset_name}_20k.out
#SBATCH --error=out_err_files/run_partis_${dataset_name}_20k.err
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=4000
#SBATCH --time=200:00:00

# Loading mafft
module load mafft/7.310 

# Defining input and output files
input_file=${mouse_dir}/${dataset_name}_processed_reads.fasta
output_file=../results/partis/seq_data_Greiff2017/${dataset_name}_20k.yaml" > $sbatch_file

# Running partis, randomly picking a subset of the sequences
echo '/project2/cobey/partis/bin/partis partition --n-procs 16 --species mouse --n-random-queries 20000 --infname $input_file --outfname $output_file --extra-annotation-columns regional_bounds:cdr3_seqs:seqs_aa:naive_seq_aa:consensus_seq:consensus_seq_aa
#' >> $sbatch_file

# After running partis, run R script to process output

echo 'module load R/3.6.1

Rscript process_partis_output_seq_data_Greiff2017.R $output_file' >> $sbatch_file

sbatch $sbatch_file
rm $sbatch_file
    
done

