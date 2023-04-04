#!/bin/bash
time_point=$1

# Create partis results directory (if it doesn't exist)
mkdir -p ../results/partis/

# Directory for keeping cluster out and err files
mkdir -p out_err_files

sbatch_file=run_partis_day_${time_point}.sbatch
# --------------- SBATCH FILE ---------
echo '#!/bin/bash' > $sbatch_file
echo "#SBATCH --job-name=run_partis_ogrdb_day_${time_point}" >> $sbatch_file
echo '#SBATCH --output=out_err_files/run_partis_ogrdb_%A_%a.out' >> $sbatch_file
echo '#SBATCH --error=out_err_files/run_partis_ogrdb_%A_%a.err' >> $sbatch_file
echo "#SBATCH --time=200:00:00" >> $sbatch_file
echo "#SBATCH --partition=cobey" >> $sbatch_file
echo "#SBATCH --nodes=1" >> $sbatch_file
echo "#SBATCH --ntasks-per-node=16" >> $sbatch_file
echo "#SBATCH --mem-per-cpu=3000" >> $sbatch_file
    

echo input_file=../processed_data/mouse_specific_data_files/$time_point-'${SLURM_ARRAY_TASK_ID}'.csv >> $sbatch_file

echo output_file=../results/partis/${time_point}-'${SLURM_ARRAY_TASK_ID}_partis_ogrdb'.yaml >> $sbatch_file

# If output file already exists
echo if [ -e '$output_file' ] >> $sbatch_file 
echo then >> $sbatch_file
echo    echo "Output file already exists. Skipping" >> $sbatch_file
# If not...
echo else >> $sbatch_file
    # Export sequences to temporary fasta file 
    echo python seqs_to_fasta.py '$input_file' >> $sbatch_file

    echo temp_fasta_file=../processed_data/mouse_specific_data_files/$time_point-'${SLURM_ARRAY_TASK_ID}'.fasta >> $sbatch_file

    # Run partis
    echo /project2/cobey/partis/bin/partis partition --n-procs 15 --species c57bl --leave-default-germline --infname '$temp_fasta_file' --outfname '$output_file' --extra-annotation-columns regional_bounds:cdr3_seqs:seqs_aa:naive_seq_aa:consensus_seq:consensus_seq_aa >> $sbatch_file

    # Remove temporary fasta file
    echo rm '$temp_fasta_file' >> $sbatch_file
    
    # Add new line to partis output to avoid R's read_yaml complaining
    echo 'echo "" >> $output_file' >> $sbatch_file

echo fi >> $sbatch_file



#echo "rm -r _output" >> $sbatch_file
# ---------------------------------------

# Find all mouse numbers for the specified time point
ids=$(ls ../processed_data/mouse_specific_data_files/$time_point-*.csv | grep -o "\-[0-9]*" | tr '\n' ',' | tr -d '-')
ids=${ids::-1}


sbatch --array=$ids $sbatch_file
rm $sbatch_file
