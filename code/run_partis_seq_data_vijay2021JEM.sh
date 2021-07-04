#!/bin/bash


sbatch_file=run_partis_seq_data_vijay2021JEM.sbatch

# --------------- SBATCH FILE ---------
echo '#!/bin/bash' > $sbatch_file
echo "#SBATCH --job-name=run_partis_seq_data_vijay2021JEM" >> $sbatch_file
echo '#SBATCH --output=out_err_files/run_partis_seq_data_vijay2021JEM.out' >> $sbatch_file
echo '#SBATCH --error=out_err_files/run_partis_seq_data_vijay2021JEM.err' >> $sbatch_file
echo "#SBATCH --time=200:00:00" >> $sbatch_file
echo "#SBATCH --partition=cobey" >> $sbatch_file
echo "#SBATCH --nodes=1" >> $sbatch_file
echo "#SBATCH --ntasks-per-node=6" >> $sbatch_file
echo "#SBATCH --mem-per-cpu=3000" >> $sbatch_file

echo input_file=../data/seq_data_vijay2021JEM/seq_data_vijay2021JEM.fasta >> $sbatch_file
echo output_file=../results/partis/seq_data_vijay2021JEM/seq_data_vijay2021JEM.yaml >> $sbatch_file

# Run partis
echo /project2/cobey/partis/bin/partis partition --n-procs 5 --species mouse --infname '$input_file' --outfname '$output_file' --extra-annotation-columns regional_bounds:cdr3_seqs:seqs_aa:naive_seq_aa:consensus_seq:consensus_seq_aa >> $sbatch_file

# ---------------------------------------

# Find all mouse numbers for the specified time point


sbatch $sbatch_file
rm $sbatch_file
