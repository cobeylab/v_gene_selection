#!/bin/bash

path_to_allele_info=$1 # Read path to directory containing allele info file
n_individuals=$2

sbatch_file_id=$(shuf -i 1-100000 -n 1) # Random string of integers, just to prevent file name conflicts

  
# ---------------------------- Create sbatch file for job array --------------------------
sbatch_file=run_simulation_${sbatch_file_id}.sbatch

echo '#!/bin/bash' > $sbatch_file
echo "#SBATCH --job-name=run_simulations_"${sbatch_file_id} >> $sbatch_file
echo "#SBATCH --nodes=1" >> $sbatch_file
echo "#SBATCH --ntasks-per-node=1" >> $sbatch_file
echo "#SBATCH --partition=cobey" >> $sbatch_file
echo "#SBATCH -o out_err_files/run_simulations_"${sbatch_file_id}"_%A_%a.out" >> $sbatch_file       
echo "#SBATCH -e out_err_files/run_simulations_"${sbatch_file_id}"_%A_%a.err" >> $sbatch_file         
echo "#SBATCH --time=100:00:00" >> $sbatch_file
echo "#SBATCH --mem-per-cpu=4000" >> $sbatch_file

echo module load R/3.6.1 >> $sbatch_file

echo Rscript simulation_model.R $path_to_allele_info '${SLURM_ARRAY_TASK_ID}' >> $sbatch_file
# ----------------------------------------------------------------------------------------


# Run job array
sbatch --array=1-${n_individuals} $sbatch_file     
# Remove sbatch file
rm $sbatch_file
