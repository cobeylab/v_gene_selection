#!/bin/bash
scenario_dir=$1 # Path to directory containing allele info file and GC pars file
n_individuals=$2 # Number of individuals to simulate

path_to_allele_info=${scenario_dir}/allele_info.csv
path_to_model_parameters=${scenario_dir}/model_parameters.csv


sbatch_file_id=$(basename $scenario_dir)

  
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
echo "#SBATCH --mem-per-cpu=15000" >> $sbatch_file

echo module load R/3.6.1 >> $sbatch_file

echo Rscript simulation_run.R $path_to_allele_info $path_to_model_parameters '${SLURM_ARRAY_TASK_ID}' >> $sbatch_file

# ----------------------------------------------------------------------------------------


# Run job array (1 job per individual)
# Randomly samples individuals from the 40 mice with 1000 or more naive seqs to use
# as basis of simulation

sampled_individuals=$(shuf -i 1-40 -n $n_individuals | tr '\n' ',')
sampled_individuals=${sampled_individuals::-1}


sbatch --array=$sampled_individuals $sbatch_file     
# Remove sbatch file
rm $sbatch_file

    

