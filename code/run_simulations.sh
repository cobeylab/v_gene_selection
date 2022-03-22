#!/bin/bash
scenario_dir=$1 # Path to directory containing allele info file and GC pars file
n_individuals=$2 # Number of individuals to simulate
n_GCs=$3 # Number of germinal centers per individual

path_to_allele_info=${scenario_dir}/allele_info.csv


# Randomly samples individuals from the 40 mice with 1000 or more naive seqs to use
# as basis of simulation

sampled_individuals=$(shuf -i 1-40 -n $n_individuals | tr '\n' ',')
sampled_individuals=${sampled_individuals::-1}


# For all directories representing parameter combinations...
for param_dir in $scenario_dir/raw_simulation_files/*/
do
    # Run job array for each individual, with as many jobs as GCs 

    for ind in ${sampled_individuals//,/ } # space before closing } is required.
    do
        sbatch_file_id=$(basename $scenario_dir)_$(basename $param_dir)_individual_$ind
    
        # -------------------------- Create sbatch file for job array ------------------------
        sbatch_file=run_simulation_${sbatch_file_id}.sbatch
    
        echo '#!/bin/bash' > $sbatch_file
        echo "#SBATCH --job-name=run_simulations_"${sbatch_file_id} >> $sbatch_file
        echo "#SBATCH --nodes=1" >> $sbatch_file
        echo "#SBATCH --ntasks-per-node=1" >> $sbatch_file
        echo "#SBATCH --partition=broadwl" >> $sbatch_file
        echo "#SBATCH -o out_err_files/run_simulations_"${sbatch_file_id}"_%A_%a.out" >> $sbatch_file       
        echo "#SBATCH -e out_err_files/run_simulations_"${sbatch_file_id}"_%A_%a.err" >> $sbatch_file         
        echo "#SBATCH --time=00:15:00" >> $sbatch_file
        echo "#SBATCH --mem-per-cpu=2000" >> $sbatch_file

        echo module load R/3.6.1 >> $sbatch_file
    
        # SLURM_ARRAY_TASK_ID will be an integer specifying a GC
    
        echo Rscript simulation_run.R $path_to_allele_info $param_dir $ind '${SLURM_ARRAY_TASK_ID}' >> $sbatch_file
    
        # ----------------------------- Run job array ----------------------------------------
    
        sbatch --array=1-$n_GCs $sbatch_file     
        # Remove sbatch file
        rm $sbatch_file
    done
done














    

