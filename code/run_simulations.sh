#!/bin/bash
scenario_dir=$1 # Path to directory containing allele info file and GC pars file
n_individuals=$2 # Number of individuals to simulate
n_GCs=$3 # Number of germinal centers per individual

path_to_allele_info=${scenario_dir}/allele_info.csv

# Run job array for each individual, with as many jobs as GCs   
for ind in $(seq 1 $n_individuals) 
do
    # Randomly sample a base individual from the 37 mice with 1000 or more naive seqs
    base_ind=$(shuf -i 1-37 -n 1)

    # For all directories representing parameter combinations...
    for param_dir in $scenario_dir/raw_simulation_files/*/
    do
        
        sbatch_file_id=$(basename $scenario_dir)_$(basename $param_dir)_individual_$ind
    
        # -------------------------- Create sbatch file for job array ------------------------
        sbatch_file=run_simulation_${sbatch_file_id}.sbatch
    
        echo '#!/bin/bash' > $sbatch_file
        echo "#SBATCH --job-name=run_simulations_"${sbatch_file_id} >> $sbatch_file
        echo "#SBATCH --nodes=1" >> $sbatch_file
        echo "#SBATCH --ntasks-per-node=1" >> $sbatch_file
	    echo "#SBATCH --partition=cobey" >> $sbatch_file
        echo "#SBATCH -o out_err_files/run_simulations_"${sbatch_file_id}"_%A_%a.out" >> $sbatch_file       
        echo "#SBATCH -e out_err_files/run_simulations_"${sbatch_file_id}"_%A_%a.err" >> $sbatch_file         
        echo "#SBATCH --time=100:00:00" >> $sbatch_file
        echo "#SBATCH --mem-per-cpu=2G" >> $sbatch_file

        echo module load R/3.6.1  >> $sbatch_file
    
        # SLURM_ARRAY_TASK_ID will be an integer specifying a GC
    
        echo Rscript run_simulations.R $path_to_allele_info $param_dir $ind $base_ind '${SLURM_ARRAY_TASK_ID}' >> $sbatch_file
    
        # ----------------------------- Run job array ----------------------------------------
    
        sbatch --array=1-$n_GCs $sbatch_file     
        # Remove sbatch file
        rm $sbatch_file
    done
done
