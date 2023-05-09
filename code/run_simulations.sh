#!/bin/bash
scenario_dir=$1 # Path to directory containing allele info file and GC pars file
n_individuals=$2 # Number of individuals to simulate
n_GCs=$3 # Number of germinal centers per individual
n_simultaneous_tasks=$4

path_to_allele_info=${scenario_dir}/allele_info.csv
  
  
# For all directories representing parameter combinations...
for param_dir in $scenario_dir/raw_simulation_files/*/
do

    # Run a parallel batch jobs for each parameter combination
    # Each job will run $n_individuals, ntasks at a time
    # https://rcc.uchicago.edu/docs/running-jobs/srun-parallel/index.html#parallel-batch
    
    # First, create a temporary shell script to be called by the sbatch file
    # This shell file will have a simulated individual id number as argument $1.
    # It will simulate n_GCs for that individual for given parameter combinations
        
    file_id=$(basename $scenario_dir)_$(basename $param_dir)
        
    temp_sh_file=temp_sh_files/run_$file_id.sh
    
    
        
    echo '#!/bin/bash' > $temp_sh_file
    # Skip preexinsing files (allows resuming interrupted jobs)
    echo if [ -e $param_dir/simulation_individual_'$1'.csv ] >> $temp_sh_file
    echo then >> $temp_sh_file
    echo echo 'Output file already exists. Skipping' >> $temp_sh_file
    echo else >> $temp_sh_file
    echo echo "Simulating individual" '$1' >> $temp_sh_file
    echo Rscript run_simulations.R $path_to_allele_info $param_dir/ $n_GCs '$1' >> $temp_sh_file
    echo fi >> $temp_sh_file
        
    chmod +x $temp_sh_file
        
    # --------------------- Create sbatch file for parallel batch --------------------
    # Will do parallel runs of the temporary shell file, ntasks at a time
       

    sbatch_file=run_simulation_${file_id}.sbatch
    
    echo '#!/bin/bash' > $sbatch_file
    echo "#SBATCH --job-name=run_simulations_"${file_id} >> $sbatch_file
    echo "#SBATCH -e out_err_files/run_simulations_"${file_id}.err >> $sbatch_file
    echo "#SBATCH -o out_err_files/run_simulations_"${file_id}.out >> $sbatch_file
    # RUN n_simultaneous_tasks individuals at a time
    echo "#SBATCH --ntasks="${n_simultaneous_tasks} >> $sbatch_file 
    echo "#SBATCH --partition=cobey" >> $sbatch_file
    echo "#SBATCH --time=100:00:00" >> $sbatch_file
    echo "#SBATCH --mem-per-cpu=2G" >> $sbatch_file

    echo module load R/3.6.1 >> $sbatch_file
    echo module load parallel >> $sbatch_file
    
    echo ulimit -n 10000 >> $sbatch_file

    echo srun='"srun --exclusive -N1 -n1"' >> $sbatch_file
        
    echo parallel='"'parallel --delay 1 -j '$SLURM_NTASKS' --joblog out_err_files/$file_id.log --resume'"' >> $sbatch_file
        
    echo '$parallel' '"$srun' ./${temp_sh_file}'"' ::: {1..$n_individuals} >> $sbatch_file
    
    # ----------------------------- Run job  ----------------------------------------
    
    sbatch $sbatch_file
             
    # Remove sbatch file
    rm $sbatch_file

done
