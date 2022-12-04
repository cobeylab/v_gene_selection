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
    
    for ind in ${sampled_individuals//,/ } # space before closing } is required.
    do
        
        # Run parallel batch jobs for each individual and parameter combination
        # https://rcc.uchicago.edu/docs/running-jobs/srun-parallel/index.html#parallel-batch
    
        # First, create a temporary shell script to be called by the sbatch file
        # This shell file will have a GC number as argument $1.
        
        file_id=$(basename $scenario_dir)_$(basename $param_dir)_individual_$ind
        
        temp_sh_file=temp_sh_files/run_$file_id.sh
        
        echo '#!/bin/bash' > $temp_sh_file
        echo Rscript run_simulations.R $path_to_allele_info $param_dir/ $ind '$1' >> $temp_sh_file
        
        chmod +x $temp_sh_file
    
    
        # --------------------- Create sbatch file for parallel batch --------------------
        # Will do parallel runs of the temporary shell file, ntasks at a time
       

        sbatch_file=run_simulation_${file_id}.sbatch
    
        echo '#!/bin/bash' > $sbatch_file
        echo "#SBATCH --job-name=run_simulations_"${file_id} >> $sbatch_file
        echo "#SBATCH -e out_err_files/run_simulations_"${file_id}.err >> $sbatch_file
        echo "#SBATCH -o out_err_files/run_simulations_"${file_id}.out >> $sbatch_file
        echo "#SBATCH --ntasks="${n_GCs} >> $sbatch_file
        echo "#SBATCH --partition=cobey" >> $sbatch_file
        echo "#SBATCH --time=100:00:00" >> $sbatch_file
        echo "#SBATCH --mem-per-cpu=2G" >> $sbatch_file

        echo module load R/3.6.1 >> $sbatch_file
        echo module load parallel >> $sbatch_file
    
        echo srun='"srun --exclusive -N1 -n1"' >> $sbatch_file
        
        echo parallel='"'parallel --delay 0.2 -j '$SLURM_NTASKS' --joblog out_err_files/$file_id.log --resume'"' >> $sbatch_file
        
        echo '$parallel' '"$srun' ./${temp_sh_file}'"' ::: {1..$n_GCs} >> $sbatch_file
    
        # ----------------------------- Run job  ----------------------------------------
    
        sbatch $sbatch_file
             
        # Remove sbatch file
        rm $sbatch_file
        #rm $temp_sh_file
    done
done














    

