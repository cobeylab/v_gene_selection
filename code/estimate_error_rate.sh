#!/bin/bash

# List of all mouse csv files, identified by mouse id
ids=$(ls ../processed_data/mouse_specific_data_files/*-*.csv | grep -o "[0-9]*\-[0-9]*" | tr '\n' ',' )
ids=${ids::-1}

# Create output directory
mkdir -p ../results/error_rate/

# Iterate across all mice
for mouse_id in $(echo $ids | sed "s/,/ /g")
do
    sbatch_file=error_rate_$mouse_id.sbatch
    
    #----------------- SBATCH FILE FOR RUNNING SCRIPT FOR A SINGLE CLONE -----------------
    echo '#!/bin/bash' > $sbatch_file
    echo "#SBATCH --job-name=error_rate_${mouse_id}" >> $sbatch_file
    echo "#SBATCH --output=out_err_files/error_rate_${mouse_id}.out" >> $sbatch_file
    echo "#SBATCH --error=out_err_files/error_rate_${mouse_id}.err" >> $sbatch_file
    echo "#SBATCH --time=20:00:00" >> $sbatch_file
    echo "#SBATCH --partition=cobey" >> $sbatch_file
    echo "#SBATCH --nodes=1" >> $sbatch_file
    echo "#SBATCH --ntasks-per-node=1" >> $sbatch_file
    echo "#SBATCH --mem-per-cpu=5000" >> $sbatch_file
    
    output_file=../results/error_rate/${mouse_id}_error_rate.csv

    echo mouse_csv_file=../processed_data/mouse_specific_data_files/${mouse_id}.csv >> $sbatch_file
    echo python estimate_error_rate.py '$mouse_csv_file' $output_file >> $sbatch_file
    # ---------------------------------------
    
    sbatch $sbatch_file
    rm $sbatch_file

done
