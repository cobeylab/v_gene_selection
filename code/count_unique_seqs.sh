#!/bin/bash

# List of all mouse csv files, identified by mouse id
ids=$(ls ../processed_data/mouse_specific_data_files/*csv | grep -o '[0-9]*-[0-9]*' | tr '\n' ','  | tr '\n' ',')
ids=${ids::-1}

# Iterate across all mice
for mouse_id in $(echo $ids | sed "s/,/ /g")
do
    output_file=../processed_data/unique_seq_counts_files/${mouse_id}_unique_seqs.csv
    
    if [ -f $output_file ]; then
        echo 'Output file already exists, skipping.'
    else
        echo "Processing mouse $mouse_id"

        sbatch_file=count_unique_seqs_${mouse_id}.sbatch
    
        #----------------- SBATCH FILE FOR RUNNING SCRIPT FOR A SINGLE MOUSE -----------------
        echo '#!/bin/bash' > $sbatch_file
        echo "#SBATCH --job-name=count_unique_seqs_${mouse_id}" >> $sbatch_file
        echo "#SBATCH --output=out_err_files/count_unique_seqs_${mouse_id}.out" >> $sbatch_file
        echo "#SBATCH --error=out_err_files/count_unique_seqs_${mouse_id}.err" >> $sbatch_file
        echo "#SBATCH --time=100:00:00" >> $sbatch_file
        echo "#SBATCH --partition=cobey" >> $sbatch_file
        echo "#SBATCH --nodes=1" >> $sbatch_file
        echo "#SBATCH --ntasks-per-node=9" >> $sbatch_file
        echo "#SBATCH --mem-per-cpu=5000" >> $sbatch_file
    
        echo module load mafft/7.310 >> $sbatch_file
    
    
        mouse_csv_file=../processed_data/mouse_specific_data_files/${mouse_id}.csv
        echo python count_unique_seqs.py $mouse_csv_file $output_file >> $sbatch_file
        # ---------------------------------------
    
        sbatch $sbatch_file
        rm $sbatch_file
    fi


done
