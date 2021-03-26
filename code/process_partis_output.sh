#!/bin/bash
# Annotates mouse-specific sequence data files with partis results, output annotated seq files

mice_ids=$(ls ../processed_data/mouse_specific_data_files/*csv | grep -o '[0-9]*-[0-9]*' | tr '\n' ','  | tr '\n' ',')
mice_ids=${mice_ids::-1}

for mid in ${mice_ids//,/ }
do
    sbatch_file=process_partis_yaml_$mid.sbatch
    yaml_file=../results/partis/${mid}_partis.yaml 
    mouse_data_file=../processed_data/mouse_specific_data_files/${mid}.csv
    
    echo "Processing mouse $mid"
    
    # ---------------------------- Create sbatch file for job array ----------------------
    echo '#!/bin/bash' > $sbatch_file
    echo "#SBATCH --job-name=process_partis_yaml_$mid" >> $sbatch_file
    echo "#SBATCH --nodes=1" >> $sbatch_file
    echo "#SBATCH --ntasks-per-node=1" >> $sbatch_file
    echo "#SBATCH --partition=cobey" >> $sbatch_file
    echo "#SBATCH -o out_err_files/process_partis_yaml_$mid.out" >> $sbatch_file       
    echo "#SBATCH -e out_err_files/process_partis_yaml_$mid.err" >> $sbatch_file         
    echo "#SBATCH --time=03:00:00" >> $sbatch_file
    echo "#SBATCH --mem-per-cpu=6000" >> $sbatch_file

    echo module load R/3.6.1 >> $sbatch_file
 
    echo Rscript process_partis_output.R ${yaml_file} ${mouse_data_file} >> $sbatch_file
        # ------------------------------------------------------------------------------------

    sbatch $sbatch_file   
    rm $sbatch_file

done



