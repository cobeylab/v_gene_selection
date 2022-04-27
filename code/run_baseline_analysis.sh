#!/bin/bash
# Runs baseline analysis for each mouse

mice_ids=$(ls ../processed_data/mouse_specific_data_files/*csv | grep -o '[0-9]*-[0-9]*' | tr '\n' ','  | tr '\n' ',')
mice_ids=${mice_ids::-1}

for mid in ${mice_ids//,/ }
do
    sbatch_file=run_baseline_$mid.sbatch
        
    # ---------------------------- Create sbatch file for job array ----------------------
    echo '#!/bin/bash' > $sbatch_file
    echo "#SBATCH --job-name=run_baseline_$mid" >> $sbatch_file
    echo "#SBATCH --nodes=1" >> $sbatch_file
    echo "#SBATCH --ntasks-per-node=11" >> $sbatch_file
    echo "#SBATCH --partition=cobey" >> $sbatch_file
    echo "#SBATCH --account=pi-cobey" >> $sbatch_file
    echo "#SBATCH -o out_err_files/run_baseline_$mid.out" >> $sbatch_file       
    echo "#SBATCH -e out_err_files/run_baseline_$mid.err" >> $sbatch_file         
    echo "#SBATCH --time=20:00:00" >> $sbatch_file
    echo "#SBATCH --mem-per-cpu=4000" >> $sbatch_file

    echo module load R/3.6.1 >> $sbatch_file
    
    # Run analysis with 10 cores
    echo Rscript run_baseline_analysis.R ${mid} 10 >> $sbatch_file
        #------------------------------------------------------------------------------------

    sbatch $sbatch_file   
    rm $sbatch_file

done



