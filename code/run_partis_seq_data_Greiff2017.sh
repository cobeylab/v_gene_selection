#!/bin/bash
for mouse_dir in ../data/seq_data_Greiff2017/ERR*
do

    dataset_name=$(basename $mouse_dir)
    sbatch_file=run_partis_${dataset_name}.sbatch
    
    echo "#!/bin/bash
#SBATCH --job-name=run_partis_${dataset_name}
#SBATCH --output=out_err_files/run_partis_${dataset_name}.out
#SBATCH --error=out_err_files/run_partis_${dataset_name}.err
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=3000
#SBATCH --time=200:00:00

input_file=${mouse_dir}/${dataset_name}_processed_reads.fasta
output_file=../results/partis/seq_data_Greiff2017/${dataset_name}.yaml" > $sbatch_file


    
echo '/project2/cobey/partis/bin/partis partition --n-procs 16 --species mouse --infname $input_file --outfname $output_file --extra-annotation-columns regional_bounds:cdr3_seqs:seqs_aa:naive_seq_aa:consensus_seq:consensus_seq_aa

' >> $sbatch_file

sbatch $sbatch_file
rm $sbatch_file
    
done