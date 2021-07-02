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
#SBATCH --ntasks-per-node=28
#SBATCH --mem-per-cpu=2000
#SBATCH --time=100:00:00
    
module load python/anaconda-2020.02
cd $mouse_dir

output_identifier=$dataset_name

# Iterate over .gz files containing the fastq files
for f in ./*.gz" > $sbatch_file
    
    echo '
    do
        gz_file=$(basename $f)
        
        # This will be the name of the fastq file after the gz is decompressed
        fastq_file=${gz_file%.gz}         
        
        # Gunzip into standard output, pipe into fastq_file (to keep original .gz)
        gunzip -c $gz_file > $fastq_file
        
    done
    
fastq_1=$(ls *_1.fastq) # File indexed 1 contains the forward read starting from V region
fastq_2=$(ls *_2.fastq) # File 2 contains the reverse read starting from constant region
        
# Paired-end assembly
AssemblePairs.py align -1 $fastq_1 -2 $fastq_2 --coord sra --rc tail --outname $output_identifier --log AP.log

ParseLog.py -l AP.log -f ID LENGTH OVERLAP ERROR PVALUE

# Filter low-quality seqs
FilterSeq.py quality -s ${output_identifier}_assemble-pass.fastq -q 20 --outname $output_identifier --log FS.log

ParseLog.py -l FS.log -f ID QUALITY

# ---- Mask primers

# Forward (V regoin) primers

MaskPrimers.py score -s ${output_identifier}_quality-pass.fastq -p ../V_primers.fasta \
    --start 4 --mode cut --pf VPRIMER --outname ${output_identifier}-FWD --log MPV.log

# Reverse (C region) primers. The --fasta flag forces output to be a fasta file

MaskPrimers.py score -s ${output_identifier}-FWD_primers-pass.fastq -p ../C_primers.fasta \
    --start 4 --mode cut --revpr --pf CPRIMER --outname ${output_identifier}-REV --log MPC.log --fasta

ParseLog.py -l MPV.log MPC.log -f ID PRIMER ERROR

# The fasta file with the processed sequences:
mv ${output_identifier}-REV_primers-pass.fasta ${output_identifier}_processed_reads.fasta 

# Create annotation table
#ParseHeaders.py table -s ${output_identifier}_processed_reads.fasta  -f ID CPRIMER VPRIMER

# Remove intermediate files
rm *pass.fastq*
rm *.log
  
# Remove decompressed fastq files with raw reads (leave only .gz files to save space)
rm $fastq_1
rm $fastq_2' >> $sbatch_file

sbatch $sbatch_file
rm $sbatch_file
    
done