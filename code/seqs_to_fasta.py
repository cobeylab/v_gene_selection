#!/usr/bin/python

"""Given csv file with a mouse's sequences (e.g. 8-1.csv), exports seqs. to temporary fasta file in same directory (8-1.fasta)
"""
import sys
import csv
import os

def main(argv):
    csv_path = str(argv[1])
    parent_dir = os.path.dirname(csv_path)
    mouse_id = os.path.basename(csv_path).replace('.csv','')
    
    fasta_file_path = parent_dir + '/' + mouse_id + '.fasta'
    
    with open(csv_path, 'rU') as f, open(fasta_file_path,'w') as fasta_file:
        original_file = csv.DictReader(f)
        for row in original_file:
            fasta_file.write('>' + row['participant_alt_label'] + '_' + row['trimmed_read_id'])
            fasta_file.write('\n' + row['trimmed_sequence'])
            fasta_file.write('\n')

if __name__ == "__main__":
    status = main(sys.argv)
    sys.exit(status)



