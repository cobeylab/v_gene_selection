#!/usr/bin/python
"""
"""
import csv
import sys
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import pairwise2

# Read primer sequences
primers = {}
with open('../data/sequence_data/primer_sequences.txt') as f:
    f.next()
    for row in f:
        primers[row.split('\t')[0]] = row.split('\t')[1].replace('\n','')
        
# Read reference constant region sequences
reference_C_regions = {}
with open('../data/sequence_data/reference_C_regions.txt') as f:
    for row in f:
        reference_C_regions[row.split(':')[1].replace('\n','')] = row.split(':')[0].upper()
        
#data_csv_path = '../data/sequence_data/mouse_specific_data_files/8-1.csv'


def match_reference_C(obs_C_region, reference_C_regions):
    """"
    Returns the reference C region that best matches observed region
    """
    match = {}
    for ref in reference_C_regions:
        ref_seq = reference_C_regions[ref]   
        #print ref
        # Match score = 1, mismatch score = -3, gap open = -10, gap ext = -1
        alignment = pairwise2.align.globalms(obs_C_region, ref_seq,1,-3,-10,-1)
        #print(format_alignment(*alignment[0]))

        aln_obs = alignment[0][0]
        aln_ref = alignment[0][1]
        non_gap_sites = [i for i in range(len(aln_obs)) if 
                   aln_obs[i]!='-' and aln_ref[i]!='-']
        diffs = [i for i in non_gap_sites if aln_obs[i] != aln_ref[i]]
        
        match[ref] = {'score': alignment[0][2], 'diffs':len(diffs),
             'n_bases':len(non_gap_sites)}
    
    maxscore = max([match[ref]['score'] for ref in match.keys()])
    
    # First, find best match by minimizing differences
    best_match = [key for key in match.keys() if match[key]['score'] == maxscore]
    
    # If more than 1 best match, pick the ref with most aligned sites without gaps
    if len(best_match) > 1:
        largest_n_bases = max([match[ref]['n_bases'] for ref in best_match])
        best_match = [ref for ref in best_match if 
                      match[ref]['n_bases'] == largest_n_bases]
        best_match = best_match[0]
    else:
        best_match = best_match[0]
        
    return {'best_match' : best_match, 'errors' : match[best_match]['diffs'],
            'n_bases' : match[best_match]['n_bases']}


def get_errors_in_dataset(mouse_csv_path, output_path):
    n_errors = 0
    n_bases = 0
    
    # Get mouse id from mouse csv path
    mouse_id = mouse_csv_path.split('/')[-1].replace('.csv','')
    

    # To compare with isotype assignments that Scott, Tho et al. did
    isotype_disagreements = 0
    seqs_with_preassigned_isotypes = 0
    

    with open(mouse_csv_path, 'rU') as f:
        data_csv = csv.DictReader(f)
        for row in data_csv:
                f_primer_id = row["forward_primer"]
                r_primer_id = row["reverse_primer"]
    
                f_primer_seq = primers[f_primer_id]
                r_primer_seq = Seq(primers[r_primer_id], generic_dna)
                r_primer_seq = str(r_primer_seq.reverse_complement())
                
                
                # Get reverse complement of demultiplexed sequence
                demux_seq = Seq(row['demuxed_sequence'], generic_dna)
                demux_seq = demux_seq.reverse_complement()
                
                # Trim off forward and reverse primers
                seq = str(demux_seq).replace(f_primer_seq, '')
                seq = seq.replace(r_primer_seq,'')
                
                # To find non-primer encoded constant region, take seq. past J gene
                j_seq = row['j_sequence'].upper()
                j_start_position = seq.find(j_seq)
                j_end_position = j_start_position + len(j_seq) - 1
                
                # Observed constant region (including last nucleotide of J)
                obs_C_region = seq[j_end_position:len(seq)]
                match = match_reference_C(obs_C_region, reference_C_regions)
                
                n_errors = n_errors + match['errors']
                n_bases = n_bases + match['n_bases']
                
                if len(row['isotype']) > 0:
                    seqs_with_preassigned_isotypes += 1
                    if match['best_match'].split('_')[0] != row['isotype']:
                        #1/0
                        isotype_disagreements += 1
                        #print match['best_match'].split('_')[0] + ',' + row['isotype']
                        
    error_rate = float(n_errors) / n_bases
    isotype_disagreement_rate = float(isotype_disagreements)/seqs_with_preassigned_isotypes                     
    
    with open(output_path, 'w') as f:
        f.write('mouse_id,errors,bases_analyzed,error_rate,seqs_with_preassigned_isotypes,isotype_disagreements,isotype_disagreement_rate\n')
        results = [mouse_id,n_errors,n_bases,error_rate,seqs_with_preassigned_isotypes,
                   isotype_disagreements,isotype_disagreement_rate]
        results = [str(result) for result in results]
        f.write(','.join(results))
        f.write('\n')
                    
def main(argv):
    mouse_csv_path = str(argv[1])
    output_path = str(argv[2])
    
    get_errors_in_dataset(mouse_csv_path, output_path)
   
if __name__ == "__main__":
    status = main(sys.argv)
    sys.exit(status)
         
            
            
            
            
            
            
            
            
            
            
            
            





