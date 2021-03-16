#!/usr/bin/python
"""
"""
import sys
import csv
import os
import re
from copy import deepcopy
from dendropy import Tree
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Align.Applications import MafftCommandline
from Bio.Phylo.Applications import RaxmlCommandline
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import *
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import Phylo

sys.setrecursionlimit(10000)

def build_alignment(seqs, naive_seq, temp_dir, temp_file_id):
    aln_seqs = [naive_seq] + seqs
    # Arbitrarily assign sequence ids
    seq_ids = [temp_file_id + '_' + str(i) for i in range(len(seqs)) ]
    aln_ids = ['naive'] + seq_ids
    
    temp_fasta_file =  temp_dir + temp_file_id + '.fasta'
    temp_aln_file = temp_dir + temp_file_id + '_alignment.fasta'

    
    fasta_seqs = [SeqRecord(Seq(aln_seqs[i], generic_dna), id = aln_ids[i]) for
                  i in range(len(aln_seqs))]
        
    # Write temporary fasta file to align
    with open(temp_fasta_file, "w") as output_handle:
        SeqIO.write(fasta_seqs, output_handle, "fasta")
    
        
    # Align with Mafft, output to temp_aln_file
    mafft_command = MafftCommandline(input=temp_fasta_file)
    stdout, stderr = mafft_command()
    
    with open(temp_aln_file, "w") as handle:
        handle.write(stdout)
    
    # Read alignment
    aln = AlignIO.read(temp_aln_file, "fasta")
    os.remove(temp_fasta_file)
    os.remove(temp_aln_file)
    
    return(aln)
    
# Aligns seqs from a subpopulation from the same clone, builds parsimony tree
def build_tree(seqs, naive_seq, temp_dir, temp_file_id):
    #temp_tree_file =  temp_dir + temp_file_id + '_tree.nex'
    
    aln = build_alignment(seqs, naive_seq, temp_dir, temp_file_id)
    
    
    # ----- Get RAXML tree ----- 
    # Write alignment to temporary file
    temp_aln_file = temp_dir + temp_file_id + '_alignment.fasta'
    SeqIO.write(aln, temp_aln_file, "fasta")
    
    raxml_cline = RaxmlCommandline(sequences=temp_aln_file,
                                   model="GTRGAMMA",
                                   name= temp_file_id)
    raxml_cline.program_name = '/project2/cobey/RaxML/raxmlHPC-PTHREADS-AVX'
    raxml_cline.outgroup = 'naive'
    raxml_cline.threads =8
    raxml_cline()
    
    temp_tree_file = 'RAxML_bestTree.' + temp_file_id
    os.remove('RAxML_info.' + temp_file_id)
    os.remove('RAxML_log.' + temp_file_id)
    os.remove('RAxML_parsimonyTree.' + temp_file_id)
    os.remove('RAxML_result.' + temp_file_id)
    tree = Tree.get_from_path(temp_tree_file,'newick')
    
    # Get NJ tree
    #calculator = DistanceCalculator('identity')
    #constructor = DistanceTreeConstructor(calculator, 'nj')
    #nj_tree = constructor.build_tree(aln)
    #Phylo.write(nj_tree, temp_tree_file, 'nexus')
    
    # Get max. parsimony tree
    #scorer = ParsimonyScorer()
    #searcher = NNITreeSearcher(scorer)
    #constructor = ParsimonyTreeConstructor(searcher,)
    #pars_tree = constructor.build_tree(aln)
    # Save tree as newick string to read with dendropy
    #Phylo.write(pars_tree, temp_tree_file, 'nexus')
    
    # Read tree with dendropy
    #tree = Tree.get_from_path(temp_tree_file,'nexus')

    # Re-reoot at naive sequence
    tree.reroot_at_node(tree.find_node_with_taxon_label('naive'))    
       
    os.remove(temp_tree_file)
    os.remove(temp_aln_file)
    
    return [tree, aln]


def collapse_polytomies(tree):
    
    collapsed_tree = deepcopy(tree)
    
    uncolapsed_polys = True
    while uncolapsed_polys is True:

        # Find all internal nodes with a child internal node with br len 0
        focal_nodes = [node for node in collapsed_tree.internal_nodes() if 
                       any([child.is_internal() and child.edge_length < 1e-5 
                            for child in node.child_nodes()])]
        
        if len(focal_nodes) == 0:
            uncolapsed_polys = False
        else:
            # Resolve first node on list, update list before choosing next node
            node  = focal_nodes[0]
            for child in node.child_nodes():
                if child.is_internal() and child.edge_length < 1e-5:
                    # Set all children of child node to node
                    for gdchild in child.child_nodes():
                        node.add_child(gdchild)
                        # Remove spurious child node
                    node.remove_child(child)
    
    return collapsed_tree
            
        
# Counts differences between pair of seqs, excluding sites with gaps in either
def count_seq_diffs(pair_seqs):
    sites = range(len(pair_seqs[0]))
    nongap_sites = [i for i in sites if all([s[i] != '-' for 
                                                         s in pair_seqs])]
    diff_sites = [i for i in nongap_sites if 
                  pair_seqs[0][i] != pair_seqs[1][i]]
    n_diffs = len(diff_sites)
    
    return n_diffs  

# Counts terminal branches after collapsing sister tips that differ <= threshold
# For polytomies, tests seqs. against parent node.
def cluster_tree(tree, aln, threshold):
    # First, collapse any polytomies
    tree = collapse_polytomies(tree)

    effective_n_seqs = 0
        
    for int_node in tree.internal_nodes():
        children = int_node.child_nodes()
        tips_children = [n for n in children if n.is_leaf()]
        n_desc_tips = len(tips_children)
        
        if n_desc_tips == 1:
            effective_n_seqs += 1     
        # If there are two tips as descendants, see if distance > threshold
        elif n_desc_tips == 2:
            pair_ids = [n.taxon.label for n in tips_children] 
            pair_ids = [pid.replace(' ','_') for pid in pair_ids]
                   
            pair_seqs = [record.seq for record in aln if record.id in pair_ids]
            pair_seqs = [str(s) for s in pair_seqs]
                
            n_diffs = count_seq_diffs(pair_seqs)
            if n_diffs <= threshold:
                effective_n_seqs += 1
            else:
                effective_n_seqs +=2
        # If there are more than two tip descendants        
        elif n_desc_tips > 2:
            # Count all descendants with br length 0 as 1
            if any([n.edge_length < 1e-5 for n in tips_children]):
                effective_n_seqs +=1
            
            # For all descendants with br length > 0
            for child in [n for n in tips_children if n.edge_length > 1e-5]:
                n_diffs = child.edge_length * len(aln[1,:])
                if n_diffs > threshold:
                    effective_n_seqs += 1

    return effective_n_seqs
            
                   
# Counts effective number of sequences, calling cluster_tree if 3 or more seqs.
def estimate_uniq_seqs(seqs, naive_seq, threshold, temp_dir, 
                              temp_file_id):
    # Analyze only unique sequences
    seqs = list(set(seqs))
    
    if len(seqs) == 1 or len(seqs) == 0:
        effective_n_seqs = len(seqs)
    elif len(seqs) == 2:
        # If two seqs., align but don't build tree
        aln = build_alignment(seqs, naive_seq, temp_dir, temp_file_id)
        aligned_seqs = [str(rec.seq) for rec in aln if rec.id != 'naive']
        
        n_diffs = count_seq_diffs(aligned_seqs)
        if n_diffs <= threshold:
            effective_n_seqs =  1
        else:
            effective_n_seqs = 2
    else:
        tree, aln = build_tree(seqs, naive_seq, temp_dir, temp_file_id)
        effective_n_seqs = cluster_tree(tree, aln, threshold)
    
    return effective_n_seqs
        
        
def master_function(mouse_csv_path, output_path, germline_V, germline_J, threshold,
                   clonal_method):

    print "Collecting clone info from" + mouse_csv_path
    mouse_id = re.search(r'[0-9]*-[0-9]*', mouse_csv_path).group()

    clone_info = {}

    empty_seq_dic = {'prod_seqs':[], 'unprod_seqs':[], 'prod_seq_ids':[], 
                     'unprod_seq_ids':[]}

    empty_isotype_dic = {"IGA":deepcopy(empty_seq_dic), "IGD":deepcopy(empty_seq_dic),
                         "IGM":deepcopy(empty_seq_dic), "IGG2B":deepcopy(empty_seq_dic),
                         "IGG2A.C":deepcopy(empty_seq_dic), "IGG1":deepcopy(empty_seq_dic),
                         "IGG3":deepcopy(empty_seq_dic), "IGE":deepcopy(empty_seq_dic),
                         "NA": deepcopy(empty_seq_dic)}
    
    
    empty_cell_type_dic = {'GC':deepcopy(empty_isotype_dic),'mem':deepcopy(empty_isotype_dic),
                           'PC':deepcopy(empty_isotype_dic),'naive':deepcopy(empty_isotype_dic)}
    
    empty_clone_dic = {'LN':deepcopy(empty_cell_type_dic), 'BM':deepcopy(empty_cell_type_dic),
                        'spleen': deepcopy(empty_cell_type_dic)}

    fieldnames = ["clone_id", "mouse_id", "tissue", "cell_type", "isotype", "uniq_prod_seqs",
                  "uniq_unprod_seqs", "total_prod_seqs", "total_unprod_seqs"]
    
 
    with open(mouse_csv_path, 'rU') as f:
        mouse_csv = csv.DictReader(f)
        for row in mouse_csv:
            # Get cell type, tissue, clone id, V, D, J genes
            cell_type = row['specimen_cell_subset']
            if cell_type == 'na\xc3\xafve':
                cell_type = 'naive'
            tissue = row['specimen_tissue']
            clone_id = row['clone_id_' + clonal_method]
            
            assert clone_id != ''

            # If 1st time clone is represented, initialize empty dictionary:
            if clone_id not in clone_info.keys():
                clone_info[clone_id] = deepcopy(empty_clone_dic)
                clone_info[clone_id]['v_gene'] = row['v_segment_' + clonal_method]
                clone_info[clone_id]['j_gene'] = row['j_segment_' + clonal_method]
                clone_info[clone_id]['d_gene'] = row['d_segment_' + clonal_method]


            seq = row['trimmed_sequence']
            seq_is_productive = row['productive_' + clonal_method] == 't'
            read_id = row['trimmed_read_id']
            
            # Retain sequence up to end of J gene
            j_seq = row['j_sequence'].upper()
            j_start_position = seq.find(j_seq)
            j_end_position = j_start_position + len(j_seq) - 1
            seq = seq[0:j_end_position+1]
            
            isotype = row['isotype']
            if isotype == '':
                isotype = 'NA'
            
            if seq_is_productive:
                clone_info[clone_id][tissue][cell_type][isotype]['prod_seqs'].append(seq)
                clone_info[clone_id][tissue][cell_type][isotype]['prod_seq_ids'].append(read_id)
            else:
                clone_info[clone_id][tissue][cell_type][isotype]['unprod_seqs'].append(seq)
                clone_info[clone_id][tissue][cell_type][isotype]['unprod_seq_ids'].append(read_id)
    #return clone_info

    # Append info to output file
    with open(output_path, 'a') as f:
        f.write(','.join(fieldnames) + '\n')
        output_csv = csv.DictWriter(f, fieldnames=fieldnames)
        for clone in clone_info:
            print 'Processing clone ' + clone
            for tissue in empty_clone_dic.keys():
                for cell_type in empty_cell_type_dic.keys():
                    for isotype in empty_isotype_dic.keys():
                        
                        prod_seqs = clone_info[clone][tissue][cell_type][isotype]['prod_seqs']
                        unprod_seqs = clone_info[clone][tissue][cell_type][isotype]['unprod_seqs']
                        
                        if len(prod_seqs) > 0 or len(unprod_seqs) > 0:
                            
                            row = {}
                            
                            clone_v = clone_info[clone]['v_gene'].replace('mm','')
                            clone_j = clone_info[clone]['j_gene'].replace('mm','')
                            
                            row['clone_id'], row['mouse_id'] = clone, mouse_id
                            row['tissue'], row['cell_type'] = tissue, cell_type
                            row['isotype'] = isotype
                            
                            row['total_prod_seqs'], row['total_unprod_seqs'] = len(prod_seqs), len(unprod_seqs)
                            
                            # Estimate number of unique seqs (productive only for now)
                            #prod_seq_ids = clone_info[clone][tissue][cell_type]['prod_seq_ids']
                            if clone_v in germline_V.keys() and clone_j in germline_J.keys():
                                naive_V = germline_V[clone_v]
                                naive_J = germline_J[clone_j]
                                naive_seq = (str(naive_V) + str(naive_J)).upper()
                                temp_dir = 'out_err_files/'
                                temp_file_id = mouse_id + '_' + clone
                                n_uniq_prod_seqs = estimate_uniq_seqs(prod_seqs, naive_seq,
                                                                      threshold,temp_dir,
                                                                      temp_file_id)
                            else:
                                if len(prod_seqs) == 1:
                                    n_uniq_prod_seqs = 1
                                else:
                                    n_uniq_prod_seqs = 'NA'
                            
                            row['uniq_prod_seqs'] = n_uniq_prod_seqs
                            # Not computing unique unproductive sequences for now.
                            row['uniq_unprod_seqs'] = 'NA'
                            output_csv.writerow(row)

def main(argv):
    mouse_csv_path = str(argv[1])
    output_path = str(argv[2])
    clonal_method = 'partis'
    threshold = 2
    
    if clonal_method == 'igblast':
        germline_V_path = "../data/sequence_data/reference_V_genes.fasta"
        germline_J_path = "../data/sequence_data/reference_J_genes.fasta"
    else:
        assert clonal_method == 'partis'
        mouse_id = re.search(r'[0-9]*-[0-9]*', mouse_csv_path).group()
        germline_V_path = '../results/partis/partis_germline_genes/v_genes_' + mouse_id + '.fasta'
        germline_J_path = '../results/partis/partis_germline_genes/j_genes_' + mouse_id + '.fasta'
    
    germline_V = {}
    with open(germline_V_path, "rU") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if clonal_method == 'igblast':
                gene_name = record.id.split('|')[1]
            else:
                gene_name = record.id
            germline_V[gene_name] = record.seq
        
    germline_J = {}
    with open(germline_J_path, "rU") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if clonal_method == 'igblast':
                gene_name = record.id.split('|')[1]
            else:
                gene_name = record.id
            germline_J[gene_name] = record.seq
   
    
    master_function(mouse_csv_path, output_path, germline_V, germline_J,
                   threshold, clonal_method)


if __name__ == "__main__":
    status = main(sys.argv)
    sys.exit(status)



