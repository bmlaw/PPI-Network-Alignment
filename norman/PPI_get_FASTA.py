#!/python
import sys
import operator
import Bio.Entrez as BE
import Bio.SeqIO as SeqIO


'''
reads a .el file and return a set of all nodes

@param: filename of .el file

@return: a set (complete and non-repeating) of all nodes (gene ids)
'''

'''
def readPPI(fileName):
    ofile = open(fileName, "r")

    vertices = set()
    
    for line in ofile:
        nodes = line.split('\t')
        for node in nodes:
            node_n = node.split('_')[1]
            vertices.add(node_n)
    
    ofile.close()
    return vertices
'''



'''
def readPPI(fileName):
    ofile = open(fileName, "r")

    vertices = set()
    id_dict = dict()
    return_me = []
    
    for line in ofile:
        nodes = line.strip().split('\t')
        for node in nodes:
            node_n = node.split('_')[1]
            vertices.add(node_n)
            id_dict[node_n] = node.split('_')[0]

    return_me.append(vertices)
    return_me.append(id_dict)
    
    ofile.close()
    return return_me
'''



'''
Take all the gene ids read from a .el file and output proper protein seq corresponds
to each gene. Longer, ref seq has higher priority to be chosen

@param: filename of the output fasta file
@param: set of all entrez_ids (aka gene ids)

@output: prints
'''

'''
# version that does not work for AThaliana
# [vertices_set, name_dict, new_el_line_list]
def readPPI(fileName):
    ofile = open(fileName, "r")

    vertices = set()
    id_dict = dict()
    new_el_line_list = []
    
    return_me = []
    
    for line in ofile:
        nodes = line.strip().split('\t')
        for node in nodes:
            node_n = node.split('_')[1]
            vertices.add(node_n)
            id_dict[node_n] = node.split('_')[0]

        new_el_line_list.append([nodes[0].split('_')[1], nodes[1].split('_')[1]])

    return_me.append(vertices) # <-- avoid vertex duplication
    return_me.append(id_dict)   # <-- naming
    return_me.append(new_el_line_list) # <-- potential lines to be printed
    
    ofile.close()
    return return_me
'''

def readPPI(fileName):
    ofile = open(fileName, "r")

    vertices = set()
    id_dict = dict()
    new_el_line_list = []
    
    return_me = []
    
    for line in ofile:
        nodes = line.strip().split('\t')
        for node in nodes:
            node_n = node.split('_')[-1]
            vertices.add(node_n)
            #id_dict[node_n] = node.split('_')[0: len(node.split('_')) - 1]
            
            before_id = node.split('_')[0: len(node.split('_')) - 1]
            name = ""
            for part in before_id:
                name += part
                name += '_'
            name = name[:-1]
            id_dict[node_n] = name
            

        new_el_line_list.append([nodes[0].split('_')[-1], nodes[1].split('_')[-1]])

    return_me.append(vertices) # <-- avoid vertex duplication
    return_me.append(id_dict)   # <-- naming
    return_me.append(new_el_line_list) # <-- potential lines to be printed
    
    ofile.close()
    return return_me

def fasta_from_entrezid(fileName, new_el_filename, new_el_with_names_filename, entrez_ids):
    ofile = open(fileName, "w")
    new_el = open(new_el_filename, "w")
    new_el_with_names = open(new_el_with_names_filename, "w")
    count = 1

    presence_set = set()
    
    for interaction in entrez_ids[2]:
        valid = True
        for eid in interaction:
            if eid not in presence_set:
                try:
                    
                    print("processing ", count, " / ", str(len(entrez_ids[0])) + ':', eid)

                    # find all the protein ids of the current gene
                    res = BE.elink(db='protein', dbfrom='gene', id=eid)
                    record = BE.read(res)
                    res.close()
                    
                    # Get a list of dict that contains 'Link', 'DbTo', 'LinkName' as keys
                    linkset = record[0]['LinkSetDb']
                    
                    # locates the dict with LinkName gene_protein, initialized to -1
                    normal_pos = -1
                    # locates the dict with LinkName gene_protein_ref_seq, initialized to -1
                    ref_pos = -1
                    
                    # find the dict of gene proteins seqs and refseq(s)
                    for i in range(0, len(linkset)):
                        if linkset[i]['LinkName'] == "gene_protein":
                            normal_pos = i
                        elif linkset[i]['LinkName'] == "gene_protein_refseq":
                            ref_pos = i
                    
                    # holder for the seq we choose to output
                    chosen_seq = ''
                    
                    # if the ref seq(s) exists
                    if ref_pos != -1:
                        # find and keep the longest ref seq
                        for id_dict in linkset[ref_pos]['Link']:
                            pid = id_dict['Id']
                            handle = BE.efetch(db='protein', id=pid, rettype='fasta')
                            seq = str(SeqIO.read(handle, 'fasta').seq)
                            
                            if len(seq) > len(chosen_seq):
                                chosen_seq = seq
                    # otherwise, find and keep the longest seq
                    else:
                        for id_dict in linkset[normal_pos]['Link']:
                            pid = id_dict['Id']
                            handle = BE.efetch(db='protein', id=pid, rettype='fasta')
                            seq = str(SeqIO.read(handle, 'fasta').seq)
                            
                            if len(seq) > len(chosen_seq):
                                chosen_seq = seq
                    

                    # output the chosen protein seq
                    #ofile.write('>' + eid + '_' + entrez_ids[1][eid] + '\n')
                    ofile.write('>' + entrez_ids[1][eid] + '_' + eid + '\n')
                    ofile.write("%s\n"%(chosen_seq))

                    

 
                
                except Exception as e:
                    print(e)
                    print("seq: ", eid, "not found")

                    valid = False

                presence_set.add(eid)
            
                count = count + 1

        if valid:
            new_el.write(interaction[0] + '\t' + interaction[1] + '\n')
            new_el_with_names.write(entrez_ids[1][interaction[0]] + '_' + interaction[0] + '\t' + entrez_ids[1][interaction[1]] + '_' + interaction[1] + '\n')

                
        
    ofile.close()
    new_el.close()
    new_el_with_names.close()

BE.email = 'mluo@iwu.edu'

PPI_File = sys.argv[1] #  <-- input .el PPI network
Fasta_File = sys.argv[2] #  <-- output .fasta file
NEW_el_File = sys.argv[3] # <-- output new .el file (only number ids)
NEW_el_File_with_names = sys.argv[4] # <-- output new .el file (name_id)


print("Loading entrez-ids from PPI file", PPI_File)

# get a set of gene ids read from a PPI file
entrez_list = readPPI(PPI_File)


print(len(entrez_list[0]))

print("Retrieving FASTA sequences", "\n")

# output all protein seqs correspond to the gene ids
nodes = fasta_from_entrezid(Fasta_File, NEW_el_File, NEW_el_File_with_names, entrez_list)




# TODO:
# ADD COMMENTS, RENAME VARIABLES MEANINGFULLY
# Run on all species and check that all get a fetched protein sequence
# Spot-check protein FASTA files, make sure they downloaded properly
# Maybe TODO:
# Always pull the protein ID from the gene_protein_refseq link if it exists
# If there is more than 1 protein link, use the LONGEST protein sequence
# Can this be sped up? I.e. fetch all the elinks at once rather than submitting as separate requests
# Similarly for fetching multiple FASTAs, if we need to, and/or determining which protein is longest without fetching all of them
# Next:
# Change script to submit jobs to run in parallel
# Add all network species from networks
# Run scripts
# 


'''
# [vertices_set, name_dict, new_el_line_list]
def readPPI(fileName):
    ofile = open(fileName, "r")

    vertices = set()
    id_dict = dict()
    new_el_line_list = []
    added_el = []
    
    return_me = []
    
    for line in ofile:
        nodes = line.strip().split('\t')
        for node in nodes:
            node_n = node.split('_')[1]
            vertices.add(node_n)
            id_dict[node_n] = node.split('_')[0]

        new_edge = [nodes[0].split('_')[1], nodes[1].split('_')[1]]
        new_edge_reversed = [nodes[1].split('_')[1], nodes[0].split('_')[1]]
        if new_edge not in added_el:
            new_el_line_list.append(new_edge)
            added_el.append(new_edge_reversed)
            added_el.append(new_edge)

    return_me.append(vertices) # <-- avoid vertex duplication
    return_me.append(id_dict)   # <-- naming
    return_me.append(new_el_line_list) # <-- potential lines to be printed
    
    ofile.close()
    return return_me

def fasta_from_entrezid(fileName, new_el_filename, new_el_with_names_filename, entrez_ids):
    ofile = open(fileName, "w")
    new_el = open(new_el_filename, "w")
    new_el_with_names = open(new_el_with_names_filename, "w")
    count = 1

    presence_set = set()
    
    for interaction in entrez_ids[2]:
        valid = True
        for eid in interaction:
            if eid not in presence_set:
                try:
                    
                    print("processing ", count, " / ", str(len(entrez_ids[0])) + ':', eid)

                    # find all the protein ids of the current gene
                    res = BE.elink(db='protein', dbfrom='gene', id=eid)
                    record = BE.read(res)
                    res.close()
                    
                    # Get a list of dict that contains 'Link', 'DbTo', 'LinkName' as keys
                    linkset = record[0]['LinkSetDb']
                    
                    # locates the dict with LinkName gene_protein, initialized to -1
                    normal_pos = -1
                    # locates the dict with LinkName gene_protein_ref_seq, initialized to -1
                    ref_pos = -1
                    
                    # find the dict of gene proteins seqs and refseq(s)
                    for i in range(0, len(linkset)):
                        if linkset[i]['LinkName'] == "gene_protein":
                            normal_pos = i
                        elif linkset[i]['LinkName'] == "gene_protein_refseq":
                            ref_pos = i
                    
                    # holder for the seq we choose to output
                    chosen_seq = ''
                    
                    # if the ref seq(s) exists
                    if ref_pos != -1:
                        # find and keep the longest ref seq
                        for id_dict in linkset[ref_pos]['Link']:
                            pid = id_dict['Id']
                            handle = BE.efetch(db='protein', id=pid, rettype='fasta')
                            seq = str(SeqIO.read(handle, 'fasta').seq)
                            
                            if len(seq) > len(chosen_seq):
                                chosen_seq = seq
                    # otherwise, find and keep the longest seq
                    else:
                        for id_dict in linkset[normal_pos]['Link']:
                            pid = id_dict['Id']
                            handle = BE.efetch(db='protein', id=pid, rettype='fasta')
                            seq = str(SeqIO.read(handle, 'fasta').seq)
                            
                            if len(seq) > len(chosen_seq):
                                chosen_seq = seq
                    

                    # output the chosen protein seq
                    #ofile.write('>' + eid + '_' + entrez_ids[1][eid] + '\n')
                    ofile.write('>' + entrez_ids[1][eid] + '_' + eid + '\n')
                    ofile.write("%s\n"%(chosen_seq))

                    

 
                
                except Exception as e:
                    print(e)
                    print("seq: ", eid, "not found")

                    valid = False

                presence_set.add(eid)
            
                count = count + 1

        if valid:
            new_el.write(interaction[0] + '\t' + interaction[1] + '\n')
            new_el_with_names.write(entrez_ids[1][interaction[0]] + '_' + interaction[0] + '\t' + entrez_ids[1][interaction[1]] + '_' + interaction[1] + '\n')

                
        
    ofile.close()
    new_el.close()
    new_el_with_names.close()

BE.email = 'mluo@iwu.edu'

PPI_File = sys.argv[1] #  <-- input .el PPI network
Fasta_File = sys.argv[2] #  <-- output .fasta file
NEW_el_File = sys.argv[3] # <-- output new .el file (only number ids)
NEW_el_File_with_names = sys.argv[4] # <-- output new .el file (name_id)


print("Loading entrez-ids from PPI file", PPI_File)

# get a set of gene ids read from a PPI file
entrez_list = readPPI(PPI_File)


print(len(entrez_list[0]))


print("Retrieving FASTA sequences", "\n")

# output all protein seqs correspond to the gene ids
nodes = fasta_from_entrezid(Fasta_File, NEW_el_File, NEW_el_File_with_names, entrez_list)
'''
