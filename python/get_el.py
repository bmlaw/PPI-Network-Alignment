"""
get_el.py extracts data from physical_interaction_codes.txt,
experimental_detected_codes.txt, and ensembl gene names.
"""

from typing import List
import sys
#from Bio import Entrez

# *Always* tell NCBI who you are
#Entrez.email = "mluo@iwu.edu"

"""
Runs the functions defined below
"""
def main():
    species1 = []

    # read ensembl file
    with open("worm_protein_ids106.txt", 'r') as f1:
        content = f1.readlines()

        for line in content:
            species1.append(get_gene_ids(line))

    print(species1)

"""
Extracts data from physical and experimental interaction files 
and places them into two respective sets.

:return: list with two values [phys, exp]
"""
def retrieve_code():
    # open both phys and exp files to read
    physical_codes = open("physical_interaction_codes.txt", 'r')
    experimental_codes = open("experimental_detected_codes.txt", 'r')

    # sets (separate from files) where data from files will be stored
    physical_set = set()
    experimental_set = set()

    # add all the physical codes to the physical code set
    for line in physical_codes:
        code = line.split(':')[1].strip()
        physical_set.add(code)

    # add all the experimental codes to the experimental code set
    for line in experimental_codes:
        code = line.split(':')[1].strip()
        experimental_set.add(code)

    physical_codes.close()
    experimental_codes.close()

    # return the data from each set
    return [physical_set, experimental_set]

"""
Extracts data from ensembl to check if the interaction is phys/exp

:param line: str line from ensembl file
:return: 2D list, each element having the name and gene id of genes in interaction
"""
def get_gene_ids(line: str) -> List:
    # create a list contains strs after line splitted by tabs
    by_tab = line.split('\t')

    name0 = by_tab[2].split(':')[2].split('|')[0]
    name1 = by_tab[3].split(':')[2].split('|')[0]

    # self-loop filter
    if name0 == name1:
        return []

    taxid0 = by_tab[9].split(':')[1].strip()
    taxid1 = by_tab[10].split(':')[1].strip()

    # interspecies interaction filter (check taxids)
    if taxid0 != taxid1:
        return []

    interaction_type = by_tab[11].split(':')[2][0:4]

    # direct interaction only
    if interaction_type != "0407":
        return []

    gene_id = []
    # gene_id = [id0, id1]
    for i in range(0, 2):
        by_colon = by_tab[i].split(':')
        gene_id.append(by_colon[1])

    if gene_id[0] == gene_id[1]:
        return []

    return [[name0, gene_id[0]], [name1, gene_id[1]]]

'''
Parallelization code

read_me = sys.argv[1] #  <-- input file
el_file = sys.argv[2] #  <-- output .el
el_n_file = sys.argv[3] # <-- output .el file with only number ids



interaction_set = set()

#{id0 - id1, id0 - id2}
#{id0 -> {id1, id2}}

# suggested: put the interactions in alphabetical order, both vertically and horizontally.

with open(read_me) as r:
    with open(el_file, 'w') as e:
        with open(el_n_file, 'w') as k:
            r.readline()
            for line in r:
                ids = get_gene_ids(line)

                if (len(ids) > 0) and (ids[0][0] + '\t' + ids[1][0] not in interaction_set) and (ids[1][0] + '\t' + ids[0][0] not in interaction_set):


                    interaction_set.add(ids[0][0] + '\t' + ids[1][0])
                    interaction_set.add(ids[1][0] + '\t' + ids[0][0])

                    k.write(ids[0][1] + '\t' + ids[1][1] + '\n')
                    e.write(ids[0][0] + '_' + ids[0][1] + '\t' + ids[1][0] + '_' + ids[1][1] + '\n')
'''

#print(retrieve_code())
#print([len(retrieve_code()[0]), len(retrieve_code()[1])]) what is the reason??


