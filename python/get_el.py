"""
get_el.py extracts data from physical_interaction_codes.txt,
experimental_detected_codes.txt, and ensembl gene names.
"""

from typing import List
import sys
from classes import Protein
#from Bio import Entrez
import json

# *Always* tell NCBI who you are
#Entrez.email = "mluo@iwu.edu"

"""
Runs the functions defined below
"""
def main():

    # read ensembl file
    # relative directory
    with open("C:\\Users\\annsb\\OneDrive\\Documents\\PPI-Network-Alignment\\ensembl\\worm_protein_ids106.txt", 'r') as f1:
        content = f1.readlines() # i shouldn't be doing this :))

        # build a list of Protein objects
        objList = []
        for line in content[1:]:
            line = line.split('\t')
            id_map = Protein.Protein(line)
            objList.append(id_map)

    # receive the physical and experimental codes
    phys_exp = retrieve_code()

    # read the BioGRID file
    with open("C:\\Users\\annsb\\OneDrive\\Documents\\PPI-Network-Alignment\\BioGRID\\BIOGRID-ORGANISM-Caenorhabditis_elegans-4.4.219.mitab.txt", 'r') as f2:

        # build a list of interactions
        inter_list = []
        for line in f2:
            # skip the first line
            if line[0] != '#':
                line = line.rstrip()
                line = line.split("\n")

                interaction = get_gene_ids(line, phys_exp)

                if interaction != []:
                    inter_list.append(interaction)

    map_list = id_to_protein(objList, inter_list)

    prots_dict = list_to_dict(map_list)

    retList = dict_to_nodes(prots_dict, objList)

    retList2 = dict_to_edges(prots_dict, retList[1])

    # print(retList[1])
    # print("...")
    # print(retList2)

    retList[0].extend(retList2)

    json_str = json.dumps(retList[0])

    file = open("test_json.json", "w")
    file.write(json_str)

"""
Converts python dictionary to JSON format with established conventions.
Writes to a JSON file, returns nothing.

:param p_dict: dictionary of proteins and their connections as def list
:param objList: list of Protein objects to get alternate ids/names
:return: final list of JSON elements without edges, dict of edge ids and e_ids
"""
def dict_to_nodes(p_dict, objList):
    # take a dictionary with the entrezgene as the key
    d = {x.get_p_sid(): x for x in objList}

    # the list that will be ultimately converted to JSON
    final_list = []

    # dictionary that will be converted to json
    json_dict = {}

    # dictionary that will hold node ids and edges
    edge_dict = {}

    # string that will be added to with every loop iteration
    json_str = ""

    count = 0

    for key in p_dict:
        json_dict = {"data": {}}
        json_dict["data"]["id"] = "n" + str(count)

        # if the id wasn't 'mappable'
        if key[0] == '!':
            # remove the ! from id
            json_dict["data"]["ncbi"] = key[2:]

        else:
            json_dict["data"]["e_id"] = key
            json_dict["data"]["name"] = d[key].get_name()
            json_dict["data"]["parent"] = "unaligned"
            json_dict["data"]["ncbi"] = d[key].get_ncbi()
            json_dict["data"]["uniprot"] = d[key].get_swissprot()## will need to alter this
                                                                  # there is also trembl
            json_dict["classes"] = "species1 unaligned protein"

        # add data to edge dictionary
        edge_dict[key] = json_dict["data"]["id"]

        final_list.append(json_dict)

        count += 1

    return [final_list, edge_dict]

"""
Takes the converted nodes and writes edges to JSON based on established conventions
and data in the established python dictionary

:param p_dict: dictionary of each protein in the network and their connections
:param edge_dict: dictionary of id and its corresponding ensembl id
:return: list of JSON edges
"""
def dict_to_edges(p_dict, edge_dict):
    # final list of JSON edges to be returned
    final_list = []

    temp = {}

    # loop through p_dict
    for ele in p_dict:
        for inter in p_dict[ele]:
            temp = {"data": {}}

            # set the source to the current protein
            temp["data"]["source"] = edge_dict[ele]

            # set the target to the interacting protein
            temp["data"]["target"] = edge_dict[inter]

            #print(temp)

            # add to the final list
            final_list.append(temp)

    return final_list

"""
Constructs a dictionary with a list of each protein it interacts with 
as def

:param map_list: list of mapped protein ids
:return: dict of each protein id and it's interactors
"""
def list_to_dict(map_list):

    prot_hash = {}
    for interaction in map_list:
        current = 0
        for prot in interaction:
            # check if the protein already exists in the dictionary
            if prot not in prot_hash:
                prot_hash[prot] = []    # create a new list if not in dict

            prot2 = interaction[current - 1]    # receive interacting protein (either 0 or -1)

            # check if the interacting protein is already in the current def list
            # and also not in the dict already so that there aren't duplicate interactions
            if prot2 not in prot_hash[prot] and prot2 not in prot_hash:
                # add the interacting protein to the list of interactors for the current
                prot_hash[prot].append(interaction[current - 1])

            current += 1

    return prot_hash

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
Take the id entry and convert it to a form that can be made into a 
Protein object

:param objList: list of protein objects from ensembl
:param id_list: list of str alternate ids from BioGRID
:return: Protein object with corresponding ids as attributes
"""
def id_to_protein(objList, ids):

    # count the number of interactions processed
    count = 0

    # count the amount of unmapped id's
    count_none = 0

    # dict of unmapped id's with name as key and count as def
    none_dict = {}

    # count the valid interactions that weren't filtered
    count_invalid = 0

    # other trackers
    non_phys = 0
    non_exp = 0
    self_loop = 0
    interspecies = 0

    # keeps track of whether the interaction has been filtered or not
    filtered = False

    # take a dictionary with the entrezgene as the key
    d = {x.get_ncbi(): x for x in objList}  # print as a repr

    # also take a dict with the name as the key
    # in case the entrezgene fails
    d2 = {x.get_name(): x for x in objList}

    # take a dict with the swissprot id as the key
    d3 = {x.get_swissprot(): x for x in objList}

    map_list = []
    for id_list in ids:
        mapped = []
        current = ""

        for ele in id_list:
            # key is the first ele in the list
            # treated like a dict but allows multiples
            key = ele[0]

            # count filtered and filter type
            if ele == 'phys':
                non_phys += 1
                count_invalid += 1
                break
            elif ele == 'exp':
                non_exp += 1
                count_invalid += 1
                break
            elif ele == 'intersp':
                interspecies += 1
                count_invalid += 1
                break
            elif ele == 'self':
                self_loop += 1
                count_invalid += 1
                break
            else:

                try:
                    # check if there is an sid for each key
                    current = d[key].get_p_sid()

                except Exception as ex:

                    # if the key does not have an associated sid, loop alternate keys
                    for alts in ele[1]:

                        # check if any of the alts have the gene name
                        try:
                            current = d2[alts].get_p_sid()
                            break

                        # if there is an error, check if the alt is the swissprot id
                        except Exception as ex2:

                            try:
                                current = d3[alts].get_p_sid()
                                break

                            except Exception as ex3:
                                current = "! " + str(key)

                                # check if current isn't already in the dict
                                if current not in none_dict:
                                    none_dict[current] = 0

                                count_none += 1

                                # increase instance of current by 1 in dict
                                none_dict[current] += 1

                # if str(current)[0] == "!":   # for checking how many BioGrid proteins were actually paired
                #     count_none += 1

                # only add actual interactions to the list
                if current != "":
                    filtered = True
                    mapped.append(current)

        # only add unfiltered interactions
        if filtered:
            map_list.append(mapped)
            filtered = False

        count += 1

    ## places output in a form that can be pasted into a spreadsheet ###
    file = open('test_mapping.txt', 'w')
    try:
        file.write("# Total interactions processed: " + str(count) + "\n")
        file.write("# Total filtered interactions: " + str(count_invalid) + "\n")
        file.write("# Percent filtered: " + str(round((count_invalid/count), 3)*100) + "%\n")
        file.write("# Total unmappable id's in file (including duplicates): " + str(count_none) + "\n")
        file.write("# Percent unmappable: " + str(round((count_none/count), 3)*100) + "%\n")
        file.write("# Non-physical interactions filtered out: " + str(non_phys) + "\n")
        file.write("# Non-experimental interactions filtered out: " + str(non_exp) + "\n")
        file.write("# Self-loops filtered out: " + str(self_loop) + "\n")
        file.write("# Interspecies interactions filtered out: " + str(interspecies) + "\n")
        file.write("\n")
        for id in none_dict:
            file.write(str(id) + " occurs " + str(none_dict[id]) + " times " + "\t" + str(round((none_dict[id]/count_none), 3)*100) + "%\n")
        file.write("\n")
        for inter in map_list:
            for id in inter:
                file.write(str(id) + '\t')
            file.write("\n")
    finally:
        file.close()

    return map_list

"""
Extracts data from BioGRID and filters 'bad' interactions

:param line: str line from BioGRID file
:param phys_exp: list of sets of physical and experimental codes
:return: 2 ele dict, protein name as key and alt id list as definition
"""
def get_gene_ids(line, phys_exp) -> List:
    # create a list contains strs after line splitted by tabs
    by_tab = line[0].split('\t')

    # set of physical codes
    phys_codes = phys_exp[0]

    # set of experimental codes
    exp_codes = phys_exp[1]

    # protein names in interaction
    name0 = by_tab[0].split(':')[1]
    name1 = by_tab[1].split(':')[1]

    # self-loop filter
    if name0 == name1:
        return ['self']

    taxid0 = by_tab[9].split(':')[1].strip()
    taxid1 = by_tab[10].split(':')[1].strip()

    # interspecies interaction filter (check taxids)
    if taxid0 != taxid1:
        return ['intersp']

    interaction_type = by_tab[11].split(':')[2][0:4]

    interaction_detection = by_tab[6].split(':')[2][0:4]

    # direct interaction only
    if interaction_type not in phys_codes:
        return ['phys']

    if interaction_detection not in exp_codes:
        return ['exp']

    gene_id = []
    # gene_id = [id0, id1]
    for i in range(0, 2):
        by_colon = by_tab[i].split(':')
        gene_id.append(by_colon[1])

    # if gene_id[0] == gene_id[1]:  is this necessary?
    #     return ['self']

    id_list0 = build_id_list(by_tab[2])
    id_list1 = build_id_list(by_tab[3])

    return [[name0, id_list0], [name1, id_list1]]

"""
Helper function for building a list of ids from the biogrid file 
under the alt ids for interactors A and B

:param line: line of the BioGRID file that includes alternate ids
:return: list of alternate ids in the order that they appear in the file
"""
def build_id_list(line):
    id = line.split("|")

    id_list = []
    for type in id:
        type = type.split(":")
        id_list.append(type[1])

    return id_list

'''
Parallelization code

read_me = sys.argv[1] #  <-- input file   #biogrid file
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

main()
#print(retrieve_code())
#print([len(retrieve_code()[0]), len(retrieve_code()[1])]) what is the reason??


