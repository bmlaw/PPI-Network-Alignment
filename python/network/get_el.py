"""
get_el.py extracts data from physical_interaction_codes.txt,
experimental_detected_codes.txt, and ensembl gene names.
"""
__author__ = "Anna Sheaffer"
__email__ = "asheaffe@iwu.edu"
__credits__ = ["Norman Luo", "Brian Law"]

import json
import sys
from classes import Protein


def main(mitab_file1=None):
	"""


    @param mitab_file1: filepath to biogrid download to be processed, passed from command line
                        defaults to hard-coded path in file for manual runs from an IDE
    @return:
    """

	# If run from PyCharm, filenames can be set here manually
	if mitab_file1 is None:
		# hold filepaths as strings
		mitab_file1 = "../../biogrid/BIOGRID-ORGANISM-Caenorhabditis_elegans-4.4.222.mitab.txt"

	species1 = extract_name(mitab_file1)

	# read ensembl file for worm
	with open("../../ensembl/" + species1 + "_ensembl-109.txt", 'r') as f:
		# initialize a list of proteins
		protein_list = []

		# skip first line
		f.readline()

		for line in f:
			# List of Protein objects
			line = line.split('\t')
			new_protein = Protein.Protein(line)
			protein_list.append(new_protein)

	# receive the physical and experimental codes
	good_codes = retrieve_code()

	# read the BioGRID file
	with open(mitab_file1, 'r') as f:
		# build a list of interactions
		inter_list1 = []

		# skip first line
		f.readline()

		for line in f:
			line = line.rstrip()
			line = line.split("\n")

			interaction = get_gene_ids(line, good_codes)

			if interaction:
				inter_list1.append(interaction)

	map_list1 = id_to_protein(protein_list, inter_list1)

	prots_dict1 = list_to_dict(map_list1)

	# query_subnetwork(prots_dict1, "R05F9.1d.1", prots_dict2, "ENSMUSP00000090649", protein_list, objList2, species1,
	#                 species2)

	print(map_list1)

	# with open(el_file, 'w') as e:
	#    e.write(ids[0][0] + '_' + ids[0][1] + '\t' + ids[1][0] + '_' + ids[1][1] + '\n')

	##### TESTING THE ENTIRE NETWORK #####
	# retList = dict_to_nodes(prots_dict, objList)
	#
	# retList2 = dict_to_edges(prots_dict, retList[1])
	#
	# retList[0].extend(retList2)
	#
	# # combines the protein and edge data for the entire network


def extract_name(filepath):
	"""
    Takes the BioGRID filepath and extracts the species name from it
    :param filepath: str BioGRID filepath
    :return: species name as str
    """
	species = filepath.split("/")
	species = species[-1].split("-")
	species = species[2]
	species = species[0] + species[species.find('_'):]
	species = species.lower()

	return species


def get_gene_ids(line, good_codes) -> ((str, [str]), (str, [str])):
	"""
    Extracts data from BioGRID and filters 'bad' interactions

    :param line: str line from BioGRID file
    :param good_codes: tuple of sets of physical and experimental codes
    :return: a tuple, representing a PPI, with 2 inner tuples, each with the EntrezGene ID and a list of
             the other listed BioGRID IDs for each interactor
    """

	# create a list contains strs after line splitted by tabs
	by_tab = line[0].split('\t')

	# set of physical codes
	phys_codes = good_codes[0]

	# set of experimental codes
	exp_codes = good_codes[1]

	# protein names in interaction
	name0 = by_tab[0].split(':')[1]
	name1 = by_tab[1].split(':')[1]

	# self-loop filter
	if name0 == name1:
		return 'self',

	taxid0 = by_tab[9].split(':')[1].strip()
	taxid1 = by_tab[10].split(':')[1].strip()

	# interspecies interaction filter (check taxids)
	if taxid0 != taxid1:
		return 'intersp',

	interaction_type = by_tab[11].split(':')[2][0:4]

	interaction_detection = by_tab[6].split(':')[2][0:4]

	# direct interaction only
	if interaction_type not in phys_codes:
		return 'phys',

	if interaction_detection not in exp_codes:
		return 'exp',

	gene_id = []
	# gene_id = [id0, id1]
	for i in range(0, 2):
		by_colon = by_tab[i].split(':')
		gene_id.append(by_colon[1])

	if gene_id[0] == gene_id[1]: # is this necessary?
	  return 'self',

	id_list0 = build_id_list(by_tab[2])
	id_list1 = build_id_list(by_tab[3])

	return (name0, id_list0), (name1, id_list1)


def build_id_list(line):
	"""
    Helper function for building a list of ids from the biogrid file
    under the alt ids for interactors A and B

    :param line: line of the BioGRID file that includes alternate ids
    :return: list of alternate ids in the order that they appear in the file
    """
	id = line.split("|")

	id_list = []
	for type in id:
		type = type.split(":")
		id_list.append(type[1])

	return id_list


def id_to_protein(protein_list: list[Protein], interactions: list[(str, list[str]), (str, list[str])]) -> list[(Protein, Protein)]:
	"""
    Takes a list of Protein objectss (with their ID aliases from Ensembl) and a list of protein interactions (with their
    ID alises from BioGRID), and merges the data into a list of protein interactions, stored as a list of tuples of 2
    Protein objects.

    :param protein_list: list of protein objects from ensembl
    :param interactions: list of PPIs, each a tuple of 2 interacting proteins, each a tuple with a string name and a
                         list of alternate ID strings from BioGRID
    :return: Protein object with corresponding ids as attributes
    """

	# count the number of interactions processed
	count = 0

	# count the amount of unmapped id's
	count_none = 0

	# dict of unmapped id's with name as key and count as value
	none_dict = {}

	# count the valid interactions that weren't filtered
	count_invalid = 0

	# other trackers
	non_phys = 0
	non_exp = 0
	self_loop = 0
	interspecies = 0

	# take a dictionary with the entrezgene as the key
	ncbi_lookup = {x.get_ncbi(): x for x in protein_list}

	# also take a dict with the name as the key
	# in case the entrezgene fails
	name_lookup = {x.get_name(): x for x in protein_list}

	# take a dict with the swissprot id as the key
	swissprot_lookup = {x.get_swissprot(): x for x in protein_list}
	trembl_lookup = {x.get_trembl(): x for x in protein_list}

	remapped_interactions = []
	unmapped_interactions = []
	for interaction in interactions:
		# keeps track of whether the interaction has been filtered or not
		failure = False


		remapped_interaction = []

		# count filtered out interactions
		if len(interaction) == 1:
			if interaction[0] == 'phys':
				non_phys += 1
			elif interaction[0] == 'exp':
				non_exp += 1
			elif interaction[0] == 'intersp':
				interspecies += 1
			elif interaction[0] == 'self':
				self_loop += 1
			continue

		# process "real" interactions
		for interactor in interaction:

			# Try to remap the BioGRID interactor to an Ensembl ID in current
			current = None

			# The BioGRID ID is stored as the first element of the interactor
			biogrid_id = interactor[0]

			if biogrid_id in ncbi_lookup:
				current = ncbi_lookup[biogrid_id]
			# If the current protein's NCBI ID from BioGRID does not match one in Ensembl, try matching its other IDs
			else:
				# check if any of the alternate IDs match any of the Ensembl IDs remaining
				for alt_id in interactor[1]:
					current = swissprot_lookup.get(alt_id, trembl_lookup.get(alt_id, name_lookup.get(alt_id, None)))

			# If the protein failed to be mapped, note it
			if current is None:
				current = "! " + str(biogrid_id)

				# increase count of current by 1 in tracking dictionary
				none_dict[current] = none_dict.get(current, 0) + 1
				count_none += 1

				failure = True

			remapped_interaction.append(current)

		# only add unfiltered interactions
		if failure:
			count_invalid += 1
			unmapped_interactions.append(remapped_interaction)
		else:
			count += 1
			remapped_interactions.append(remapped_interaction)



	## places output in a form that can be pasted into a spreadsheet ###
	with open('test_mapping.txt', 'w') as f:
		f.write("# Total interactions processed: " + str(count) + "\n")
		f.write("# Total filtered interactions: " + str(count_invalid) + "\n")
		f.write("# Percent filtered: " + str(round((count_invalid / count), 3) * 100) + "%\n")
		f.write("# Total unmappable id's in file (including duplicates): " + str(count_none) + "\n")
		f.write("# Percent unmappable: " + str(round((count_none / count), 3) * 100) + "%\n")
		f.write("# Non-physical interactions filtered out: " + str(non_phys) + "\n")
		f.write("# Non-experimental interactions filtered out: " + str(non_exp) + "\n")
		f.write("# Self-loops filtered out: " + str(self_loop) + "\n")
		f.write("# Interspecies interactions filtered out: " + str(interspecies) + "\n")
		f.write("\n")
		for id in none_dict:
			f.write(str(id) + " occurs " + str(none_dict[id]) + " times " + "\t" + str(
				round((none_dict[id] / count_none), 3) * 100) + "%\n")
		f.write("\n")
		for inter in remapped_interactions:
			for id in inter:
				f.write(str(id) + '\t')
			f.write("\n")

	return remapped_interactions


def query_subnetwork(p_dict_s1, prot_s1, p_dict_s2, prot_s2, objList_s1, objList_s2, s1, s2):
	"""
    Takes in a query protein and query species and returns a subnetwork in JSON format

    :param p_dict: dictionary of proteins and their connections as def list
    :param prot: desired protein as str
    :param species: desired species as str
    :return: list of dict values that corresponds with the desired subnetwork
    """
	prot_def1 = p_dict_s1[prot_s1]
	prot_def2 = p_dict_s2[prot_s2]

	# add the query protein to the list of its interacting proteins
	temp1 = prot_def1
	temp1.append(prot_s1)

	temp2 = prot_def2
	temp2.append(prot_s2)

	retList1 = list_to_nodes(temp1, objList_s1, 1)
	nodes1 = retList1[0]  # protein nodes to become JSON
	edge_dict1 = retList1[1]  # dict of protein name as key and corresponding id as def

	retList2 = list_to_nodes(temp2, objList_s2, 2)
	nodes2 = retList2[0]  # protein nodes to become JSON
	edge_dict2 = retList2[1]  # dict of protein name as key and corresponding id as def

	edges1 = list_to_edges(prot_s1, prot_def1, edge_dict1)
	edges2 = list_to_edges(prot_s2, prot_def2, edge_dict2)

	json_header = [
		{"data":
			 {"id": "species1", "name": s1},
		 "_comment": "Test output for a JSON file -- contains a subnetwork of C. elegans and a subnetwork of M. musculus",
		 "classes": "container s1"},
		{"data": {"id": "species2", "name": s2}, "classes": "container s2"},
		{"data": {"id": "aligned non-ortho", "name": "aligned non-orthology"}, "classes": "container"},
		{"data": {"id": "aligned ortho", "name": "aligned orthology"}, "classes": "container"}]

	# add all lists of dict objects together to form the whole JSON file
	json_header.extend(nodes1)
	json_header.extend(edges1)
	json_header.extend(nodes2)
	json_header.extend(edges2)

	json_str = json.dumps(json_header)

	file = open("test_json.json", "w")
	file.write(json_str)


def list_to_nodes(p_inters, objList, species_num):
	"""
    Converts python dictionary to JSON format with established conventions.
    Writes to a JSON file, returns nothing.

    :param p_dict: dictionary of proteins and their connections as def list
    :param objList: list of Protein objects to get alternate ids/names
    :return: final list of JSON elements without edges, dict of edge ids and e_ids
    """
	# take a dictionary with the entrezgene as the key
	d = {x.get_p_sid(): x for x in objList}

	# the list that will be ultimately converted to JSON
	# start with comment and protein category boxes
	final_list = []

	# dictionary that will hold node ids and edges
	edge_dict = {}

	count = 0
	for key in p_inters:
		json_dict = {"data": {}}
		json_dict["data"]["id"] = str(species_num) + "." + str(count)

		# if the id wasn't 'mappable'
		if key[0] == '!':
			# remove the ! from id
			json_dict["data"]["name"] = key[2:]
			json_dict["data"]["ncbi"] = key[2:]

		else:
			json_dict["data"]["e_id"] = key
			json_dict["data"]["name"] = d[key].get_name()
			json_dict["data"]["parent"] = "unaligned"
			json_dict["data"]["ncbi"] = d[key].get_ncbi()
			json_dict["data"]["uniprot"] = d[key].get_swissprot()  ## will need to alter this
			# there is also trembl

		json_dict["classes"] = "species" + str(species_num) + " unaligned protein"

		# add data to edge dictionary
		edge_dict[key] = json_dict["data"]["id"]

		final_list.append(json_dict)

		count += 1

	# the last entry in the dict will be the query
	json_dict["classes"] = json_dict["classes"] + " query"

	return [final_list, edge_dict]


def list_to_edges(prot, e_list, e_dict):
	"""
    Takes the converted nodes and writes edges to JSON based on established conventions
    and data in the established python dictionary

    :param p_dict: dictionary of each protein in the network and their connections
    :param edge_dict: dictionary of id and its corresponding ensembl id
    :return: list of JSON edges
    """

	# final list of JSON edges to be returned
	final_list = []

	temp = {}

	# loop through p_dict
	for ele in e_list:
		temp = {"data": {}}

		# set the source to the current protein
		temp["data"]["source"] = e_dict[prot]

		# set the target to the interacting protein
		temp["data"]["target"] = e_dict[ele]

		temp["classes"] = "species" + str(e_dict[prot][0]) + " edge"

		# print(temp)

		# add to the final list
		final_list.append(temp)

	return final_list


def list_to_dict(map_list):
	"""
    Constructs a dictionary with a list of each protein it interacts with
    as def

    :param map_list: list of mapped protein ids
    :return: dict of each protein id and it's interactors
    """

	prot_hash = {}
	for interaction in map_list:
		current = 0
		for prot in interaction:
			# check if the protein already exists in the dictionary
			if prot not in prot_hash:
				prot_hash[prot] = []  # create a new list if not in dict

			prot2 = interaction[current - 1]  # receive interacting protein (either 0 or -1)

			# check if the interacting protein is already in the current def list
			# and also not in the dict already so that there aren't duplicate interactions
			if prot2 not in prot_hash[prot]:
				# add the interacting protein to the list of interactors for the current
				prot_hash[prot].append(interaction[current - 1])

			current += 1

	return prot_hash


def retrieve_code() -> (set[str], set[str]):
	"""
    Extracts data from physical and experimental interaction files
    and places them into two respective sets.

    :return: 2 sets of strings, with MI codes for physical interactions and experimentally-detected interaction
    """
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
	return physical_set, experimental_set








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

if __name__ == '__main__':
	from argparse import ArgumentParser

	ap = ArgumentParser()

	args, unknown_args = ap.parse_known_args()
	print(args)
	print(unknown_args)

	if len(unknown_args) == 1:
		main(unknown_args[0])
	else:
		main()

	# print(retrieve_code())
	# print([len(retrieve_code()[0]), len(retrieve_code()[1])]) what is the reason??
