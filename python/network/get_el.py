"""
get_el.py extracts data from physical_interaction_codes.txt,
experimental_detected_codes.txt, and ensembl gene names.
"""
__author__ = "Anna Sheaffer"
__email__ = "asheaffe@iwu.edu"
__credits__ = ["Norman Luo", "Brian Law"]


import os
import datetime

from python.classes import Species
from python.classes import Protein
from python.classes import Network
import python.utility as utility
import python.network.interaction_codes as interaction_codes


def main(species_name=None, ensembl_version=None, biogrid_version=None):
  """

  :param species_name: filepath to biogrid download to be processed, passed from command line
                       defaults to hard-coded path in file for manual runs from an IDE
  :param ensembl_version: version of Ensembl to use; if not specified, function will search the biogrid subfolder for
                          the latest version
  :param biogrid_version: version of Ensembl to use; if not specified, function will search the biogrid subfolder for
                          the latest version
  :return:
  """

  # If run from PyCharm, filenames can be set here manually
  if species_name is None:
    # hold filepaths as strings
    species_name = "yeast"

  species = utility.get_species(species_name)

  if ensembl_version is None:
    ensembl_species_name = species.short_name.lower().replace(' ', '_')
    filenames = os.listdir(f'{utility.get_project_root()}/ensembl/')
    filenames = [f for f in filenames if ensembl_species_name in f]
    file_versions = [int(filename[filename.rfind('-') + 1:filename.rfind('.')]) for filename in filenames if
                     os.path.isfile(os.path.join(f'{utility.get_project_root()}/ensembl', filename))]

    ensembl_version = max(file_versions)

  if biogrid_version is None:
    filenames = os.listdir(f'{utility.get_project_root()}/biogrid/')
    file_versions = [filename[filename.rfind('-') + 1:filename.rfind('.mitab')] for filename in filenames if
                     os.path.isfile(os.path.join(f'{utility.get_project_root()}/biogrid', filename))]
    file_versions = [tuple([int(x) for x in version.split('.')]) for version in file_versions]

    biogrid_version = '.'.join(map(str, max(file_versions, key=lambda x: (x[0], x[1], x[2]))))

  biomart_filepath = f'{os.path.join(utility.get_project_root())}/biogrid/BIOGRID-ORGANISM-{species.long_name.replace(" ", "_")}-{biogrid_version}.mitab.txt'
  print(biomart_filepath)

  # initialize a dictionary of Proteins, with key as their Ensembl gene IDs for easy lookup
  protein_dict = utility.get_ensembl_data(species_name, ensembl_version)

  # get the physical and experimental MI codes
  good_codes = (interaction_codes.get_experimental_codes(), interaction_codes.get_physical_codes)

  # read the BioGRID file
  with open(biomart_filepath, 'r', encoding='UTF-8') as f:
    # build a list of interactions
    inter_list1 = []

    # skip first line
    f.readline()

    for line in f:
      line = line.rstrip()
      line = line.split("\n")

      interaction = Network.get_gene_ids(line, good_codes)

      if interaction:
        inter_list1.append(interaction)


  # print(map_list1)


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
  # print(args)
  # print(unknown_args)

  if len(unknown_args) == 1:
    main(unknown_args[0])
  else:
    # for species in Species.species_list:
    #   print(species)
    #   main(species.short_name)
    main('zebrafish')


# print(retrieve_code())
# print([len(retrieve_code()[0]), len(retrieve_code()[1])]) what is the reason??
