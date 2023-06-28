"""
get_el.py extracts data from physical_interaction_codes.txt,
experimental_detected_codes.txt, and ensembl gene names.
"""
__author__ = "Anna Sheaffer"
__email__ = "asheaffe@iwu.edu"
__credits__ = ["Norman Luo", "Brian Law"]

import json
import os

from python.classes.Species import Species, species_dict, species_list
from python.classes import Protein


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

  species = species_dict[species_name]

  if ensembl_version is None:
    filenames = os.listdir('../../ensembl/')
    file_versions = [int(filename[filename.rfind('-') + 1:filename.rfind('.')]) for filename in filenames if
                     os.path.isfile(os.path.join('../../ensembl', filename))]

    ensembl_version = max(file_versions)

  if biogrid_version is None:
    filenames = os.listdir('../../biogrid/')
    file_versions = [filename[filename.rfind('-') + 1:filename.rfind('.mitab')] for filename in filenames if
                     os.path.isfile(os.path.join('../../biogrid', filename))]
    file_versions = [tuple([int(x) for x in version.split('.')]) for version in file_versions]

    biogrid_version = '.'.join(map(str, max(file_versions, key=lambda x: (x[0], x[1], x[2]))))

  ensembl_filepath_others = f'../../ensembl/{species.short_name.replace(" ", "_").lower()}_ensembl_others-{ensembl_version}.txt'
  ensembl_filepath_ncbi = f'../../ensembl/{species.short_name.replace(" ", "_").lower()}_ensembl_ncbi-{ensembl_version}.txt'
  biomart_filepath = f'../../biogrid/BIOGRID-ORGANISM-{species.long_name.replace(" ", "_")}-{biogrid_version}.mitab.txt'

  # initialize a dictionary of Proteins, with key as their Ensembl gene IDs for easy lookup
  protein_dict = {}

  # read ensembl file with swissprot, trembl and refseq
  with open(ensembl_filepath_others, 'r') as f:
    f.readline()  # skips the first line
    for line in f:

      # List of Protein objects
      line = line.strip().split('\t')
      gene_id = line[0]
      if gene_id not in protein_dict:
        protein = Protein.Protein(gene_id)
        protein_dict[gene_id] = protein
      else:
        protein = protein_dict[gene_id]
      protein.add_other_ids(line)

  # read ensembl file with ncbi ids
  with open(ensembl_filepath_ncbi, 'r') as f:
    f.readline() # skips the first line
    for line in f:

      # List of Protein objects
      line = line.strip().split('\t')
      gene_id = line[0]
      if gene_id not in protein_dict:
        protein = Protein.Protein(gene_id)
        protein_dict[gene_id] = protein
      else:
        protein = protein_dict[gene_id]
      protein.add_ncbi_ids(line)

  # receive the physical and experimental codes
  good_codes = retrieve_code()

  # read the BioGRID file
  with open(biomart_filepath, 'r') as f:
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

  map_list1 = id_to_protein(protein_dict, inter_list1, species, ensembl_version, biogrid_version)

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

  if gene_id[0] == gene_id[1]:  # is this necessary?
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
  ids = line.split("|")

  id_list = []
  for id_type in ids:
    id_type = id_type.split(":")
    id_list.append(id_type[1])

  return id_list


def id_to_protein(protein_dict: dict[str, Protein], interactions: list[(str, list[str]), (str, list[str])],
                  species: Species, ensembl_version: str, biogrid_version: str) -> list[(Protein, Protein)]:
  """
  Takes Proteins (with their ID aliases from Ensembl) and a list of protein interactions (with their
  ID aliases from BioGRID), and merges the data into a list of protein interactions, stored as a list of tuples of 2
  Protein objects.

  :param protein_dict: dictionary mapping Ensembl Gene IDs to Protein objects built from Ensembl data
  :param interactions: list of PPIs, each a tuple of 2 interacting proteins, each a tuple with a string name and a
                       list of alternate ID strings from BioGRID
  :param ensembl_version:
  :param biogrid_version:
  :return:
  """
  # containers for final interactions
  remapped_interactions = set()
  unmapped_interactions = set()

  # count the number of interactions processed
  count = 0

  # count the number of unmapped id's
  count_none = 0
  none_dict = {}

  # count the valid interactions that weren't filtered
  count_invalid = 0

  # count the # interations mapped to noncoding genes
  count_noncoding = 0
  noncoding_dict = {}

  # other trackers
  non_phys = 0
  non_exp = 0
  self_loop = 0
  interspecies = 0

  # take a dictionary with the entrezgene as the key
  ncbi_lookup = {y: protein_dict[x] for x in protein_dict for y in protein_dict[x].ncbis}

  # also take a dict with the name as the key
  # in case the entrezgene fails
  name_lookup = {y: protein_dict[x] for x in protein_dict for y in protein_dict[x].names}

  # take a dict with the swissprot id as the key
  swissprot_lookup = {y: protein_dict[x] for x in protein_dict for y in protein_dict[x].swissprots}
  trembl_lookup = {y: protein_dict[x] for x in protein_dict for y in protein_dict[x].trembls}
  refseq_lookup = {y: protein_dict[x] for x in protein_dict for y in protein_dict[x].refseqs}

  # Debugging counters for checking protein mapping source efficacy
  # ncbi_count = 0
  # ensembl_count = 0
  # swissprot_count = 0
  # trembl_count = 0
  # refseq_count = 0
  # name_count = 0

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
        # ncbi_count += 1
      # If the current protein's NCBI ID from BioGRID does not match one in Ensembl, try matching its other IDs
      else:
        # check if any of the alternate IDs match any of the Ensembl IDs remaining
        for alt_id in interactor[1]:
          current = protein_dict.get(alt_id, swissprot_lookup.get(alt_id, trembl_lookup.get(alt_id,
                    refseq_lookup.get(alt_id, name_lookup.get(alt_id, None)))))

          # Debugging code for checking protein mapping source efficacy
          # if alt_id in ensembl_lookup:
          #   current = ensembl_lookup[alt_id]
          #   ensembl_count += 1
          #   break
          # elif alt_id in swissprot_lookup:
          #   current = swissprot_lookup
          #   swissprot_count += 1
          #   break
          # elif alt_id in trembl_lookup:
          #   current = swissprot_lookup
          #   trembl_count += 1
          #   break
          # elif alt_id in refseq_lookup:
          #   current = swissprot_lookup
          #   refseq_count += 1
          #   break
          # elif alt_id in name_lookup:
          #   current = swissprot_lookup
          #   name_count += 1
          #   break

      # If the protein failed to be mapped, note it
      if current is None:
        current = "! " + str(biogrid_id)

        # increase count of current by 1 in tracking dictionary
        none_dict[current] = none_dict.get(current, 0) + 1
        count_none += 1

        failure = True
      elif len(current.p_ids) == 0:
        noncoding_dict[current] = noncoding_dict.get(current, 0) + 1
        count_noncoding += 1
        failure = True

      remapped_interaction.append(current)

    # only add unfiltered interactions
    if failure:
      count_invalid += 1
      unmapped_interactions.add(tuple(remapped_interaction))
    else:
      count += 1
      # If protein #2 is alphabetically before protein #1, swap the two so they're alphabetically ordered.
      if remapped_interaction[0].gene_id < remapped_interaction[1].gene_id:
        remapped_interaction[0], remapped_interaction[1] = remapped_interaction[1], remapped_interaction[0]

      remapped_interactions.add(tuple(remapped_interaction))

  # Debugging code for checking protein mapping source efficacy
  # print(ncbi_count, ensembl_count, swissprot_count, trembl_count, refseq_count, name_count)

  remapped_interactions = list(remapped_interactions)
  remapped_interactions.sort(key=lambda x: (x[0].gene_id.upper(), x[1].gene_id.upper()))

  # for i in range(len(remapped_interactions) - 1):
  #   if remapped_interactions[i][0].get_ncbi().upper() == remapped_interactions[i+1][0].get_ncbi().upper() and \
  #      remapped_interactions[i][1].get_ncbi().upper() == remapped_interactions[i + 1][1].get_ncbi().upper() and \
  #      remapped_interactions[i][0] != remapped_interactions[i+1][0] and remapped_interactions[i][1] != remapped_interactions[i+1][1]:
  #     print(remapped_interactions[i])
  #     print(remapped_interactions[i+1])
  #     print('')

  total_processed = count + count_invalid + non_phys + non_exp + self_loop + interspecies

  ## places output in a form that can be pasted into a spreadsheet ###
  with open(f'../../networks/{species.short_name.replace(" ", "_").lower()}.network-{ensembl_version}-{biogrid_version}.txt', 'w') as f1, \
    open(f'../../networks/SANA/{species.short_name.replace(" ", "_").lower()}.network-{ensembl_version}-{biogrid_version}.el', 'w') as f2:

    f1.write(f'! Generated using Ensembl {ensembl_version} and BioGRID {biogrid_version}\n')
    f1.write("! Total interactions processed: " + str(total_processed) + "\n")
    f1.write(f'! {count} total interactions mapped successfully into {len(remapped_interactions)} unique interactions\n')
    f1.write("! non-physical interactions filtered out: " + str(non_phys) + "\n")
    f1.write("! non-experimental interactions filtered out: " + str(non_exp) + "\n")
    f1.write("! self-loops filtered out: " + str(self_loop) + "\n")
    f1.write("! interspecies interactions filtered out: " + str(interspecies) + "\n")
    f1.write(f'! interactions mapped to invalid "proteins": {count_invalid}\n')
    f1.write(f'! % interactions mapped to invalid "proteins": {round((count_invalid / total_processed) * 100, 3)}%\n')
    f1.write("! total interactions with unmappable IDs (including duplicates): " + str(count_none) + "\n")
    f1.write("! % IDs unmappable (including duplicates): " + str(round((count_none / (count * 2 + count_none)) * 100, 3)) + "%\n")
    f1.write("! total unmappable proteins (excluding duplicates): " + str(len(none_dict)) + "\n")
    f1.write("! total interactions with noncoding genes (including duplicates): " + str(count_noncoding) + "\n")
    f1.write(f'! % IDs mapped to noncoding genes (including duplicates): {round((count_noncoding / (count * 2 + count_noncoding)) * 100, 3)}%\n')
    f1.write("! total mapped noncoding proteins (excluding duplicates): " + str(len(noncoding_dict)) + "\n")
    f1.write("!\n")
    f1.write('! Non-coding genes:\n')
    for protein in sorted(noncoding_dict, key=lambda x: x.gene_id.upper()):
      f1.write(f'! {protein.gene_id} occurs {noncoding_dict[protein]} times\t{round((noncoding_dict[protein] / total_processed) * 100, 3)}%\n')
    f1.write("!\n")
    f1.write('! Non-existent proteins:\n')
    for gene_id in sorted(none_dict, key=lambda x: none_dict[x], reverse=True):
      f1.write(str(gene_id) + " occurs " + str(none_dict[gene_id]) + " times " + "\t" + str(
        round((none_dict[gene_id] / total_processed) * 100, 3)) + "%\n")
    f1.write("!\n")

    for interaction in remapped_interactions:
      f1.write(interaction[0].gene_id + '\t' + interaction[1].gene_id + '\n')
      f2.write(interaction[0].gene_id + '\t' + interaction[1].gene_id + '\n')


  return remapped_interactions


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
  # print(args)
  # print(unknown_args)

  if len(unknown_args) == 1:
    main(unknown_args[0])
  else:
    for species in species_list:
      print(species)
      main(species.short_name)
    #main()


# print(retrieve_code())
# print([len(retrieve_code()[0]), len(retrieve_code()[1])]) what is the reason??