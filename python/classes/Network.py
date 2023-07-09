__author__ = "Brian Law"
__email__ = "blaw@iwu.edu"

import os

from python.classes.Species import Species, species_dict, species_list
from python.classes import Protein
import python.utility as utility

class Network:

  def __init__(self, species_name=None, ensembl_version=None, biogrid_version=None):
    self.forward = {}
    self.backward = {}
    self.nodes = []

    if species_name is None:
      return

    species = utility.get_species(species_name)

    protein_dict = utility.get_ensembl_data(species.name, ensembl_version)

    if biogrid_version is None:
      filenames = os.listdir('../../biogrid/')
      file_versions = [filename[filename.rfind('-') + 1:filename.rfind('.mitab')] for filename in filenames if
                       os.path.isfile(os.path.join('../../biogrid', filename))]
      file_versions = [tuple([int(x) for x in version.split('.')]) for version in file_versions]

      biogrid_version = '.'.join(map(str, max(file_versions, key=lambda x: (x[0], x[1], x[2]))))

    ensembl_filepath_others = f'{utility.get_project_root()}/ensembl/{species.short_name.replace(" ", "_").lower()}_ensembl_others-{ensembl_version}.txt'
    ensembl_filepath_ncbi = f'{utility.get_project_root()}/ensembl/{species.short_name.replace(" ", "_").lower()}_ensembl_ncbi-{ensembl_version}.txt'
    network_filepath = f'{utility.get_project_root()}/networks/{species.short_name.lower().replace(" ", "_")}-{ensembl_version}-{biogrid_version}.txt'

    with open(network_filepath) as f:
      # build a list of interactions
      interaction_list = []

      # skip first line
      f.readline()

      for line in f:
        line = line.rstrip()
        line = line.split("\n")

        interaction = get_gene_ids(line, good_codes)

        if interaction:
          interaction_list.append(interaction)


def get_gene_ids(line, good_codes) -> tuple[tuple[str, dict[str: list[str]]], tuple[str, dict[str: list[str]]]]:
  """
  Processes 1 line of the BioGRID file, filtering out filters "bad" interactions, and returning a representation of
  the remaining "good" interactions.

  :param line: str line from BioGRID file
  :param good_codes: tuple of sets of physical and experimental codes
  :return: - if "bad," a single string tuple with a string code indicating why it was bad ("self"-loops, "intersp"-ecies
           interactions, non-"phys"-ical interactions, non-"exp"-erimentally detected interactions)
           - if "good", a tuple, representing a PPI, with 2 inner tuples, each with the EntrezGene ID and a dictionary
           of the alt IDs for that protein, read from the BioGRID file, as per the helper function build_id_list(...)
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

  # interspecies interaction filter (check taxids)
  taxid0 = by_tab[9].split(':')[1].strip()
  taxid1 = by_tab[10].split(':')[1].strip()
  if taxid0 != taxid1:
    return 'intersp',

  # direct PPI physical interactions only
  interaction_type = by_tab[11].split(':')[2][0:4]
  if interaction_type not in phys_codes:
    return 'phys',

  # experimentally detected interactions only
  interaction_detection = by_tab[6].split(':')[2][0:4]
  if interaction_detection not in exp_codes:
    return 'exp',

  # interaction is "good" - passed all filters - and ready to be processed
  gene_id = []
  # gene_id = [id0, id1]
  for i in range(0, 2):
    by_colon = by_tab[i].split(':')
    gene_id.append(by_colon[1])

  if gene_id[0] == gene_id[1]:  # is this necessary?
    return 'self',

  id_dict1 = build_id_list(by_tab[2])
  id_dict2 = build_id_list(by_tab[3])

  return (name0, id_dict1), (name1, id_dict2)


def build_id_list(line: str) -> dict[str, list[str]]:
  """
  Helper function for building a list of ids from the biogrid file under the alt ids for interactors A and B

  :param line: line of the BioGRID file that includes alternate ids
  :return: a dictionary mapping ID types, as entered in the BioGRID Alt IDs column, to a list of the alternate IDs for
           that type. E.g. {'refseq':['NP_571578', 'NP_439203'], 'uniprot/swiss-prot':['P57094']}
  """
  id_text_split = line.split('|')
  id_dict = {}

  for id_type in id_text_split:
    id_type = id_type.split(':')
    id_dict.setdefault(id_type[0], list())
    id_dict[id_type[0]].append(id_type[1])

  return id_dict
