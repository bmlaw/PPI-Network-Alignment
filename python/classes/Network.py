__author__ = "Brian Law"
__email__ = "blaw@iwu.edu"

import os
import datetime
import networkx

from python.classes.Species import Species, species_dict, species_list
from python.classes import Protein
import python.network.interaction_codes as interaction_codes
import python.utility as utility


class Network:

  def __init__(self, species_name=None, ensembl_version=None, biogrid_version=None):
    self.graph = networkx.Graph()

    if species_name is None:
      return

    species = utility.get_species(species_name)

    # If no Ensembl version is specified, search the directory for the latest version and use that.
    if ensembl_version is None:
      ensembl_species_name = species.short_name.lower().replace(' ', '_')
      filenames = os.listdir(f'{utility.get_project_root()}/ensembl/')
      filenames = [f for f in filenames if ensembl_species_name in f]
      file_versions = [int(filename[filename.rfind('-') + 1:filename.rfind('.')]) for filename in filenames if
                      os.path.isfile(os.path.join(f'{utility.get_project_root()}/ensembl', filename))]

      ensembl_version = max(file_versions)

    # initialize a dictionary of Proteins, with key as their Ensembl gene IDs for easy lookup
    protein_dict = utility.get_ensembl_data(species.name, ensembl_version)

    if biogrid_version is None:
      filenames = os.listdir(f'{utility.get_project_root()}/biogrid/')
      file_versions = [filename[filename.rfind('-') + 1:filename.rfind('.mitab')] for filename in filenames if
                       os.path.isfile(os.path.join(f'{utility.get_project_root()}/biogrid', filename))]
      file_versions = [tuple([int(x) for x in version.split('.')]) for version in file_versions]

      biogrid_version = '.'.join(map(str, max(file_versions, key=lambda x: (x[0], x[1], x[2]))))
    network_filepath = f'{utility.get_project_root()}/networks/{species.short_name.replace(" ", "_").lower()}-{ensembl_version}-{biogrid_version}.txt'

    # get the physical and experimental MI codes
    good_codes = (interaction_codes.get_experimental_codes(), interaction_codes.get_physical_codes)

    s = f'{utility.get_project_root()}/networks/{species.short_name.replace(" ", "_").lower()}.network-{ensembl_version}-{biogrid_version}.txt'
    print(s)
    f = open(s)
    print(f.readline())
    f.close()

    with open(network_filepath, encoding='UTF-8') as f:
      # build a list of interactions
      interaction_list = []

      # skip first line
      f.readline()

      for line in f:
        line = line.rstrip()
        line = line.split("\n")

        interaction = get_gene_ids(line, good_codes)

        interaction_list.append(interaction)

    interaction_list = id_to_protein(protein_dict, interaction_list, species, ensembl_version, biogrid_version)

    self.graph.add_edge(interaction_list[0], interaction_list[1])  

    print(self.graph.nodes)


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
    return ('self',)

  # interspecies interaction filter (check taxids)
  taxid0 = by_tab[9].split(':')[1].strip()
  taxid1 = by_tab[10].split(':')[1].strip()
  if taxid0 != taxid1:
    return ('intersp',)

  # direct PPI physical interactions only
  interaction_type = by_tab[11].split('"')[1]
  if interaction_type not in phys_codes:
    return ('phys',)

  # experimentally detected interactions only
  interaction_detection = by_tab[6].split('"')[1]
  if interaction_detection not in exp_codes:
    return ('exp',)

  # interaction is "good" - passed all filters - and ready to be processed
  gene_id = []
  # gene_id = [id0, id1]
  for i in range(0, 2):
    by_colon = by_tab[i].split(':')
    gene_id.append(by_colon[1])

  if gene_id[0] == gene_id[1]:  # is this necessary?
    return ('self',)

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
    id_dict.setdefault(id_type[0], [])
    id_dict[id_type[0]].append(id_type[1])

  return id_dict


def id_to_protein(protein_dict: dict[str, Protein.Protein],
                  interactions: list[tuple[tuple[str, dict[str: list[str]]], tuple[str, dict[str: list[str]]]]],
                  species: Species,
                  ensembl_version: str,
                  biogrid_version: str) -> list[(Protein.Protein, Protein.Protein)]:
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

  ensembl_lookup = {x: protein_dict[x] for x in protein_dict} | \
                   {y: protein_dict[x] for x in protein_dict for y in protein_dict[x].t_ids} | \
                   {y: protein_dict[x] for x in protein_dict for y in protein_dict[x].p_ids}

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
  ncbi_count = 0
  ensembl_count = 0
  swissprot_count = 0
  trembl_count = 0
  refseq_count = 0
  name_count = 0

  for interaction in interactions:
    # Keeps track of whether the interaction has been filtered or not. Each interactor may flip this to True if
    # mapping fails.
    failure = False

    # After mapping to proteins, store new interaction in this list.
    remapped_interaction = []

    # Count filtered out interactions.
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

    # Process "real" interactions by trying to map each interactor
    for interactor in interaction:

      # Try to remap the BioGRID interactor to an Ensembl ID in current
      current = None

      # The BioGRID ID is stored as the first element of the interactor
      biogrid_id = interactor[0]

      # First look-up is BioGRID's NCBI IDs
      if biogrid_id in ncbi_lookup:
        current = ncbi_lookup[biogrid_id]
        ncbi_count += 1

      # If the current protein's NCBI ID from BioGRID does not match one in Ensembl, try matching its other IDs
      if current is None:
        for alt_id in interactor[1]['entrez gene/locuslink']:
          current = ensembl_lookup.get(alt_id)
          if current is not None:
            ensembl_count += 1
            break

      # if current is None:
      #   for alt_id in interactor[1]['entrez gene/locuslink']:
      #     current = ensembl_lookup.get(alt_id.replace('CELE_', ''))

      if current is None and 'uniprot/swiss-prot' in interactor[1]:
        for alt_id in interactor[1]['uniprot/swiss-prot']:
          current = swissprot_lookup.get(alt_id)
          if current is not None:
            swissprot_count += 1
            break

      if current is None and 'uniprot/trembl' in interactor[1]:
        for alt_id in interactor[1]['uniprot/trembl']:
          current = trembl_lookup.get(alt_id)
          if current is not None:
            trembl_count += 1
            break

      if current is None and 'refseq' in interactor[1]:
        for alt_id in interactor[1]['refseq']:
          current = refseq_lookup.get(alt_id)
          if current is not None:
            refseq_count += 1
            break

      if current is None:
        for alt_id in interactor[1]['entrez gene/locuslink']:
          current = name_lookup.get(alt_id)
          if current is not None:
            name_count += 1
            break

      # If the protein failed to be mapped, note it
      if current is None:
        # If no protein exists, just use a string mapping instead.
        current = "! " + str(biogrid_id)

        # Increase count of current by 1 in tracking dictionary
        none_dict[current] = none_dict.get(current, 0) + 1
        count_none += 1

        # Flag as an invalid interaction
        failure = True

      # If the gene is a non-coding gene, it shouldn't be included in the final network. This is a possible problem
      # with BioGRID data; we are assuming Ensembl is correctly identifying coding/non-coding genes using a protein
      # id (or the absence of one)
      elif len(current.p_ids) == 0:
        noncoding_dict[current] = noncoding_dict.get(current, 0) + 1
        count_noncoding += 1

        # Flag as an invalid interaction
        failure = True

      # Add the current interactor, whether successful or failed, to the remapped interaction.
      remapped_interaction.append(current)

    # Only add interactions where both interactors were successfully mapped to coding genes
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
  # print('Lookup counts: ', ncbi_count, ensembl_count, swissprot_count, trembl_count, refseq_count, name_count, count_none)

  # Sort all interactions alphabetically for output
  remapped_interactions = list(remapped_interactions)
  remapped_interactions.sort(key=lambda x: (x[0].gene_id.upper(), x[1].gene_id.upper()))

  # Count number of interactions processed (different from # interactions in final network) as some may be duplicates
  total_processed = count + count_invalid + non_phys + non_exp + self_loop + interspecies

  ## places output in a form that can be pasted into a spreadsheet ###
  with open(f'{utility.get_project_root()}/networks/{species.short_name.replace(" ", "_").lower()}.network-{ensembl_version}-{biogrid_version}.txt', 'w', encoding='UTF-8') as f1, \
       open(f'{utility.get_project_root()}/networks/SANA/{species.short_name.replace(" ", "_").lower()}.network-{ensembl_version}-{biogrid_version}.el', 'w', encoding='UTF-8') as f2:

    f1.write(f'! Generated using Ensembl {ensembl_version} and BioGRID {biogrid_version} on {datetime.datetime.now()}\n')
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
    for gene_id in sorted(none_dict, key=lambda x: (none_dict[x], x), reverse=True):
      f1.write(str(gene_id) + " occurs " + str(none_dict[gene_id]) + " times " + "\t" + str(round((none_dict[gene_id] / total_processed) * 100, 3)) + "%\n")
    f1.write("!\n")

    # Write each interaction to file
    for interaction in remapped_interactions:
      f1.write(interaction[0].gene_id + '\t' + interaction[1].gene_id + '\n')
      f2.write(interaction[0].gene_id + '\t' + interaction[1].gene_id + '\n')

  return remapped_interactions

if __name__ == '__main__':
  n = Network('zebrafish')