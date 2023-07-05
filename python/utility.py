__author__ = "Brian Law"
__email__ = "blaw@iwu.edu"

import os

from python.classes.Species import Species, species_dict, species_list
from python.classes import Protein
from python import utility
from pathlib import Path

g_ensembl_data = {}

def get_project_root() -> Path:
  return Path(__file__).parent.parent


def get_species(species_name: str) -> Species:
  try:
    species = species_dict[species_name]
  except KeyError as ke:
    try:
      species_name = species.name.capitalize()
      species = species_dict[species_name]
    except KeyError as ke:
      species_name = species.name.replace('.', '').replace('_', '')
      species = species_dict[species_name]

  return species


def get_ensembl_data(species_name: str, ensembl_version: str=None) -> dict[str, Protein.Protein]:
  species = get_species(species_name)

  if ensembl_version is None:
    filenames = os.listdir(f'{get_project_root()}/ensembl/')
    file_versions = [int(filename[filename.rfind('-') + 1:filename.rfind('.')]) for filename in filenames if
                     os.path.isfile(os.path.join('../../ensembl', filename))]

    ensembl_version = max(file_versions)

  if species.name not in g_ensembl_data:
    g_ensembl_data[species.name] = read_ensembl_data(species_name, ensembl_version)

  return g_ensembl_data[species.name]


def read_ensembl_data(species_name: str, ensembl_version: str) -> dict[str, Protein.Protein]:
  """
  Utility function to read in Ensembl data for a specific species. Will generate a list (dictionary) of Protein objects
  and their IDs from downloaded Ensembl files.

  Accesses Ensembl files in the ensembl subdirectory.

  :param species_name: the name of the species to load Ensembl data from. Can be colloquial (human) or scientific (homo sapiens).
  :param ensembl_version:
  :return:
  """

  # Fetch species based on provided species name
  species = get_species(species_name)

  ensembl_filepath_others = f'{utility.get_project_root()}/ensembl/{species.short_name.replace(" ", "_").lower()}_ensembl_others-{ensembl_version}.txt'
  ensembl_filepath_ncbi = f'{utility.get_project_root()}/ensembl/{species.short_name.replace(" ", "_").lower()}_ensembl_ncbi-{ensembl_version}.txt'

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
      protein.add_ncbi_ids(line)

  return protein_dict


if __name__ == '__main__':
  print(f'{get_project_root()}/ensembl/')
