__author__ = "Brian Law"
__email__ = "blaw@iwu.edu"

import os

from pathlib import Path
from python.classes.Species import Species, species_dict, species_list
from python.classes import Protein

g_ensembl_data = {}


def get_project_root() -> Path:
  """ Function to return a filepath to the project root. This function exists to simplify issues with different IDEs and
  / or command-line usage of these Python files, which may be stored in different subdirectory (levels).

  :return: A filepath to the project's root directory.
  """
  return Path(__file__).parent.parent


def get_species(species_name: str) -> Species:
  """ Global getter method to get the Species object for a species name, in colloquial or scientific format. Mostly for
  convenient translation from one name-type to another.

  :param species_name: the name of the species
  :return: the global Species object for that species, should be treated as immutable
  """
  try:
    species = species_dict[species_name]
  except KeyError:
    try:
      species_name = species_name.capitalize()
      species = species_dict[species_name]
    except KeyError:
      species_name = species_name.replace('.', '').replace('_', '')
      species = species_dict[species_name]

  return species


def get_ensembl_data(species_name: str, ensembl_version: str=None) -> dict[str, Protein.Protein]:
  """ Gets the Ensembl data for a species from the downloaded Ensembl BioMart files, loaded into a dictionary for easy
      lookup. Accessor for a global variable / singleton pattern.

      Warning: the version number is ONLY for reading the data off disk. If another version of Ensembl data for this 
      species has already been loaded off the disk, Ensembl data for that version will be returned instead of the 
      specified version.

  :param species_name: the name of the species, can be in colloquial (human), short (H sapiens), or long form (Homo 
                       sapiens)
  :param ensembl_version: the specific version of Ensembl to load up, defaults to the latest version that can be found
                          in the directory if not provided
  :return: dictionary mapping Ensembl gene IDs to Protein objects
  """
  
  # Transform the provides species name into the global Species object for translational purposes
  species = get_species(species_name)

  # If species data is already loaded into the global, return it. Note this is UNVERSIONED.
  if species.name in g_ensembl_data:
    return g_ensembl_data[species.name]

  # If no Ensembl version is specified, search the directory for the latest version and use that.
  if ensembl_version is None:
    ensembl_species_name = species.short_name.lower().replace(' ', '_')
    filenames = os.listdir(f'{get_project_root()}/ensembl/')
    filenames = [f for f in filenames if ensembl_species_name in f]
    file_versions = [int(filename[filename.rfind('-') + 1:filename.rfind('.')]) for filename in filenames if
                     os.path.isfile(os.path.join(f'{get_project_root()}/ensembl', filename))]

    ensembl_version = max(file_versions)

  # Read data from disk, store it as global, then return it.
  if species.name not in g_ensembl_data:
    g_ensembl_data[species.name] = read_ensembl_data(species_name, ensembl_version)

  return g_ensembl_data[species.name]


def read_ensembl_data(species_name: str, ensembl_version: str) -> dict[str, Protein.Protein]:
  """ Utility function to read in Ensembl data for a specific species. Will generate a list (dictionary) of Protein 
  objects and their IDs from downloaded Ensembl files, saves it into a global variable, and returns the data for
  immediate access.

  Accesses Ensembl files in the ensembl subdirectory.

  :param species_name: the name of the species, can be in colloquial (human), short (H sapiens), or long form (Homo 
                       sapiens)
  :param ensembl_version: the specific version of Ensembl to load up, defaults to the latest version that can be found
                          in the directory if not provided
  :return: dictionary mapping Ensembl gene IDs to Protein objects
  """

  # Fetch species based on provided species name
  species = get_species(species_name)

  ensembl_filepath_others = f'{get_project_root()}/ensembl/{species.short_name.replace(" ", "_").lower()}_ensembl_others-{ensembl_version}.txt'
  ensembl_filepath_ncbi = f'{get_project_root()}/ensembl/{species.short_name.replace(" ", "_").lower()}_ensembl_ncbi-{ensembl_version}.txt'

  # initialize a dictionary of Proteins, with key as their Ensembl gene IDs for easy lookup
  protein_dict = {}

  # read ensembl file with swissprot, trembl and refseq
  with open(ensembl_filepath_others, 'r', encoding='UTF-8') as f:
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
  with open(ensembl_filepath_ncbi, 'r', encoding='UTF-8') as f:
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
