'''
By Brian Law
'''
import os.path
import python.utility as utility

G_PHYSICAL_CODES = None
G_EXPERIMENTAL_CODES = None


def get_physical_codes():
  """_summary_

  :return: _description_
  """
  global G_PHYSICAL_CODES
  if G_PHYSICAL_CODES is None:
    G_PHYSICAL_CODES, G_EXPERIMENTAL_CODES = retrieve_codes()
  return G_PHYSICAL_CODES

def get_experimental_codes():
  """_summary_

  :return: _description_
  """
  global G_EXPERIMENTAL_CODES
  if G_EXPERIMENTAL_CODES is None:
    G_PHYSICAL_CODES, G_EXPERIMENTAL_CODES = retrieve_codes()
  return G_EXPERIMENTAL_CODES


def read_mi_tree(mi_file: str) -> dict[str, set[str]]:
  """ Read in the provided molecular interactions vocabulary file (i.e. mi.owl, downloaded from EBI.

  :param mi_file: the filepath for the vocabulary file.
  :return: a dictionary mapping each MI term to its immediate parent.
  """
  with open(mi_file, encoding='utf-8') as f:
    mi_tree = {}
    test = {}

    # while line := f.readline():
    #   print(line)

    new_term = None
    for line in f:
      if line[:3] == 'id:':
        new_term = line[4:].strip()
      elif line[:5] == 'is_a:':
        if new_term not in mi_tree:
          mi_tree[new_term] = set()
        parent = line[6:line.find('!')].strip()
        mi_tree[new_term].add(parent)
        if parent not in test:
          test[parent] = [new_term]
        else:
          test[parent].append(new_term)

  # print(test)
  return mi_tree


def process_mi_data(input_obo_filename: str):
  """_summary_

  :param input_obo_filename: _description_
  """
  mi_tree = read_mi_tree(input_obo_filename)

  with open(f'{utility.get_project_root()}/mi/physical_interaction_codes.txt', 'w', encoding='UTF-8') as physical_interaction_file,\
          open(f'{utility.get_project_root()}/mi/experimental_detected_codes.txt', 'w', encoding='UTF-8') as experimental_detected_file:

    physical_code_set = search_physical_codes(mi_tree)
    experimental_code_set = search_experimental_codes(mi_tree)

    for code in physical_code_set:
      physical_interaction_file.write(code + '\n')

    for code in experimental_code_set:
      experimental_detected_file.write(code + '\n')


def search_mi_descendants(mi_tree: dict[str,  set[str]], codes: set[str]) -> set[str]:
  """ Find all the descendants of specified MI codes in a tree of MI codes.

  :param mi_tree: a dictionary mapping MI codes to their parents
  :param codes: the codes for which to find all descendants
  :return: all the codes in the MI-tree that are descendants of the specified codes
  """
  size = 0

  while size != len(codes):
    remove_me = set()
    size = len(codes)
    for term in mi_tree:
      for parent in mi_tree[term]:
        if parent in codes:
          codes.add(term)
          remove_me.add(term)

    for code in remove_me:
      del mi_tree[code]

  return codes


def search_physical_codes(mi_tree: dict[str,  set[str]]):
  """ Get all the MI codes that represent physical molecular interactions.

  :param mi_tree: a dictionary mapping MI codes to their parents
  :return: all the codes in the MI-tree that are descendants of physical interaction codes
  """
  physical_codes = {'MI:0407'}
  # print(sorted(get_mi_descendants(mi_tree.copy(), physical_codes)))

  return search_mi_descendants(mi_tree.copy(), physical_codes)


def search_experimental_codes(mi_tree):
  """ Get all the MI codes that represent experimentally-detected molecular interactions.

  :param mi_tree: a dictionary mapping MI codes to their parents
  :return: all the codes in the MI-tree that are descendants of experimentally-detected interaction codes
  """

  experimental_codes = {'MI:0045'}
  return search_mi_descendants(mi_tree.copy(), experimental_codes)


def retrieve_codes() -> tuple[set[str], set[str]]:
  """
  Extracts data from physical and experimental interaction files
  and places them into two respective sets.

  :return: 2 sets of strings, with MI codes for physical interactions and experimentally-detected interactions
  """

  if not os.path.isfile(f'{utility.get_project_root()}/mi/physical_interaction_codes.txt') or \
     not os.path.isfile(f'{utility.get_project_root()}/mi/experimental_detected_codes.txt'):
    process_mi_data(f'{utility.get_project_root()}/mi/psi-mi.obo')

  # open both phys and exp files to read
  with open(f'{utility.get_project_root()}/mi/physical_interaction_codes.txt', 'r', encoding='UTF-8') as physical_codes, \
          open(f'{utility.get_project_root()}/mi/experimental_detected_codes.txt', 'r', encoding='UTF-8') as experimental_codes:

    # sets (separate from files) where data from files will be stored
    physical_set = set()
    experimental_set = set()

    # add all the physical codes to the physical code set
    for line in physical_codes:
      code = line.strip()
      physical_set.add(code)

    # add all the experimental codes to the experimental code set
    for line in experimental_codes:
      code = line.strip()
      experimental_set.add(code)

    # return the data from each set
    return physical_set, experimental_set


if __name__ == '__main__':
  # psi-mi.obo is downloaded from the EBI Molecular Interactions Ontology
  # GitHub is located: https://github.com/HUPO-PSI/psi-mi-CV
  # process_mi_data(f'{utility.get_project_root()}/mi/psi-mi.obo')

  print(get_physical_codes())
  print(get_experimental_codes())
