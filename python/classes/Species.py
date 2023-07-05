class Species:

  def __init__(self, name: str, long_name: str, short_name: str, taxid: int):
    self.name = name
    self.long_name = long_name
    self.short_name = short_name
    self.taxid = taxid

  def __str__(self) -> str:
    return self.name

  def __repr__(self) -> str:
    return f'Species({self.name}, {self.long_name}, {self.short_name}, {self.taxid})'


species_list = [Species('yeast', 'Saccharomyces cerevisiae S288c', 'S cerevisiae', 4932),
                Species('worm', 'Caenorhabditis elegans', 'C elegans', 6239),
                #Species('zebrafish', 'Danio rerio', 'D rerio', 7955),
                Species('rat', 'Rattus norvegicus', 'R norvegicus', 10116),
                Species('mouse', 'Mus musculus', 'M musculus', 10090),
                Species('human', 'Homo sapiens', 'H sapiens', 9606)]

species_dict = {}
for species in species_list:
  species_dict[species.name] = species
  species_dict[species.short_name] = species
  species_dict[species.long_name] = species
  species_dict[species.taxid] = species

