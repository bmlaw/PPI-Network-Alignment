"""
Class that maps all protein id's to one another
"""


class Protein:

  def __init__(self, gene_id: str):
    self.gene_id = gene_id
    self.t_ids = set()
    self.p_sids = set()
    self.ncbis = set()
    self.swissprots = set()
    self.trembls = set()
    self.names = set()
    self.refseqs = set()

  def add_ids(self, list):
    self.t_ids.add(list[1])
    self.p_sids.add(list[2])
    if len(list) == 5:
      self.names.add(list[3])
      self.ncbis.add(list[4])
    elif len(list) == 7:
      self.t_ids.add(list[1])
      self.p_sids.add(list[2])
      self.names.add(list[3])
      self.swissprots.add(list[4])
      self.trembls.add(list[5])
      self.refseqs.add(list[6])
    return self

  def __eq__(self, o: object) -> bool:
    if not isinstance(o, Protein):
      return False
    return self.gene_id == o.gene_id

  def __repr__(self) -> str:
    return f'{self.__class__.__name__}({",".join([repr(v) for v in self.__dict__.values()])})'

  def __hash__(self) -> int:
    return hash(self.gene_id)

  # d = {x.get_name(): x for x in list_of_mappings}
  # d['ITSN-1'].get_swissprot()

# if __name__ == '__main__':
#     p = Protein(['1', '2', '3', '4', '5', '6', '7'])
#     print(p)
#     print(p.__dict__)
