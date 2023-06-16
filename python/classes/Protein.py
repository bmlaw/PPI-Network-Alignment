"""
Class that maps all protein id's to one another
"""


class Protein:

  def __init__(self, ids):
    self.name = ids[6]
    self.gene_sid = ids[0]
    self.t_sid = ids[1]
    self.p_sid = ids[2]
    self.ncbi = ids[3]
    self.swissprot = ids[4]
    self.trembl = ids[5]
    self.refseq = ids[7]


  def __eq__(self, o: object) -> bool:
    if not isinstance(o, Protein):
      return False
    return self.gene_sid == o.gene_sid and self.t_sid == o.t_sid and self.p_sid == o.t_sid and self.ncbi == o.ncbi and \
           self.swissprot == o.swissprot and self.trembl == o.trembl and self.refseq == o.refseq and self.name == o.name

  def __str__(self) -> str:
    return self.get_p_sid()

  def __repr__(self) -> str:
    return f'{self.__class__.__name__}({",".join([repr(v) for v in self.__dict__.values()])})'

  def __hash__(self) -> int:
    return sum(hash(x) for x in self.__dict__.values())

  def get_gene_sid(self):
    if self.gene_sid == "":
      return None
    return self.gene_sid

  def get_t_sid(self):
    if self.t_sid == "":
      return None
    return self.t_sid

  def get_p_sid(self):
    if self.p_sid == "":
      return None
    return self.p_sid

  def get_ncbi(self):
    if self.ncbi == "":
      return None
    return self.ncbi

  def get_swissprot(self):
    if self.swissprot == "":
      return None
    return self.swissprot

  def get_trembl(self):
    if self.trembl == "":
      return None
    return self.trembl

  def get_name(self):
    if self.name == "":
      return None
    return self.name

  def get_refseq(self):
    if self.refseq == '':
      return None
    return self.refseq

  # d = {x.get_name(): x for x in list_of_mappings}
  # d['ITSN-1'].get_swissprot()

# if __name__ == '__main__':
#     p = Protein(['1', '2', '3', '4', '5', '6', '7'])
#     print(p)
#     print(p.__dict__)
