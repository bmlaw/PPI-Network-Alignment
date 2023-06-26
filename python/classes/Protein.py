"""
Class that maps all protein id's to one another
"""
from __future__ import annotations

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

  def add_ncbi_ids(self, ids: list[str]) -> Protein:
    '''
    Adds IDs to this protein from an Ensembl BioMart download. The list of IDs should be a split of a line from a
    BioMart download. First four elements should be Ensembl gene ID, Ensembl transcript ID, Ensembl gene name.
    Then if the file is an NCBI file, then the fifth element should be the associated NCBI GenBank ID. Otherwise,
    here should be 3 more elements with Swissport, Trembl, or Refseq ID.

    :param ids: list of external IDs from Ensembl BioMart download to be added to this protein
    :return: a reference to this Protein
    '''
    # Adding additional transcript / protein IDs associated with this gene
    if len(ids) >= 2 and len(ids[1].strip()) != 0:
      self.t_ids.add(ids[1])
    if len(ids) >= 3 and len(ids[2].strip()) != 0:
      self.p_sids.add(ids[2])

    # Capitalizing all gene names
    if len(ids) >= 4 and len(ids[3].strip()) != 0:
      self.names.add(ids[3].upper())

    # NCBI file should result in 5-column file, last column is NCBI IDs.
    if len(ids) == 5:
      if len(ids[4].strip()) != 0:
        self.ncbis.add(ids[4])

    return self

  def add_other_ids(self, ids: list[str]):
    '''
    Adds IDs to this protein from an Ensembl BioMart download. The list of IDs should be a split of a line from a
    BioMart download. First four elements should be Ensembl gene ID, Ensembl transcript ID, Ensembl gene name.
    Then if the file is an NCBI file, then the fifth element should be the associated NCBI GenBank ID. Otherwise,
    here should be 3 more elements with Swissport, Trembl, or Refseq ID.

    :param ids: list of external IDs from Ensembl BioMart download to be added to this protein
    :return: a reference to this Protein
    '''
    # Adding additional transcript / protein IDs associated with this gene
    if len(ids) >= 2 and len(ids[1].strip()) != 0:
      self.t_ids.add(ids[1])
    if len(ids) >= 3 and len(ids[2].strip()) != 0:
      self.p_sids.add(ids[2])

    # Capitalizing all gene names
    if len(ids) >= 4 and len(ids[3].strip()) != 0:
      self.names.add(ids[3].upper())

    for i in range(4, len(ids)):
      if len(ids[i].strip()) != 0:
        if i == 4:
          self.swissprots.add(ids[i])
        if i == 5:
          self.trembls.add(ids[i])
        if i == 6:
          self.refseqs.add(ids[i])

    # Otherwise, external ID file should have 3 more columns, a SwissProt column, a Trembl column, and a Refseq column.
    # if len(ids) == 7:
    #   if len(ids[1].strip()) != 0:
    #     self.swissprots.add(ids[4])
    #   if len(ids[1].strip()) != 0:
    #     self.trembls.add(ids[5])
    #   if len(ids[1].strip()) != 0:
    #     self.refseqs.add(ids[6])
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
