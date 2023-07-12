"""
    combine_ensembl_processing.py: Adds refseq id column to an existing ensembl_processing file
"""
__author__ = "Anna Sheaffer"

class RefSeq:
    """Takes a line of ensembl_processing file with the first ele as the ensembl_processing id and the second ele as the refseq id and
    adds pair to a dict"""
    _instance = None
    _currentID = None   # keeps track of the ID of an entry
    _count = 0  # counts the number of entries with a given id

    def __init__(self):
        raise RuntimeError("Call instance() instead")

    @classmethod
    def instance(cls, entry, id):
        if entry is None:
          return cls._instance
        if cls._instance is None:
          cls._instance = {}
        # checks for a new entry
        if len(entry) > 1:
          while len(entry) < 5:
            entry.append(None)
          if id in cls._instance:
            cls._instance[id].append([None if x == '' else x for x in entry])
          else:
            cls._instance[id] = [entry]
            cls._instance[id][0] = [None if x == '' else x for x in cls._instance[id][0]]

        return cls._instance

    @classmethod
    def check_id(cls, ensembl, id):
      """Takes the full ensembl file with the ensembl_processing id as the first element and checks if the id is in _instance,
      returns updated line if it is"""
      if cls._currentID is None or id != cls._currentID:
        cls._currentID = id
        cls._count = 0

      # adds a new entry if one file has more entries for an id than the first
      if cls._count >= len(cls._instance[cls._currentID]):
        ensembl.insert(4, '')
        cls._instance[cls._currentID].append(ensembl)

      else:
        del cls._instance[cls._currentID][cls._count][5:]
        cls._instance[cls._currentID][cls._count] = cls._instance[cls._currentID][cls._count] + ensembl[4:]

      while len(cls._instance[cls._currentID][cls._count]) < 8:
        cls._instance[cls._currentID][cls._count].append(' ')

      cls._count += 1
      return cls._instance

def main():

  with open("r_norvegicus_ensembl_ncbi.txt", 'r') as f:
    # skip first line
    f.readline()

    for line in f:
      line = line.replace('\n', '')
      line = line.split('\t')
      RefSeq.instance(line, line[0])

  with open("r_norvegicus_ensembl_external.txt", 'r+') as f:
    # skip first line
    f.readline()

    result = []
    for line in f:
      line = line.replace('\n', '')
      line = line.split('\t')
      RefSeq.instance(line, line[0])

  result = RefSeq.instance(None, None)
  #print(result)

  with open("../../ensembl/r_norvegicus_ensembl-109.txt", "w") as f:
    f.write("# Combined BioMart data for r_norvegicus version 109\n" +
            "# Elements 1-4 on every line are the Gene stable ID,	Transcript stable ID,	Protein stable ID and	Gene name\n"
            "# On the lines with 5 elements the fifth is the NCBI gene (formerly Entrezgene) ID\n"
            "# On the lines with 7 elements the elements 5-7 are UniProtKB/Swiss-Prot ID,	UniProtKB/TrEMBL ID and	RefSeq peptide ID")
    for id in result:
      for group in result[id]:
        for ele in group:
          f.write(str(ele))
          f.write("\t")
        f.write("\n")

main()

