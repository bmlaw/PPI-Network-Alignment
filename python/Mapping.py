

"""
Class that maps all protein id's to one another
"""

class Mapping:
    def __init__(self, gene_sid, t_sid, p_sid, ncbi, swissprot, trembl, name):
        self.gene_sid = gene_sid
        self.t_sid = t_sid
        self.p_sid = p_sid
        self.ncbi = ncbi
        self.swissprot = swissprot
        self.trembl = trembl
        self.name = name

    def get_gene_sid(self):
        return self.gene_sid

    def get_t_sid(self):
        return self.t_sid

    def get_p_sid(self):
        return self.p_sid

    def get_ncbi(self):
        return self.ncbi

    def get_swissprot(self):
        return self.swissprot

    def get_trembl(self):
        return self.trembl

    def get_name(self):
        return self.name