__author__ = "Brian Law"
__email__ = "blaw@iwu.edu"

import os

from python.classes.Species import Species, species_dict, species_list
from python.classes import Protein

class Network:

  def __init__(self, species, ensembl=None, biogrid=None):
    self.edges = {}

  