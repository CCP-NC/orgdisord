"""Main module.
Parses the cif file using the parse_cif_file module and
generates the enumerated ordered confugurations.
"""
import numpy as np
import warnings
from parse_cif_file import CifParser
from enumerate import OrderedfromDisordered
from merge import merge_structures

cif_file = 'tests/EROHEA_modified.cif'

cif = CifParser(cif_file)
od = OrderedfromDisordered(cif)
