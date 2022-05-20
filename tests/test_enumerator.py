#!/usr/bin/env python

"""Tests for the enumerate.py module."""
import pytest
from sodorg_renewal.enumerate import OrderedfromDisordered
from sodorg_renewal.parse_cif_file import CifParser
import numpy as np

cif = CifParser('tests/EROHEA_modified.cif')
od = OrderedfromDisordered(cif)




# exhaustive enumeration in 111 supercell:
supercell = [1,1,1]
images = od.get_supercell_configs(
                supercell = supercell, 
                maxiters = 100, # should be more than enough 
                exclude_ordered = False, 
                random_configs=False)

# there should be 16 configurations generated
#(binary options and four possible sites, 2^4 = 16)
assert len(images) == 16
# check that the number of sites is correct (including ordered sites)
natoms = 288*supercell[0]*supercell[1]*supercell[2]
assert all([len(atoms) == natoms for atoms in images])
# check that the number of ordered sites is correct
# (ordered sites are tagged 0)
tag = 0
nordered = 252*supercell[0]*supercell[1]*supercell[2]
assert len(np.where(images[0].get_tags() == tag)[0]) == nordered







# Try a supercell of [2,1,1]
supercell = [2,1,1]
images = od.get_supercell_configs(
                supercell = supercell, 
                maxiters = 500, # should be more than enough 
                exclude_ordered = False, 
                random_configs=False)

# there should be 100 configurations generated (< 2^8 = 256)
assert len(images) == 256
# check that the number of sites is correct (including ordered sites)
natoms = 288*supercell[0]*supercell[1]*supercell[2]
assert all([len(atoms) == natoms for atoms in images])
# check that the number of ordered sites is correct
# (ordered sites are tagged 0)
tag = 0
nordered = 252*supercell[0]*supercell[1]*supercell[2]
assert len(np.where(images[0].get_tags() == tag)[0]) == nordered





# Random configurations
supercell = [1,1,1]
images = od.get_supercell_configs(
                supercell = supercell, 
                maxiters = 34, # should be more than enough 
                exclude_ordered = False, 
                random_configs=True)

# requested 34 random configurations
assert len(images) == 34
# check that the number of sites is correct (including ordered sites)
natoms = 288*supercell[0]*supercell[1]*supercell[2]
assert all([len(atoms) == natoms for atoms in images])
# check that the number of ordered sites is correct
# (ordered sites are tagged 0)
tag = 0
nordered = 252*supercell[0]*supercell[1]*supercell[2]
assert len(np.where(images[0].get_tags() == tag)[0]) == nordered



# from ase.io import write
# supercell = [1,1,1]
# images_explicit = od.get_supercell_configs(
#                         supercell = supercell, 
#                         maxiters = 100, # should be more than enough 
#                         exclude_ordered = False, 
#                         random_configs=False)
# temp = images + images_explicit

# write('tests/images_50-5_groups.traj', temp)




cif = CifParser('tests/VIFQIL01.cif')
od = OrderedfromDisordered(cif)




# exhaustive enumeration in 111 supercell:
supercell = [1,1,1]
images = od.get_supercell_configs(
                supercell = supercell, 
                maxiters = 100, # should be more than enough 
                exclude_ordered = False, 
                random_configs=False)

# there should be 16 configurations generated
# (binary options and four possible sites, 2^4 = 16)
assert len(images) == 16

# check that the number of sites is correct (including ordered sites)
natoms = 112*supercell[0]*supercell[1]*supercell[2]
assert all([len(atoms) == natoms for atoms in images])

# check that the number of ordered sites is correct
# (ordered sites are tagged 0)
tag = 0
nordered = 100*supercell[0]*supercell[1]*supercell[2]
