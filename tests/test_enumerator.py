#!/usr/bin/env python

"""Tests for the enumerate.py module."""
import pytest
from orgdisord.enumerate import OrderedfromDisordered
from orgdisord.parse_cif_file import CifParser
from orgdisord.disordered_structure import from_disorder_components
import numpy as np
from ase.io import read
from ase.spacegroup import Spacegroup

cif = CifParser("tests/EROHEA_modified.cif")
disordered_structure = cif.get_disordered_structure()
od = OrderedfromDisordered(disordered_structure)


# exhaustive enumeration in 111 supercell:
supercell = [1, 1, 1]
images = od.get_supercell_configs(
    supercell=supercell,
    maxiters=100,  # should be more than enough
    exclude_ordered=False,
    random_configs=False,
)

# there should be 16 configurations generated
# (binary options and four possible sites, 2^4 = 16)
assert len(images) == 16
# check that the number of sites is correct (including ordered sites)
natoms = 288 * supercell[0] * supercell[1] * supercell[2]
assert all([len(atoms) == natoms for atoms in images])
# check that the number of ordered sites is correct
# (ordered sites are tagged 0)
tag = 0
nordered = 252 * supercell[0] * supercell[1] * supercell[2]
assert len(np.where(images[0].get_tags() == tag)[0]) == nordered


# Try a supercell of [2,1,1]
supercell = [2, 1, 1]
images = od.get_supercell_configs(
    supercell=supercell,
    maxiters=500,  # should be more than enough
    exclude_ordered=False,
    random_configs=False,
)

# there should be 100 configurations generated (< 2^8 = 256)
assert len(images) == 256
# check that the number of sites is correct (including ordered sites)
natoms = 288 * supercell[0] * supercell[1] * supercell[2]
assert all([len(atoms) == natoms for atoms in images])
# check that the number of ordered sites is correct
# (ordered sites are tagged 0)
tag = 0
nordered = 252 * supercell[0] * supercell[1] * supercell[2]
assert len(np.where(images[0].get_tags() == tag)[0]) == nordered


# Random configurations
supercell = [1, 1, 1]
images = od.get_supercell_configs(
    supercell=supercell,
    maxiters=34,  # should be more than enough
    exclude_ordered=False,
    random_configs=True,
)

# requested 34 random configurations
assert len(images) == 34
# check that the number of sites is correct (including ordered sites)
natoms = 288 * supercell[0] * supercell[1] * supercell[2]
assert all([len(atoms) == natoms for atoms in images])
# check that the number of ordered sites is correct
# (ordered sites are tagged 0)
tag = 0
nordered = 252 * supercell[0] * supercell[1] * supercell[2]
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


cif = CifParser("tests/VIFQIL01.cif")
disordered_structure = cif.get_disordered_structure()
od = OrderedfromDisordered(disordered_structure)


# exhaustive enumeration in 111 supercell:
supercell = [1, 1, 1]
images = od.get_supercell_configs(
    supercell=supercell,
    maxiters=100,  # should be more than enough
    exclude_ordered=False,
    random_configs=False,
)

# there should be 16 configurations generated
# (binary options and four possible sites, 2^4 = 16)
assert len(images) == 16

# check that the number of sites is correct (including ordered sites)
natoms = 112 * supercell[0] * supercell[1] * supercell[2]
assert all([len(atoms) == natoms for atoms in images])

# check that the number of ordered sites is correct
# (ordered sites are tagged 0)
tag = 0
nordered = 100 * supercell[0] * supercell[1] * supercell[2]


########################################################################################
# Test enumeration from two ordered structures
atoms_maj = read("tests/ABABUB_maj.xyz")
atoms_min = read("tests/ABABUB_min.xyz")
disordered_structure = from_disorder_components(atoms_maj, atoms_min)

# make sure the disordered structure is correct:
# check number of assemblies and groups is correct:
assert disordered_structure.get_number_of_disorder_groups_per_assembly() == [2]
# check that the correct Z and number of symops are assigned
assert disordered_structure.Z == 4
assert disordered_structure.spacegroup == Spacegroup(14)
# check that each group contains the correct number of sites
assert len(disordered_structure.disorder_assemblies[0].disorder_groups[0].atoms) == 12
assert len(disordered_structure.disorder_assemblies[0].disorder_groups[1].atoms) == 12

od = OrderedfromDisordered(disordered_structure)
# exhaustive enumeration in 111 supercell:
supercell = [1, 1, 1]
images = od.get_supercell_configs(
    supercell=supercell,
    maxiters=100,  # should be more than enough
    exclude_ordered=False,
    random_configs=False,
)

# there should be 16 configurations generated
# (binary options and four possible sites, 2^4 = 16)
assert len(images) == 16
########################################################################################


#  Explicitly create DisorderedAtoms object
from ase.build import bulk
from ase.build import molecule

# create a disordered structure
atoms = bulk("Cu", cubic=True)
atoms *= (2, 2, 2)

# add a disordered molecule
atoms += molecule("H2O")
