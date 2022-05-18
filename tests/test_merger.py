#!/usr/bin/env python

"""Tests for the merge.py module.

I used the enumerate module to 
generate 50 configs containing 
34 random configs and 16 exhaustive configs.

There should merge to just 5 groups of symmetry-inequivalent configs.
"""
import pytest
from sodorg_renewal.merge import merge_structures
from ase.spacegroup import get_spacegroup
import numpy as np
from ase.io import read
import spglib

original_atoms = read('tests/EROHEA_modified.cif')
images = read('tests/images_50-5_groups.traj', index=':')
# more copies?
# images += images
# images += images


symops = get_spacegroup(original_atoms).get_symop()
groups = merge_structures(
                        images,
                        algo='symm',
                        symops=symops,
                        use_disordered_only = True,
                        symprec=1e-4,
                        verbose=False)


print(f'Found {len(groups)} groups')
for g in groups:
    print(f"Spacegroup: {spglib.get_spacegroup(g[0])}\t multiplicity: {g[1]}")

# Sanity check!
# note that I have run single point energy calculations on the 16 exhaustive configs and 
# verified that the energy grouping is exactly the same as the grouping by symmetry (phew!).
assert len(groups) == 5



print('REmatch approach')
groups = merge_structures(
                        images,
                        algo='rematch',
                        symops=symops,
                        use_disordered_only = True,
                        symprec=1e-2,
                        verbose=False)

print(f'Found {len(groups)} groups')
for g in groups:
    print(f"Spacegroup: {spglib.get_spacegroup(g[0])}\t multiplicity: {g[1]}")
assert len(groups) == 5



print('EWALD approach')
groups = merge_structures(
                        images,
                        oxidation_states = {'N': -3.0, 'H': 1.0, 'C': 1.0, 'O': -2.0},
                        algo='ewald',
                        symops=symops,
                        use_disordered_only = True,
                        symprec=1e-2,
                        verbose=False)

print(f'Found {len(groups)} groups')
for g in groups:
    print(f"Spacegroup: {spglib.get_spacegroup(g[0])}\t multiplicity: {g[1]}")
assert len(groups) == 5
