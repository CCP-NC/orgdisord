#!/usr/bin/env python
"""Regression tests for supercell generation (PR #4 fixes).

Checks that supercell configs keep molecules whole, use the correct cell
matrix, and keep molecules wrapped inside the cell.
"""
from pathlib import Path

import numpy as np

from orgdisord.enumerate import (
    OrderedfromDisordered,
    connected_components_indices,
    select_configs,
)
from orgdisord.parse_cif_file import CifParser

TESTDIR = Path(__file__).parent


def test_supercell_generation():
    supercell = [2, 1, 1]

    cif = CifParser(TESTDIR / "EROHEA_modified.cif", molecular_crystal=False)
    disordered_structure = cif.get_disordered_structure()
    od = OrderedfromDisordered(disordered_structure, quiet=True)

    images, configs = od.get_supercell_configs(
        supercell=supercell, maxiters=8, return_configs=True
    )

    # supercell cell matrix must scale the cell rows, not the columns
    expected_cell = np.diag(supercell) @ np.asarray(disordered_structure.cell)
    for atoms in images:
        assert np.allclose(np.asarray(atoms.cell), expected_cell)

    # molecules must not be split across periodic boundaries: every molecule
    # of each primitive image must end up in a single connected component of
    # the supercell (the supercell is built tile by tile, so atom indices map
    # directly). NOTE: distinct molecules may merge in clashing configs -- that
    # is the known supercell overlap limitation, so only splitting is checked.
    natoms_prim = len(od.images[0])
    for atoms, config in zip(images, configs):
        comp_of = np.zeros(len(atoms), dtype=int)
        for icomp, comp in enumerate(connected_components_indices(atoms)):
            comp_of[comp] = icomp
        for itile, img_idx in enumerate(
            select_configs(config, supercell, od.ndisordergroups)
        ):
            offset = itile * natoms_prim
            for mol in connected_components_indices(od.images[img_idx]):
                assert len(set(comp_of[np.array(mol) + offset])) == 1

    # each molecule must be wrapped into the cell (centroid in [0, 1)).
    # Molecules that straddle the cell boundary keep atoms outside [0, 1) --
    # that is unavoidable while keeping them whole, so only the centroid is
    # guaranteed.
    inv_cell = np.linalg.inv(expected_cell)
    for atoms in images:
        pos = atoms.get_positions()
        for comp in connected_components_indices(atoms):
            frac = pos[comp].mean(axis=0) @ inv_cell
            assert np.all(frac >= 0.0) and np.all(frac < 1.0)
