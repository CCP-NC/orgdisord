#!/usr/bin/env python
"""Regression tests for supercell generation (PR #4 fixes).

Checks that supercell configs keep molecules whole, use the correct cell
matrix, keep every molecular site in the same lattice image across configs
(so configs can be mixed tile by tile without artificial overlaps), and that
molecules are unwrapped correctly however far they extend.
"""

import contextlib
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from ase import Atoms
from ase.data import covalent_radii

from orgdisord.enumerate import (
    OrderedfromDisordered,
    component_centroids,
    connected_components_indices,
    make_whole_molecule,
    select_configs,
)
from orgdisord.parse_cif_file import CifParser

TESTDIR = Path(__file__).parent
MATRIX_ARGS = ["1", "0", "-1", "0", "2", "0", "1", "0", "1"]


def _od(cif, **kwargs):
    parser = CifParser(TESTDIR / cif, **kwargs)
    disordered_structure = parser.get_disordered_structure()
    return disordered_structure, OrderedfromDisordered(disordered_structure, quiet=True)


def _min_heavy_distance(atoms):
    """Smallest periodic distance between two non-hydrogen atoms."""
    atoms = atoms.copy()
    atoms.pbc = True
    d = atoms.get_all_distances(mic=True)
    np.fill_diagonal(d, np.inf)
    heavy = atoms.get_atomic_numbers() > 1
    return d[np.ix_(heavy, heavy)].min()


@pytest.mark.parametrize("molecular_crystal", [False, True])
def test_supercell_generation(molecular_crystal):
    supercell = [2, 1, 1]
    disordered_structure, od = _od(
        "EROHEA_modified.cif", molecular_crystal=molecular_crystal
    )

    images, configs = od.get_supercell_configs(
        supercell=supercell, maxiters=256, return_configs=True
    )
    assert len(images) == 256

    # supercell cell matrix must scale the cell rows, not the columns
    expected_cell = np.diag(supercell) @ np.asarray(disordered_structure.cell)
    for atoms in images:
        assert np.allclose(np.asarray(atoms.cell), expected_cell)

    # no artificial overlaps when mixing configs between tiles: a disordered
    # molecule must not be placed one lattice vector away in different configs
    # (EROHEA has molecules with centroids ~0.001 from a cell face).
    for atoms in images:
        assert _min_heavy_distance(atoms) > 1.0

    # molecules must not be split across periodic boundaries: every molecule
    # of each primitive image must end up in a single connected component of
    # the supercell (the supercell is built tile by tile, so atom indices map
    # directly when the primitive images keep their atom order).
    if not molecular_crystal:
        natoms_prim = len(od.images[0])
        for atoms, config in list(zip(images, configs))[:32]:
            comp_of = np.zeros(len(atoms), dtype=int)
            for icomp, comp in enumerate(connected_components_indices(atoms)):
                comp_of[comp] = icomp
            for itile, img_idx in enumerate(
                select_configs(config, int(np.prod(supercell)), od.ndisordergroups)
            ):
                offset = itile * natoms_prim
                for mol in connected_components_indices(od.images[img_idx]):
                    assert len(set(comp_of[np.array(mol) + offset])) == 1


def test_primitive_images_share_lattice_image():
    """Each molecular site must sit in the same lattice image in every config."""
    _, od = _od("EROHEA_modified.cif", molecular_crystal=False)
    images = od.get_all_configs()
    cell = np.asarray(images[0].cell)
    inv_cell = np.linalg.inv(cell)

    def centroids(atoms):
        return component_centroids(atoms, connected_components_indices(atoms))

    ref = centroids(images[0])
    for atoms in images[1:]:
        cents = centroids(atoms)
        assert len(cents) == len(ref)
        for c in cents:
            # nearest reference centroid in plain Cartesian space must also be
            # the nearest one under the minimum image convention
            cart = np.linalg.norm(ref - c, axis=1)
            frac = (ref - c) @ inv_cell
            mic = np.linalg.norm((frac - np.round(frac)) @ cell, axis=1)
            assert np.isclose(cart.min(), mic.min(), atol=1e-6)
            assert cart.min() < 0.5 * np.linalg.norm(cell, axis=1).min()


@pytest.mark.parametrize("cif", ["704722.cif", "VAGKUM.cif"])
def test_supercell_cell_matches_primitive_angles(cif):
    """Non-orthogonal cells: supercell must be diag(sc) @ cell for any supercell."""
    disordered_structure, od = _od(cif, molecular_crystal=False)
    prim = np.asarray(disordered_structure.cell)
    prim_angles = disordered_structure.cell.cellpar()[3:]
    for supercell in ([2, 1, 1], [1, 1, 2], [1, 2, 3]):
        images = od.get_supercell_configs(
            supercell=supercell, maxiters=2, random_configs=True
        )
        for atoms in images:
            assert np.allclose(np.asarray(atoms.cell), np.diag(supercell) @ prim)
            assert np.allclose(atoms.cell.cellpar()[3:], prim_angles)


def test_tile_mapping_is_bijective():
    """Every tile of the supercell must receive exactly one primitive image."""
    _, od = _od("704722.cif", molecular_crystal=False)
    supercell = [2, 3, 2]
    od.get_supercell_configs(supercell=supercell, maxiters=1, random_configs=True)
    nsuper = int(np.prod(supercell))
    # all tiles different: tile k gets image k % nimages -> distinct translations
    config = [k % od.ndisordergroups for k in range(od.disordered_structure.Z * nsuper)]
    atoms = od.get_supercell_config(config)
    natoms_prim = len(od.images[0])
    # the first atom of every tile (before translation) is different; compare
    # tile origins via the centroid shift relative to the unshifted image
    cell = np.asarray(od.cell)
    shifts = set()
    for itile, img_idx in enumerate(select_configs(config, nsuper, od.ndisordergroups)):
        tile = atoms[itile * natoms_prim : (itile + 1) * natoms_prim]
        shift = tile.get_positions() - od.images[img_idx].get_positions()
        frac = np.round(shift[0] @ np.linalg.inv(cell)).astype(int)
        assert np.allclose(shift, shift[0])
        shifts.add(tuple(frac))
    assert shifts == {(a, b, c) for a in range(2) for b in range(3) for c in range(2)}


def _zigzag_alkane(n=16):
    """All-trans CnH2n+2 along x (C-C 1.54, C-C-C 112 deg, C-H 1.09)."""
    cc, ang, ch = 1.54, np.radians(112.0), 1.09
    dx, dy = cc * np.sin(ang / 2), cc * np.cos(ang / 2)
    theta = np.radians(54.75)
    pos, sym = [], []
    for i in range(n):
        s = 1 if i % 2 == 0 else -1
        c = np.array([i * dx, s * dy / 2, 0.0])
        pos.append(c)
        sym.append("C")
        for sz in (-1, 1):
            pos.append(c + ch * np.array([0.0, s * np.cos(theta), sz * np.sin(theta)]))
            sym.append("H")
    # terminal hydrogens along the chain axis
    pos.append(pos[0] + np.array([-ch, 0.0, 0.0]))
    sym.append("H")
    pos.append(pos[3 * (n - 1)] + np.array([ch, 0.0, 0.0]))
    sym.append("H")
    return Atoms(sym, pos)


@pytest.mark.parametrize("ax", [5.0, 8.0, 12.0])
@pytest.mark.parametrize("cx", [0.0, 5.0])
def test_make_whole_molecule_long_molecule(ax, cx):
    """A molecule spanning several cells must be unwrapped without tearing bonds.

    The a vector has a lateral component so that the periodic images of the chain
    pack side by side (as in a real crystal) rather than overlapping.
    """
    original = _zigzag_alkane(16)
    original.set_cell([[ax, 5.5, 0.0], [0.0, 30.0, 0.0], [cx, 0.0, 20.0]])
    original.pbc = True

    wrapped = original.copy()
    wrapped.wrap()
    whole = make_whole_molecule(wrapped)

    radii = covalent_radii[original.get_atomic_numbers()]
    cutoff = 1.25 * (radii[:, None] + radii[None, :])
    d0 = np.linalg.norm(
        original.positions[:, None] - original.positions[None, :], axis=-1
    )
    d1 = np.linalg.norm(whole.positions[:, None] - whole.positions[None, :], axis=-1)
    z = original.get_atomic_numbers()
    bonded = (d0 < cutoff) & ~np.eye(len(original), dtype=bool)
    bonded &= ~((z[:, None] == 1) & (z[None, :] == 1))
    assert bonded.sum() > 0
    assert np.allclose(d0[bonded], d1[bonded], atol=1e-6)
    assert len(connected_components_indices(whole)) == 1


def test_optimal_supercell_generation():
    """Target size integer generates the supercell closest to the target shape."""
    disordered_structure, od = _od("EROHEA_modified.cif", molecular_crystal=False)
    prim_cell = np.asarray(disordered_structure.cell)

    images, configs = od.get_supercell_configs(
        supercell=4, maxiters=16, return_configs=True
    )
    assert len(images) == 16
    assert od.nsuper == 4
    # transformation matrix should be 3x3 with det == 4
    assert np.shape(od.transformation_matrix) == (3, 3)
    assert int(round(abs(np.linalg.det(od.transformation_matrix)))) == 4

    expected_cell = od.transformation_matrix @ prim_cell
    for atoms in images:
        assert np.allclose(np.asarray(atoms.cell), expected_cell)
        assert len(atoms) == len(od.images[0]) * 4
        assert _min_heavy_distance(atoms) > 1.0


def test_explicit_non_diagonal_matrix():
    """Explicit 3x3 integer matrix generates expected non-diagonal supercell."""
    disordered_structure, od = _od("EROHEA_modified.cif", molecular_crystal=False)
    prim_cell = np.asarray(disordered_structure.cell)
    P = np.array([[1, 0, -1], [0, 2, 0], [1, 0, 1]], dtype=int)

    images, configs = od.get_supercell_configs(
        supercell=P, maxiters=8, return_configs=True
    )
    assert len(images) == 8
    assert od.nsuper == 4
    assert np.array_equal(od.transformation_matrix, P)

    expected_cell = P @ prim_cell
    for atoms in images:
        assert np.allclose(np.asarray(atoms.cell), expected_cell)
        assert len(atoms) == len(od.images[0]) * 4
        assert _min_heavy_distance(atoms) > 1.0


def _parse_supercell(text):
    """Parse the 'Supercell' column: either 'axbxc' or 'r00,r01,r02;r10,...'."""
    if ";" in text:
        return np.array([[int(v) for v in row.split(",")] for row in text.split(";")])
    return np.diag([int(v) for v in text.split("x")])


def test_cli_non_diagonal_options(tmp_path):
    """CLI supports --target-size and --matrix and records the supercell matrix."""
    from click.testing import CliRunner

    from orgdisord.cli import cli

    runner = CliRunner()
    cif_path = str(TESTDIR / "EROHEA_modified.cif")
    with contextlib.chdir(tmp_path):
        res = runner.invoke(
            cli,
            [
                "enumerate",
                cif_path,
                "--target-size",
                "2",
                "-N",
                "2",
                "--no_write",
                "-q",
            ],
        )
        assert res.exception is None, res.output
        assert res.exit_code == 0, res.output
        P = _parse_supercell(pd.read_csv("orgdisord.csv")["Supercell"].iloc[0])
        assert int(round(abs(np.linalg.det(P)))) == 2

        res = runner.invoke(
            cli,
            [
                "enumerate",
                cif_path,
                "--matrix",
                *MATRIX_ARGS,
                "-N",
                "2",
                "--no_write",
                "-q",
            ],
        )
        assert res.exception is None, res.output
        assert res.exit_code == 0, res.output
        supercell = pd.read_csv("orgdisord.csv")["Supercell"].iloc[0]
        assert supercell == "1,0,-1;0,2,0;1,0,1"


@pytest.mark.parametrize(
    "extra",
    [
        ["-s", "2", "1", "1", "--target-size", "2"],
        ["--matrix", *MATRIX_ARGS, "--target-size", "2"],
    ],
)
def test_cli_conflicting_size_options(extra):
    from click.testing import CliRunner

    from orgdisord.cli import cli

    runner = CliRunner()
    cif_path = str(TESTDIR / "EROHEA_modified.cif")
    res = runner.invoke(cli, ["enumerate", cif_path, *extra, "--no_write", "-q"])
    assert res.exit_code != 0
    assert "Use only one of" in res.output


def test_supercell_argument_validation():
    _, od = _od("EROHEA_modified.cif", molecular_crystal=False)
    with pytest.raises(ValueError):
        od.get_supercell_configs(supercell=[1, 1, 1], target_size=2)
    with pytest.raises(ValueError):
        od.get_supercell_configs(supercell=np.eye(3) * 1.5)
    with pytest.raises(ValueError):
        od.get_supercell_configs(supercell=True)


@pytest.mark.parametrize("result", [ValueError("matmul failed"), None])
def test_no_optimal_supercell_found(monkeypatch, result):
    import orgdisord.enumerate as enumerate_mod

    def fake_find(*args, **kwargs):
        if isinstance(result, Exception):
            raise result
        return result

    monkeypatch.setattr(enumerate_mod, "find_optimal_cell_shape", fake_find)
    _, od = _od("EROHEA_modified.cif", molecular_crystal=False)
    with pytest.raises(ValueError, match="Could not find a supercell"):
        od.get_supercell_configs(supercell=4, maxiters=8)
