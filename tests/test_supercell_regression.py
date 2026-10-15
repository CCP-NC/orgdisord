#!/usr/bin/env python
"""Regression tests for the failure modes of supercell generation seen on `main`.

The five systems (named after the refcodes used in the PR discussion: ESEQEB,
ERITEE, FUWRAR, FADFIZ01, DOVCEX) cover orthorhombic, monoclinic and triclinic
cells, whole-molecule and partial-molecule disorder, and Z = 2, 4 and 8. The
fixtures in tests/disorder_cifs are reduced to the cell, the symmetry operations
and the atomic sites (with disorder assembly/group); they carry no other data.

Failure modes captured (all reproduced on `main` before the fix):

1. supercell cell built as ``cell * sc`` (scales columns) instead of
   ``diag(sc) @ cell``: sheared lattice for any non-orthogonal cell
2. molecules split across periodic boundaries, so mixed-config supercells contain
   fragmented and/or merged (overlaid) molecules
3. the same molecular site placed in different lattice images in different
   configs (disordered molecule near a cell face): atoms land on top of each other
4. unwrapping that tears bonds for molecules spanning more than one cell
   (see test_supercell_fix.py)
5. fixed-ratio sampling returning supercells with the above problems
6. numpy 2 incompatibilities (``np.NaN``, ``np.product``) in the CLI/enumerator
7. ``--prefix`` containing a directory breaking structure writing in the CLI

In the FADFIZ01 and DOVCEX fixtures the two disorder groups have been marked up
by hand: the source data tagged all four alternative O atoms of FADFIZ01 as one
group (-> a 8**8 enumeration) and DOVCEX had no disorder markup at all. Here
O15/O16 vs O1/O2 (FADFIZ01) and H1N vs H2N (DOVCEX) are groups 1 and 2.
"""

import ast
import random
from functools import cache
from pathlib import Path

import numpy as np
import pytest
from ase.data import covalent_radii
from ase.io import read
from ase.neighborlist import neighbor_list
from click.testing import CliRunner
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

import orgdisord
from orgdisord import cli
from orgdisord.enumerate import (
    OrderedfromDisordered,
    component_centroids,
    connected_components_indices,
)
from orgdisord.parse_cif_file import CifParser

TESTDIR = Path(__file__).parent
CIFDIR = TESTDIR / "disorder_cifs"

CASES = {
    "ESEQEB": "mono_cc_z4.cif",  # Cc, beta=102
    "ERITEE": "mono_p21c_z4.cif",  # P21/c, beta=112
    "FUWRAR": "tricl_p1bar_z2.cif",  # P-1 triclinic
    "FADFIZ01": "orth_pbca_z8.cif",  # Pbca, Z=8, ring spans a short axis
    "DOVCEX": "mono_c2c_z8.cif",  # C2/c, Z=8, beta=112.5
}
SMALL = ["ESEQEB", "ERITEE", "FUWRAR"]  # cheap enough for a 2x2x2 supercell
NCONFIGS = 8  # random supercells per test
SEED = 1234


def cif_path(name):
    return CIFDIR / CASES[name]


@cache
def disordered(name, molecular_crystal=True):
    parser = CifParser(cif_path(name), molecular_crystal=molecular_crystal)
    return parser.get_disordered_structure()


def enumerator(name, molecular_crystal=True):
    return OrderedfromDisordered(disordered(name, molecular_crystal), quiet=True)


# ------------------------------------------------------------------ helpers
def _bond_pairs(atoms):
    """Bonded (i, j) pairs under the minimum image convention, H-H excluded."""
    atoms = atoms.copy()
    atoms.pbc = True
    radii = covalent_radii[atoms.numbers] * 1.25
    i, j = neighbor_list("ij", atoms, radii)
    keep = ~((atoms.numbers[i] == 1) & (atoms.numbers[j] == 1))
    return i[keep], j[keep], radii


def molecule_labels(atoms):
    """Periodic connected-component label of every atom."""
    i, j, _ = _bond_pairs(atoms)
    n = len(atoms)
    graph = coo_matrix((np.ones(len(i)), (i, j)), shape=(n, n))
    return connected_components(graph, directed=False)


def molecule_formulas(atoms, labels):
    symbols = np.array(atoms.get_chemical_symbols())
    return {"".join(sorted(symbols[labels == label])) for label in np.unique(labels)}


def min_intermolecular_distance(atoms, labels):
    atoms = atoms.copy()
    atoms.pbc = True
    dist = atoms.get_all_distances(mic=True)
    np.fill_diagonal(dist, np.inf)
    dist[labels[:, None] == labels[None, :]] = np.inf
    return dist.min()


def torn_bonds(atoms):
    """Bonds that exist under the MIC but are not short in plain Cartesian space."""
    i, j, radii = _bond_pairs(atoms)
    d = np.linalg.norm(atoms.positions[i] - atoms.positions[j], axis=1)
    return int((d > radii[i] + radii[j]).sum() // 2)


def supercell_problems(supercell_atoms, images, supercell):
    """Everything that must hold for a supercell built from primitive `images`.

    Returns a list of problem descriptions (empty if the supercell is sound).
    `images` are the primitive configs the supercell was tiled from.
    """
    problems = []
    nsuper = int(np.prod(supercell))
    image0 = images[0]
    prim_cell = np.asarray(image0.cell)

    # 1. cell scales the lattice vectors (rows), keeping the angles
    expected_cell = np.diag(supercell) @ prim_cell
    if not np.allclose(np.asarray(supercell_atoms.cell), expected_cell, atol=1e-6):
        problems.append("cell is not diag(sc) @ cell")

    # composition conserved
    if len(supercell_atoms) != nsuper * len(image0):
        problems.append("wrong number of atoms")
    expected = sorted(image0.get_chemical_symbols() * nsuper)
    if sorted(supercell_atoms.get_chemical_symbols()) != expected:
        problems.append("wrong composition")

    # 2. molecules whole and distinct: right number of connected components,
    # and each one is a complete molecule
    ncomp0, labels0 = molecule_labels(image0)
    ncomp, labels = molecule_labels(supercell_atoms)
    if ncomp != nsuper * ncomp0:
        problems.append(f"{ncomp} molecules, expected {nsuper * ncomp0}")
    if not molecule_formulas(supercell_atoms, labels) <= molecule_formulas(
        image0, labels0
    ):
        problems.append("fragmented molecules")

    # 3. no artificial overlaps: closest approach between different molecules is
    # comparable to the one already present in the primitive images
    baseline = min(
        min_intermolecular_distance(img, molecule_labels(img)[1]) for img in images[:16]
    )
    dmin = min_intermolecular_distance(supercell_atoms, labels)
    if dmin < 0.8 * baseline:
        problems.append(f"overlap: {dmin:.2f} A vs primitive {baseline:.2f} A")
    return problems


def make_supercells(name, supercell, molecular_crystal=True, **kwargs):
    random.seed(SEED)
    od = enumerator(name, molecular_crystal)
    kwargs.setdefault("random_configs", True)
    kwargs.setdefault("maxiters", NCONFIGS)
    supercells, configs = od.get_supercell_configs(
        supercell=supercell, return_configs=True, **kwargs
    )
    return od, supercells, configs


# -------------------------------------------------------------------- tests
@pytest.mark.parametrize("name", CASES)
def test_primitive_images(name):
    """Primitive images: whole molecules, none torn, same lattice image everywhere."""
    od = enumerator(name, molecular_crystal=False)
    images = od.get_all_configs()
    assert len(images) > 1
    natoms = len(images[0])
    ref_comps = connected_components_indices(images[0])
    ref_centroids = component_centroids(images[0], ref_comps)
    cell = np.asarray(images[0].cell)
    inv_cell = np.linalg.inv(cell)
    half_shortest = 0.5 * np.linalg.norm(cell, axis=1).min()

    for atoms in images:
        assert len(atoms) == natoms
        assert torn_bonds(atoms) == 0
        assert len(connected_components_indices(atoms)) == len(ref_comps)
        for centroid in component_centroids(atoms, connected_components_indices(atoms)):
            # nearest reference molecule in plain Cartesian space must be the
            # nearest one under the minimum image convention too, i.e. the same
            # molecular site is in the same lattice image in every config
            cart = np.linalg.norm(ref_centroids - centroid, axis=1)
            frac = (ref_centroids - centroid) @ inv_cell
            mic = np.linalg.norm((frac - np.round(frac)) @ cell, axis=1)
            assert np.isclose(cart.min(), mic.min(), atol=1e-6)
            assert cart.min() < half_shortest


@pytest.mark.parametrize("name", CASES)
@pytest.mark.parametrize("molecular_crystal", [True, False])
def test_exhaustive_1x1x1(name, molecular_crystal):
    od, supercells, _ = make_supercells(
        name, [1, 1, 1], molecular_crystal, random_configs=False, maxiters=16
    )
    assert len(supercells) == min(16, len(od.images))
    for atoms in supercells:
        assert supercell_problems(atoms, od.images, [1, 1, 1]) == []


SUPERCELLS = [
    pytest.param(name, sc, mc, id=f"{name}-{'x'.join(map(str, sc))}-mc{int(mc)}")
    for name in CASES
    for sc, mc in (([2, 1, 1], True), ([2, 1, 1], False), ([1, 1, 4], True))
] + [pytest.param(name, [2, 2, 2], True, id=f"{name}-2x2x2-mc1") for name in SMALL]


@pytest.mark.parametrize("name,supercell,molecular_crystal", SUPERCELLS)
def test_mixed_supercells(name, supercell, molecular_crystal):
    """Mixed-config supercells: right cell, whole and non-overlapping molecules."""
    od, supercells, configs = make_supercells(name, supercell, molecular_crystal)
    assert len(supercells) == NCONFIGS
    for atoms, _config in zip(supercells, configs):
        assert supercell_problems(atoms, od.images, supercell) == []
    # the sample must really mix disorder groups between tiles
    assert any(len(set(config)) > 1 for config in configs)


@pytest.mark.parametrize(
    "name,supercell,ratio",
    [(name, [2, 1, 1], 0.5) for name in CASES]
    + [(name, [4, 1, 1], 0.25) for name in SMALL if name != "FUWRAR"],
)
def test_fixed_ratio_sampling(name, supercell, ratio):
    """Fixed mole-fraction sampling returns the ratio asked for, with sound supercells."""
    od, supercells, configs = make_supercells(
        name, supercell, fix_ratio=True, ratio=ratio, ratio_tolerance=0.01
    )
    assert len(supercells) > 0
    for atoms, config in zip(supercells, configs):
        counts = np.bincount(config, minlength=2)
        assert np.allclose(counts / counts.sum(), [ratio, 1 - ratio], atol=0.01)
        assert supercell_problems(atoms, od.images, supercell) == []


@pytest.mark.parametrize("name", ["ESEQEB", "FUWRAR"])
@pytest.mark.parametrize("merge", [False, True])
def test_cli_enumerate_with_directory_prefix(tmp_path, name, merge):
    """CLI end to end (numpy 2 safe) with a prefix that includes a directory."""
    prefix = tmp_path / "does" / "not" / "exist" / "run"
    args = [
        "enumerate",
        str(cif_path(name)),
        "-s", "2", "1", "1",
        "-N", "4",
        "--random",
        "--prefix", str(prefix),
        "--log", str(tmp_path / "orgdisord.log"),
        "--quiet",
    ]  # fmt: skip
    if merge:
        args.append("--merge")
    result = CliRunner().invoke(cli.cli, args)
    assert result.exit_code == 0, result.output + repr(result.exception)

    assert Path(f"{prefix}.csv").is_file()
    written = sorted(Path(f"{prefix}-results").glob("*.xyz"))
    assert 1 <= len(written) <= 4
    assert all(path.name.startswith("run-") for path in written)
    # nothing may be written to a doubled path
    assert not (Path(f"{prefix}-results") / str(prefix).lstrip("/")).exists()
    natoms_prim = len(enumerator(name).get_all_configs()[0])
    assert len(read(written[0])) == 2 * natoms_prim


def test_no_numpy2_removed_aliases():
    """No use of numpy aliases removed in numpy 1.24/2.0 (np.NaN, np.product, ...)."""
    removed = {
        "NaN", "NAN", "Inf", "Infinity", "infty", "PINF", "NINF", "PZERO", "NZERO",
        "product", "cumproduct", "sometrue", "alltrue", "in1d", "trapz", "row_stack",
        "asscalar", "msort", "find_common_type", "cast", "unicode_", "string_",
        "float_", "complex_", "mat", "float", "int", "bool", "object", "str",
        "complex", "long", "unicode",
    }  # fmt: skip
    offenders = []
    for path in sorted(Path(orgdisord.__file__).parent.glob("*.py")):
        for node in ast.walk(ast.parse(path.read_text())):
            if (
                isinstance(node, ast.Attribute)
                and isinstance(node.value, ast.Name)
                and node.value.id in ("np", "numpy")
                and node.attr in removed
            ):
                offenders.append(f"{path.name}:{node.lineno} np.{node.attr}")
    assert offenders == []
