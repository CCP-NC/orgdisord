#!/usr/bin/env python
"""Tests for the merge.py module."""

import logging

import numpy as np
import pytest
from ase import Atoms
from ase.io import read
from ase.spacegroup import Spacegroup, get_spacegroup
from click.testing import CliRunner
from scipy.spatial import cKDTree

from orgdisord import cli
from orgdisord.enumerate import OrderedfromDisordered
from orgdisord.merge import (
    NonDiscreteSitesError,
    _wrap,
    merge_structures,
    supercell_symops,
)
from orgdisord.parse_cif_file import CifParser


def _partition(groups):
    return sorted(sorted(int(i) for i in g) for g in groups)


def _is_closed_group(ops, tol=1e-8):
    """True if the composition of any two operations is again in ``ops``."""

    def known(rot, trans):
        for r, t in ops:
            d = np.abs(trans - t)
            d = np.minimum(d, 1.0 - d)
            if np.allclose(rot, r, atol=tol) and np.all(d < 1e-6):
                return True
        return False

    for r1, t1 in ops:
        for r2, t2 in ops:
            if not known(r2 @ r1, np.mod(t1 @ r2.T + t2, 1.0)):
                return False
    return True


def _equivalent_copies(originals, symops, P, rng):
    """Copy each supercell structure by a random symmetry operation and a random
    unit-cell lattice translation (expressed in the supercell frame P)."""
    Pinv = np.linalg.inv(P)

    def transformed(atoms, r, t, n):
        frac = atoms.get_scaled_positions(wrap=False)
        new = (frac @ P @ r.T + t + n) @ Pinv
        return Atoms(
            symbols=atoms.get_chemical_symbols(),
            scaled_positions=np.mod(new, 1.0),
            cell=atoms.cell,
            pbc=True,
            tags=atoms.get_tags(),
        )

    return [
        transformed(
            a,
            *symops[rng.integers(0, len(symops))],
            rng.integers(0, 3, size=3),
        )
        for a in originals
    ]


def _erohea():
    ds = CifParser("examples/EROHEA_modified.cif").get_disordered_structure()
    return ds, ds.spacegroup.get_symop(), OrderedfromDisordered(ds, quiet=True)


def test_wrap_stays_inside_periodic_box():
    """np.mod(-1e-17, 1.0) is exactly 1.0, which cKDTree(boxsize=1) rejects."""
    x = np.array([[-1e-17, 0.2, 0.3], [1.0, -1.0, 0.5]])
    assert np.mod(x, 1.0).max() == 1.0
    wrapped = _wrap(x)
    assert wrapped.min() >= 0.0 and wrapped.max() < 1.0
    cKDTree(wrapped, boxsize=1.0)


def test_merge_symm_50_images():
    """Verify 50 configs in 1x1x1 reduce to 5 symmetry-inequivalent groups."""
    original_atoms = read("tests/EROHEA_modified.cif")
    images = read("tests/images_50-5_groups.traj", index=":")
    symops = get_spacegroup(original_atoms).get_symop()

    for algo in ["symm", "canonical", "pairwise"]:
        groups = merge_structures(
            images,
            algo=algo,
            symops=symops,
            use_disordered_only=True,
            symprec=1e-4,
            quiet=True,
        )
        assert len(groups) == 5, f"Expected 5 groups with algo={algo}, got {len(groups)}"


def test_canonical_is_strict_and_symm_falls_back(caplog):
    """'canonical' raises on non-discrete input; 'symm' warns and falls back."""
    original_atoms = read("tests/EROHEA_modified.cif")
    images = read("tests/images_50-5_groups.traj", index=":")
    symops = get_spacegroup(original_atoms).get_symop()

    # two disordered atoms on top of each other: not a discrete site set
    bad = images[0].copy()
    idx = np.where(bad.get_tags() != 0)[0]
    pos = bad.get_scaled_positions()
    pos[idx[1]] = pos[idx[0]]
    bad.set_scaled_positions(pos)
    images = images + [bad]

    with pytest.raises(NonDiscreteSitesError):
        merge_structures(images, algo="canonical", symops=symops, quiet=True)

    with caplog.at_level(logging.WARNING, logger="orgdisord.merge"):
        groups = merge_structures(
            images, algo="symm", symops=symops, quiet=True, return_group_indices=True
        )[1]
    assert "falling back" in caplog.text
    pairwise = merge_structures(
        images, algo="pairwise", symops=symops, quiet=True, return_group_indices=True
    )[1]
    assert _partition(groups) == _partition(pairwise)
    assert len(groups) == 6


@pytest.mark.parametrize(
    "supercell, n_ops",
    [
        ([1, 1, 1], 3),
        # 3-fold maps a -> b, which a 2x1x1 cell does not respect: only identity survives
        ([2, 1, 1], 2),
        ([2, 2, 1], 12),
        # sqrt(3) x sqrt(3) cell: 3-fold symmetric, not diagonal
        ([[1, -1, 0], [1, 2, 0], [0, 0, 1]], 9),
    ],
)
def test_supercell_symops_is_closed_group(supercell, n_ops):
    """P3 operations in various supercells: only lattice-preserving operations
    are kept, and together with the lattice translations they form a group."""
    symops = Spacegroup(143).get_symop()
    ops = supercell_symops(symops, supercell)
    assert len(ops) == n_ops
    assert _is_closed_group(ops)


def test_supercell_symops_rejects_bad_input():
    symops = Spacegroup(143).get_symop()
    with pytest.raises(ValueError):
        supercell_symops(symops, [2, 0, 1])
    with pytest.raises(ValueError):
        supercell_symops(symops, [1.5, 1, 1])
    with pytest.raises(ValueError):
        supercell_symops(symops, [[1, 2, 3], [2, 4, 6], [0, 0, 1]])


def test_supercell_synthetic_duplicates():
    """Verify that copies generated by true supercell spacegroup operations
    (including primitive-cell translations) merge into the expected groups.
    """
    _, symops, od = _erohea()
    supercell = [2, 2, 2]
    originals = od.get_supercell_configs(supercell, maxiters=20, random_configs=True)
    n = len(originals)
    copies = _equivalent_copies(
        originals, symops, np.diag(supercell).astype(float), np.random.default_rng(42)
    )
    images = originals + copies

    for algo in ["canonical", "pairwise"]:
        _, group_indices = merge_structures(
            images,
            algo=algo,
            symops=symops,
            supercell=supercell,
            symprec=1e-3,
            quiet=True,
            return_group_indices=True,
        )
        # Every copy k + n must belong to the same group as original k
        group_lookup = {idx: gi for gi, grp in enumerate(group_indices) for idx in grp}
        for k in range(n):
            assert (
                group_lookup[k] == group_lookup[k + n]
            ), f"Copy {k} did not merge with original in algo={algo}"


def test_supercell_exhaustive_orbits():
    """2x1x1 exhaustive enumeration (256 configs) collapses to 26 orbits under
    both algorithms.

    26 was cross-checked once, outside the test suite, in two independent ways:
    - pymatgen's StructureMatcher on the disordered atoms: every member matches
      its group representative and no two representatives match.
    - MACE-OFF24-medium single points on all 256 configs (CPU, float64): energy
      and sorted force magnitudes are constant within every group (spread
      5e-8 eV, noise floor 1.4e-7 eV). Different groups have different energies,
      except one pair (groups 7 and 15) whose force spectra still differ by
      1e-2 eV/A.
    Before the supercell-cell fix in the enumeration (it sheared non-orthogonal
    cells) this gave 70 groups, of which 27 pairs were energy- and
    force-degenerate, i.e. unmerged duplicates.
    """
    _, symops, od = _erohea()
    supercell = [2, 1, 1]
    images = od.get_supercell_configs(supercell, maxiters=5000)
    assert len(images) == 256

    results = {}
    for algo in ["canonical", "pairwise"]:
        results[algo] = merge_structures(
            images,
            algo=algo,
            symops=symops,
            supercell=supercell,
            symprec=1e-3,
            quiet=True,
            return_group_indices=True,
        )[1]

    assert len(results["canonical"]) == 26
    assert _partition(results["canonical"]) == _partition(results["pairwise"])


def test_non_diagonal_supercell():
    """A non-diagonal supercell is merged in its own frame: no operation is lost
    for a lattice-preserving matrix, copies merge, and algorithms agree."""
    _, symops, od = _erohea()
    P = np.array([[1, 1, 0], [-1, 1, 0], [0, 0, 1]])
    images = od.get_supercell_configs(P, maxiters=5000)
    assert len(images) == 256
    assert len(supercell_symops(symops, P)) == 2 * len(symops)

    merged = {
        algo: merge_structures(
            images,
            algo=algo,
            symops=symops,
            supercell=P,
            symprec=1e-3,
            quiet=True,
            return_group_indices=True,
        )[1]
        for algo in ["canonical", "pairwise"]
    }
    assert _partition(merged["canonical"]) == _partition(merged["pairwise"])
    assert len(merged["canonical"]) < len(images)

    originals = images[:40]
    copies = _equivalent_copies(originals, symops, P.astype(float), np.random.default_rng(7))
    groups = merge_structures(
        originals + copies,
        algo="canonical",
        symops=symops,
        supercell=P,
        symprec=1e-3,
        quiet=True,
        return_group_indices=True,
    )[1]
    lookup = {idx: gi for gi, grp in enumerate(groups) for idx in grp}
    assert all(lookup[k] == lookup[k + 40] for k in range(40))


def test_cli_passes_supercell_matrix_to_merge(tmp_path, monkeypatch):
    """--matrix must reach merge_structures, not the default --supercell 1 1 1."""
    seen = {}
    real = cli.merge_structures

    def spy(*args, **kwargs):
        seen["supercell"] = kwargs["supercell"]
        return real(*args, **kwargs)

    monkeypatch.setattr(cli, "merge_structures", spy)
    result = CliRunner().invoke(
        cli.cli,
        [
            "enumerate", "examples/EROHEA_modified.cif",
            "--matrix", "1", "1", "0", "-1", "1", "0", "0", "0", "1",
            "-N", "8", "--random", "--merge",
            "--prefix", str(tmp_path / "run"),
            "--log", str(tmp_path / "orgdisord.log"),
            "--quiet",
        ],
    )  # fmt: skip
    assert result.exit_code == 0, result.output + repr(result.exception)
    assert np.array_equal(seen["supercell"], [[1, 1, 0], [-1, 1, 0], [0, 0, 1]])


def test_merge_ewald():
    """Verify Ewald merging on test trajectory."""
    pytest.importorskip("pymatgen")
    original_atoms = read("tests/EROHEA_modified.cif")
    images = read("tests/images_50-5_groups.traj", index=":")
    symops = get_spacegroup(original_atoms).get_symop()

    groups = merge_structures(
        images,
        oxidation_states={"N": -3.0, "H": 1.0, "C": 1.0, "O": -2.0},
        algo="ewald",
        symops=symops,
        use_disordered_only=True,
        symprec=1e-2,
        quiet=True,
    )
    assert len(groups) == 5
