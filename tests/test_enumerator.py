#!/usr/bin/env python

"""Tests for the enumerate module.

Tests enumeration of ordered configurations from disordered structures:
1. Exhaustive enumeration in different supercell sizes
2. Random configuration generation
3. Generation from major/minor component structures
"""

from pathlib import Path
import shutil
import pytest
from orgdisord.enumerate import OrderedfromDisordered
from orgdisord.parse_cif_file import CifParser
from orgdisord.disordered_structure import from_disorder_components
import numpy as np
from ase.io import read
from ase.spacegroup import Spacegroup

@pytest.fixture
def test_files(tmp_path):
    """Copy test files to a temporary directory."""
    test_dir = tmp_path / "tests"
    # Get the absolute path to the tests directory
    src_dir = Path(__file__).parent
    shutil.copytree(src_dir, test_dir)
    return test_dir

@pytest.fixture
def disordered_structure(test_files):
    """Load test disordered structure."""
    cif = CifParser(test_files / "EROHEA_modified.cif")
    return cif.get_disordered_structure()

@pytest.fixture
def disordered_structure_ababub(test_files):
    """Load test disordered structure."""
    cif = CifParser(test_files / "ABABUB.cif")
    return cif.get_disordered_structure()

@pytest.mark.parametrize("supercell,expected_configs", [
    ([1, 1, 1], 16),  # Unit cell should have 16 configs (2^4)
    ([2, 1, 1], 256)  # 2x1x1 cell should have 256 configs
])
def test_exhaustive_enumeration(disordered_structure, supercell, expected_configs):
    """Test exhaustive enumeration of configurations."""
    od = OrderedfromDisordered(disordered_structure)
    
    images = od.get_supercell_configs(
        supercell=supercell,
        maxiters=500,  # Should be enough for test cases
        exclude_ordered=False,
        random_configs=False
    )
    
    # Check number of configurations
    assert len(images) == expected_configs, f"Wrong number of configurations for {supercell} supercell"
    
    # Check atom counts
    natoms = 288 * np.prod(supercell)  # Total number of atoms
    assert all([len(atoms) == natoms for atoms in images]), "Wrong number of atoms in configurations"
    
    # Check ordered sites (tagged as 0)
    nordered = 252 * np.prod(supercell)  # Number of ordered sites
    assert len(np.where(images[0].get_tags() == 0)[0]) == nordered, "Wrong number of ordered sites"

@pytest.mark.parametrize("n_configs", [10, 34, 100])
def test_random_enumeration(disordered_structure, n_configs):
    """Test random configuration generation."""
    od = OrderedfromDisordered(disordered_structure)
    
    images = od.get_supercell_configs(
        supercell=[1, 1, 1],
        maxiters=n_configs,
        exclude_ordered=False,
        random_configs=True
    )
    
    assert len(images) == n_configs, "Wrong number of random configurations generated"
    
    # Check that configurations are valid (proper atom counts and tags)
    natoms = 288
    nordered = 252
    for atoms in images:
        assert len(atoms) == natoms, "Wrong number of atoms in random configuration"
        n_ordered = len(np.where(atoms.get_tags() == 0)[0])
        assert n_ordered == nordered, "Wrong number of ordered sites in random configuration"

def test_component_enumeration(test_files):
    """Test enumeration from major/minor component structures."""
    atoms_maj = read(test_files / "ABABUB_maj.xyz")
    atoms_min = read(test_files / "ABABUB_min.xyz")
    
    disordered_structure = from_disorder_components(
        atoms_maj, atoms_min,
        tolerance=0.01,
        symprec=1e-3
    )
    
    # Validate disorder pattern
    assert disordered_structure.get_number_of_disorder_groups_per_assembly() == [2]
    assert disordered_structure.Z == 4
    assert disordered_structure.spacegroup == Spacegroup(14)
    
    # Check disorder group sizes
    assembly = disordered_structure.disorder_assemblies[0]
    assert len(assembly.disorder_groups[0].atoms) == 12
    assert len(assembly.disorder_groups[1].atoms) == 12
    
    # Test configuration generation
    od = OrderedfromDisordered(disordered_structure)
    images = od.get_supercell_configs(
        supercell=[1, 1, 1],
        maxiters=100,
        exclude_ordered=False,
        random_configs=False
    )
    
    assert len(images) == 16, "Wrong number of configurations from component structures"

def test_component_enumeration_erohea(test_files):
    """Test enumeration from erohea component structures.
    
    TODO: these tests pass but the generated structures are wrong! 
    """
    atoms_maj = read(test_files / "EROHEA_maj_P1.cif")
    atoms_min = read(test_files / "EROHEA_min_P1.cif")

    # Ensure arrays are set
    if not atoms_maj.has('labels'):
        atoms_maj.set_array('labels', np.array([f"{sym}{i}" for i, sym in enumerate(atoms_maj.get_chemical_symbols())]))
    if not atoms_maj.has('tags'):
        atoms_maj.set_array('tags', np.zeros(len(atoms_maj), dtype=int))


    if not atoms_min.has('labels'):
        atoms_min.set_array('labels', np.array([f"{sym}{i}" for i, sym in enumerate(atoms_min.get_chemical_symbols())]))
    if not atoms_min.has('tags'):
        atoms_min.set_array('tags', np.zeros(len(atoms_min), dtype=int))
    
    disordered_structure = from_disorder_components(
        atoms_maj, atoms_min,
        tolerance=0.01,
        symprec=1e-4
    )
    
    # Validate disorder pattern
    assert disordered_structure.get_number_of_disorder_groups_per_assembly() == [2]
    assert disordered_structure.Z == 4
    # assert disordered_structure.spacegroup == Spacegroup(14)
    
    # Check disorder group sizes
    assembly = disordered_structure.disorder_assemblies[0]
    # assert len(assembly.disorder_groups[0].atoms) == 12
    # assert len(assembly.disorder_groups[1].atoms) == 12
    
    # Test configuration generation
    od = OrderedfromDisordered(disordered_structure)
    images = od.get_supercell_configs(
        supercell=[1, 1, 1],
        maxiters=100,
        exclude_ordered=False,
        random_configs=False
    )
    
    assert len(images) == 16, "Wrong number of configurations from tricky component structures"

def test_split_major_minor(disordered_structure_ababub):
    """Test splitting structure into major/minor components.
    
    TODO: rethink these tests - what should we check for in the case of EROHEA - i.e. special cases?
    
    """
    od = OrderedfromDisordered(disordered_structure_ababub)
    
    # Get fully occupied configurations
    configs = od.split_major_minor(exclude_ordered=False)
    
    # Should get 2 configurations (one per disorder group)
    assert len(configs) == 2, "Wrong number of major/minor configurations"
    
    # Each configuration should have same total number of atoms
    natoms = 120  # Total atoms in structure
    assert all(len(atoms) == natoms for atoms in configs), "Wrong number of atoms in configurations"
    
    # Check disordered site occupancy
    disordered_tags = [tag for tag in configs[0].get_tags() if tag != 0]
    unique_tags = np.unique(disordered_tags)
    
    # Should only have one unique tag value for disordered sites in each config
    assert len(unique_tags) == 1, "Multiple disorder groups present in single configuration"
    
    # Test with exclude_ordered=True
    configs_no_ordered = od.split_major_minor(exclude_ordered=True)
    
    # Verify configs have different tags (different disorder groups)
    tags_1 = np.unique([tag for tag in configs_no_ordered[0].get_tags() if tag != 0])
    tags_2 = np.unique([tag for tag in configs_no_ordered[1].get_tags() if tag != 0])
    assert len(tags_1) == 1, "Multiple disorder groups present in first configuration"
    assert len(tags_2) == 1, "Multiple disorder groups present in second configuration"
    assert tags_1[0] != tags_2[0], "Both configurations have same disorder group"
    
    # Should still get 2 configurations
    assert len(configs_no_ordered) == 2, "Wrong number of major/minor configurations"
    
    # But with fewer atoms (only disordered sites)
    n_disordered = natoms - 72  # 72 is number of ordered atoms
    assert all(len(atoms) == n_disordered for atoms in configs_no_ordered), \
        "Wrong number of atoms in configurations with exclude_ordered=True"
