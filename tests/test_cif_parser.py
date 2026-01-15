#!/usr/bin/env python

"""Tests for the parse_cif_file module.

Tests parsing of CIF files with different disorder patterns and validates:
- Assembly/group parsing
- Space group determination 
- Site occupancy handling
- Ordered/disordered site identification
"""
import pytest
from orgdisord.parse_cif_file import CifParser
import numpy as np
from pathlib import Path
import shutil

# Suppress warnings for these tests since we're testing error handling
import warnings
warnings.filterwarnings("ignore")

@pytest.fixture
def test_files(tmp_path):
    """Copy test files to a temporary directory."""
    test_dir = tmp_path / "tests"
    # Get the absolute path to the tests directory
    src_dir = Path(__file__).parent
    shutil.copytree(src_dir, test_dir)
    return test_dir

@pytest.mark.parametrize("cif_name,expected", [
    ("EROHEA_modified.cif", {
        "assemblies": 1,
        "spacegroup": 15,
        "nops": 8,
        "Z": 4,
        "ordered_atoms": 33,
        "disorder_groups": [5, 4]  # 5 sites in group -2, 4 sites in group -1
    }),
    ("704722.cif", {
        "assemblies": 1,
        "spacegroup": 14,
        "nops": 4,
        "Z": 4,
        "ordered_atoms": 0,
        "disorder_groups": [14, 14]  # Two groups with 14 sites each
    }),
    ("VAGKUM.cif", {
        "assemblies": 1,
        "spacegroup": 15,
        "nops": 8,
        "Z": 8,
        "ordered_atoms": 37,  # Including partially occupied water
        "disorder_groups": [6, 6]  # Two SO2NH2 groups with 6 sites each
    }),
    ("VIFQIL01.cif", {
        "assemblies": 1,
        "spacegroup": 14,
        "nops": 4,
        "Z": 4,
        "ordered_atoms": 25,
        "disorder_groups": [3, 3]  # Two CF3 groups with 3 sites each
    })
])
def test_cif_parser(test_files, cif_name, expected):
    """Test parsing of CIF files with different disorder patterns."""
    cif_file = test_files / cif_name
    cif = CifParser(cif_file)
    struct = cif.get_disordered_structure()

    # Check basic structure properties
    assert len(struct.disorder_assemblies) == expected["assemblies"]
    assert struct.spacegroup.no == expected["spacegroup"]
    assert len(struct.spacegroup.get_symop()) == expected["nops"]
    assert struct.Z == expected["Z"]
    
    # Check disorder pattern
    assembly = struct.disorder_assemblies[0]
    assert sum(cif.ordered_mask) == expected["ordered_atoms"]  # Check ordered sites in asymmetric unit
    assert [len(g.atoms) for g in assembly.disorder_groups] == expected["disorder_groups"]

@pytest.mark.parametrize("cif_name", [
    "EROHEA_modified.cif",
    "704722.cif",
    "VAGKUM.cif",
    "VIFQIL01.cif"
])
def test_disorder_assemblies(test_files, cif_name):
    """Test disorder assembly identification."""
    cif_file = test_files / cif_name
    cif = CifParser(cif_file)
    struct = cif.get_disordered_structure()
    
    # Check assembly and group properties
    assembly = struct.disorder_assemblies[0]
    assert len(assembly.disorder_groups) > 0
    assert all(len(g.atoms) > 0 for g in assembly.disorder_groups)
    assert all(g.occupancy > 0 for g in assembly.disorder_groups)
    assert np.isclose(sum(g.occupancy for g in assembly.disorder_groups), 1.0)

def test_invalid_cif(test_files):
    """Test handling of invalid CIF files."""
    with pytest.raises(FileNotFoundError):
        CifParser(test_files / "nonexistent.cif")
