#!/usr/bin/env python

"""Tests for the merge.py module.

Tests the different structure merging algorithms:
1. Symmetry-based merging
2. REMatch kernel merging 
3. Ewald energy-based merging

The test uses a set of 50 configurations (34 random + 16 exhaustive)
that should merge into 5 symmetry-inequivalent groups.
"""
import pytest
from orgdisord.merge import merge_structures
from ase.spacegroup.symmetrize import check_symmetry
import numpy as np
from ase.io import read
import spglib
from pathlib import Path
import shutil

try:
    from dscribe.descriptors import SOAP
    from dscribe.kernels import REMatchKernel
    HAS_DSCRIBE = True
except ImportError:
    HAS_DSCRIBE = False

@pytest.fixture
def test_files(tmp_path):
    """Copy test files to a temporary directory."""
    test_dir = tmp_path / "tests"
    # Get the absolute path to the tests directory
    src_dir = Path(__file__).parent
    shutil.copytree(src_dir, test_dir)
    return test_dir

@pytest.fixture
def test_structures(test_files):
    """Load test structures and their symmetry operations."""
    original_atoms = read(test_files / "EROHEA_modified.cif")
    images = read(test_files / "images_50-5_groups.traj", index=":")
    dataset = check_symmetry(original_atoms)
    symops = [(rot, trans) for rot, trans in zip(dataset['rotations'], dataset['translations'])]
    return original_atoms, images, symops

def test_symmetry_merging(test_structures):
    """Test merging structures using symmetry operations."""
    original_atoms, images, symops = test_structures
    
    # Run symmetry merging
    groups = merge_structures(
        images,
        algo="symm",
        symops=symops,
        use_disordered_only=True,
        symprec=1e-4,
        quiet=True
    )
    
    # Should find 5 unique configurations
    assert len(groups) == 5, "Wrong number of unique configurations"

@pytest.mark.skipif(not HAS_DSCRIBE, reason="dscribe package not available")
def test_rematch_merging(test_structures):
    """Test merging structures using REMatch kernel."""
    original_atoms, images, symops = test_structures
    
    # Run REMatch merging
    unique_indices = merge_structures(
        images,
        algo="rematch",
        symprec=0.3,
        quiet=True
    )
    
    # Should find approximately 5 unique configurations
    assert 4 <= len(unique_indices) <= 6, "Wrong number of unique configurations"

def test_ewald_merging(test_structures):
    """Test merging structures using Ewald energy."""
    original_atoms, images, symops = test_structures
    
    oxidation_states = {
        'C': 0,
        'H': 1,
        'N': -3,
        'O': -2
    }
    
    # Run Ewald merging
    unique_indices = merge_structures(
        images,
        algo="ewald",
        oxidation_states=oxidation_states,
        symprec=0.3,
        quiet=True
    )
    
    # Should find approximately 5 unique configurations
    assert 4 <= len(unique_indices) <= 6, "Wrong number of unique configurations"
