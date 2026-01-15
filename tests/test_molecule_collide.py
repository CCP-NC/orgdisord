#!/usr/bin/env python

"""Tests for molecule collision detection functionality.

Tests the molecule_collide function which detects overlapping/colliding
molecular fragments in crystal structures.
"""

import pytest
from ase import Atoms, Atom
from orgdisord.utils import molecule_collide
import numpy as np

@pytest.mark.skip(reason="Molecule collide function under revision")
def test_molecule_collide():
    """Test the molecule_collide function."""
    h2o = Atoms("O", positions=[[0, 0, 0]])
    h2 = Atoms("H2", positions=[[0, 0, 0.1], [0, 0, -0.1]])
    assert molecule_collide(h2o, h2) is True
    assert molecule_collide(h2o, h2, tolerance=0.05) is False

@pytest.mark.skip(reason="Molecule collide function under revision") 
def test_molecule_collide_tolerance():
    """Test the tolerance parameter of molecule_collide function."""
    h2o = Atoms("O", positions=[[0, 0, 0]])
    h2 = Atoms("H2", positions=[[0, 0, 0.1], [0, 0, -0.1]])
    assert molecule_collide(h2o, h2, tolerance=0.15) is True
    assert molecule_collide(h2o, h2, tolerance=0.05) is False
