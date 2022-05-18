#!/usr/bin/env python

"""Tests for `sodorg_renewal` package."""

import pytest

from click.testing import CliRunner

from sodorg_renewal import sodorg_renewal
from sodorg_renewal import cli
import numpy as np


@pytest.fixture
def erohea():
    """Read in the EROHEA_modified.cif file.

    return a CifParser object
    
    """
    import sodorg_renewal.parse_cif_file as CifParser
    return CifParser.CifParser('tests/EROHEA_modified.cif')

@pytest.fixture
def test_enumerator(cif):
    """Use the pytest fixture cif to test the enumerator module."""
    from sodorg_renewal.enumerate import OrderedfromDisordered
    od = OrderedfromDisordered(cif)
    # exhaustive enumeration in 111 supercell:
    supercell = [1,1,1]
    images = od.get_supercell_configs(
                    supercell = supercell,
                    maxiters = 100, # should be more than enough
                    include_ordered = True,
                    random_configs=False)
    # there should be 16 configurations generated
    #(binary options and four possible sites, 2^4 = 16)
    assert len(images) == 16
    # check that the number of sites is correct (including ordered sites)
    natoms = 288*supercell[0]*supercell[1]*supercell[2]
    assert all([len(atoms) == natoms for atoms in images])
    # check that the number of ordered sites is correct
    # (ordered sites are tagged 0)
    tag = 0
    nordered = 252*supercell[0]*supercell[1]*supercell[2]
    assert len(np.where(images[0].get_tags() == tag)[0]) == nordered

    return images




def test_command_line_interface():
    """Test the CLI."""
    runner = CliRunner()
    result = runner.invoke(cli.main)
    assert result.exit_code == 0
    assert 'sodorg_renewal.cli.main' in result.output
    help_result = runner.invoke(cli.main, ['--help'])
    assert help_result.exit_code == 0
    assert '--help  Show this message and exit.' in help_result.output
