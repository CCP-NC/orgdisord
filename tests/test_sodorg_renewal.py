#!/usr/bin/env python

"""Tests for `orgdisord` package."""

import contextlib

import numpy as np
import pytest
from click.testing import CliRunner

# from orgdisord import orgdisord
from orgdisord import cli


@pytest.fixture
def erohea():
    """Read in the EROHEA_modified.cif file.

    return a CifParser object

    """
    import orgdisord.parse_cif_file as CifParser

    return CifParser.CifParser("tests/EROHEA_modified.cif")


@pytest.fixture
def test_enumerator(cif):
    """Use the pytest fixture cif to test the enumerator module."""
    from orgdisord.enumerate import OrderedfromDisordered

    od = OrderedfromDisordered(cif)
    # exhaustive enumeration in 111 supercell:
    supercell = [1, 1, 1]
    images = od.get_supercell_configs(
        supercell=supercell,
        maxiters=100,  # should be more than enough
        include_ordered=True,
        random_configs=False,
    )
    # there should be 16 configurations generated
    # (binary options and four possible sites, 2^4 = 16)
    assert len(images) == 16
    # check that the number of sites is correct (including ordered sites)
    natoms = 288 * supercell[0] * supercell[1] * supercell[2]
    assert all([len(atoms) == natoms for atoms in images])
    # check that the number of ordered sites is correct
    # (ordered sites are tagged 0)
    tag = 0
    nordered = 252 * supercell[0] * supercell[1] * supercell[2]
    assert len(np.where(images[0].get_tags() == tag)[0]) == nordered

    return images


def test_command_line_interface(tmp_path):
    """Test the CLI."""
    runner = CliRunner()
    with contextlib.chdir(tmp_path):
        help_result = runner.invoke(cli.cli, ["--help"])
        assert help_result.exit_code == 0
        assert "enumerate" in help_result.output
        assert "analyse" in help_result.output
        assert "--help  Show this message and exit." in help_result.output

        enum_help = runner.invoke(cli.cli, ["enumerate", "--help"])
        assert enum_help.exit_code == 0
        for option in ["--supercell", "--matrix", "--target-size"]:
            assert option in enum_help.output

        result = runner.invoke(cli.cli, ["enumerate"])
        assert result.exit_code == 2
        assert "Expected 1 or 2 CIF files" in result.output
