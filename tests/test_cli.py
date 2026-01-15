#!/usr/bin/env python

"""Tests for the sodorg renewal CLI functionality."""

import pytest
from typer.testing import CliRunner
from orgdisord.cli import app
import shutil
import os
from pathlib import Path

@pytest.fixture
def test_files(tmp_path):
    """Copy test files to a temporary directory."""
    test_dir = tmp_path / "tests"
    # Get the absolute path to the tests directory
    src_dir = Path(__file__).parent
    shutil.copytree(src_dir, test_dir)
    return test_dir

@pytest.mark.parametrize("options", [
    [],  # Basic enumeration
    ["-s", "2", "1", "1"],  # With supercell
    ["--maxiters", "10"],  # Limited iterations
    ["--random"],  # Random configurations
    ["--merge"],  # With merging
])
def test_enumerate_configs_command(test_files, options):
    """Test enumeration command with various options."""
    runner = CliRunner()
    try:
        os.chdir(test_files)  # Change to test directory
        cmd = ["enumerate", str("EROHEA_modified.cif"), "--quiet"] + options
        result = runner.invoke(app, cmd, catch_exceptions=False)
        assert result.exit_code == 0, f"Command failed with output: {result.stdout}"
    finally:
        os.chdir(test_files)  # Change back to tests directory

def test_enumerate_configs_two_inputs(test_files):
    """Test enumeration from two input structures."""
    runner = CliRunner()
    try:
        os.chdir(test_files)  # Change to temp directory
        cmd = ["enumerate",
              str("EROHEA_maj_P1.cif"),
              str("EROHEA_min_P1.cif"),
              "--quiet"]
        result = runner.invoke(app, cmd, catch_exceptions=False)
        print(result.stdout)
        assert result.exit_code == 0, f"Command failed with output: {result.stdout}"
    finally:
        os.chdir(test_files)
