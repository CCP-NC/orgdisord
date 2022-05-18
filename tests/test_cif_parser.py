#!/usr/bin/env python

"""Tests for the parse_cif_file module."""
import pytest
from sodorg_renewal.parse_cif_file import CifParser
import numpy as np


def test_cifparser():
    cif = CifParser('tests/EROHEA_modified.cif')
    assert cif.nassemblies == 1
    group1 = [26, 31, 32, 33, 34]
    assert all([a==b for a,b in zip(group1, cif.disorder_groups[0][0])])
    assert cif.ndisordergroups == 2
    # # test that the spacegroup is correct
    assert cif.sg.no == 15
    assert cif.nops == 8
    assert cif.Z == 4
    # check number of ordered sites
    assert len(np.where(cif.ordered_mask)[0]) == 33

if __name__ == '__main__':
    test_cifparser()