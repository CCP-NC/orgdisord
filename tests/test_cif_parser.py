#!/usr/bin/env python

"""Tests for the parse_cif_file module."""
import pytest
from sodorg_renewal.parse_cif_file import CifParser
import numpy as np
# suppress warnings for these tests
import warnings
warnings.filterwarnings("ignore")


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



    cif = CifParser('tests/704722.cif')
    assert cif.nassemblies == 1
    assert cif.ndisordergroups == 2
    # # test that the spacegroup is correct
    assert cif.sg.no == 14
    assert cif.nops == 4
    assert cif.Z == 4
    # check number of ordered sites (none in this case!)
    assert len(np.where(cif.ordered_mask)[0]) == 0

    cif = CifParser('tests/VIFQIL01.cif')
    assert cif.nassemblies == 1
    assert cif.ndisordergroups == 2
    # # test that the spacegroup is correct
    assert cif.sg.no == 14
    assert cif.nops == 4
    assert cif.Z == 4
    # check number of ordered sites
    assert len(np.where(cif.ordered_mask)[0]) == 25
    

    cif = CifParser('tests/VAGKUM.cif')
    from ase.visualize import view
    import spglib
    from ase.spacegroup import get_spacegroup
    # print('my final sg is: ',cif.sg, cif.sg.nsymop)
    
    assert cif.nassemblies == 1
    assert cif.ndisordergroups == 2
    # # test that the spacegroup is correct
    assert cif.sg.no == 15
    assert cif.nops == 8
    assert cif.Z == 8
    # check number of ordered sites
    assert len(np.where(cif.ordered_mask)[0]) == 34
    

if __name__ == '__main__':
    test_cifparser()