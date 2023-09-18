#!/usr/bin/env python

"""Tests for the parse_cif_file module."""
import pytest
from orgdisord.parse_cif_file import CifParser
import numpy as np
# suppress warnings for these tests
import warnings
warnings.filterwarnings("ignore")


def test_cifparser():
    cif = CifParser('tests/EROHEA_modified.cif')
    assert cif.nassemblies == 1
    group1 = [26, 31, 32, 33, 34]
    disorder_groups = cif.disorder_groups
    assert len(disorder_groups) == 1
    assert all([a==b for a,b in zip(group1, disorder_groups['A'][-2])])
    # # test that the spacegroup is correct
    assert cif.sg.no == 15
    assert cif.nops == 8
    assert cif.Z == 4
    # check number of ordered sites
    assert len(np.where(cif.ordered_mask)[0]) == 33



    cif = CifParser('tests/704722.cif')
    disorder_groups = cif.disorder_groups
    ngroupsperass = [len(x) for x in disorder_groups.values()]
    assert cif.nassemblies == 1
    assert len(ngroupsperass) == 1
    assert ngroupsperass[0] == 2
    # # test that the spacegroup is correct
    assert cif.sg.no == 14
    assert cif.nops == 4
    assert cif.Z == 4
    # check number of ordered sites (none in this case!)
    assert len(np.where(cif.ordered_mask)[0]) == 0

    cif = CifParser('tests/VIFQIL01.cif')
    disorder_groups = cif.disorder_groups
    ngroupsperass = [len(x) for x in disorder_groups.values()]
    assert cif.nassemblies == 1
    assert len(ngroupsperass) == 1
    assert ngroupsperass[0] == 2
    # # test that the spacegroup is correct
    assert cif.sg.no == 14
    assert cif.nops == 4
    assert cif.Z == 4
    # check number of ordered sites
    assert len(np.where(cif.ordered_mask)[0]) == 25
    

    cif = CifParser('tests/VAGKUM.cif')
    disorder_groups = cif.disorder_groups
    ngroupsperass = [len(x) for x in disorder_groups.values()]
    assert cif.nassemblies == 1
    assert ngroupsperass[0] == 2
    # # test that the spacegroup is correct
    assert cif.sg.no == 15
    assert cif.nops == 8
    assert cif.Z == 8
    # check number of ordered sites
    # TODO: what is the correct number of ordered sites here?
    ## The H2O molecules are not treated correctly at present, I think...
    # assert len(np.where(cif.ordered_mask)[0]) == 34


    # cif = CifParser('tests/MICKEP.cif')
    # assert cif.nassemblies == 7
    # # assert cif.ndisordergroups == 2
    # # # test that the spacegroup is correct
    # assert cif.sg.no == 14
    # assert cif.nops == 4
    # assert cif.Z == 4
    # # check number of ordered sites
    # # (note that this is based on the occupancy of the sites, not the
    # #  assembly/disorder group labels!)
    # assert len(np.where(cif.ordered_mask)[0]) == 22


if __name__ == '__main__':
    test_cifparser()