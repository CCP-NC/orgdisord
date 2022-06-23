==============
SODORG Renewal
==============


.. image:: https://img.shields.io/pypi/v/sodorg_renewal.svg
        :target: https://pypi.python.org/pypi/sodorg_renewal

.. image:: https://img.shields.io/travis/jkshenton/sodorg_renewal.svg
        :target: https://travis-ci.com/jkshenton/sodorg_renewal

.. image:: https://readthedocs.org/projects/sodorg-renewal/badge/?version=latest
        :target: https://sodorg-renewal.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status




Python code to handle disordered molecular crystals.


* Free software: MIT license
* Documentation: https://sodorg-renewal.readthedocs.io.


Features
--------


Example use cases
-----------------
* Generate 10 000 randomly ordered structures in a specified supercell
`python sodorg_renewal/cli.py tests/EROHEA_modified.cif -v --view --supercell 2 1 2 --random --maxiters 10000`

* EROHEA merged using structure comparison
`python sodorg_renewal/cli.py tests/EROHEA_modified.cif -v -m --view --algo symm`

* EROHEA merged using local descriptors
`python sodorg_renewal/cli.py tests/EROHEA_modified.cif -v -m --view --algo rematch --symprec 1e-2`

* EROHEA merged using Electrostatic energies
`python sodorg_renewal/cli.py tests/EROHEA_modified.cif -v -m --view --algo ewald --ox C 1 --ox H 1 --ox N -3 --ox O -2`

* VAGKUM merge, only look at disordered bit (makes it easier to see what is going on sometimes.)
`python sodorg_renewal/cli.py tests/VAGKUM.cif -v --view -m --exclude_ordered`




Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
