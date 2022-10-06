.. highlight:: shell


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

.. warning::

    This is a work in progress. The code is not yet ready for general use.



* Free software: MIT license
* Documentation: https://sodorg-renewal.readthedocs.io.


Features
--------

* Parse CIF files marked up with disorder.
* Enumerate all possible ordered configurations from a disordered CIF file, given a supercell size.
* Merge symmetry equivalent ordered configurations using a variety of strategies: 

   * Using symmetry operations defined in the CIF file
   * Computing and comparing Ewald energies
   * Computing a SOAP descriptor and comparing using the REMatch kernel

* Visualise the merged configurations using ASE, with different disorder assemblies/group tagged.
* Output the generated (merged) configurations in a number of formats (cif, xyz, poscar, cell).


Limitations
-----------

* Disorder in the CIF file must be marked up using the ``_atom_site_disorder_group`` and ``_atom_site_disorder_assembly`` tags.
* The handling of Z' < 1 is not currently robust and will not work in many cases. 
In addition, the generation of supercells in Z' < 1 cases currently generates structures with overlapping atoms. 



Example use cases
-----------------

Note you have to first install the package using ``pip install``. See the :doc:`installation` instructions for full details.

See the :doc:`command line interface documentation <cli>` for full details of the available commands and options. Alternatively you can run ``sodorg_renewal --help`` to see the available commands and e.g. ``sodorg_renewal enumerate --help`` to see the available options for the ``enumerate`` command.

To run the specific examples below, you can download the example files (you may have to right-click and save-as):

* :download:`EROHEA <../examples/EROHEA_modified.cif>`
* :download:`VAGKUM <../examples/VAGKUM.cif>`

The examples below assume you have the example CIF files in the current directory.


* Generate all possible ordered configurations in the primitive cell: ::

        $ sodorg_renewal EROHEA_modified.cif --no_write --view

* Generate all possible ordered configurations in a 2x1x1 supercell: ::

        $ sodorg_renewal EROHEA_modified.cif --supercell 2 1 1 --no_write --view

* Merge fully enumerated configurations using *structure comparison*: ::

        $ sodorg_renewal enumerate EROHEA_modified.cif -m --algo symm --no_write --view

* Merge fully enumerated configurations using *local descriptors*: ::

        $ sodorg_renewal enumerate EROHEA_modified.cif -m --algo rematch --symprec 1e-2 --no_write --view

* Merge fully enumerated configurations using *electrostatic energies*: ::
        
        $ sodorg_renewal enumerate EROHEA_modified.cif -m --algo ewald --ox C 1 --ox H 1 --ox N -3 --ox O -2 --no_write --view

* Generate 10 000 randomly ordered structures in a specified supercell: ::

        $ sodorg_renewal enumerate EROHEA_modified.cif --supercell 2 1 2 --random --maxiters 10000 --no_write --view

Notice that a `sodorg.log` file is generated in the current directory. This file contains the full command line options used to generate the output, as well as information about how the CIF file was parsed, the enumeration was done and, if applicable, how the merging was done. This can be useful for debugging and reproducing results.

.. note::

        For each of these examples we included the ``--no_write`` option to prevent the configurations from being written to file. Without this flag, the code will write the configurations, in extended xyz format, into a directory named `sodorg-results`, and it will also generate a .csv file containing a summary of the generated structures. 
        
        We also included the ``--view`` option to visualise the merged configurations using ASE. In the ASE GUI, you can click View→Colors and select color By Tag to easily see the different disorder assemblies and groups. You can also click View→Show Tags to see the tags associated with each atom.



Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
