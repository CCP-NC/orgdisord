.. highlight:: shell


==============
orgdisord
==============


.. image:: https://github.com/CCP-NC/orgdisord/actions/workflows/ci.yml/badge.svg
        :target: https://github.com/CCP-NC/orgdisord/actions/workflows/ci.yml

.. image:: https://img.shields.io/pypi/v/orgdisord.svg
        :target: https://pypi.python.org/pypi/orgdisord

.. image:: https://img.shields.io/pypi/pyversions/orgdisord.svg
        :target: https://pypi.python.org/pypi/orgdisord

.. image:: https://img.shields.io/badge/docs-GitHub%20Pages-blue.svg
        :target: https://CCP-NC.github.io/orgdisord

.. image:: https://img.shields.io/github/license/CCP-NC/orgdisord.svg
        :target: https://github.com/CCP-NC/orgdisord/blob/main/LICENSE





Python code to handle disordered molecular crystals.

.. warning::

    This is a work in progress. The code is not yet ready for general use.



* Free software: MIT license
* Documentation: https://CCP-NC.github.io/orgdisord


Installation
------------

orgdisord requires Python 3.11 or newer. Install the latest release from PyPI with pip::

    $ pip install orgdisord

or with `uv <https://docs.astral.sh/uv/>`_, either into a project environment::

    $ uv add orgdisord

or as a command-line tool in its own isolated environment::

    $ uv tool install orgdisord

Some machine-learning features need the optional ``ML`` extra (``pip install "orgdisord[ML]"``). See the :doc:`installation` page for full details, including installing from source.


Features
--------

* Parse CIF files marked up with disorder.
* Enumerate all possible ordered configurations from a disordered CIF file, given a supercell size.
* Merge symmetry equivalent ordered configurations using a variety of strategies:

   * Using crystal symmetry operations (canonical orbit labelling in the supercell frame, with a pairwise fallback when sites are not discrete)
   * Computing and comparing Ewald energies
   * Computing a SOAP descriptor and comparing using the REMatch kernel

* Visualise the merged configurations using ASE, with different disorder assemblies/group tagged.
* Output the generated (merged) configurations in a number of formats (cif, xyz, poscar, cell).


Limitations
-----------

* Disorder in the CIF file must be marked up using the ``_atom_site_disorder_group`` and ``_atom_site_disorder_assembly`` tags.
* The handling of Z' < 1 is not currently robust and will not work in many cases. In addition, the generation of supercells in Z' < 1 cases currently generates structures with overlapping atoms.



Quick start examples
---------------------

The examples below assume you have installed the package as per the :doc:`installation` instructions and have the example CIF file in the current directory. You can download the example CIF file here:
:download:`ABABUB.cif <../examples/ABABUB.cif>`

#. Generate all possible ordered configurations in the primitive cell: ::

        $ orgdisord enumerate ABABUB.cif --no_write --view

#. Generate all possible ordered configurations in a 2x1x1 supercell: ::

        $ orgdisord enumerate ABABUB.cif --supercell 2 1 1 --no_write --view

#. Merge fully enumerated configurations: ::

        $ orgdisord enumerate ABABUB.cif -m --no_write --view

#. Merge randomly generated configurations in the primitive cell. If you have enough of them, you will end up with the same groups as the fully enumerated case (though this is much less efficient): ::

        $ orgdisord enumerate ABABUB.cif -m --no_write --view --random -N 1000

#. Generate 5000 randomly ordered structures in a specified supercell: ::

        $ orgdisord enumerate ABABUB.cif --supercell 2 1 2 --random -N 5000 --no_write --view

#. Generate 3 randomly ordered structures in a very large supercell: ::

        $ orgdisord enumerate ABABUB.cif --supercell 4 6 4 --random -N 3 --no_write --view

#. Generate all possible ordered configurations in a non-diagonal supercell given by a 3x3 integer matrix (row-major, the supercell cell is the matrix times the primitive cell): ::

        $ orgdisord enumerate ABABUB.cif --matrix 1 0 -1 0 2 0 1 0 1 --no_write --view

#. Generate all possible ordered configurations in the supercell of 4 primitive cells that is closest to a simple cubic shape: ::

        $ orgdisord enumerate ABABUB.cif --target-size 4 --target-shape sc --no_write --view

The supercell is set with one of ``--supercell``, ``--matrix`` or ``--target-size``, and only one of them may be given. ``--target-size N`` searches for the supercell of N primitive cells whose shape is closest to ``--target-shape`` (``sc``, the default, or ``fcc``). ``--score-key`` sets how candidates are ranked. ``length`` (the default) compares only the three cell-vector lengths with the ideal length for the target shape and ignores angles. ``metric`` also compares the cell angles, using the full metric tensor.

Passing ``-m`` (``--merge``) groups symmetry-equivalent configurations in the supercell frame, including primitive-cell lattice translations. The default algorithm (``--algo symm``) labels each configuration on a shared discrete site template and falls back to pairwise coordinate matching if the sites are not discrete. ``--algo canonical`` disables the fallback, and ``--algo pairwise`` always compares pairs. By default only the disordered atoms are compared; pass ``--use_all_atoms`` to include the ordered framework.


Notice that a `orgdisord.log` file is generated in the current directory. This file contains the full command line options used to generate the output, as well as information about how the CIF file was parsed, the enumeration was done and, if applicable, how the merging was done. This can be useful for debugging and reproducing results.


The :doc:`command line interface documentation <cli>` has full details of the available commands and options. Alternatively you can run ``orgdisord --help`` to see the available commands and e.g. ``orgdisord enumerate --help`` to see the available options for the ``enumerate`` command.


.. note::

        For each of these examples we included the ``--no_write`` option to prevent the configurations from being written to file. Without this flag, the code will write the configurations, in extended xyz format, into a directory named `orgdisord-results`, and it will also generate a .csv file containing a summary of the generated structures.

        We also included the ``--view`` option to visualise the merged configurations using ASE. In the ASE GUI, you can click View→Colors and select color By Tag to easily see the different disorder assemblies and groups. You can also click View→Show Tags to see the tags associated with each atom.


For more examples and case studies, see the :doc:`case_studies` page.
