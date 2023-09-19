================
Case Studies
================


Here are some case studies of specific systems that have been explored using orgdisord.

The names of the structures refer to the CSD code/CIF file in the ./examples/ directory.


Straightforward case (ABABUB)
--------------------------------


This is a relatively straightforward system consisting of one disorder assembly and two disorder groups. 

* :download:`ABABUB CIF file <../examples/ABABUB.cif>`

When we run the code with (mostly) default options: ::

    orgdisord enumerate ABABUB.cif --prefix ababub

Parsing the CIF file
^^^^^^^^^^^^^^^^^^^^

We get the following output:

.. code-block:: console
    
    Parsing disorder in cif file: ABABUB.cif
    Warning: Site 0  (O1)   in assembly A has no disorder group. This site will be considered *ordered*. 
    Warning: Site 1  (O2)   in assembly A has no disorder group. This site will be considered *ordered*. 
    Warning: Site 2  (N1)   in assembly A has no disorder group. This site will be considered *ordered*. 
    Warning: Site 29 (C6)   in assembly A has no disorder group. This site will be considered *ordered*. 
    Warning: Site 32 (H8A)  in assembly A has no disorder group. This site will be considered *ordered*. 
    Warning: Site 34 (C9)   in assembly A has no disorder group. This site will be considered *ordered*. 
    Warning: Site 38 (H10A) in assembly A has no disorder group. This site will be considered *ordered*. 
    Warning: Site 40 (C11)  in assembly A has no disorder group. This site will be considered *ordered*. 
    
    Disordered structure:
    Disorder assembly: A
    Contains the following groups:
    Disorder group: 1 contains 4 symm. ops. and 12 sites:
            label:       C2, species:    C, occupancy:  0.59
            label:      H2A, species:    H, occupancy:  0.59
            label:      H2B, species:    H, occupancy:  0.59
            label:       C3, species:    C, occupancy:  0.59
            label:      H3A, species:    H, occupancy:  0.59
            label:      H3B, species:    H, occupancy:  0.59
            label:       C4, species:    C, occupancy:  0.59
            label:      H4A, species:    H, occupancy:  0.59
            label:      H4B, species:    H, occupancy:  0.59
            label:       C5, species:    C, occupancy:  0.59
            label:      H5A, species:    H, occupancy:  0.59
            label:      H5B, species:    H, occupancy:  0.59

    Disorder group: 2 contains 4 symm. ops. and 12 sites:
            label:      C2', species:    C, occupancy:  0.41
            label:     H2'1, species:    H, occupancy:  0.41
            label:     H2'2, species:    H, occupancy:  0.41
            label:      C3', species:    C, occupancy:  0.41
            label:     H3'1, species:    H, occupancy:  0.41
            label:     H3'2, species:    H, occupancy:  0.41
            label:      C4', species:    C, occupancy:  0.41
            label:     H4'1, species:    H, occupancy:  0.41
            label:     H4'2, species:    H, occupancy:  0.41
            label:      C5', species:    C, occupancy:  0.41
            label:     H5'1, species:    H, occupancy:  0.41
            label:     H5'2, species:    H, occupancy:  0.41


As we can see, the code warns about certain sites that are part of disorder assembly A yet have no disorder group (they also have an occupancy of 1).
The sites are therefore considered ordered by the code.


Enumerating ordered structures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The code correctly generates a total of 16 (2\ :sup:`4` ) ordered structures and outputs the following table:

.. csv-table:: ABABUB enumerated configurations
   :file: ../examples/ababub.csv
   :header-rows: 1


Where columns for the (free) energies have been added in anticipation of next steps. 
Once these energies have been and added to the table, the code can then be used to analyse the results.
See the :doc:`command line interface documentation <cli>` for the analysis options.

The code will also make a directory (<prefix>-results) in which the structure files are written.

.. tip::
    By default the code will output the structures in the `extended .xyz format <https://github.com/libAtoms/extxyz>`_.
    This has the advantage of being able to store artibrary additional information site-specific and structure-specific information.
    For example, tags specfying the disorder group/assembly and the original CIF labels are preserved. 

    For CIF and CASTEP .cell files, the code will also preserve the CIF-labels, but these may be harder to post-process using `ASE <https://wiki.fysik.dtu.dk/ase/>`_ should you wish to do so.

    You can specify the output format using the ``-f`` option. See the :doc:`command line interface documentation <cli>` for more details.


Merging equivalent structures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If we instead re-run the code with the ``-m`` flag, the code will merge symmetry-equivalent structures:

.. code-block:: console

    Enumerating ordered configurations.
    Generating 16 out of the 16 possible configurations in the (1, 1, 1) supercell:
    100%|██████████████████████████████| 16/16 [00:00<00:00, 3912.83it/s]
    Merging structures...
    Checking symmetry-equivalence: 100%|█| 16/16 [00:00<00:00, 1967.25it/
    Merging took     0.01 s and found 7 groups
    Spacegroup: P2_1/c (14)          multiplicity: 1
    Spacegroup: P1 (1)               multiplicity: 4
    Spacegroup: P2_1 (4)             multiplicity: 2
    Spacegroup: P-1 (2)              multiplicity: 2
    Spacegroup: Pc (7)               multiplicity: 2
    Spacegroup: P1 (1)               multiplicity: 4
    Spacegroup: P2_1/c (14)          multiplicity: 1

and the corresponding table is:

.. csv-table:: ABABUB enumerated configurations (merged)
   :file: ../examples/ababub_merged.csv
   :header-rows: 1

Supercells and random structures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We would typically want to then generate ordered structures for a larger supercell.
Because of the combinatorial explosion, caution is advised when doing this and a conservative number of maximum structures is generated: 512 by default.
To generate more than this, the ``-N`` flag can be used. For example, to generate 1024 structures:

.. code-block:: console

    orgdisord enumerate ABABUB.cif --prefix ababub -N 1024 --supercell 2 2 1

This will enumerate, in order, all possible configurations in a 2x2x1 supercell up to a maximum of 1024 structures.


Sometimes, however, we may want to generate a sample of random configurations for a given supercell size.
This can be done with the ``-r`` flag:

.. code-block:: console

    orgdisord enumerate ABABUB.cif --prefix ababub -r -N 1024 --supercell 2 2 1

.. tip:: 
    To visualise the generated structures, you can add the ``--view`` flag. 
    This will open the structures in `ASE's GUI viewer <https://wiki.fysik.dtu.dk/ase/ase/gui/gui.html>`_.
    
    The ordered sites are tagged with 0 and the disordered sites are tagged according to the disorder group/assembly.
    It's useful to color-by tag to see the disorder groups/assembly in the GUI viewer. 
    You can do this by clicking on ``View->Colors`` button and selecting ``By tag``.
    You can then choose the range of tags you're interested in and also change the color map (cmap).
    
Constrained ratios of disorder components
-----------------------------------------

.. note::
    This feature is currently limited to simple cases where only 1 disorder assembly and 2 disorder groups are present.

The code can also be used to generate structures with a specified ratio of disorder components.
If the ``--fix_ratio`` flag is present, the code will only generate structures with a fixed ratio of disorder components.
By default this ratio comes from the occupancies specified in the CIF file, but they can be overridden using the ``--ratios`` flag.
For example, to generate only structures with 0.75 disorder group 1 and 0.25 disorder group 2:


.. code-block:: console

    orgdisord enumerate ABABUB.cif --prefix ababub --fix_ratio --ratios 0.75 0.25


.. tip::
    The ``--ratio-tol`` flag is useful when the occupancies in the CIF file are not exact.
    For example, if the occupancies are 0.75 and 0.25, the code will generate structures with a ratio of 0.75 and 0.25, but also structures with a ratio of 0.76 and 0.24.
    This is because the code will generate structures with occupancies that are within the tolerance of the specified ratio.
    The default tolerance is 0.01, but this can be changed using the ``--ratio-tol`` flag.





Another straightforward case (AXURIX)
-------------------------------------

This is another relatively straightforward system consisting of one disorder assembly and two disorder groups. 
The main difference is that here Z=8 and so the code generates 256 (2\ :sup:`8` ) ordered structures.
This takes significantly longer than the previous example, but the code still completes in a reasonable time. 
The main additional cost is in reloading this configurations as molecular crystals (i.e. making sure molecular units are connected together in a reasonable way.). 
You can disable this check using the ``--not_molecular_crystal`` flag, though the results of this are not well-tested!

* :download:`AXURIX CIF file <../examples/AXURIX.cif>`

We can run the code and merge the structures as before: ::
        orgdisord enumerate AXURIX.cif --prefix axurix_merged -m

The CIF file is parsed as:

.. code-block:: console

    Parsing disorder in cif file: AXURIX.cif
    Disordered structure:
    Disorder assembly: A
    Contains the following groups:
    Disorder group: 1 contains 8 symm. ops. and 12 sites:
            label:     C28A, species:    C, occupancy:  0.75
            label:     H28A, species:    H, occupancy:  0.75
            label:     H28B, species:    H, occupancy:  0.75
            label:     H28C, species:    H, occupancy:  0.75
            label:     C29A, species:    C, occupancy:  0.75
            label:     H29A, species:    H, occupancy:  0.75
            label:     H29B, species:    H, occupancy:  0.75
            label:     H29C, species:    H, occupancy:  0.75
            label:     C30A, species:    C, occupancy:  0.75
            label:     H30A, species:    H, occupancy:  0.75
            label:     H30B, species:    H, occupancy:  0.75
            label:     H30C, species:    H, occupancy:  0.75

    Disorder group: 2 contains 8 symm. ops. and 12 sites:
            label:     C28B, species:    C, occupancy:  0.25
            label:     H28D, species:    H, occupancy:  0.25
            label:     H28E, species:    H, occupancy:  0.25
            label:     H28F, species:    H, occupancy:  0.25
            label:     C29B, species:    C, occupancy:  0.25
            label:     H29D, species:    H, occupancy:  0.25
            label:     H29E, species:    H, occupancy:  0.25
            label:     H29F, species:    H, occupancy:  0.25
            label:     C30B, species:    C, occupancy:  0.25
            label:     H30D, species:    H, occupancy:  0.25
            label:     H30E, species:    H, occupancy:  0.25
            label:     H30F, species:    H, occupancy:  0.25

The 256 structures are generated and merged as follows:

.. code-block:: console

    Enumerating ordered configurations.
    Generating 256 out of the 256 possible configurations in the (1, 1, 1) supercell:
    100%|████████████████████████████| 256/256 [00:00<00:00, 1939.00it/s]
    Merging structures...
    Checking symmetry-equivalence: 100%|█| 256/256 [00:04<00:00, 57.51it/
    Merging took     4.48 s and found 46 groups
    Spacegroup: Pbca (61)            multiplicity: 1
    Spacegroup: P1 (1)               multiplicity: 8
    Spacegroup: P2_1 (4)             multiplicity: 4
    Spacegroup: P2_1 (4)             multiplicity: 4
    Spacegroup: P2_1 (4)             multiplicity: 4
    Spacegroup: P1 (1)               multiplicity: 8
    Spacegroup: P2_12_12_1 (19)      multiplicity: 2
    Spacegroup: P-1 (2)              multiplicity: 4
    Spacegroup: Pc (7)               multiplicity: 4
    Spacegroup: P1 (1)               multiplicity: 8
    Spacegroup: Pc (7)               multiplicity: 4
    Spacegroup: P1 (1)               multiplicity: 8
    Spacegroup: P1 (1)               multiplicity: 8
    Spacegroup: P1 (1)               multiplicity: 8
    Spacegroup: Pc (7)               multiplicity: 4
    Spacegroup: P1 (1)               multiplicity: 8
    Spacegroup: P1 (1)               multiplicity: 8
    Spacegroup: P1 (1)               multiplicity: 8
    Spacegroup: P1 (1)               multiplicity: 8
    Spacegroup: P1 (1)               multiplicity: 8
    Spacegroup: P1 (1)               multiplicity: 8
    Spacegroup: P1 (1)               multiplicity: 8
    Spacegroup: P2_1/c (14)          multiplicity: 2
    Spacegroup: P1 (1)               multiplicity: 8
    Spacegroup: P1 (1)               multiplicity: 8
    Spacegroup: P1 (1)               multiplicity: 8
    Spacegroup: Pca2_1 (29)          multiplicity: 2
    Spacegroup: P1 (1)               multiplicity: 8
    Spacegroup: P2_1 (4)             multiplicity: 4
    Spacegroup: P2_1/c (14)          multiplicity: 2
    Spacegroup: P1 (1)               multiplicity: 8
    Spacegroup: P1 (1)               multiplicity: 8
    Spacegroup: Pca2_1 (29)          multiplicity: 2
    Spacegroup: P1 (1)               multiplicity: 8
    Spacegroup: P2_1 (4)             multiplicity: 4
    Spacegroup: P2_1/c (14)          multiplicity: 2
    Spacegroup: P1 (1)               multiplicity: 8
    Spacegroup: Pca2_1 (29)          multiplicity: 2
    Spacegroup: P1 (1)               multiplicity: 8
    Spacegroup: P2_1 (4)             multiplicity: 4
    Spacegroup: P-1 (2)              multiplicity: 4
    Spacegroup: Pc (7)               multiplicity: 4
    Spacegroup: Pc (7)               multiplicity: 4
    Spacegroup: Pc (7)               multiplicity: 4
    Spacegroup: P1 (1)               multiplicity: 8
    Spacegroup: Pbca (61)            multiplicity: 1

These timings don't include the (in this case significant) time to 'reload as molecular crystal'. On my machine it took about 1 minute!






Z' < 1 and multiple disorder groups (EROHEA)
---------------------------------------------



This system contains a caffeine molecule that is disordered at a special symmetry site.  
i.e. it is "disordered by symmetry" where the structure is described by a small unit cell in which the caffeine is disordered over a symmetry axis. 
Unpicking this into two separate orientations is a difficult challenge for the current version of orgdisord.
It currently requires some manual intervention on the CIF file as described below.

This system has one disorder assembly (A) with two groups ("-1" and "-2"). 
The "-" sign in the group label indicates that these sites are at special symmetry positions.
Indeed, although there are 8 symmetry operations, Z = 4.
To generate all of the ordered structures, the code must therefore apply a subset of symmetry operations to 
each group, with the configuration (e.g. (0,0,1,0)) indicating which subgroup of symmetry operations to apply.

Note that the original EROHEA CIF file had to be modified by manually moving the O10 site to another symmetry equivalent position in the CIF file.


.. code-block:: diff

    -O10 O 0.3754(5) 0.4007(5) 0.6805(5) 0.0223(11)
    +O10 O 0.62460(5) 0.4007(5) 0.81950(5) 0.0223(11)

In addition, the C12, H12A, H12B and H12C sites were manually moved to disorder group -2. 

In the future, we might be able to deal with original CIF file without these manual interventions.

* :download:`(modified) EROHEA CIF file <../examples/EROHEA_modified.cif>`

When we run the code with (mostly) default options: ::

    orgdisord enumerate EROHEA_modified.cif --prefix erohea

We get the following output:



Parsing the CIF file
^^^^^^^^^^^^^^^^^^^^


.. code-block:: console
    :caption: parse_cif output for EROHEA

    Parsing disorder in cif file: EROHEA_modified.cif
    Warning: Site 21 (C9) in assembly A has no disorder group. This site will be considered *ordered*. 
    Double check the CIF file.
    Disordered structure:
    Disorder assembly: A
    Contains the following groups:
    Disorder group: -2 has special disorder symmetry and contains [4, 4] symm. ops. and 5 sites:
            label:      O10, species:    O, occupancy:  0.50
            label:      C12, species:    C, occupancy:  0.50
            label:     H12A, species:    H, occupancy:  0.50
            label:     H12B, species:    H, occupancy:  0.50
            label:     H12C, species:    H, occupancy:  0.50

    Disorder group: -1 has special disorder symmetry and contains [4, 4] symm. ops. and 4 sites:
            label:       N2, species:    N, occupancy:  0.50
            label:       N3, species:    N, occupancy:  0.50
            label:      C11, species:    C, occupancy:  0.50
            label:      H11, species:    H, occupancy:  0.50

Enumerating ordered structures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The code correctly generates a total of 16 (2\ :sup:`4` ) ordered structures and outputs the following table:

.. csv-table:: EROHEA enumerated configurations
   :file: ../examples/erohea.csv
   :header-rows: 1

Merging equivalent structures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If we instead re-run the code with the ``-m`` flag, the code will merge symmetry-equivalent structures:

.. code-block:: console

    Enumerating ordered configurations.
    Generating 16 out of the 16 possible configurations in the (1, 1, 1) supercell:
    100%|██████████████████████████████| 16/16 [00:00<00:00, 3958.99it/s]
    Merging structures...
    Checking symmetry-equivalence: 100%|█| 16/16 [00:00<00:00, 2059.06it/
    Merging took     0.01 s and found 5 groups
    Spacegroup: P-1 (2)              multiplicity: 2
    Spacegroup: P1 (1)               multiplicity: 8
    Spacegroup: Cc (9)               multiplicity: 2
    Spacegroup: P2_1/c (14)          multiplicity: 2
    Spacegroup: P2_1/c (14)          multiplicity: 2

and the corresponding table is:

.. csv-table:: EROHEA enumerated configurations (merged)
   :file: ../examples/erohea_merged.csv
   :header-rows: 1


.. warning::

    For systems with special symmetry groups (i.e. with disorder group labels starting with "-"),
    The code may not correctly partition the symmetry operations into subgroups. So please check your output carefully!

    Another serious limitation in these cases is that generated supercells may have overlapping sites.

    We're working to make the code more robust for these cases. 





Z'< 1 and only one disorder group (DASRAU)
-------------------------------------------


This is a Ruddlesden-Popper phase with a disordered butylammonium cation. 
The cation is disordered at a special symmetry position such that each of the 32 symmetry operations generates another configuration.
Assuming only one cation per unit cell, the code will generate, by default, 32 ordered structures.

* :download:`DASRAU CIF file <../examples/DASRAU.cif>`

When we run the code with (mostly) default options: ::

    orgdisord enumerate DASRAU.cif --prefix dasrau

We get the following output:

.. code-block:: console
    :caption: parse_cif output for DASRAU

    Parsing disorder in cif file: DASRAU.cif
    Disordered structure:
    Disorder assembly: A
    Contains the following groups:
    Disorder group: -1 has special disorder symmetry and contains [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] symm. ops. and 17 sites:
            label:       N1, species:    N, occupancy:  0.12
            label:      H1A, species:    H, occupancy:  0.12
            label:      H1B, species:    H, occupancy:  0.12
            label:      H1C, species:    H, occupancy:  0.12
            label:       C1, species:    C, occupancy:  0.12
            label:      H1D, species:    H, occupancy:  0.12
            label:      H1E, species:    H, occupancy:  0.12
            label:       C2, species:    C, occupancy:  0.12
            label:      H2A, species:    H, occupancy:  0.12
            label:      H2B, species:    H, occupancy:  0.12
            label:       C3, species:    C, occupancy:  0.12
            label:      H3A, species:    H, occupancy:  0.12
            label:      H3B, species:    H, occupancy:  0.12
            label:       C4, species:    C, occupancy:  0.12
            label:      H4A, species:    H, occupancy:  0.12
            label:      H4B, species:    H, occupancy:  0.12
            label:      H4C, species:    H, occupancy:  0.12


The code correctly generates a total of 32 (32\ :sup:`1` ) ordered structures and outputs the following table: :download:`DASRAU results table <../examples/dasrau.csv>`



Specifying the disorder components using two ordered structures
----------------------------------------------------------------

Rather than working with CIF files tagged with disorder information, 
it's sometimes easier to work with two ordered structures, one for the major disorder component and one for the minor.
The code can use the difference between the two structures to determine the disorder structure and enumerate based on that.

This functionality is currently limited to the case of one disorder assembly and two disorder groups.

Taking the previous ABABUB example, I manually split the structure into two P1 ordered structures, one for the major disorder component and one for the minor. You can download these here:

* :download:`ABABUB major <../examples/ABABUB_maj.xyz>`
* :download:`ABABUB minor <../examples/ABABUB_min.xyz>`

Note that you can provide the P1 structures in any format that ASE can read; here I've used the extended .xyz format.

If you pass two files to the command-line-interface, it will assume these are the two disorder components ::
    
        orgdisord enumerate ABABUB_maj.xyz ABABUB_min.xyz --no_write







Multiple assemblies, each with same number of groups (HAXPIH)
--------------------------------------------------------------


Correlated assemblies
^^^^^^^^^^^^^^^^^^^^^

This structure (Z=6) has three disorder assemblies, each with two disorder groups. 
If we consider these assemblies to be independent, we can generate 2\ :sup:`6` x 2\ :sup:`6` x 2\ :sup:`6` = 264 144 ordered structures.
This would crash the code and is outside our scope anyway. 
However, there may be cases in which we might want to consider the assemblies as correlated, i.e. we pick from the same disorder group index in each assembly.

* :download:`HAXPIH CIF file <../examples/HAXPIH.cif>`

You can do this by specifying the ``-c`` flag. ::
    
        orgdisord enumerate HAXPIH.cif --prefix haxpih -c

Parsing the CIF file
^^^^^^^^^^^^^^^^^^^^

The code correctly parses the disorder groups and outputs the following:

.. code-block:: console

    Parsing disorder in cif file: HAXPIH.cif
    Disordered structure:
    Disorder assembly: 2
    Contains the following groups:
    Disorder group: 1 contains 6 symm. ops. and 15 sites:
            label:     C230, species:    C, occupancy:  0.26
            label:     C240, species:    C, occupancy:  0.26
            label:     C250, species:    C, occupancy:  0.26
            label:     C260, species:    C, occupancy:  0.26
            label:     H221, species:    H, occupancy:  0.26
            label:     H222, species:    H, occupancy:  0.26
            label:    H2301, species:    H, occupancy:  0.26
            label:    H2302, species:    H, occupancy:  0.26
            label:    H2401, species:    H, occupancy:  0.26
            label:    H2501, species:    H, occupancy:  0.26
            label:    H2502, species:    H, occupancy:  0.26
            label:    H2503, species:    H, occupancy:  0.26
            label:    H2601, species:    H, occupancy:  0.26
            label:    H2602, species:    H, occupancy:  0.26
            label:    H2603, species:    H, occupancy:  0.26

    Disorder group: 2 contains 6 symm. ops. and 15 sites:
            label:     C231, species:    C, occupancy:  0.74
            label:     C241, species:    C, occupancy:  0.74
            label:     C251, species:    C, occupancy:  0.74
            label:     C261, species:    C, occupancy:  0.74
            label:     H223, species:    H, occupancy:  0.74
            label:     H224, species:    H, occupancy:  0.74
            label:    H2311, species:    H, occupancy:  0.74
            label:    H2312, species:    H, occupancy:  0.74
            label:    H2411, species:    H, occupancy:  0.74
            label:    H2511, species:    H, occupancy:  0.74
            label:    H2512, species:    H, occupancy:  0.74
            label:    H2513, species:    H, occupancy:  0.74
            label:    H2611, species:    H, occupancy:  0.74
            label:    H2612, species:    H, occupancy:  0.74
            label:    H2613, species:    H, occupancy:  0.74


    Disorder assembly: 3
    Contains the following groups:
    Disorder group: 1 contains 6 symm. ops. and 12 sites:
            label:     C340, species:    C, occupancy:  0.24
            label:     C350, species:    C, occupancy:  0.24
            label:     C360, species:    C, occupancy:  0.24
            label:     H331, species:    H, occupancy:  0.24
            label:     H332, species:    H, occupancy:  0.24
            label:    H3401, species:    H, occupancy:  0.24
            label:    H3501, species:    H, occupancy:  0.24
            label:    H3502, species:    H, occupancy:  0.24
            label:    H3503, species:    H, occupancy:  0.24
            label:    H3601, species:    H, occupancy:  0.24
            label:    H3602, species:    H, occupancy:  0.24
            label:    H3603, species:    H, occupancy:  0.24

    Disorder group: 2 contains 6 symm. ops. and 12 sites:
            label:     C341, species:    C, occupancy:  0.76
            label:     C351, species:    C, occupancy:  0.76
            label:     C361, species:    C, occupancy:  0.76
            label:     H333, species:    H, occupancy:  0.76
            label:     H334, species:    H, occupancy:  0.76
            label:    H3411, species:    H, occupancy:  0.76
            label:    H3511, species:    H, occupancy:  0.76
            label:    H3512, species:    H, occupancy:  0.76
            label:    H3513, species:    H, occupancy:  0.76
            label:    H3611, species:    H, occupancy:  0.76
            label:    H3612, species:    H, occupancy:  0.76
            label:    H3613, species:    H, occupancy:  0.76


    Disorder assembly: 4
    Contains the following groups:
    Disorder group: 1 contains 6 symm. ops. and 18 sites:
            label:     C421, species:    C, occupancy:  0.46
            label:     C431, species:    C, occupancy:  0.46
            label:     C441, species:    C, occupancy:  0.46
            label:     C451, species:    C, occupancy:  0.46
            label:     C461, species:    C, occupancy:  0.46
            label:     H411, species:    H, occupancy:  0.46
            label:     H412, species:    H, occupancy:  0.46
            label:    H4211, species:    H, occupancy:  0.46
            label:    H4212, species:    H, occupancy:  0.46
            label:    H4311, species:    H, occupancy:  0.46
            label:    H4312, species:    H, occupancy:  0.46
            label:    H4411, species:    H, occupancy:  0.46
            label:    H4511, species:    H, occupancy:  0.46
            label:    H4512, species:    H, occupancy:  0.46
            label:    H4513, species:    H, occupancy:  0.46
            label:    H4611, species:    H, occupancy:  0.46
            label:    H4612, species:    H, occupancy:  0.46
            label:    H4613, species:    H, occupancy:  0.46

    Disorder group: 2 contains 6 symm. ops. and 18 sites:
            label:     C420, species:    C, occupancy:  0.54
            label:     C430, species:    C, occupancy:  0.54
            label:     C440, species:    C, occupancy:  0.54
            label:     C450, species:    C, occupancy:  0.54
            label:     C460, species:    C, occupancy:  0.54
            label:     H413, species:    H, occupancy:  0.54
            label:     H414, species:    H, occupancy:  0.54
            label:    H4201, species:    H, occupancy:  0.54
            label:    H4202, species:    H, occupancy:  0.54
            label:    H4301, species:    H, occupancy:  0.54
            label:    H4302, species:    H, occupancy:  0.54
            label:    H4401, species:    H, occupancy:  0.54
            label:    H4501, species:    H, occupancy:  0.54
            label:    H4502, species:    H, occupancy:  0.54
            label:    H4503, species:    H, occupancy:  0.54
            label:    H4601, species:    H, occupancy:  0.54
            label:    H4602, species:    H, occupancy:  0.54
            label:    H4603, species:    H, occupancy:  0.54


The code will generate 2\ :sup:`6` = 64 ordered structures and outputs the following table: :download:`HAXPIH results table <../examples/haxpih.csv>`



Merging the equivalent structures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Running the code again with the ``-m`` flag: ::

    orgdisord enumerate HAXPIH.cif --prefix haxpih_merged -c -m

We get the following table: :download:`HAXPIH results table <../examples/haxpih_merged.csv>` and this information printed to the terminal:

.. code-block:: console

    Generating 64 out of the 64 possible configurations in the (1, 1, 1) supercell:
    100%|██████████████████████████████| 64/64 [00:00<00:00, 3721.86it/s]
    Merging structures...
    Checking symmetry-equivalence: 100%|█| 64/64 [00:03<00:00, 18.17it/s]
    Merging took     3.53 s and found 16 groups
    Spacegroup: P3_221 (154)         multiplicity: 1
    Spacegroup: P1 (1)               multiplicity: 6
    Spacegroup: P1 (1)               multiplicity: 6
    Spacegroup: P3_2 (145)           multiplicity: 2
    Spacegroup: C2 (5)               multiplicity: 3
    Spacegroup: C2 (5)               multiplicity: 3
    Spacegroup: P1 (1)               multiplicity: 6
    Spacegroup: C2 (5)               multiplicity: 3
    Spacegroup: P1 (1)               multiplicity: 6
    Spacegroup: P1 (1)               multiplicity: 6
    Spacegroup: P1 (1)               multiplicity: 6
    Spacegroup: C2 (5)               multiplicity: 3
    Spacegroup: C2 (5)               multiplicity: 3
    Spacegroup: C2 (5)               multiplicity: 3
    Spacegroup: P1 (1)               multiplicity: 6
    Spacegroup: P3_221 (154)         multiplicity: 1


We can see that the 64 ordered configurations have been merged into 16 groups of equivalent structures.
The first group has a spacegroup of P3_221 (154) and a multiplicity of 1, 
the second group has a spacegroup of P1 (1) and a multiplicity of 6, and so on. 




Manually/programmatically specifying disorder
----------------------------------------------

You can manually create a DisorderedStructure object and pass that to the enumerator to generate the ordered configurations.

As a basic example, let's put a thiophene (C4H4S) molecule in a cubic cell with no symmetry. 
We can then generate several disorder groups that each contains a different configuration of the molecule -- in this case let's just rotate the molecule by 360/5 degrees to move the S atom around the ring. 

Using ASE we can do this as follows:

.. code-block:: python

    # ASE already knows about some molecules :)
    from ase.build import molecule
    from orgdisord.utils import get_new_labels

    # Spacegroup 1 == P1
    sg_number = 1
    # Unit cell. If 3 numbers are given, 
    # they are interpreted as the lengths of the unit cell vectors.
    # If a 3x3 list or array is given, 
    # they are interpreted as the three vectors.
    cell = [6.5, 6.5, 6.5]
    # Create the main atoms object
    atoms = molecule('C4H4S', cell = cell, pbc = True)
    # Good to always have useful labels to help keep track of things
    atoms.set_array('labels', np.array(get_new_labels(atoms)))
    
    # where is the centre of ring:
    ring_c = atoms.positions[:5].mean(axis=0)
    
    # place the ring in the centre of the cell
    atoms.translate(np.array(cell)*0.5 - ring_c)

    # Now we generate 5 structures that will comprise our 5 disorder groups
    # in assembly A, by rotating the ring around.
    images = [atoms]
    for i in range(1,5):
        atoms_temp = atoms.copy()
        # rotate the ring by 360/5 degrees
        atoms_temp.rotate(i*360/5.0, 'x', center='COU')
        images.append(atoms_temp)
    


Now we can create a DisorderedStructure object to later be passed to the enumerator:


.. code-block:: python

    from orgdisord.disordered_structure import DisorderAssembly,
                                                    DisorderGroup,
                                                    DisorderedStructure
    from orgdisord.enumerate import OrderedfromDisordered

    # symmetry operations -- in this case just the identity
    sg = Spacegroup(sg_number)
    symops = sg.get_symop()

    # For the ordered part, we just have an empty Atoms object.
    ordered_atoms = Atoms(cell=cell)

    # make list of DisorderGroup objects to pass into DisorderAssembly
    disorder_groups = []
    for i, rot in enumerate(images):
        disorder_group = DisorderGroup(
            label = str(i+1),
            atoms = rot,
            symmetry_operations=symops,
            tag = i+1,
            occupancy = 1/len(images),
            )
        disorder_groups.append(disorder_group)

    # create disorder assembly
    da = DisorderAssembly(
        label = 'A',
        disorder_groups = disorder_groups,
        tag = 0,
        )
        
    # create disordered structure
    ds = DisorderedStructure(
        ordered_atoms= ordered_atoms,
        Z = 1,
        spacegroup =  sg,
        disorder_assemblies = [da],
        molecular_crystal=True,
        )

    # now we can pass this to the enumerator
    od = OrderedfromDisordered(ds)

    # Let's generate 20 random configurations in a 1x5x5 supercell:
    supercell = [1, 5, 5]
    images = od.get_supercell_configs(
                    supercell = supercell, 
                    maxiters = 20,
                    exclude_ordered = False, 
                    random_configs=True)
    
    # Now we could write them to disk
    # in lots of different formats. E.g.:
    for i, image in enumerate(images):
        image.write(f'C4H4S_example_{i+1:02d}.cif')
        #image.write(f'C4H4S_example_{i+1:02d}.cell')
        #image.write(f'C4H4S_example_{i+1:02d}.xyz')

    # Or we could use the ASE GUI to view them:
    from ase.visualize import view
    view(images)


I used the ASE POVRAY writer to create nice images of each configuration in a 1x6x3 supercell and stitched them into a gif using ImageMagick:

.. image:: images/C4H4S_6x3.gif
    :width: 800px
    :align: center
