=======
History
=======

Unreleased
----------

* Support non-diagonal supercell matrices (``--matrix``) and automatic supercell selection by target size and shape (``--target-size``, ``--target-shape``, ``--score-key``).
* Fix supercell cell construction and molecule unwrapping for non-orthogonal cells.
* Convert unit-cell symmetry operations to the supercell frame when merging, including primitive-cell lattice translations modulo the supercell. Translation-equivalent configurations in supercells larger than 1x1x1 are now merged into the same orbit, which changes their reported multiplicities.
* Use canonical orbit labelling on a discrete site template for ``--algo symm``, falling back to pairwise periodic KD-tree matching when sites are not discrete. Add ``--algo canonical`` (no fallback) and ``--algo pairwise``.
* Allow ``--use_disordered_only`` to be turned off with ``--use_all_atoms``.
* Remove the unused ``coords_match``, ``coords_match_symmops`` and ``compare_ref_unmatched`` helpers from ``orgdisord.merge``.

0.1.0 (2022-05-11)
------------------

* First release on PyPI.
