=======
History
=======

Unreleased
----------

* Packaging modernisation:

  * Use `uv <https://docs.astral.sh/uv/>`_ for environments, dependencies and the lockfile (``uv sync``, ``uv run pytest``), with dependency groups for test, docs, lint, typing and dev. The ``ML`` extra is unchanged.
  * Move the package to a ``src/orgdisord`` layout and switch the build backend to ``uv_build``.
  * Require Python 3.11 or newer; Python 3.11, 3.12 and 3.13 are tested on Linux, macOS and Windows.
  * Lint and format with ruff, add pre-commit hooks, and run a non-blocking ``ty`` type check.
  * Replace Travis CI with GitHub Actions for CI, and publish releases to PyPI from version tags via trusted publishing.
  * Build the documentation with Sphinx and deploy it to GitHub Pages.
  * Remove ``tox``, ``setup.py``, ``requirements*.txt``, ``MANIFEST.in``, ``Makefile`` and ``bumpversion``. Versions are now bumped with ``uv version``.

* Support non-diagonal supercell matrices (``--matrix``) and automatic supercell selection by target size and shape (``--target-size``, ``--target-shape``, ``--score-key``).
* Fix supercell cell construction and molecule unwrapping for non-orthogonal cells.
* Convert unit-cell symmetry operations to the supercell frame when merging, including primitive-cell lattice translations modulo the supercell. Translation-equivalent configurations in supercells larger than 1x1x1 are now merged into the same orbit, which changes their reported multiplicities.
* Use canonical orbit labelling on a discrete site template for ``--algo symm``, falling back to pairwise periodic KD-tree matching when sites are not discrete. Add ``--algo canonical`` (no fallback) and ``--algo pairwise``.
* Allow ``--use_disordered_only`` to be turned off with ``--use_all_atoms``.
* Remove the unused ``coords_match``, ``coords_match_symmops`` and ``compare_ref_unmatched`` helpers from ``orgdisord.merge``.

0.1.0 (2022-05-11)
------------------

* First release on PyPI.
