.. highlight:: shell

============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every little bit
helps, and credit will always be given.

You can contribute in many ways:

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at https://github.com/CCP-NC/orgdisord/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs
~~~~~~~~

Look through the GitHub issues for bugs. Anything tagged with "bug" and "help
wanted" is open to whoever wants to implement it.

Implement Features
~~~~~~~~~~~~~~~~~~

Look through the GitHub issues for features. Anything tagged with "enhancement"
and "help wanted" is open to whoever wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

orgdisord could always use more documentation, whether as part of the
official orgdisord docs, in docstrings, or even on the web in blog posts,
articles, and such.

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://github.com/CCP-NC/orgdisord/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

Get Started!
------------

Ready to contribute? Here's how to set up `orgdisord` for local development.

1. Fork the `orgdisord` repo on GitHub.
2. Clone your fork locally::

    $ git clone git@github.com:your_name_here/orgdisord.git

3. Install `uv <https://docs.astral.sh/uv/>`_ if you don't have it, then install
   your local copy into a project-local virtual environment. This installs the
   package in editable mode together with the ``test``, ``lint``, ``typing``,
   ``docs`` and ``dev`` dependency groups::

    $ cd orgdisord/
    $ uv sync

4. Install the git hooks, so that ruff runs automatically on each commit::

    $ uv run pre-commit install

5. Create a branch for local development::

    $ git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally.

6. When you're done making changes, check that your changes pass the linter,
   the formatter, the type checker and the tests::

    $ uv run ruff check
    $ uv run ruff format
    $ uv run ty check
    $ uv run pytest

   ``ty`` is currently advisory and does not block CI, but please do not add
   new type errors. ``ruff format`` rewrites files in place; CI runs it in
   check mode, so commit its output.

   To build the documentation locally::

    $ uv run --group docs sphinx-build -W --keep-going -b html docs docs/_build/html

   Then open ``docs/_build/html/index.html`` in your browser.

7. Commit your changes and push your branch to GitHub::

    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature

8. Submit a pull request through the GitHub website.

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring, and add the
   feature to the list in README.rst.
3. The pull request should pass the ``Lint (ruff)`` and ``Test`` checks on
   GitHub Actions. Tests run on Python 3.11, 3.12 and 3.13 on Linux, macOS and
   Windows. You can follow the results at
   https://github.com/CCP-NC/orgdisord/actions.
4. If you change the public behaviour of the package, add an entry to the
   "Unreleased" section of HISTORY.rst.

Tips
----

To run a subset of tests::

    $ uv run pytest tests/test_cif_parser.py

To run a single test by name::

    $ uv run pytest -k "name_of_test"


Deploying
---------

A reminder for the maintainers on how to release a new version.
Make sure all your changes are merged to ``main`` and CI is green, and that
HISTORY.rst has an entry for the release. Then bump the version (this updates
the ``version`` field in ``pyproject.toml`` and ``uv.lock``) and commit it::

    $ uv version --bump patch # possible: major / minor / patch
    $ git commit -am "Release vX.Y.Z"
    $ git tag vX.Y.Z
    $ git push
    $ git push --tags

Pushing the ``v*`` tag triggers the release workflow
(``.github/workflows/release.yml``). It checks that the tag matches the
project version, builds the distributions, publishes them to PyPI using
trusted publishing and creates a GitHub release.
