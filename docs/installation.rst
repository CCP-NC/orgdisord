.. highlight:: shell

============
Installation
============

orgdisord requires Python 3.11 or newer.

Stable release
--------------

Install the latest release from `PyPI <https://pypi.org/project/orgdisord/>`_ with pip:

.. code-block:: console

    $ pip install orgdisord

If you use `uv <https://docs.astral.sh/uv/>`_, you can add it to a project:

.. code-block:: console

    $ uv add orgdisord

or install the command-line tool into its own isolated environment:

.. code-block:: console

    $ uv tool install orgdisord

Optional extras
^^^^^^^^^^^^^^^

Machine-learning related features need the optional ``ML`` extra:

.. code-block:: console

    $ pip install "orgdisord[ML]"

From source
-----------

The sources for orgdisord can be found on the `GitHub repository`_. Clone it with:

.. code-block:: console

    $ git clone https://github.com/CCP-NC/orgdisord.git
    $ cd orgdisord

For development, we use `uv <https://docs.astral.sh/uv/>`_ to manage the environment.
From the root of the repository, running:

.. code-block:: console

    $ uv sync

creates a virtual environment in ``.venv`` with orgdisord installed in editable mode,
together with the development dependency groups (tests, linting, type checking and documentation).
To work on the documentation only, use:

.. code-block:: console

    $ uv sync --group docs

and then build the HTML documentation with:

.. code-block:: console

    $ uv run --group docs sphinx-build -W --keep-going -b html docs docs/_build/html

The output is written to ``docs/_build/html``; open ``docs/_build/html/index.html`` in a browser.

.. tip::

    Using a virtual environment or conda environment rather than installing orgdisord system wide is recommended.
    If you prefer plain pip, the equivalent is ``python -m venv .venv``, activating it, and then ``pip install -e .``.

The package has not been tested on the Windows operating system. It's recommended to use the WSL (Windows Subsystem for Linux) if you are using Windows.
Though in that case you may need to `set up an X11 server`_ to (optionally) view the generated structures in the ASE GUI.


.. _GitHub repository: https://github.com/CCP-NC/orgdisord
.. _set up an X11 server: https://stackoverflow.com/a/61110604
