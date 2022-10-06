.. highlight:: shell

============
Installation
============



From sources
------------

The sources for SODORG Renewal can be downloaded from the `Github repo`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone git://github.com/jkshenton/sodorg_renewal

Or download the `tarball`_:

.. code-block:: console

    $ curl -OJL https://github.com/jkshenton/sodorg_renewal/tarball/main

Once you have a copy of the source, you can install it with:

.. code-block:: console

    $ python setup.py install

or 

.. code-block:: console

    $ pip install .


.. tip::
    
    It is recommended to use a new `virtual environment`_ or new `conda environment`_ rather than installing SODORG Renewal system wide.

    The package has not been tested on the Windows operating system. It's recommended to use the WSL (Windows Subsystem for Linux) if you are using Windows.
    Though in that case you may need to `set up an X11 server`_ to view the generated structures in the ASE GUI.


.. _Github repo: https://github.com/jkshenton/sodorg_renewal
.. _tarball: https://github.com/jkshenton/sodorg_renewal/tarball/main
.. _virtual environment: http://docs.python-guide.org/en/latest/dev/virtualenvs/
.. _conda environment: https://conda.io/docs/user-guide/tasks/manage-environments.html
.. _set up an X11 server: https://stackoverflow.com/a/61110604





Stable release
--------------

We will eventually publish a PyPi package such that you can install it with:

.. code-block:: console

    $ pip install sodorg_renewal


**but this is not yet available**. See above for how to build it from source. 

.. To install SODORG Renewal, run this command in your terminal:

.. .. code-block:: console

..     $ pip install sodorg_renewal

.. This is the preferred method to install SODORG Renewal, as it will always install the most recent stable release.

.. If you don't have `pip`_ installed, this `Python installation guide`_ can guide
.. you through the process.

.. .. _pip: https://pip.pypa.io
.. .. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/
