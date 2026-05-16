Installation
=========================================================

SOURSOP is a pure-Python package built on top of the excellent `MDTraj
<https://mdtraj.org/>`_, which provides the underlying trajectory reading
and representation. SOURSOP is distributed on `PyPI
<https://pypi.org/project/soursop/>`_ and can be installed with either
``pip`` or `uv <https://docs.astral.sh/uv/>`_. It can also be installed
directly from the `GitHub repository
<https://github.com/holehouse-lab/soursop>`_ if you want the latest
(unreleased) development version.

SOURSOP requires **Python 3.7 or newer**. All other dependencies
(including MDTraj) are resolved automatically by the installer, so in
most cases a single command is all that is needed.

.. contents:: On this page
   :local:
   :depth: 2


Quick start
----------------------

If you just want the latest stable release into your current
environment:

.. code-block:: bash

   pip install soursop

or, equivalently, with uv:

.. code-block:: bash

   uv pip install soursop

Then verify the install (see `Verifying the installation`_ below):

.. code-block:: bash

   python -c "import soursop; print(soursop.__version__)"


Install using pip
----------------------

From PyPI (recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~

The simplest route is to install the most recent release from PyPI. We
strongly recommend doing this inside an isolated virtual environment so
SOURSOP and its dependencies do not interfere with other projects:

.. code-block:: bash

   # create and activate a virtual environment
   python -m venv soursop-env
   source soursop-env/bin/activate        # on Windows: soursop-env\Scripts\activate

   # install the latest stable release (pulls in mdtraj, numpy, scipy, ...)
   pip install soursop

   # check everything works
   python -c "import soursop; print(soursop.__version__)"

From GitHub (development version)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To install the `current development version
<https://github.com/holehouse-lab/soursop>`_ directly from the GitHub
``master`` branch, use ``pip`` with a VCS URL.

Over HTTPS (no GitHub credentials required):

.. code-block:: bash

   pip install "git+https://github.com/holehouse-lab/soursop.git"

Over SSH (requires an SSH key configured with GitHub):

.. code-block:: bash

   pip install "git+ssh://git@github.com/holehouse-lab/soursop.git"

You can pin to a specific branch, tag, or commit by appending
``@<ref>``, e.g.:

.. code-block:: bash

   pip install "git+https://github.com/holehouse-lab/soursop.git@master"


Install using uv
----------------------

`uv <https://docs.astral.sh/uv/>`_ is a fast, modern Python package and
environment manager that is a drop-in replacement for many ``pip`` and
``virtualenv`` workflows. SOURSOP installs cleanly with uv.

From PyPI (recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~

Into an existing environment using the pip-compatible interface:

.. code-block:: bash

   uv pip install soursop

Or, to create a fresh, self-contained environment for SOURSOP:

.. code-block:: bash

   # create and activate a uv-managed virtual environment
   uv venv soursop-env
   source soursop-env/bin/activate        # on Windows: soursop-env\Scripts\activate

   uv pip install soursop

If you are managing a project with a ``pyproject.toml``, you can instead
add SOURSOP as a project dependency:

.. code-block:: bash

   uv add soursop

From GitHub (development version)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

uv accepts the same VCS URLs as ``pip``.

Over HTTPS:

.. code-block:: bash

   uv pip install "git+https://github.com/holehouse-lab/soursop.git"

Over SSH:

.. code-block:: bash

   uv pip install "git+ssh://git@github.com/holehouse-lab/soursop.git"

To add the development version as a project dependency:

.. code-block:: bash

   uv add "git+https://github.com/holehouse-lab/soursop.git"


Install using conda
----------------------

SOURSOP itself is installed from PyPI with ``pip``, but if you prefer a
`[mini]conda <https://docs.conda.io/en/latest/miniconda.html>`_
environment you can install MDTraj from ``conda-forge`` first and then
install SOURSOP with ``pip`` into the same environment:

.. code-block:: bash

   # create and activate a new conda environment (the name is arbitrary)
   conda create -n soursop python=3.10
   conda activate soursop

   # use the conda-forge channel
   conda config --add channels conda-forge
   conda config --set channel_priority strict

   # install mdtraj from conda-forge
   conda install mdtraj

   # install soursop
   pip install soursop

   # check everything has worked
   python -c "import soursop; print(soursop.__version__)"


Install from source (editable / development install)
-------------------------------------------------------

If you want to modify SOURSOP, contribute changes, or track the
bleeding-edge code, clone the repository and perform an *editable*
install. With an editable install, changes you make to the source tree
are picked up immediately without reinstalling.

With pip:

.. code-block:: bash

   git clone https://github.com/holehouse-lab/soursop.git
   cd soursop

   # editable install plus the optional test dependencies
   pip install -e ".[test]"

With uv:

.. code-block:: bash

   git clone https://github.com/holehouse-lab/soursop.git
   cd soursop

   uv venv
   source .venv/bin/activate
   uv pip install -e ".[test]"

The ``[test]`` extra additionally installs PyTest so you can run the
test suite (see below). Omit it (``pip install -e .``) if you do not
need to run the tests.


Verifying the installation
-----------------------------

After installing, confirm SOURSOP imports and reports a version:

.. code-block:: bash

   python -c "import soursop; print(soursop.__version__)"

A quick functional smoke test:

.. code-block:: python

   from soursop.sstrajectory import SSTrajectory
   help(SSTrajectory)


Running tests
----------------------

SOURSOP ships with an extensive test suite (including a large
regression suite that recomputes every public observable against stored
reference values). Running it is the most thorough way to confirm a
correct installation. The tests require PyTest, which is included in the
``[test]`` optional dependency group (or can be installed directly):

.. code-block:: bash

   pip install pytest        # or: uv pip install pytest

From a source checkout, run the full battery of tests from the
repository root:

.. code-block:: bash

   pytest -v

or run just a single module, e.g.:

.. code-block:: bash

   pytest soursop/tests/test_sstrajectory.py -v

The full suite can take several minutes because of the regression
tests.


Dependencies
----------------------

All of the following are installed automatically when you install
SOURSOP with ``pip`` or ``uv``. They are listed here for completeness:

* ``mdtraj`` (>= 1.9.5) - underlying trajectory reading and
  representation.

* ``numpy`` (>= 1.20.0) - core numerical computing.

* ``scipy`` (>= 1.5.0) - core scientific computing routines.

* ``pandas`` (>= 0.23.0) - data structures (also required by MDTraj).

* ``threadpoolctl`` (>= 2.2.0) - control over native thread pools.

* ``natsort`` - natural sorting of trajectory/replicate file paths.

* ``matplotlib`` - plotting support used by some analysis helpers.

* ``cython`` - build-time/optional acceleration dependency.

Optional (only needed to run the test suite):

* ``pytest`` (>= 6.1.2) - test runner; installed via the ``[test]``
  extra (``pip install "soursop[test]"``).
