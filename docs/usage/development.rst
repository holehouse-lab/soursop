Development
=============

SOURSOP was designed to be highly extendable. In particular, third
parties (that's YOU!) are invited to contribute plugins in the
``soursop/plugins`` directory. These can be simple stand-alone functions
or entire classes that implement stateless or stateful functionality.
The bar to developing and sharing a plugin is deliberately kept as low
as possible. If you have an analysis routine you think would be useful
to include, please consider making a pull request, or contact Alex or
Jared about integrating new code into SOURSOP.

.. contents:: On this page
   :local:
   :depth: 2


Setting up a development environment
----------------------------------------

To work on SOURSOP itself, clone the repository and perform an
*editable* install with the test dependencies. With an editable install,
changes you make to the source tree take effect immediately without
reinstalling.

Using pip:

.. code-block:: bash

   git clone https://github.com/holehouse-lab/soursop.git
   cd soursop

   python -m venv .venv
   source .venv/bin/activate
   pip install -e ".[test]"

Using `uv <https://docs.astral.sh/uv/>`_:

.. code-block:: bash

   git clone https://github.com/holehouse-lab/soursop.git
   cd soursop

   uv venv
   source .venv/bin/activate
   uv pip install -e ".[test]"

See the :doc:`installation` page for more installation options.


Repository layout
----------------------

The most relevant parts of the repository are:

* ``soursop/`` - the package source.

  * ``sstrajectory.py`` - the ``SSTrajectory`` class (system-level /
    multi-chain analysis and trajectory loading).
  * ``ssprotein.py`` - the ``SSProtein`` class (single-chain analysis).
  * ``ssnmr.py`` - random-coil chemical-shift prediction.
  * ``sspre.py`` - synthetic PRE profile calculations.
  * ``sssampling.py`` - sampling-quality / PENGUIN tools.
  * ``sstools.py`` - shared numerical helpers.
  * ``plugins/`` - user-contributed analysis plugins.
  * ``tests/`` - the PyTest suite, including the parametrized
    regression suite that recomputes every public observable against
    stored reference values.

* ``docs/`` - the Sphinx documentation sources.


Running the tests
----------------------

SOURSOP ships with an extensive test suite. Any change to the codebase
should be validated against it before opening a pull request. From the
repository root:

.. code-block:: bash

   # run the full battery of tests (this can take several minutes)
   pytest -v

   # run a single module while iterating
   pytest soursop/tests/test_sstrajectory.py -v

   # quietly run the regression and trajectory suites
   pytest soursop/tests/test_reference_observables.py soursop/tests/test_sstrajectory.py -q

The regression suite (``test_reference_observables.py``) recomputes
every public observable against stored reference pickles, so it will
immediately surface any numerical drift introduced by a change. If you
intentionally change the numerical behaviour of an observable, the
reference pickles must be regenerated (see
``soursop/tests/build_reference/``).


Coding and documentation conventions
----------------------------------------

To keep the codebase consistent and the API documentation usable:

* Public functions and methods use **numpy-style docstrings**: a
  one-line summary, an extended description, typed ``Parameters``,
  ``Returns`` (with array shapes where relevant), ``Raises`` where
  applicable, and a short ``Example``. See ``sstrajectory.py`` for the
  canonical style.

* Prefer ``numpy`` over ``pandas`` for large-scale numerical data.

* Add tests for new functionality. New observables should ideally be
  added to the regression suite so future changes cannot silently
  perturb them.

* Keep changes focused; behaviour-preserving cleanup should not change
  numerical results unless that is the explicit intent of the change.


General workflow for adding plugins
----------------------------------------

1. Fork the SOURSOP repository and add your code into the
   ``soursop/plugins`` directory.

2. Once your code is complete, write a set of stand-alone tests using
   PyTest. If this is challenging, please don't hesitate to reach out to
   Alex and Jared about how best to do this.

3. Once your code is working and the tests pass, open a pull request to
   merge your fork back into the main SOURSOP branch.


Plugin example
----------------------

By way of example, there is a simple demo plugin in the `soursop/plugins
<https://github.com/holehouse-lab/soursop/tree/master/soursop/plugins/>`_
directory. To use it:

.. code-block:: python

   from soursop.sstrajectory import SSTrajectory
   T = SSTrajectory('ntl9.xtc', 'ntl9.pdb')

   NTL9_CP = T.proteinTrajectoryList[0]

   # import the sparrow_plugin module
   from soursop.plugins import sparrow_plugin

   print(sparrow_plugin.get_protein_net_charge(NTL9_CP))

   spobj = sparrow_plugin.SparrowProtein(NTL9_CP)
   print(spobj.NCPR)


Reporting bugs and requesting features
------------------------------------------

If you find a bug, typo, or error, please `raise an issue on GitHub
<https://github.com/holehouse-lab/soursop/issues>`_. Feature requests
and pull requests are welcome - for larger contributions it is worth
contacting Alex or Jared first so we can help align the design with the
rest of the package.
