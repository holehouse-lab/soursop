Installation
=========================================================

To install camparitraj during development run.

``pip install .``

Dependencies
************************

``soursop`` has a number of standard dependencies that we suggest should be installed in the following order.

``numpy`` - core numerical computing algorithms
``scipy`` - core scientific computing algorithms
``pandas`` - core data structure (required by ``mdtraj`` but not always included as a dependency)
``mdtraj`` - underlying trajectory reading and representation is done by ``mdtraj``

If installed from source, these are included in the dependency tree defined in ``setup.py`` so in principle should not cause any issues.


