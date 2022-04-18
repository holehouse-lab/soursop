Installation
=========================================================

SOURSOP is built on top of the incredible `MDTraj
<https://mdtraj.org/>`_. As such, once MDTraj is installed, SOURSOP can then be easily installed using ``pip``.


Install using conda
----------------------

We HIGHLY recommend installing ``mdtraj`` using ``conda`` via the ``conda-forge`` channel. Specifically, this can be done using:

.. code-block:: 

   conda install -c conda-forge mdtraj
	
Once ``mdtraj`` is installed, SOURSOP can be installed via

.. code-block:: 
   
   pip install soursop

This will install the current stable version. 

To install the `current development version
<https://github.com/holehouse-lab/soursop>`_ of SOURSOP you can use `pip` to install directly from the git repository using: 

.. code-block::

   pip install git+ssh://git@github.com/holehouselab/soursop.git


Compile from source
----------------------
To build from source run:

.. code-block::
   
   git clone git@github.com:holehouse-lab/soursop.git

   cd soursop

   pip install -e .


Dependencies
----------------------

``soursop`` has a number of standard dependencies that we suggest should be installed in the following order.

``numpy`` - core numerical computing algorithms
``scipy`` - core scientific computing algorithms
``pandas`` - core data structure (required by ``mdtraj`` but not always included as a dependency)
``mdtraj`` - underlying trajectory reading and representation is done by ``mdtraj``

