Installation
=========================================================

SOURSOP is built on top of the incredible `MDTraj
<https://mdtraj.org/>`_. As such, once MDTraj is installed, SOURSOP can then be easily installed using ``pip``.

Install using conda
----------------------
We recommend using SOURSOP with `[mini]conda
<https://docs.conda.io/en/latest/miniconda.html>`_, and the link here can get get miniconda set up and installed. 

Once conda is set up, you can open a terminal and install MDTraj followed by SOURSOP, and you should be good to go. Specifically, to start from scratch and install soursop in a clean environment:

.. code-block ::

	# create a new conda environment called soursop. Note the
	# name (after -n flag) can be anything
	conda create -n soursop python=3.9 
	
	# activate your new conda environment
	conda activate soursop
	
	# set conda forge (if this is already done, re-running it
	# is no issue
	conda config --add channels conda-forge
	conda config --set channel_priority strict
	
	# install mdtraj
	conda install mdtraj
	
	# install soursop
	pip install soursop
	
	# check everything has worked OK
	python -c "import soursop; soursop.version()"


Alternatively, if you already have an environment set up, you can install MDTraj and then SOURSOP into that environment, although there may be some dependency issues/clashes with MDTraj.

To install the `current development version
<https://github.com/holehouse-lab/soursop>`_ of SOURSOP you can use `pip` to install directly from the git repository using: 

.. code-block::

   pip install git+ssh://git@github.com/holehouselab/soursop.git

This assumes you already have MDTraj installed.

Compile from source
----------------------
To build from source run:

.. code-block::

	git clone git@github.com:holehouse-lab/soursop.git
	
	# move into the soursop root directory where the setup.py file is
	cd soursop
	
	# the -e flag means soursop is linked against a live version of 
	# the code here, so any changes to the files made will be reflected
	# in your system-wide imports.
	pip install -e .
	

This also assumes you already have MDTraj installed.	
	
	
Running tests
----------------------
To ensure SOURSOP has been installed correctly you can run our complete battery of tests using `PyTest
<https://docs.pytest.org/en/7.2.x/>`_. If PyTest is not installed you can install using

.. code-block::

	pip install pytest
	
and from the main SOURSOP directory (where `setup.py` is) run

.. code-block::

	# move to the tests directory
	cd soursop/tests
	
	# run ALL tests (takes some time...)
	pytest -v
	

Dependencies
----------------------
If you install MDTraj first via conda, all these dependencies will be met. However, for completeness, we provide the set of 


* ``mdtraj`` - underlying trajectory reading and representation is done by ``mdtraj``

* ``threadpoolctl`` - package for dealing with multiple threads 

* ``numpy`` - core numerical computing algorithms

* ``scipy`` - core scientific computing algorithms

* ``pandas`` - core data structure (required by ``mdtraj`` but not always included as a dependency)

* ``pytest`` - required for running test


