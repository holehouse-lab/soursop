.. camparitraj documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the CAMPARITraj docs!
=========================================================
*Last updated: June 5th 2019*

**CAMPARITraj** is a Python package for the analysis of all-atom simulations of unfolded/disordered proteins. It provides a wide range of functionalities that may not be of use for folded proteins, but provide important insight into simulations of intrinsically disordered proteins. *CAMPARITraj* was formerly known as *CTraj*, and includes all the original functionality therein. 

**CAMPARITraj** was built for the analysis of simulations from the CAMPARI simulation engine in mind. However, it has been succesfully tested on a wide-range of trajectories generated from different software packages. It utilizes the ``mdtraj`` (`http://mdtraj.org <http://mdtraj.org/>`_.) backend for trajectory reading and representation, providing effectively an additional layer of analysis utilities on top of ``mdtraj``.

This documentation is currently being generated as we finalize **CAMPARITraj**. 


.. toctree::
   :maxdepth: 1
   :caption: Contents:
   
   usage/quickstart
   usage/examples
   modules/cttrajectory
   modules/ctprotein
   modules/ctnmr


Quickstart!
==============
*Update: May 2019*
**CAMPARITraj** is currently still being finalized, so we do not currently recommend installation directly from the GitHub repository as tests, final code tweaks, and even major changes to the code-base are occurring constantly. However, in principle, assuming **camparitraj** is building, once downloaded from GitHub it can be installed via

``pip install .``

We plan to release **CAMPARITraj** over ``conda-forge`` once it has finalized and all tests are building. For more information see the installation page.


About
========
**CAMPARITraj** is built and maintained by `Alex Holehouse <https://holehouselab.com/>`_. It's development was supported financially and intelluctually by the  `Molecular Sciences Software Institute (MOLSSI) <https://molssi.org/>`_.
