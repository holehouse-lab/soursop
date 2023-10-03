.. soursop documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the SOURSOP Documentation!
=========================================================

*Last updated: July 2023*

**SOURSOP** (**S**\imulation analysis **O**\f **U**\nfolded **R**\egion\ **S** **O**\f **P**\roteins) is a Python package for the analysis of all-atom and coarse-grained simulations of unfolded and disordered proteins. It provides a wide range of functionalities that may not be of use for folded proteins, but provide important insight into simulations of intrinsically disordered proteins. **SOURSOP** was formerly *CAMPARITraj*, which was formerly *CTraj*, and includes all the original functionality therein. 

The **SOURSOP** GitHub page can be accessed here: `https://github.com/holehouse-lab/soursop <https://github.com/holehouse-lab/soursop/>`_.

**SOURSOP** was built for the analysis of simulations with the CAMPARI simulation engine in mind. However, it has been successfully tested on a wide range of trajectories generated from different software packages. It utilizes the ``mdtraj`` (`http://mdtraj.org <http://mdtraj.org/>`_.) backend for trajectory reading and representation, but focusses on providing analysis routines relevant for characterizing ensembles of disordered and unfolded proteins through the lens of polymer physics.

This documentation is currently being generated as we finalize **SOURSOP**. 

.. toctree::
   :maxdepth: 1
   :caption: Contents:
   
   usage/overview
   usage/installation
   usage/examples
   modules/sstrajectory
   modules/ssprotein
   modules/ssnmr
   modules/sspre
   usage/development

Changelog
==========
*Update: July. 2023*
Added ``explicit_residue_checking`` option to SSTrajectory constructor to make parsing solvated .gro files or files where non-protein molecules are included in the same chain possible and easy.

*Update: Jan. 2023*
Finalization of code and documentation ahead of preprint deposition.

*Update: Sept. 2022*
Additional tests and docs updates ahead of preprint.

*Update: April 2022*
Moved SOURSOP onto PyPI in anticipation of the final release. Additional tests, code clean up etc. 

*Update: March 2022*
Numerous updates to internal code documentation, removal of `soursop_cli`, update to pip install over git to use https.

*Update: July 2021*
The anticipated release of `soursop` is August 2021! We have nearly finished all testing and are finalizing the associated manuscript.

*Update: Jan 2021*
**CAMPARITraj** was recently restructured to homogenize a number of functions, as well as ensure future and backwards compatibility with `mdtraj` version 1.9.5. 
 
*Update: May 2019*
**CAMPARITraj** is currently still being finalized, so we do not currently recommend installation directly from the GitHub repository as tests, final code tweaks, and even major changes to the code base are constantly occurring. However, in principle, assuming **camparitraj** is building, once downloaded from GitHub, it can be installed via


About
========
**SOURSOP** is built and maintained by Jared Lalmansingh (Pappu lab) and `Alex Holehouse <https://holehouselab.com/>`_. Its development was supported financially and intellectually by the  `Molecular Sciences Software Institute (MOLSSI) <https://molssi.org/>`_. It was also supported by NSF grant no. 2128068 to Alex, and we thank members of the `Water and Life Interface Institute (WALII) <https://www.walii.science/>`_, supported by NSF DBI grant #2213983, for helpful discussions. 
