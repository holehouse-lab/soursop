.. soursop documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the SOURSOP docs!
=========================================================

*Last updated: April 2022*

**SOURSOP** ( **S**\imulation analysis **O**\f **U**\nstructured and disordered **R**\egion\ **S** **O**\rchestrated in **P**\ython) is a Python package for the analysis of all-atom simulations of unfolded/disordered proteins. It provides a wide range of functionalities that may not be of use for folded proteins, but provide important insight into simulations of intrinsically disordered proteins. **SOURSOP** was formerly *CAMPARITraj*, which was formerly *CTraj*, and includes all the original functionality therein. 

The **SOURSOP** GitHub page can be accessed here: `https://github.com/holehouse-lab/soursop <https://github.com/holehouse-lab/soursop/>`_.

**SOURSOP** was built for the analysis of simulations from the CAMPARI simulation engine in mind. However, it has been succesfully tested on a wide-range of trajectories generated from different software packages. It utilizes the ``mdtraj`` (`http://mdtraj.org <http://mdtraj.org/>`_.) backend for trajectory reading and representation, providing effectively an additional layer of analysis utilities on top of ``mdtraj``.

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

*Update: April 2022*
Moved SOURSOP onto PyPI in anticipation of final release. Additional tests, code clean up etc. Also added 

*Update: March 2022*
Numerous updates to internal code documentation, removal of `soursop_cli`, update to pip install over git to use https.

*Update: July 2021*
The anticipated release of `soursop` is August 2021! We have nearly finished all testing and are finalizing the associated manuscript.

*Update: Jan 2021*
**CAMPARITraj** was recently restructured to homogenize a number of functions, as well as ensure future and backwards compatibility with `mdtraj` version 1.9.5. 
 
*Update: May 2019*
**CAMPARITraj** is currently still being finalized, so we do not currently recommend installation directly from the GitHub repository as tests, final code tweaks, and even major changes to the code-base are occurring constantly. However, in principle, assuming **camparitraj** is building, once downloaded from GitHub it can be installed via




About
========
**SOURSOP** is built and maintained by Jared Lalmansingh (Pappu lab) and `Alex Holehouse <https://holehouselab.com/>`_. It's development was supported financially and intelluctually by the  `Molecular Sciences Software Institute (MOLSSI) <https://molssi.org/>`_.
