SOURSOP
==============================
[![Build Status](https://github.com/jaredl7/soursop/actions/workflows/soursop-ci.yml/badge.svg?branch=master)](https://github.com/jaredl7/soursop/actions)
[![codecov](https://codecov.io/gh/jaredl7/soursop/branch/master/graph/badge.svg?token=RHGII0235L)](https://codecov.io/gh/jaredl7/soursop)
[![Documentation Status](https://readthedocs.org/projects/soursop/badge/?version=latest)](https://soursop.readthedocs.io/en/latest/?badge=latest)
## ABOUT
SOURSOP is a Python-based simulation analysis package for working with intrinsically disordered and unfolded proteins. It is built on top of [mdtraj](https://mdtraj.org/), and was developed by Jared Lalmansingh and Alex Holehouse. 

The current stable release candidate on PyPI is 0.2.3 (Feb 2023).

## DOCUMENTATION
All documentation, including installation information [can be found here](https://soursop.readthedocs.io/). 

## ERRORS, FEATURES, REQUESTS
If you find a bug, typo, or error [please raise an issue or GitHub](https://github.com/holehouse-lab/soursop/issues).

If you wish to add a new feature, please see our Development information in the docs (especially for adding plugins).

## MISCELLANEOUS
* As of right now, the continuous integration fails because of a specific mismatch in how an edge-case error is handled between different versions of mdtraj. All other tests are passing and SOURSOP is ready for production.

## PREPRINT
To read about SOURSOP please see our preprint:

**SOURSOP: A Python package for the analysis of simulations of intrinsically disordered proteins**

Jared Lalmansingh, Alex Keeley, Kiersten Ruff, Rohit Pappu, Alex Holehouse


#### Copyright
Copyright (c) 2015-2023 under the GNU LESSER GENERAL PUBLIC LICENSE 

# Changelog
#### Update February 2023
* Added plugins example
* Added additional tests and finalized documentation

#### Update July 2022
For version 0.2.1 we introduce potentially breaking changes into how COM distances are reported. 

##### Details:
All center-of-mass (COM)-based functions now return distances (as before) and relative positions in x/y/z (this is new) in Angstroms, not nanometers. Previously, the various center-of-mass based functions that returned absolute positions (i.e. 3xn matrices of x/y/z vs frame number) returned them where x/y/z was in units nm. This is fine if you know this, but means if you manually calculate distance between two COM vectors you'd get a distance in nanometers and not Angstroms. This is an unhelpful and unexpected behavior given all other distances are in Angstroms, so we have made the decision to fully update to Angstroms even for vector positions. This DOES NOT break any code internal to soursop, but if you were using COM positions to manually calculate distances these may need to be recalculated.

#### Update April 2022
For version 0.1.9 the documentation has been extensively extended

#### Update July 2021
CAMPARITraj is SOURSOP! For the final release we have re-named and re-branded CAMPARITraj as SOURSOP. This change in the name serves two important purposes.

Firstly, CAMPARITraj was borne out of a collection of scripts and code built to work with CAMPARI. However, it has evolved into a stand-alone package for the analysis of all-atom simulations of IDPs and IDRs, and importantly, much of the analysis it performs can also be done in CAMPARI. As such, we felt it was important to decoupled SOURSOP from CAMPARI, both to avoid the implication that CAMPARI cannot perform analyses itself, and to avoid a scenario in which it may appear that this package only works with CAMPARI simulations.

Secondly, there is a long, rich tradition of naming software tools after drinks in the Pappu lab (CAMPARI, ABSINTH, CIDER, LASSI etc.). As such SOURSOP is first-author Jared's Caribbean twist on this theme!

### WARNING
We are currently and systematically updating all of the CAMPARITraj codebase to SOURSOP. As such, for now, we recommend using a previous version. Note the SOURSOP change breaks all backwards compatibility. Sorry about that.

#### Update December 2020
Release 0.1.2 includes updated support to ensure CAMPARITraj will continue to work with MDTraj 1.9.5, as well as numerous additional updated.

#### Update May 2019
This is the *development* repository of CAMPARITraj and SHOULD NOT be used for production. Seriously, it is being modified constantly and with no building requirements during code pushes. If you want a building copy PLEASE contact Alex directly! [last touched June 24th 2019].

#### Acknowledgements
Project based on the
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.0.
