SOURSOP
==============================
[![Build Status](https://github.com/holehouse-lab/soursop/actions/workflows/soursop-ci.yml/badge.svg?branch=master)](https://github.com/holehouse-lab/soursop/actions)
[![codecov](https://codecov.io/gh/holehouse-lab/soursop/branch/master/graph/badge.svg)](https://codecov.io/gh/holehouse-lab/soursop)
[![Documentation Status](https://readthedocs.org/projects/soursop/badge/?version=latest)](https://soursop.readthedocs.io/en/latest/?badge=latest)

## ABOUT
SOURSOP is a Python-based simulation analysis package for working with intrinsically disordered and unfolded proteins. It is built on top of [mdtraj](https://mdtraj.org/), and was developed by Jared Lalmansingh and Alex Holehouse. 

The current stable release on PyPI is 0.2.7 (May 2026).

## INSTALLATION
SOURSOP can be installed from PyPI with either `pip` or [uv](https://docs.astral.sh/uv/):

```bash
pip install soursop
# or
uv pip install soursop
```

To install the latest development version directly from GitHub:

```bash
pip install "git+https://github.com/holehouse-lab/soursop.git"
```

Verify the install with:

```bash
python -c "import soursop; print(soursop.__version__)"
```

Full installation instructions (conda, editable/source installs, running the tests) are in the [documentation](https://soursop.readthedocs.io/).

## DOCUMENTATION
All documentation, including installation information [can be found here](https://soursop.readthedocs.io/). 

## ERRORS, FEATURES, REQUESTS
If you find a bug, typo, or error [please raise an issue on GitHub](https://github.com/holehouse-lab/soursop/issues).

If you wish to add a new feature, please see our Development information in the docs (especially for adding plugins).

## PUBLICATION
To read about SOURSOP please see our paper:

Lalmansingh, J. M., Keeley, A. T., Ruff, K. M., Pappu, R. V. & Holehouse, A. S. SOURSOP: A Python Package for the Analysis of Simulations of Intrinsically Disordered Proteins. J. Chem. Theory Comput. (2023). doi:10.1021/acs.jctc.3c00190

* [Journal link](https://pubs.acs.org/doi/full/10.1021/acs.jctc.3c00190)
* [Paper PDF](https://www.dropbox.com/s/bd5szapvxpn83r6/soursop_jctc.pdf?dl=0)

#### Copyright
Copyright (c) 2014-2026 under the GNU LESSER GENERAL PUBLIC LICENSE 

# Changelog

#### Update May 2026 (0.2.7)

**Bug fixes**
* `ssnmr`: Corrected three copy-paste errors in the propagation of Ser/Thr/Tyr neighbour-correction values to pSer/pThr/pTyr slots in the glycine-specific random-coil correction tables (`gly_ca_c`, `gly_co_c`, `gly_n_c`). Phospho-residues at the i+1 neighbour position of glycine now correctly apply the intended Kjaergaard correction; previously these slots were silently zeroed or used the wrong row, causing errors of up to ~2.7 ppm in rare Gly–phospho-Gly contexts.
* `ssprotein.get_Q`: Replaced manual sigmoid with `scipy.special.expit` to avoid `np.exp` overflow for large positive arguments.
* `ssprotein.get_Q`: Fixed a division-by-zero in residue–residue contact matrix normalisation when a pair had zero atomic contacts in all frames; affected entries are now safely set to zero rather than producing NaN.
* `ssprotein.get_contact_map`: Added an explicit check for glycine residues before computing `sidechain-heavy` contact maps; an `SSException` is now raised consistently rather than producing mdtraj-version-dependent behaviour (this resolves the historical edge case that previously caused the continuous-integration failure noted above).
* `ssprotein.get_clusters`: Fixed a division-by-zero in RMSD centroid selection when all pairwise distances within a cluster are identical (std == 0); the first frame is returned as the centroid in that degenerate case.
* `ssutils.set_numpy_threads`: Raises a descriptive `SSException` when neither MKL nor OpenBLAS is detected (e.g. macOS Apple Accelerate) rather than silently failing.
* `sstrajectory.get_interchain_contact_map` / `get_interchain_distance`: Fixed IndexError raised for terminal cap residues that lack a CA atom; added eager protein-ID and mode validation; wrapped per-residue computations in try/except so a single problematic residue no longer aborts the entire map.
* `sstrajectory.parallel_load_trjs`: Fixed a broken topology-broadcast condition that prevented a single shared topology (or a one-element list) from being applied to multiple trajectories; the argument now accepts either a `str` or a per-trajectory list and raises a descriptive `SSException` on a topology/trajectory count mismatch. This previously broke PENGUIN / `sssampling` parallel loading.
* `sstrajectory` (`__readTrajectory`): Fixed the PDB-vs-trajectory unit-cell comparison, which compared the z-dimension three times and never checked the y-dimension; the warning is now an f-string that reports the actual differing cell dimensions (and several message typos were corrected).
* `ssprotein` (initialization, `__check_cg_onebead`): No longer builds the one-letter sequence purely to obtain a residue count - that path raised `KeyError` on coarse-grained or non-standard residue names and could crash construction. It now compares the atom and residue counts directly.
* `ssprotein` (`reset_cache` / `__init__`): Consolidated the duplicated cache-initialisation logic into a single private helper so a freshly constructed object and one after `reset_cache()` are identical. Previously `reset_cache()` forced `__cg_onechain = False`, silently disabling coarse-grained handling and forcing the slow per-residue path after a reset.

**New features and improvements**
* `sstrajectory.get_interchain_contact_map` and `sstrajectory.get_interchain_distance`: Added a `stride` parameter that subsamples frames *before* the mdtraj distance computation, giving substantial speed-ups on large trajectories without any post-hoc loss of information.
* `sstrajectory`: Applied `functools.wraps` to the `lazy_loading_single_protein_trajectory` decorator so that `get_overall_radius_of_gyration`, `get_overall_hydrodynamic_radius`, and `get_overall_asphericity` preserve their docstrings and signatures at runtime (this also fixes Sphinx autodoc rendering for these three methods).
* `ssprotein`: Converted all legacy `%`-format strings to f-strings for consistency and readability.
* `ssprotein` (initialization): Replaced the per-residue `topology.select()` loop in `__get_resid_with_CA` (the documented initialization bottleneck) with a single topology pass, faithfully reproducing the original memoisation side effects so `resid_with_CA`, cap flags and `get_CA_index` results (and their dtypes) are unchanged.
* `ssprotein`: Replaced `== None` checks with `is None`; `__get_subtrajectory` now uses native mdtraj slicing (`traj[::stride]`).
* `sstrajectory`: Removed dead code and unused imports (`time`, `copy`) and the unused `single_chain_sim` value threaded through `__get_proteins`; simplified per-chain atom collection in `__get_all_proteins`.
* `sstrajectory.get_interchain_distance_map`: Hoisted the per-residue centre-of-mass computation out of the O(n1*n2) nested loop (chain-2 COMs were previously recomputed once per chain-1 residue). The non-periodic distance is now vectorised per row; the periodic path is preserved per-pair. Byte-for-byte identical output, with the redundant O(n1*n2) COM computations reduced to O(n1+n2).
* `sstrajectory.get_interchain_contact_map`: For the default `mode='atom'` (non-periodic), each residue's named-atom position is now computed once instead of via an O(n1*n2) loop of `get_interchain_distance` calls. Every other mode and the periodic path are unchanged. Output is byte-for-byte identical (verified against a pre-refactor baseline across CA/COM/atom and the untouched ca/closest-heavy/sidechain modes).
* `ssprotein` (`get_internal_scaling`, `get_scaling_exponent`, `get_local_to_global_correlation`): For `mode='CA'`, hoisted the per-residue CA centre-of-mass out of the O(n^2) sequence-separation loops via a new shared `__ca_position_cache` helper. Previously `get_inter_residue_atomic_distance` recomputed each residue's CA position from scratch on every pair (this path was *not* memoised, unlike the COM path). The helper reproduces the exact `__residue_atom_lookup -> atom_slice -> __get_subtrajectory -> 10*compute_center_of_mass` sequence, so results - including the seeded error-bootstrap / pair-resampling - are byte-for-byte identical, while the CA-position work drops from O(n^2) to O(n) (≈20-40x faster on the affected calls).

**Testing**
* Renamed all test-data trajectory files to carry an explicit resolution suffix (`_AA` for all-atom, `_CG` for coarse-grained), e.g. `ctl9.pdb` → `ctl9_AA.pdb`. This makes the resolution unambiguous for both tests and build scripts.
* Added two coarse-grained test trajectories (`sigA_CG`, `synth_1_CG`) to extend coverage of CG loading paths.
* Added a pickle-based regression test suite (`test_reference_observables.py`, ~600+ parametrized tests) that recomputes every public `SSProtein` observable against stored reference values and asserts numerical agreement. Reference pickles are built by `tests/build_reference/build_references.py`.
* Extended `test_sstrajectory.py` with systematic parametrized tests covering properties, consistency checks, contact-map correctness, and stride behaviour across all test trajectories.
* Updated `conftest.py` with new module-scoped fixtures for the additional trajectories.
* Fixed `test_trajectory_repr_string` which was inadvertently returning a tuple rather than asserting, so the test never validated the repr string.
* `test_set_numpy_threads` now skips gracefully on platforms without a controllable BLAS backend rather than failing.

**Documentation**
* Rewrote all public docstrings in `ssprotein.py`, `sstrajectory.py`, `sssampling.py`, `ssutils.py`, `ssmutualinformation.py`, `sspre.py`, `ssnmr.py`, and `sstools.py` to follow the numpy docstring standard (one-line summary, extended description, typed Parameters, Returns/Yields with shapes, Raises, short Example), and completed the previously empty `SSPRE.__repr__` docstring.
* Added detailed overview preambles to all five module RST pages (`ssprotein`, `sstrajectory`, `ssnmr`, `sspre`, `sssampling`), describing what each module does, how its classes relate, and what categories of analysis are available.
* Rewrote `examples.rst` with nine worked IDP analysis examples covering trajectory loading, global dimensions, polymer scaling, secondary structure, contact maps, solvent accessibility, multi-chain systems, NMR comparison (PRE and random-coil chemical shifts), and sampling quality assessment with PENGUIN.
* Expanded the narrative documentation pages: `installation.rst` now has separate `pip` and `uv` instructions (PyPI and direct-from-GitHub, plus conda, editable/source installs and test instructions); `overview.rst` now covers the core `SSTrajectory`/`SSProtein` concepts and a module map; and `development.rst` now documents dev-environment setup, repository layout, running the tests and coding/docstring conventions.
* Corrected the version-check snippet (`soursop.version()` does not exist) to `python -c "import soursop; print(soursop.__version__)"` across the docs and README, refreshed the README (CI/codecov badge links now point to `holehouse-lab/soursop`, updated version/copyright lines, added a concise installation section), and fixed Sphinx build errors caused by duplicate `SSProtein` object descriptions and RST citation-label collisions (numeric `[1]_` / `.. [1]` footnotes shared across methods; resolved by adding `:no-index:` to the second autoclass block and converting all numeric footnotes to inline prose citations).

**Packaging and repository maintenance**
* Bumped the GitHub Actions CI Python matrix from 3.7-3.9 to 3.9-3.12 (3.7 and 3.8 are end-of-life).
* Removed dead CI configuration: `.travis.yml` (Travis CI retired; still referenced the old `camparitraj` name) and `.lgtm.yml` (LGTM.com was shut down).
* Reconciled `requirements.txt` and `anaconda_requirements.txt` with `pyproject.toml`: removed non-dependencies (`cx_Freeze`, `PyYAML`, `ruamel_yaml`), added the missing real dependencies (`natsort`, `matplotlib`, `cython`), and unified the `mdtraj>=1.9.5` specifier.
* Fixed a malformed `MANIFEST.in` data-file line and a stale `metapredict/_version.py` path in the `setup.cfg` coverage-omit list (now `soursop/_version.py`).

#### Update November 2024 (0.2.6)
* Major update to SOURSOP. 
* Explicit support for single-chain one-bead-per-residue trajectory loading that dramatically improves load time (~30x improvement)
* Officially added the `SSSampling` class and support for PENGUIN (see [Lotthammer & Holehouse, bioRxiv](https://www.biorxiv.org/content/10.1101/2024.11.06.622270v1.abstract)
* Added support for NH3 and FOR residue types in get_amino_acid_sequence()
* Switched over packaging to use `pyproject.toml` instead of `setup.py` and switched versioning to use `versioningit`


#### Update March 2024 (0.2.5 [patch])
* Added `return_instantaneous_maps=False` keyword to `get_distance_maps()` function so we can return a [t,n,n] matrix where t = number of frames and n=number of residues for the instantaneous conformer-specific distance maps.


#### Update Feb 2024 (0.2.5)
* Removed periodic correction in SSProtein functions. We could provide periodic boundary condition checks for some but not all functions, and ensuring this flag was possible everywhere is not feasable. With this in mind, we opted to make a design decision to remove the `periodic` flag from the small number of functions that had it, such that there's no risk of a user forgetting and analyzing two distinct properties for a trajectory that requires PBC fixing where one analysis used `periodic=True` whereas another did not have this option. Now, the user must ensure their protein trajectories are corrected ahead of time. Note that we also added support for periodic=True for some of the SSTrajectory analysis functions, because in cases where multiple chains are present reconstructing a non-PBC trajectory becomes much more difficult. As such, intermolecular analysis does provide intrinsic PBC correction, whereas intramolecular analysis does not.

#### Update July 2023
* Added in `explicit_residue_checking` flag into SSTrajectory constructor, which makes it possible to use a solvated `.gro` file as an input file.

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
