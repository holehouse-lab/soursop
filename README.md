# SOURSOP

*A Python package for the analysis of simulations of intrinsically disordered and unfolded proteins.*

[![Build Status](https://github.com/holehouse-lab/soursop/actions/workflows/soursop-ci.yml/badge.svg?branch=master)](https://github.com/holehouse-lab/soursop/actions)
[![codecov](https://codecov.io/gh/holehouse-lab/soursop/branch/master/graph/badge.svg)](https://codecov.io/gh/holehouse-lab/soursop)
[![Documentation Status](https://readthedocs.org/projects/soursop/badge/?version=latest)](https://soursop.readthedocs.io/en/latest/?badge=latest)
[![PyPI version](https://img.shields.io/pypi/v/soursop.svg)](https://pypi.org/project/soursop/)
[![Python versions](https://img.shields.io/badge/python-3.9%20%7C%203.10%20%7C%203.11%20%7C%203.12%20%7C%203.13%20%7C%203.14-blue.svg)](https://pypi.org/project/soursop/)
[![License: LGPL v3](https://img.shields.io/badge/license-LGPLv3-blue.svg)](LICENSE)
[![DOI](https://img.shields.io/badge/DOI-10.1021%2Facs.jctc.3c00190-blue.svg)](https://doi.org/10.1021/acs.jctc.3c00190)
[![Last commit](https://img.shields.io/github/last-commit/holehouse-lab/soursop.svg)](https://github.com/holehouse-lab/soursop/commits/master)

---

## Overview

**SOURSOP** is a Python-based simulation analysis package built for the
conformational analysis of intrinsically disordered regions (IDRs), unfolded
states, and other flexible biopolymers. It is built on top of
[MDTraj](https://mdtraj.org/), which handles trajectory I/O and the low-level
atomic representation, and adds an analysis layer of polymer-physics-aware
observables designed specifically for disordered ensembles.

SOURSOP was originally developed by Jared Lalmansingh in the [Pappu lab](https://pappulab.wustl.edu/) and Alex Holehouse in the [Holehouse Lab](https://www.holehouselab.com/) at Washington University in St. Louis.

## Features

- **Global dimensions & shape** — radius of gyration, hydrodynamic radius,
  end-to-end distance, asphericity, acylindricity, prolateness, and the full
  gyration tensor.
- **Polymer scaling** — internal scaling profiles, the apparent scaling
  exponent ν (with bootstrap confidence intervals), and local scaling
  heterogeneity.
- **Distances & contacts** — inter-residue and inter-atomic distance maps,
  polymer-scaled distance maps, contact maps, and fraction of native contacts
  (*Q*).
- **Secondary structure** — per-frame DSSP assignments and BBSEG
  backbone-torsion classification.
- **Solvent accessibility** — per-residue, regional, and site-level SASA.
- **NMR observables** — random-coil chemical shifts, ³J(HN, Hα) scalar
  couplings, NOE distances (`ssnmr`), and paramagnetic relaxation enhancement
  (`sspre`).
- **Ensemble reweighting** — a consistent, deterministic per-frame `weights`
  system across the package, with Bayesian Maximum Entropy (`ssbme`: BME / iBME
  / BMECustom) and Convex Optimization for Ensemble Reweighting (`sscoper`:
  COPER / iCOPER) for reweighting ensembles against experimental data.
- **HDX protection factors** — Best–Vendruscolo ln(P) predictions (`sshdx`).
- **Sampling diagnostics** — convergence assessment of disordered-protein
  ensembles via PENGUIN (`sssampling`).
- **Multiple resolutions** — all-atom and one-bead-per-residue coarse-grained trajectories, with automatic detection.

## Installation

SOURSOP can be installed from PyPI with either `pip` or
[uv](https://docs.astral.sh/uv/):

```bash
pip install soursop
# or
uv pip install soursop
```

To install the latest development version directly from GitHub:

```bash
pip install "git+https://github.com/holehouse-lab/soursop.git"
```

Verify the installation with:

```bash
python -c "import soursop; print(soursop.__version__)"
```

Full installation instructions (conda, editable/source installs, and running the
tests) are in the [documentation](https://soursop.readthedocs.io/).

## Quickstart

```python
import numpy as np
from soursop.sstrajectory import SSTrajectory

# read in a trajectory (trajectory file + topology file)
traj = SSTrajectory('traj.xtc', 'start.pdb')

# extract the first protein chain (an SSProtein object)
protein = traj.proteinTrajectoryList[0]

# ensemble-average radius of gyration and end-to-end distance
rg  = np.mean(protein.get_radius_of_gyration())
ree = np.mean(protein.get_end_to_end_distance())

# ensemble-average inter-residue distance map
dmap = protein.get_distance_map()

print(f"Rg  = {rg:.2f} Å")
print(f"Ree = {ree:.2f} Å")
```

See the [worked examples](https://soursop.readthedocs.io/en/latest/usage/examples.html)
in the documentation for end-to-end analyses.

## Documentation

Full documentation, including installation, tutorials, worked examples, and the
complete API reference, is hosted at
**[soursop.readthedocs.io](https://soursop.readthedocs.io/)**.

## Versioning and changelog

The current PyPI release is **0.2.7**. The upcoming **2.0.0** release (in
development on this repository) is a large maintenance, performance,
documentation, and feature release that adds a consistent ensemble-reweighting
(`weights`) system across the package, two new modules for deriving frame
weights from experimental data (`ssbme`: BME / iBME / BMECustom, and `sscoper`:
COPER / iCOPER), new experimental forward-model observables (scalar
`³J(HN, Hα)` couplings and NOE distances in `ssnmr`, plus HDX protection factors
in the new `sshdx` module), and first-class support for SWAN two-bead
coarse-grained models — alongside wide-ranging bug fixes and
behaviour-preserving speed-ups.

The full, versioned changelog is in [CHANGELOG.md](CHANGELOG.md).

## Reporting bugs and requesting features

If you find a bug, typo, or error,
[please raise an issue on GitHub](https://github.com/holehouse-lab/soursop/issues).

If you wish to add a new feature or contribute a plugin, please see the
[development documentation](https://soursop.readthedocs.io/en/latest/usage/development.html).

## Citing SOURSOP

If you use SOURSOP in your work, please cite:

> Lalmansingh, J. M., Keeley, A. T., Ruff, K. M., Pappu, R. V. & Holehouse, A. S.
> **SOURSOP: A Python Package for the Analysis of Simulations of Intrinsically
> Disordered Proteins.** *J. Chem. Theory Comput.* **19**, 5609–5620 (2023).
> doi:[10.1021/acs.jctc.3c00190](https://doi.org/10.1021/acs.jctc.3c00190)

- [Journal link](https://pubs.acs.org/doi/full/10.1021/acs.jctc.3c00190)
- [Paper PDF](https://www.dropbox.com/s/bd5szapvxpn83r6/soursop_jctc.pdf?dl=0)

## License

SOURSOP is distributed under the GNU Lesser General Public License v3.0
(LGPLv3). See [LICENSE](LICENSE) for the full text.

Copyright © 2014–2026 Alex Holehouse and contributors.

## Acknowledgements

Project structure based on the
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms)
version 1.0.
