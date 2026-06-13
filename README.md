SOURSOP
==============================
[![Build Status](https://github.com/holehouse-lab/soursop/actions/workflows/soursop-ci.yml/badge.svg?branch=master)](https://github.com/holehouse-lab/soursop/actions)
[![codecov](https://codecov.io/gh/holehouse-lab/soursop/branch/master/graph/badge.svg)](https://codecov.io/gh/holehouse-lab/soursop)
[![Documentation Status](https://readthedocs.org/projects/soursop/badge/?version=latest)](https://soursop.readthedocs.io/en/latest/?badge=latest)

## ABOUT
SOURSOP is a Python-based simulation analysis package for working with intrinsically disordered and unfolded proteins. It is built on top of [mdtraj](https://mdtraj.org/), and was developed by Jared Lalmansingh and Alex Holehouse. 

The current stable release on PyPI is 2.0.0 (June 2026).

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

The full, versioned changelog has moved to [CHANGELOG.md](CHANGELOG.md).

The most recent release is **2.0.0 (June 2026)** - a large maintenance,
performance, documentation, and feature release that adds a consistent
ensemble-reweighting (`weights`) system across the package, two new modules for
deriving frame weights from experimental data (`ssbme`: BME / iBME / BMECustom,
and `sscoper`: COPER / iCOPER), and new experimental forward-model observables
(scalar `³J(HN, Hα)` couplings and NOE distances in `ssnmr`, plus HDX protection
factors in the new `sshdx` module), alongside wide-ranging bug fixes and
behaviour-preserving speed-ups. See [CHANGELOG.md](CHANGELOG.md) for the complete
list and all earlier releases.

#### Acknowledgements
Project based on the
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.0.
