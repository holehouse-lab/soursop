"""
Unit and regression test for the cttrajectory module.
"""
# Import package, test suite, and other packages as needed
import camparitraj
import hashlib
from camparitraj import cttrajectory
from camparitraj.ctexceptions import CTException
from pathlib import Path
import pytest
import sys
import random
import itertools
import platform
import tempfile
import numpy as np
import os
from camparitraj.configs import TMP_DIR


from . import conftest
#import conftest

# -------------------------------------------------------------------------------------------------
# CTTrajectory `__init__` tests.
#
# Also implicitly tests:
#       `__readTrajectory`
#       `__get_proteins`
#       `__get_proteins_by_residue`




def test_get_overal_rh(GS6_CO):
    assert np.allclose(GS6_CO.get_overall_hydrodynamic_radius(), GS6_CO.proteinTrajectoryList[0].get_hydrodynamic_radius())

def test_get_overal_asphericity(GS6_CO):
    assert np.allclose(GS6_CO.get_overall_asphericity(), GS6_CO.proteinTrajectoryList[0].get_asphericity())

def test_get_overal_rg(GS6_CO):
    assert np.allclose(GS6_CO.get_overall_radius_of_gyration(), GS6_CO.proteinTrajectoryList[0].get_radius_of_gyration())

def test_get_overal_rg(GMX_2CHAINS):
    # compares rg of 2 chains vs. same values calculated by VMD

    # rg calculated in vmd
    vmd_rg = np.array([24.497053146362305, 24.938467025756836, 25.010461807250977, 27.36289405822754, 26.659685134887695, 24.56108283996582, 24.705211639404297, 27.144346237182617, 25.043672561645508, 27.67377471923828, 48.036869049072266, 23.930212020874023, 24.449384689331055, 24.491701126098633, 25.248193740844727, 27.898420333862305, 47.11823272705078, 47.213069915771484, 49.21323776245117, 27.429519653320312])

    assert np.allclose(GMX_2CHAINS.get_overall_radius_of_gyration(), vmd_rg)


def test_get_overal_asphericity(GMX_2CHAINS):
    # compares rg of 2 chains vs. same values calculated by VMD

    asph = np.array([0.36274249, 0.4079101, 0.5837932, 0.3413606, 0.30022731, 0.25667767, 0.24720234, 0.32076792, 0.26974903, 0.21184996, 0.62855842, 0.23021391, 0.16662495, 0.45061943, 0.54406795, 0.4482196, 0.53423135, 0.73352191, 0.75562579, 0.32473148])

    assert np.allclose(GMX_2CHAINS.get_overall_asphericity(), asph)

def test_get_overal_rh(GMX_2CHAINS):
    # compares rg of 2 chains vs. same values calculated by VMD

    asph = np.array([26.99809304, 27.25590015, 27.29754026, 28.59875282, 28.22148641, 27.03576021, 27.12021072, 28.48253016, 27.31671616, 28.76251356, 36.57647347, 26.66061797, 26.96999394, 26.99493729, 27.43425109, 28.87969746, 36.31506845, 36.3423472, 36.90213999, 28.63399989])

    assert np.allclose(GMX_2CHAINS.get_overall_hydrodynamic_radius(), asph)

