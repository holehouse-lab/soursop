"""
Unit and regression test for the soursop package.
"""

# Import package, test suite, and other packages as needed
import numpy as np
import soursop
import pytest
import sys
import random
import itertools
from copy import deepcopy
from soursop import sstrajectory, ssprotein
from soursop.sssampling import SamplingQuality, glob_traj_paths
from soursop.ssexceptions import SSException
from soursop.configs import DEBUGGING
from contextlib import contextmanager


def test_compute_pdf(qual_assess):
    bins = qual_assess.get_degree_bins()
    test_traj = qual_assess.trajs[0]
    test_pol_traj = qual_assess.polymer_trajs[0]

    test_psi = test_traj.proteinTrajectoryList[0].get_angles("psi")[1]
    test_pol_psi = test_pol_traj.proteinTrajectoryList[0].get_angles("psi")[1]

    test_pdf = np.apply_along_axis(lambda col: np.histogram(col, bins=bins, density=True)[0], axis=1, arr=test_psi)
    test_pol_pdf = np.apply_along_axis(lambda col: np.histogram(col, bins=bins, density=True)[0], axis=1, arr=test_pol_psi)
    assert np.all(test_pdf == qual_assess.compute_pdf(qual_assess.psi_angles[0],bins=bins))
    assert np.all(test_pol_pdf == qual_assess.compute_pdf(qual_assess.polymer_psi_angles[0], bins=bins))

