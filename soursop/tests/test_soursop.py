"""
Unit and regression test for the soursop package.
"""
# Import package, test suite, and other packages as needed
import hashlib
import soursop
from soursop import sstrajectory
import pytest
import sys


from . import conftest


def test_soursop_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "soursop" in sys.modules


def test_validate_tests(GS6_CO, NTL9_CO):
    """
    key function that ensures we can read in the input trajectories and they _exactly_ match what
    we expect them to be - if this fails no reason to expect anything else to work. This is basically
    a validation/verification step for the test suite. This basically provides a high-resolution way to
    ensure we know what data is being fed into the tests, and avoids the need to manually populated
    this data from within the test-file (i.e. we're effectively storing the test data in a binary, 
    reading it in, and validating)

    """

    assert hashlib.sha1(str(GS6_CO.traj.xyz).encode('utf-8')).hexdigest() == '4d72177296bad848aa52a8fdc31f6ac4c6f0161d'
    assert hashlib.sha1(str(NTL9_CO.traj.xyz).encode('utf-8')).hexdigest() == 'fa3214bfa20c767c7463b4f787057a2ab7f44361'


def test_read_in_trajectory(GS6_CO, NTL9_CO):
    """
    This function tests whether we can access the underlying `mdtraj` trajectory, and that it has been properly loaded
    by checking the values of the attributes.
    """
    assert GS6_CO.traj.n_frames == 5
    assert GS6_CO.traj.n_atoms == 66
    assert GS6_CO.traj.n_residues == 8

    assert NTL9_CO.traj.n_frames == 10
    assert NTL9_CO.traj.n_atoms == 908
    assert NTL9_CO.traj.n_residues == 56
