"""
Unit and regression test for the camparitraj package.
"""
# Import package, test suite, and other packages as needed
import hashlib
import camparitraj
from camparitraj import cttrajectory
import pytest
import sys


from . import conftest


def test_camparitraj_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "camparitraj" in sys.modules


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




def test_read_in_trajectory(GS6_CO):        
    assert len(GS6_CO.traj) == 5

