"""
Unit and regression test for the camparitraj package.
"""
# Import package, test suite, and other packages as needed
import camparitraj
from camparitraj import cttrajectory
import pytest
import sys


from . import conftest


def test_camparitraj_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "camparitraj" in sys.modules


def test_read_in_trajectory(GS6_CO):        
    assert len(GS6_CO.traj) == 5

