"""
Unit and regression test for the cttrajectory module.
"""
# Import package, test suite, and other packages as needed
import camparitraj
import hashlib
from camparitraj import cttrajectory
import pytest
import sys


from . import conftest


def test_read_in_trajectory(GS6_CO):        

    # ensure 
    assert len(GS6_CO.traj) == 5
    assert len(GS6_CO.proteinTrajectoryList) == 1


    
