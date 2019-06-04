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

    # assert we've read in an identical file to what is expected'
    assert hashlib.sha1(str(GS6_CO.traj.xyz).encode('utf-8')).hexdigest() == '4d72177296bad848aa52a8fdc31f6ac4c6f0161d'
    assert len(GS6_CO.traj) == 5


    
