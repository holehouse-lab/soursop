"""
Unit and regression test for the CTraj package.
"""

# Import package, test suite, and other packages as needed
import CTraj
import pytest
import sys

def test_CTraj_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "CTraj" in sys.modules
