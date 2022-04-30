"""
Unit and regression test for the sstrajectory module.
"""
# Import package, test suite, and other packages as needed
import os
import sys
import pytest
import ctypes
import numpy as np
from soursop import ssutils
from soursop.ssexceptions import SSException
from threadpoolctl import threadpool_info, threadpool_limits


def test_set_numpy_threads():
    num_threads = 2
    set_threads, blas_library = ssutils.set_numpy_threads(num_threads)
    assert blas_library != 'unknown'
    assert set_threads == num_threads


def test_validate_keyword_option():
    allowed_modes = ['COM', 'CA']
    for mode in allowed_modes:
        ssutils.validate_keyword_option(mode, allowed_modes, 'mode')
        ssutils.validate_keyword_option(mode, allowed_modes, 'mode', 'Unsupported mode.')

    # Attempt a test which will fail - i.e. mismatched keyword.
    with pytest.raises(SSException):
        ssutils.validate_keyword_option('invalid_mode', allowed_modes, 'mode')

    # Attempt a test which will fail, but including a custom string error message.
    with pytest.raises(SSException):
        ssutils.validate_keyword_option('invalid_mode', allowed_modes, 'mode', 'Unsupported mode.')

    # Attempt a test which will fail, including a non-string error message.
    with pytest.raises(RuntimeError):
        ssutils.validate_keyword_option('invalid_mode', allowed_modes, 'mode', 123)
