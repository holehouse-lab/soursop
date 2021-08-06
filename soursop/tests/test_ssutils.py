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


def test_mkl_set_num_threads():
    # Determine the BLAS implementation in use - e.g. MKL (Intel), OpenBLAS, etc.
    info = threadpool_info()

    num_threads = 2
    mkl_path = None
    for lib in info:
        filepath = lib['filepath']
        base_filepath = os.path.basename(filepath)
        if 'mkl' in base_filepath:
            mkl_path = filepath[:]
            break

    if mkl_path is not None:
        # Traditional UNIX-like systems will have shared objects available.
        if 'bsd' in sys.platform or 'lin' in sys.platform:
            mkl_rt = ctypes.CDLL(mkl_path)
        
        # Darwin / Apple uses `*.dylib` by default for included Intel compiler libraries.
        # Traditional UNIX-like shared objects can be created (`*.so`), but are more
        # represented in third-party libraries. This is a more dynamic way of finding
        # the MKL library and using it on a Mac that has an Intel compiler installed.
        elif sys.platform == 'darwin':
            mkl_rt = ctypes.CDLL(mkl_path)
        ssutils.mkl_set_num_threads(num_threads)
        num_threads_set = mkl_rt.mkl_get_max_threads()
        assert num_threads_set == num_threads
    else:
        # While it may seem to be a good idea to raise an error at this point,
        # doing so puts the user in a catch-22, therefor no error should be
        # generated. A warning should be printed.
        print('Warning: No Intel BLAS implementation found.')


def test_set_numpy_threads():
    num_threads = 2
    blas_library, set_threads = ssutils.set_numpy_threads(num_threads)
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
