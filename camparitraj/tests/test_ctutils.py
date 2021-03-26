"""
Unit and regression test for the cttrajectory module.
"""
# Import package, test suite, and other packages as needed
import os
import sys
import ctypes
import numpy as np
import camparitraj.ctutils as ctutils
from camparitraj.ctexceptions import CTException
from threadpoolctl import threadpool_info, threadpool_limits


def test_mkl_set_num_threads():
    # Determine the BLAS implementation in use - e.g. MKL (Intel), OpenBLAS, etc.
    info = threadpool_info()

    num_threads = 2
    mkl_path = None
    openblas_path = None
    other = None
    for lib in info:
        filepath = lib['filepath']
        base_filepath = os.path.basename(filepath)
        if 'mkl' in base_filepath:
            mkl_path = filepath[:]
        elif 'openblas' in base_filepath:
            openblas_path = filepath[:]
        else:
            other = filepath[:]

    if mkl_path is not None:
        # Traditional UNIX-like systems will have shared objects available.
        if 'bsd' in sys.platform or 'lin' in sys.platform:
            mkl_rt = ctypes.CDLL('libmkl_rt.so')
        
        # Darwin / Apple uses `*.dylib` by default for included Intel compiler libraries.
        # Traditional UNIX-like shared objects can be created (`*.so`), but are more
        # represented in third-party libraries. This is a more dynamic way of finding
        # the MKL library and using it on a Mac that has an Intel compiler installed.
        elif sys.platform == 'darwin':
            mkl_rt = ctypes.CDLL(mkl_path)
        ctutils.mkl_set_num_threads(num_threads)
        num_threads_set = mkl_rt.mkl_get_max_threads
        assert num_threads_set == num_threads

    # Most BLAS implementations use OpenBLAS, this section will likely be used the most
    # in production.
    elif openblas_path is not None:
        openblas_lib = ctypes.cdll.LoadLibrary(openblas_path)
        openblas_lib.openblas_set_num_threads(num_threads)
        num_threads_set = openblas_lib.openblas_get_num_threads()
        assert num_threads_set == num_threads

    elif other is None:
        raise RuntimeError('No supported libraries found.')


def test_validate_keyword_option():
    allowed_modes = ['COM', 'CA']
    for mode in allowed_modes:
        ctutils.validate_keyword_option(mode, allowed_modes, 'mode')
        ctutils.validate_keyword_option(mode, allowed_modes, 'mode', 'Unsupported mode.')

    # Attempt a test which will fail - i.e. mismatched keyword.
    keyword_not_found = False
    try:
        ctutils.validate_keyword_option('invalid_mode', allowed_modes, 'mode')
    except CTException as e:
        keyword_not_found = True
    assert keyword_not_found == True

    # Attempt a test which will fail, but including a custom string error message.
    keyword_not_found = False
    try:
        ctutils.validate_keyword_option('invalid_mode', allowed_modes, 'mode', 'Unsupported mode.')
    except CTException as e:
        keyword_not_found = True
    assert keyword_not_found == True


    # Attempt a test which will fail, including a non-string error message.
    keyword_not_found = False
    try:
        ctutils.validate_keyword_option('invalid_mode', allowed_modes, 'mode', 123)
    except RuntimeError as e:
        keyword_not_found = True
    assert keyword_not_found == True
