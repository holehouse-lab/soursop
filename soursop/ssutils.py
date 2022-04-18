
##     _____  ____  _    _ _____   _____  ____  _____  
##   / ____|/ __ \| |  | |  __ \ / ____|/ __ \|  __ \ 
##  | (___ | |  | | |  | | |__) | (___ | |  | | |__) |
##   \___ \| |  | | |  | |  _  / \___ \| |  | |  ___/ 
##   ____) | |__| | |__| | | \ \ ____) | |__| | |     
##  |_____/ \____/ \____/|_|  \_\_____/ \____/|_|     

## Alex Holehouse (Pappu Lab and Holehouse Lab) and Jared Lalmansing (Pappu lab)
## Simulation analysis package
## Copyright 2014 - 2022
##


import os
import sys
import numpy
import ctypes
from . ssexceptions import SSException
from threadpoolctl import threadpool_info, threadpool_limits


##
## This is included the force numpy to use defined number of cores. For
## some of the linear algebra routines numpy will default to using as many
## cores as it can get its greedy little hands on - this function allows that
## thirst to be quenched...
##

def mkl_set_num_threads(cores):
    mkl_rt = ctypes.CDLL('libmkl_rt.so')
    mkl_get_max_threads = mkl_rt.mkl_get_max_threads()
    mkl_rt.mkl_set_num_threads(ctypes.byref(ctypes.c_int(cores)))


def _set_mkl_numpy_threads(mkl_path, num_threads):

    # Traditional UNIX-like systems will have shared objects available.
    if 'bsd' in sys.platform or 'lin' in sys.platform:
        mkl_rt = ctypes.CDLL(mkl_path)
    
    # Darwin / Apple uses `*.dylib` by default for included Intel compiler libraries.
    # Traditional UNIX-like shared objects can be created (`*.so`), but are more
    # represented in third-party libraries. This is a more dynamic way of finding
    # the MKL library and using it on a Mac that has an Intel compiler installed.
    elif sys.platform == 'darwin':
        mkl_rt = ctypes.CDLL(mkl_path)
    mkl_set_num_threads(num_threads)
    set_threads = mkl_rt.mkl_get_max_threads()
    return set_threads


def _set_openblas_numpy_threads(openblas_path, num_threads):
    # Most BLAS implementations use OpenBLAS, this section will likely be used the most
    # in production.
    openblas_lib = ctypes.cdll.LoadLibrary(openblas_path)
    openblas_lib.openblas_set_num_threads(num_threads)
    set_threads = openblas_lib.openblas_get_num_threads()
    return set_threads


def set_numpy_threads(num_threads):
    # Determine the BLAS implementation in use - e.g. MKL (Intel), OpenBLAS, etc.
    info = threadpool_info()

    # Set the threads based on the library available.
    set_threads = 0
    blas_implementation = None

    for lib in info:
        filepath = lib['filepath']
        base_filepath = os.path.basename(filepath)
        
        if 'mkl' in base_filepath:
            blas_implementation = 'mkl'
            set_threads = _set_mkl_numpy_threads(filepath, num_threads)
            break
        
        elif 'openblas' in base_filepath:
            blas_implementation = 'openblas'
            set_threads = _set_openblas_numpy_threads(filepath, num_threads)
            break

        else:
            blas_implementation = 'unknown'
    return blas_implementation, set_threads


def validate_keyword_option(keyword, allowed_vals, keyword_name, error_message=None):
    """
    Helper function that checks a passed keyword is only one of a set of possible
    valid keywords

    Parameters
    -----------

    keyword : str
        The actual passed keyword value

    allowed_vals : list of str
        A list of possible keywords

    keyword_name : str
        the name of the keyword as the user would select it in the function call

    error_message : str
        Allows the user to pass a custom error message. Default is None.


    Returns
    --------
    None

        No return value, but raises ssexceptions.SSException if ``keyword `` is not
        found in the allowed_vals list

    """


    if keyword not in allowed_vals:
        message = None
        if error_message is None:
            message = 'Keyword %s passed value [%s], but this is not valid.\nMust be one of :%s'  % (keyword_name, keyword, ", ".join(allowed_vals))
        else:
            error_type = type(error_message)
            if error_type is not str:
                raise RuntimeError('Invalid error message type: "{}". The error message must be a string.'.format(error_type))
            message = error_message[:]
        raise SSException(message)
