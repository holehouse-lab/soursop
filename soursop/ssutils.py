
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
import platform
import warnings
from . ssexceptions import SSException
from threadpoolctl import threadpool_info, threadpool_limits


MKL_LIBRARY = 'mkl_rt'
OPENBLAS_LIBRARY = 'openblas'


def _set_mkl_numpy_threads(mkl_path, num_threads):
    # Traditional UNIX-like systems will have shared objects available.
    #
    # Darwin / Apple uses `*.dylib` by default for included Intel compiler libraries.
    # Traditional UNIX-like shared objects can be created (`*.so`), but are more
    # represented in third-party libraries. This is a more dynamic way of finding
    # the MKL library and using it on a Mac that has an Intel compiler installed.
    mkl_rt = ctypes.CDLL(mkl_path)
    mkl_get_max_threads = mkl_rt.mkl_get_max_threads()
    mkl_rt.mkl_set_num_threads(ctypes.byref(ctypes.c_int(num_threads)))
    set_threads = mkl_rt.mkl_get_max_threads()
    return set_threads


def _set_openblas_numpy_threads(openblas_path, num_threads):
    # Most BLAS implementations use OpenBLAS, this section will likely be used the most
    # in production.
    openblas_lib = ctypes.cdll.LoadLibrary(openblas_path)
    openblas_lib.openblas_set_num_threads(num_threads)
    set_threads = openblas_lib.openblas_get_num_threads()
    return set_threads


def _locate_libraries(library_name):
    import fnmatch
    # Since `threadctl` is hit or miss on a Mac (especially for the latest versions),
    # we implement a custom finder that examines the virtual environment to find
    # library path candidates.
    os_name = platform.system().lower()
    if os_name == 'darwin':
        libname = f'*{library_name}*.dylib*' # fuzzy match for filtering with find
    elif os_name == 'linux':
        libname = f'*{library_name}*.so*'  # fuzzy match for filtering with find
    else:
        warnings.warn(f'Unsupported OS: {os_name}.')

    # Checking existing environment variables and stop on the first match. The
    # basis for this approach is that only one should be active.
    virtualized_env = None
    for env_var in 'CONDA_PREFIX,VIRTUAL_ENV'.split(','):
        env_path = os.environ.get(env_var, None)
        if env_path is not None:
            virtualized_env = env_path
            break

    # If no virtual environment is active, terminate. Support for system-wide
    # Python installations is not yet supported.
    if virtualized_env is None:
        raise SSException('No Anaconda or Python Virtual Environment found. Exiting.')

    lib_locations = list()
    include_filenames = [libname]
    for root, dirs, files in os.walk(virtualized_env, topdown=True):
        for filename_pattern in include_filenames:
            for filename in fnmatch.filter(files, filename_pattern):
                filepath = os.path.join(root, filename)
                lib_locations.append(filepath)
    return lib_locations


def _identify_library_paths():
    # Identifying the right path and library to load is somewhat tricky as Python
    # environments live across multiple OSes, and can be packaged in different
    # ways. Here we limit the results to Anaconda and regular Python environment
    # installations.

    # Identify the library paths and split them into two:
    # 1) candidates - these will be searched and attempted to be set first.
    # 2) other candidates - backup library paths to check & set.
    libraries = [MKL_LIBRARY, OPENBLAS_LIBRARY]
    candidates = list()
    other_candidates = list()
    for library in libraries:
        lib_paths = _locate_libraries(library)
        for lib_path in lib_paths:
            numpy_path_fragment = os.path.join('site-packages', 'numpy') # os-agnostic
            if numpy_path_fragment in lib_path:
                candidates.append(lib_path)
            else:
                other_candidates.append(lib_path)
    return candidates, other_candidates


def _set_numpy_threads(candidate_library_paths, num_threads):
    # This shadows the main entry point for setting the numpy threads. The libraries
    # are set automatically on most *nix OSes.
    set_threads = 0
    library = None
    for lib_path in candidate_library_paths:
        if MKL_LIBRARY in lib_path:
            set_threads = _set_mkl_numpy_threads(lib_path, num_threads)
            library = MKL_LIBRARY
        elif OPENBLAS_LIBRARY in lib_path:
            set_threads = _set_openblas_numpy_threads(lib_path, num_threads)
            library = OPENBLAS_LIBRARY
        else:
            warnings.warn('Unsupported library. Please install OPENBLAS or the Intel MKL library. No threads set.')
            library = 'unknown'
        break # stop on first set library
    return set_threads, library


##
## This is included the force numpy to use defined number of cores. For
## some of the linear algebra routines numpy will default to using as many
## cores as it can get its greedy little hands on - this function allows that
## thirst to be quenched...
##
def set_numpy_threads(num_threads):
    # Currently only MKL is supported on Windows as it's installed alongside
    # the other packages via conda. A "traditional" virtual environment requires
    # access to a compiler and other libraries for successful compilation.
    if platform.system().lower() == 'windows':
        import mkl
        mkl.set_num_threads(num_threads)
        return mkl.get_max_threads(), MKL_LIBRARY

    candidates, other_candidates = _identify_library_paths()
    if len(candidates) > 0:
        set_threads, library = _set_numpy_threads(candidates, num_threads)
    else:
        set_threads, library = _set_numpy_threads(other_candidates, num_threads)
    return set_threads, library


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
