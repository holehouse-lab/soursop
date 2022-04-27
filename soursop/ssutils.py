
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


def _locate_libraries(library_name):
    # Since `threadctl` is hit or miss on a Mac (especially for the latest versions),
    # we have to do what the package does but in an ad-hoc manner using `locate`.
    # This reasonably well in principle on both the Mac and on Linux. However, the
    # requirement for `locatedb` may be a security issue for sensitive environments.
    # But, for the vast majority of other use-cases, it should be sufficient.
    import subprocess as sp
    os_name = platform.system().lower()
    if os_name == 'darwin':
        libname = f"'*{library_name}*.dylib*'" # fuzzy match for locate
    elif os_name == 'linux':
        libname = f"'*{library_name}*.so*'"  # fuzzy match for locate
    else:
        warnings.warn(f'Unsupported OS: {os_name}.')

    proc = sp.Popen(['/usr/bin/locate', library_name], stdout=sp.PIPE, stderr=sp.PIPE)
    lib_locations = list()
    for line in proc.stdout:
        loc = line.decode('utf-8').replace('\n', '').strip()
        lib_locations.append(loc)
    return lib_locations


def _identify_library_paths():
    # Identifying the right path and library to load is somewhat tricky as Python
    # environments live across multiple OSes, and can be packaged in different
    # ways. Here we limit the results to Anaconda and regular Python environment
    # installations.

    # Check existing environment variables and stop on the first match. The basis
    # for this approach is that only one should be active.
    virtualized_env = None
    for env_var in 'CONDA_DEFAULT_ENV,VIRTUAL_ENV'.split(','):
        env_path = os.environ.get(env_var, None)
        if env_path is not None:
            virtualized_env = env_path
            break

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
            if virtualized_env is not None and virtualized_env in lib_path and numpy_path_fragment in lib_path:
                candidates.append(lib_path)
            elif env_path in lib_path:
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
    return library, set_threads


def set_numpy_threads(num_threads):
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
