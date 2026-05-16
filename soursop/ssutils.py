
##     _____  ____  _    _ _____   _____  ____  _____
##   / ____|/ __ \| |  | |  __ \ / ____|/ __ \|  __ \
##  | (___ | |  | | |  | | |__) | (___ | |  | | |__) |
##   \___ \| |  | | |  | |  _  / \___ \| |  | |  ___/
##   ____) | |__| | |__| | | \ \ ____) | |__| | |
##  |_____/ \____/ \____/|_|  \_\_____/ \____/|_|

## Alex Holehouse (Pappu Lab and Holehouse Lab) and Jared Lalmansing (Pappu lab)
## Simulation analysis package
## Copyright 2014 - 2026
##


import os
import sys
import numpy
import ctypes
import platform
import warnings
from soursop.ssexceptions import SSException
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
    """Set the number of threads NumPy's BLAS backend is allowed to use.

    Some of NumPy's linear-algebra routines will grab every CPU core they can
    find by default. That is rarely what you want when you are already
    running many SOURSOP analyses in parallel (or sharing a node), so this
    helper provides a single place to clamp BLAS thread use.

    Implementation differs across platforms:

    * **Windows**: uses the ``mkl`` Python package (installed via conda).
    * **macOS / Linux**: locates the MKL or OpenBLAS shared library inside
      the current virtual environment (or conda env) and calls the C-level
      thread-setter via ``ctypes``.

    Some BLAS backends (notably Apple's Accelerate framework) do not expose
    a thread-control API and will raise :class:`SSException` here.

    Parameters
    ----------
    num_threads : int
        Maximum number of threads to allow. Must be a positive integer.

    Returns
    -------
    tuple of (int, str)
        ``(set_threads, library)`` where ``set_threads`` is the actual
        thread count BLAS now reports (which may differ from
        ``num_threads`` if the BLAS implementation clamps it), and
        ``library`` is one of ``'mkl_rt'``, ``'openblas'``, or
        ``'unknown'``.

    Raises
    ------
    SSException
        If no MKL or OpenBLAS library can be found in the active
        environment (typical of Apple Accelerate setups).

    Example
    -------
    >>> from soursop.ssutils import set_numpy_threads
    >>> set_numpy_threads(2)
    (2, 'openblas')
    """
    # Currently only MKL is supported on Windows as it's installed alongside
    # the other packages via conda. A "traditional" virtual environment requires
    # access to a compiler and other libraries for successful compilation.
    if platform.system().lower() == 'windows':
        import mkl
        mkl.set_num_threads(num_threads)
        return mkl.get_max_threads(), MKL_LIBRARY

    candidates, other_candidates = _identify_library_paths()
    if len(candidates) == 0 and len(other_candidates) == 0:
        raise SSException(
            "No MKL or OpenBLAS library found in the current environment. "
            "Thread count control is not available for this BLAS backend (e.g., Apple Accelerate)."
        )
    if len(candidates) > 0:
        set_threads, library = _set_numpy_threads(candidates, num_threads)
    else:
        set_threads, library = _set_numpy_threads(other_candidates, num_threads)
    return set_threads, library


def validate_keyword_option(keyword, allowed_vals, keyword_name, error_message=None):
    """Raise :class:`SSException` unless ``keyword`` is one of the allowed values.

    Used throughout SOURSOP to validate ``mode`` / ``scheme`` / similar
    string options at the top of public methods. If ``error_message`` is
    not supplied a standard message naming the offending keyword and
    listing the allowed values is constructed automatically.

    Parameters
    ----------
    keyword : str
        The value the caller passed.
    allowed_vals : list of str
        The complete set of accepted values.
    keyword_name : str
        Human-readable name of the parameter (for the error message).
    error_message : str, optional
        Custom message to use instead of the default. Must be a string.

    Returns
    -------
    None

    Raises
    ------
    SSException
        If ``keyword`` is not in ``allowed_vals``.
    RuntimeError
        If ``error_message`` is supplied as a non-string type.

    Example
    -------
    >>> from soursop.ssutils import validate_keyword_option
    >>> validate_keyword_option('CA', ['CA', 'COM'], 'mode')   # passes silently
    >>> validate_keyword_option('xyz', ['CA', 'COM'], 'mode')  # raises SSException
    """


    if keyword not in allowed_vals:
        message = None
        if error_message is None:
            message = f'Keyword {keyword_name} passed value [{keyword}], but this is not valid.\nMust be one of :%s'  % (", ".join(allowed_vals))
        else:
            error_type = type(error_message)
            if error_type is not str:
                raise RuntimeError('Invalid error message type: "{}". The error message must be a string.'.format(error_type))
            message = error_message[:]
        raise SSException(message)


def validate_weights(weights, n_frames, stride=1, etol=1e-7):
    """Validate (and stride-normalise) a per-frame re-weighting vector.

    This is the single, shared entry point used by every SOURSOP function
    that accepts a ``weights`` keyword. It enforces that the weights form
    a proper per-frame probability vector so that all downstream
    deterministic weighted averages (:func:`weighted_mean`,
    :func:`weighted_std`, :func:`weighted_rms`, :func:`weighted_corr`) are
    well defined and consistent.

    The literal ``False`` (or ``None``) is the "no weighting" sentinel and
    is passed straight through, so a default of ``weights=False`` is a
    strict no-op.

    Validation (each failure raises :class:`SSException`):

    1. ``False`` / ``None`` -> returned unchanged (no-op).
    2. castable to a 1-D ``numpy.float64`` array.
    3. all elements finite (no ``nan`` / ``inf``).
    4. exactly one weight per frame (``len == n_frames``).
    5. every element in the closed interval ``[0, 1]``.
    6. if ``stride > 1``: the vector is subsampled (``weights[::stride]``)
       and **renormalised to sum to 1** so the strided weighted average is
       still a proper expectation (this fixes the historical behaviour
       where strided weights silently no longer summed to 1).
    7. ``|sum(weights) - 1| < etol`` (checked after any stride/renormalise).

    Parameters
    ----------
    weights : array_like, False or None
        Per-frame weights, or the ``False``/``None`` no-op sentinel.
    n_frames : int
        Number of frames the (unstrided) weight vector must match.
    stride : int, optional
        Frame stride that will be applied to the trajectory. Default 1.
    etol : float, optional
        Tolerance on ``|sum(weights) - 1|``. Default ``1e-7``.

    Returns
    -------
    numpy.ndarray or False
        ``False`` if the input was ``False``/``None``; otherwise a
        ``numpy.float64`` array of length ``ceil(n_frames / stride)`` that
        sums to 1 within ``etol``.

    Raises
    ------
    SSException
        If any of the validation conditions above fail.

    Example
    -------
    >>> import numpy as np
    >>> from soursop.ssutils import validate_weights
    >>> validate_weights(False, 10) is False
    True
    >>> w = validate_weights(np.full(10, 0.1), 10)
    >>> float(w.sum())
    1.0
    """

    if weights is False or weights is None:
        return False

    try:
        w = numpy.asarray(weights, dtype=numpy.float64)
    except (ValueError, TypeError) as e:
        raise SSException(
            "Unable to convert the passed weights to a numeric "
            f"numpy.float64 array (likely non-numerical input):\n{weights}\n{e}"
        )

    if w.ndim != 1:
        raise SSException(
            f"Frame weights must be a 1-D vector, got an array with shape {w.shape}"
        )

    if not numpy.all(numpy.isfinite(w)):
        raise SSException("Frame weights contain non-finite values (nan/inf)")

    if len(w) != n_frames:
        raise SSException(
            f"Passed frame weights array is {len(w)} in length, while there "
            f"are actually {n_frames} frames - these must match"
        )

    if numpy.any(w < 0.0) or numpy.any(w > 1.0):
        raise SSException(
            "Every frame weight must lie in the closed interval [0, 1] "
            f"(min={w.min():g}, max={w.max():g})"
        )

    if stride > 1:
        from soursop import ssio
        ssio.warning_message(
            "WARNING: Using stride with weights is ALMOST certainly not a good "
            "idea unless the weights are\ncalculated for every stride-th frame. "
            "The strided weights will be renormalised to sum to 1.",
            with_frills=True,
        )
        w = w[::stride]
        wsum = numpy.sum(w)
        if wsum <= 0.0:
            raise SSException(
                "After applying the frame stride the remaining weights sum to "
                f"{wsum:g} (<= 0); cannot renormalise"
            )
        w = w / wsum

    abs_diff = abs(numpy.sum(w) - 1.0)
    if abs_diff >= etol:
        raise SSException(
            "The passed weights do not sum to 1 within the specified floating "
            f"point tolerance (etol={etol:g}). | sum(weights) - 1 | = {abs_diff:g}"
        )

    return w


def weighted_mean(x, weights, axis=0):
    """Deterministic frame-weighted mean (``sum(w * x)`` along ``axis``).

    Thin wrapper around :func:`numpy.average`. Because ``weights`` is a
    validated probability vector (``sum == 1``) this is exactly the
    re-weighted ensemble expectation and is numerically identical to the
    ``numpy.average(x, axis, weights=weights)`` calls already used
    elsewhere in SOURSOP.

    Parameters
    ----------
    x : array_like
        Per-frame data; the frame axis is ``axis``.
    weights : numpy.ndarray
        Validated per-frame weights (see :func:`validate_weights`).
    axis : int, optional
        Axis to average over (the frame axis). Default 0.

    Returns
    -------
    numpy.ndarray or float
        ``x`` averaged over ``axis`` with the given weights.

    Example
    -------
    >>> import numpy as np
    >>> weighted_mean(np.array([1.0, 3.0]), np.array([0.25, 0.75]))
    2.5
    """
    return numpy.average(x, axis=axis, weights=weights)


def weighted_rms(x, weights, axis=0):
    """Deterministic frame-weighted root-mean-square (``sqrt(sum(w*x^2))``).

    The polymer-physics order parameter used by the internal-scaling and
    scaling-exponent routines is an RMS distance, so this is the weighted
    analogue of ``sqrt(mean(x**2))``.

    Parameters
    ----------
    x : array_like
        Per-frame data; the frame axis is ``axis``.
    weights : numpy.ndarray
        Validated per-frame weights (see :func:`validate_weights`).
    axis : int, optional
        Axis to reduce over. Default 0.

    Returns
    -------
    numpy.ndarray or float
        The weighted RMS of ``x`` over ``axis``.

    Example
    -------
    >>> import numpy as np
    >>> float(weighted_rms(np.array([3.0, 4.0]), np.array([0.5, 0.5])))
    3.5355339059327378
    """
    return numpy.sqrt(numpy.average(numpy.square(x), axis=axis, weights=weights))


def weighted_std(x, weights, axis=0):
    """Deterministic frame-weighted (population) standard deviation.

    Uses the reliability-weighted **population** estimator
    ``sqrt(sum(w * (x - mean)**2))``. Frame weights here are probability
    weights with no associated sample size, so there is no well-defined
    ``ddof`` (Bessel-style) correction; the population estimator is the
    unambiguous, reproducible choice and the de-facto standard for
    re-weighted molecular-dynamics ensembles.

    Parameters
    ----------
    x : array_like
        Per-frame data; the frame axis is ``axis``.
    weights : numpy.ndarray
        Validated per-frame weights (see :func:`validate_weights`).
    axis : int, optional
        Axis to reduce over. Default 0.

    Returns
    -------
    numpy.ndarray or float
        The weighted population standard deviation of ``x`` over ``axis``.

    Example
    -------
    >>> import numpy as np
    >>> float(weighted_std(np.array([1.0, 1.0]), np.array([0.5, 0.5])))
    0.0
    """
    x = numpy.asarray(x, dtype=numpy.float64)
    m = numpy.average(x, axis=axis, weights=weights)
    # keep the reduced axis so (x - m) broadcasts for any axis
    m_b = numpy.expand_dims(m, axis) if x.ndim > 1 else m
    var = numpy.average(numpy.square(x - m_b), axis=axis, weights=weights)
    return numpy.sqrt(var)


def weighted_corr(a, b, weights):
    """Deterministic frame-weighted Pearson correlation between two vectors.

    Computed from the weighted covariance matrix
    (``numpy.cov(..., ddof=0, aweights=weights)``), matching the
    ``ddof=0`` / ``aweights`` convention already used by
    ``get_local_to_global_correlation``.

    Parameters
    ----------
    a, b : array_like
        Equal-length per-frame vectors.
    weights : numpy.ndarray
        Validated per-frame weights (see :func:`validate_weights`).

    Returns
    -------
    float
        The weighted Pearson correlation coefficient of ``a`` and ``b``.

    Example
    -------
    >>> import numpy as np
    >>> w = np.full(4, 0.25)
    >>> round(float(weighted_corr(np.array([1.,2,3,4]), np.array([2.,4,6,8]), w)), 6)
    1.0
    """
    cov = numpy.cov(numpy.vstack((a, b)), ddof=0, aweights=weights)
    denom = numpy.sqrt(cov[0, 0] * cov[1, 1])
    return cov[0, 1] / denom
