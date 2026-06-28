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
import numpy
import ctypes
import platform
import warnings
from dataclasses import dataclass
from typing import Optional, Tuple
from soursop.ssexceptions import SSException


MKL_LIBRARY = "mkl_rt"
OPENBLAS_LIBRARY = "openblas"


def _set_mkl_numpy_threads(mkl_path, num_threads):
    # Traditional UNIX-like systems will have shared objects available.
    #
    # Darwin / Apple uses `*.dylib` by default for included Intel compiler libraries.
    # Traditional UNIX-like shared objects can be created (`*.so`), but are more
    # represented in third-party libraries. This is a more dynamic way of finding
    # the MKL library and using it on a Mac that has an Intel compiler installed.
    mkl_rt = ctypes.CDLL(mkl_path)
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
    if os_name == "darwin":
        libname = f"*{library_name}*.dylib*"  # fuzzy match for filtering with find
    elif os_name == "linux":
        libname = f"*{library_name}*.so*"  # fuzzy match for filtering with find
    else:
        warnings.warn(f"Unsupported OS: {os_name}.")

    # Checking existing environment variables and stop on the first match. The
    # basis for this approach is that only one should be active.
    virtualized_env = None
    for env_var in "CONDA_PREFIX,VIRTUAL_ENV".split(","):
        env_path = os.environ.get(env_var, None)
        if env_path is not None:
            virtualized_env = env_path
            break

    # If no virtual environment is active, terminate. Support for system-wide
    # Python installations is not yet supported.
    if virtualized_env is None:
        raise SSException("No Anaconda or Python Virtual Environment found. Exiting.")

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
            numpy_path_fragment = os.path.join("site-packages", "numpy")  # os-agnostic
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
            warnings.warn(
                "Unsupported library. Please install OPENBLAS or the Intel MKL library. No threads set."
            )
            library = "unknown"
        break  # stop on first set library
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
    if platform.system().lower() == "windows":
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
            message = (
                f"Keyword {keyword_name} passed value [{keyword}], but this is not valid.\nMust be one of :%s"
                % (", ".join(allowed_vals))
            )
        else:
            error_type = type(error_message)
            if error_type is not str:
                raise RuntimeError(
                    'Invalid error message type: "{}". The error message must be a string.'.format(
                        error_type
                    )
                )
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


# ======================================================================
#
# Reweighting primitives shared by soursop.ssbme (BME / iBME) and
# soursop.sscoper (COPER / iCOPER). Both modules import these and
# re-export the public names, so user code and the per-module APIs stay
# identical regardless of which reweighter is used.
#
# ======================================================================

#: Weights below this threshold are treated as zero in relative-entropy
#: sums (avoids ``log(0)`` for de-populated frames).
MIN_WEIGHT_THRESHOLD = 1e-50

#: Valid experimental constraint types for :class:`ExperimentalObservable`.
VALID_CONSTRAINTS = {"equality", "upper", "lower"}


# ........................................................................
#
@dataclass
class ExperimentalObservable:
    """Container for a single experimental observable.

    Shared by :mod:`soursop.ssbme` and :mod:`soursop.sscoper` (both
    re-export this class), so the user-facing syntax is identical for BME,
    iBME, COPER and iCOPER.

    Parameters
    ----------
    value : float
        The experimental value of the observable.
    uncertainty : float
        The experimental uncertainty (standard deviation). Must be positive.
    constraint : str, optional
        Type of constraint, one of:

        - ``"equality"`` (default): observable should match ``value`` within
          ``uncertainty``.
        - ``"upper"``: observable should not exceed ``value`` (deviations
          below ``value`` are not penalized).
        - ``"lower"``: observable should not fall below ``value`` (deviations
          above ``value`` are not penalized).
    name : str, optional
        Optional human-readable name/description.
    group : str, optional
        Optional data-type label. Used by :class:`soursop.sscoper.COPER` to
        impose a separate per-group chi-squared constraint
        (``chi2_alpha <= limit`` for each group, as in Leung et al. 2016);
        ignored by BME / iBME. Observables without a group are pooled into a
        single default group.

    Raises
    ------
    SSException
        If ``uncertainty`` is not positive or ``constraint`` is invalid.
    """

    value: float
    uncertainty: float
    constraint: str = "equality"
    name: Optional[str] = None
    group: Optional[str] = None

    def __post_init__(self):
        if self.uncertainty <= 0:
            raise SSException(f"Uncertainty must be positive, got {self.uncertainty}")

        if not isinstance(self.constraint, str):
            raise SSException(
                "constraint must be a string ('equality', 'upper', or "
                f"'lower'), got {type(self.constraint).__name__}"
            )

        constraint_lower = self.constraint.lower().strip()
        if constraint_lower not in VALID_CONSTRAINTS:
            raise SSException(
                f"Invalid constraint: '{self.constraint}'. "
                "Must be 'equality', 'upper', or 'lower'"
            )

        self.constraint = constraint_lower

    def get_bounds(self) -> Tuple[Optional[float], Optional[float]]:
        """Optimization bounds on the Lagrange multiplier for this observable.

        Returns
        -------
        tuple
            ``(None, None)`` for ``equality``, ``(0.0, None)`` for ``upper``,
            ``(None, 0.0)`` for ``lower``.
        """
        if self.constraint == "equality":
            return (None, None)
        elif self.constraint == "upper":
            return (0.0, None)
        else:  # "lower"
            return (None, 0.0)


# ........................................................................
#
def relative_entropy(w0, w1):
    """Relative entropy (Kullback-Leibler divergence) of ``w1`` from ``w0``.

    Parameters
    ----------
    w0 : numpy.ndarray
        Reference (prior) weights, normalized to sum to 1.
    w1 : numpy.ndarray
        Posterior weights, normalized to sum to 1.

    Returns
    -------
    float
        ``sum_i w1_i * log(w1_i / w0_i)`` over frames with non-negligible
        posterior weight.
    """
    idxs = numpy.where(w1 > MIN_WEIGHT_THRESHOLD)
    return float(numpy.sum(w1[idxs] * numpy.log(w1[idxs] / w0[idxs])))


# ........................................................................
#
def weighted_linear_regression(x, y, sample_weight, fit_intercept=True):
    """Closed-form weighted least-squares regression of ``y`` on ``x``.

    A small numpy replacement for ``sklearn.linear_model.LinearRegression``
    (SOURSOP does not depend on scikit-learn). Solves
    ``min_{a,b} sum_i s_i (y_i - (a x_i + b))^2``. Used by the iterative
    scale/offset reweighters (iBME, iCOPER).

    Parameters
    ----------
    x : numpy.ndarray
        Independent variable, shape ``(n,)``.
    y : numpy.ndarray
        Dependent variable, shape ``(n,)``.
    sample_weight : numpy.ndarray
        Per-sample weights, shape ``(n,)``.
    fit_intercept : bool, optional
        If True fit slope and intercept; if False force the intercept to
        zero (slope only). Default True.

    Returns
    -------
    tuple of float
        ``(slope, intercept)``. ``intercept`` is ``0.0`` when
        ``fit_intercept`` is False.
    """
    x = numpy.asarray(x, dtype=numpy.float64).ravel()
    y = numpy.asarray(y, dtype=numpy.float64).ravel()
    s = numpy.asarray(sample_weight, dtype=numpy.float64).ravel()

    if fit_intercept:
        sw = numpy.sum(s)
        x_mean = numpy.sum(s * x) / sw
        y_mean = numpy.sum(s * y) / sw
        cov_xy = numpy.sum(s * (x - x_mean) * (y - y_mean))
        var_x = numpy.sum(s * (x - x_mean) ** 2)
        slope = cov_xy / var_x
        intercept = y_mean - slope * x_mean
    else:
        slope = numpy.sum(s * x * y) / numpy.sum(s * x * x)
        intercept = 0.0

    return float(slope), float(intercept)


# ........................................................................
#
def _find_knee_perpendicular(x, y):
    """Knee index by maximum perpendicular distance to the endpoint chord."""
    x_n = (x - x.min()) / (x.max() - x.min() + 1e-10)
    y_n = (y - y.min()) / (y.max() - y.min() + 1e-10)

    p1 = numpy.array([x_n[0], y_n[0]])
    p2 = numpy.array([x_n[-1], y_n[-1]])
    line_vec = p2 - p1
    line_len = numpy.linalg.norm(line_vec)
    if line_len < 1e-10:
        return len(x) // 2
    line_unit = line_vec / line_len

    distances = []
    for i in range(len(x_n)):
        point = numpy.array([x_n[i], y_n[i]])
        vec = point - p1
        proj = p1 + numpy.dot(vec, line_unit) * line_unit
        distances.append(numpy.linalg.norm(point - proj))
    return int(numpy.argmax(distances))


# ........................................................................
#
def _find_knee_curvature(x, y):
    """Knee index by maximum Menger curvature (3-point estimate)."""
    x_n = (x - x.min()) / (x.max() - x.min() + 1e-10)
    y_n = (y - y.min()) / (y.max() - y.min() + 1e-10)
    n = len(x_n)
    curvature = numpy.zeros(n)
    for i in range(1, n - 1):
        p0 = numpy.array([x_n[i - 1], y_n[i - 1]])
        p1 = numpy.array([x_n[i], y_n[i]])
        p2 = numpy.array([x_n[i + 1], y_n[i + 1]])
        v1 = p1 - p0
        v2 = p2 - p1
        area = abs(v1[0] * v2[1] - v1[1] * v2[0]) / 2.0
        a = numpy.linalg.norm(p2 - p1)
        b = numpy.linalg.norm(p0 - p2)
        c = numpy.linalg.norm(p1 - p0)
        if a * b * c > 1e-10:
            curvature[i] = 4 * area / (a * b * c)
    if n > 2:
        curvature[0] = curvature[1]
        curvature[-1] = curvature[-2]
    return int(numpy.argmax(curvature))


# ........................................................................
#
def find_optimal_theta(x_values, y_values, method="perpendicular"):
    """Select the L-curve knee from two paired metric arrays.

    Generic knee-finder shared by BME's ``theta_scan`` (chi-squared vs.
    relative entropy across theta) and COPER's ``chi2_limit_scan``
    (chi-squared vs. relative entropy across the chi-squared limit).

    Parameters
    ----------
    x_values : numpy.ndarray
        First metric per scan point (e.g. final chi-squared).
    y_values : numpy.ndarray
        Second metric per scan point (e.g. relative entropy).
    method : str, optional
        ``"perpendicular"`` (default) or ``"curvature"``.

    Returns
    -------
    tuple
        ``(optimal_idx, method_name)``.

    Raises
    ------
    SSException
        If ``method`` is unknown.
    """
    if method == "curvature":
        return _find_knee_curvature(x_values, y_values), "Menger curvature"
    elif method == "perpendicular":
        return _find_knee_perpendicular(x_values, y_values), "Perpendicular distance"
    raise SSException(
        f"Unknown method: {method}, must be 'curvature' or 'perpendicular'"
    )


# ........................................................................
#
def validate_reweighting_inputs(observables, calculated_values, initial_weights):
    """Validate constructor arguments shared by the reweighter classes.

    Used by :class:`soursop.ssbme.BME` / ``iBME`` and
    :class:`soursop.sscoper.COPER` / ``iCOPER``.

    Parameters
    ----------
    observables : list of ExperimentalObservable
        Experimental observables (non-empty).
    calculated_values : numpy.ndarray
        Per-frame calculated values, shape ``(n_frames, n_observables)``.
    initial_weights : numpy.ndarray or None
        Optional prior weights, one per frame.

    Raises
    ------
    SSException
        If any argument is malformed or dimensions are inconsistent.
    """
    if not isinstance(observables, (list, tuple)) or len(observables) == 0:
        raise SSException("observables must be a non-empty list")
    if not all(isinstance(obs, ExperimentalObservable) for obs in observables):
        raise SSException("All observables must be ExperimentalObservable instances")
    if not isinstance(calculated_values, numpy.ndarray):
        raise SSException("calculated_values must be a numpy array")
    if calculated_values.ndim != 2:
        raise SSException("calculated_values must be 2D (n_frames, n_observables)")
    if calculated_values.shape[1] != len(observables):
        raise SSException(
            f"Number of observables ({len(observables)}) must match "
            f"calculated_values columns ({calculated_values.shape[1]})"
        )
    if initial_weights is not None:
        if not isinstance(initial_weights, numpy.ndarray):
            raise SSException("initial_weights must be a numpy array")
        if len(initial_weights) != calculated_values.shape[0]:
            raise SSException("initial_weights length must match number of frames")


# ........................................................................
#
def constraint_chi_squared(weights, calculated_values, observables, indices=None):
    """Constraint-aware reduced chi-squared for a weight vector.

    ``equality`` observables always penalize deviations; ``upper`` /
    ``lower`` only penalize the disallowed side. Shared by BME (total
    chi-squared) and COPER (per-data-type chi-squared via ``indices``).

    Parameters
    ----------
    weights : numpy.ndarray
        Frame weights, shape ``(n_frames,)``.
    calculated_values : numpy.ndarray
        Per-frame calculated values, shape ``(n_frames, n_observables)``.
    observables : list of ExperimentalObservable
        The experimental observables (columns of ``calculated_values``).
    indices : sequence of int, optional
        Restrict the chi-squared to this subset of observable columns
        (used for COPER per-group constraints). Defaults to all columns.

    Returns
    -------
    float
        Mean of ``(diff / sigma)^2`` over the selected observables.
    """
    if indices is None:
        indices = range(len(observables))
    else:
        indices = list(indices)

    chi_squared = 0.0
    count = 0
    for idx in indices:
        obs = observables[idx]
        calc_avg = numpy.sum(calculated_values[:, idx] * weights)
        diff = calc_avg - obs.value
        if obs.constraint == "equality":
            penalize = True
        elif obs.constraint == "upper":
            penalize = diff > 0
        else:  # "lower"
            penalize = diff < 0
        if penalize:
            chi_squared += (diff / obs.uncertainty) ** 2
        count += 1
    return chi_squared / count


## ------------------------------------------------------------------------
##
## Two-bead (CA/CB) coarse-grained model support
##
## The constants and ideal-helix geometry below are vendored from an internal reference implementation of the
## package (``swan/helix.py`` and ``swan/trajectory.py``) so that SOURSOP can
## detect two-bead coarse-grained trajectories and assign secondary structure from the CA trace
## *without* taking a runtime dependency on that reference implementation. If the model ever
## changes its ideal alpha-helix parameters these must be kept in sync.

# Ideal alpha-helix parameters (see swan/helix.py)
SWAN_HELIX_RISE = 1.5  # angstrom per residue along the axis
SWAN_HELIX_TWIST_DEG = 100.0  # degrees per residue (3.6 residues / turn)
SWAN_HELIX_CA_RADIUS = 2.3  # angstrom of Calpha from the helix axis

# Ideal extended (beta) strand parameters. This model does not generate beta, so there
# is no reference geometry from that model; these describe a canonical pleated extended CA
# strand: a virtual CA-CA bond of ~3.8 A with a CA(i)..CA(i+2) span of ~6.8 A
# (clearly distinct from the ~5.4 A helical value). Modelled as a planar zigzag
# with axial spacing SWAN_BETA_AXIAL and transverse amplitude SWAN_BETA_AMPLITUDE.
SWAN_BETA_AXIAL = 3.4  # angstrom along the strand axis (half of CA(i)..CA(i+2))
SWAN_BETA_AMPLITUDE = 1.70  # angstrom transverse pleat amplitude


def is_swan_topology(topology):
    """Return ``True`` if a topology is a two-bead (CA/CB) coarse-grained model.

    A two-bead topology represents every residue with a single backbone ``CA`` bead
    and (for every residue except glycine) a single sidechain ``CB`` bead. This
    is distinct from SOURSOP's existing one-bead-per-residue coarse-grained model
    (``CA`` only), which is why the presence of at least one ``CB`` bead is
    required.

    The check is intentionally strict: every atom in the topology must be named
    ``CA`` or ``CB``, every residue must contain exactly one ``CA``, and every
    residue must contain exactly one ``CB`` unless it is glycine (which must have
    none).

    Parameters
    ----------
    topology : mdtraj.Topology
        The topology to inspect.

    Returns
    -------
    bool
        ``True`` if the topology matches the two-bead model, ``False``
        otherwise.

    Example
    -------
    >>> is_swan_topology(traj.topology)
    True
    """

    n_cb_total = 0
    for residue in topology.residues:
        n_ca = 0
        n_cb = 0
        for atom in residue.atoms:
            if atom.name == "CA":
                n_ca += 1
            elif atom.name == "CB":
                n_cb += 1
            else:
                # any non-CA/CB atom immediately disqualifies the topology
                return False

        # every residue must have exactly one CA
        if n_ca != 1:
            return False

        # glycine carries no sidechain bead; everything else carries exactly one
        if residue.name == "GLY":
            if n_cb != 0:
                return False
        else:
            if n_cb != 1:
                return False

        n_cb_total += n_cb

    # a topology with no CB at all is the existing CA-only 1-bead CG model, not two-bead
    return n_cb_total > 0


def ideal_helix_ca(
    n, rise=SWAN_HELIX_RISE, twist_deg=SWAN_HELIX_TWIST_DEG, radius=SWAN_HELIX_CA_RADIUS
):
    """Ideal alpha-helix Calpha coordinates about the +z axis.

    Vendored from ``swan.helix.ideal_helix_ca``. Generates the Calpha trace of
    an idealized alpha-helix, which is used as the reference fragment when
    detecting helicity from a CA-only coarse-grained trajectory.

    Parameters
    ----------
    n : int
        Number of consecutive Calpha beads to generate.
    rise : float, optional
        Rise per residue along the helix axis, in Angstroms. Default ``1.5``.
    twist_deg : float, optional
        Twist per residue, in degrees. Default ``100.0`` (3.6 residues/turn).
    radius : float, optional
        Radius of the Calpha from the helix axis, in Angstroms. Default ``2.3``.

    Returns
    -------
    numpy.ndarray
        Array of shape ``(n, 3)`` of ideal Calpha coordinates in Angstroms.

    Example
    -------
    >>> ideal_helix_ca(4).shape
    (4, 3)
    """

    t = numpy.radians(twist_deg)
    i = numpy.arange(n, dtype=numpy.float64)
    return numpy.column_stack(
        [radius * numpy.cos(i * t), radius * numpy.sin(i * t), rise * i]
    )


def ideal_extended_ca(n, axial=SWAN_BETA_AXIAL, amplitude=SWAN_BETA_AMPLITUDE):
    """Idealized extended (beta) strand Calpha coordinates.

    Generates the Calpha trace of a canonical planar pleated extended strand,
    used as the reference fragment when detecting beta content from a CA-only
    coarse-grained trajectory. The default geometry gives a virtual CA-CA bond
    of ~3.8 Angstrom and a CA(i)..CA(i+2) span of ~6.8 Angstrom.

    Parameters
    ----------
    n : int
        Number of consecutive Calpha beads to generate.
    axial : float, optional
        Spacing along the strand axis, in Angstroms (half of the CA(i)..CA(i+2)
        span). Default ``3.4``.
    amplitude : float, optional
        Transverse pleat amplitude, in Angstroms. Default ``1.70``.

    Returns
    -------
    numpy.ndarray
        Array of shape ``(n, 3)`` of ideal extended-strand Calpha coordinates.

    Example
    -------
    >>> ideal_extended_ca(5).shape
    (5, 3)
    """

    i = numpy.arange(n, dtype=numpy.float64)
    x = i * axial
    y = (i.astype(numpy.int64) % 2) * amplitude
    z = numpy.zeros(n, dtype=numpy.float64)
    return numpy.column_stack([x, y, z])


def kabsch_rmsd(P, Q):
    """Minimal root-mean-square deviation between two point sets after optimal superposition.

    Computes the optimal rigid (rotation + translation) alignment of ``P`` onto
    ``Q`` via the Kabsch algorithm and returns the resulting RMSD. Used to score
    how closely a fragment of a CA trace matches an idealized alpha-helix.

    Parameters
    ----------
    P : numpy.ndarray
        Array of shape ``(K, 3)`` -- the mobile point set.
    Q : numpy.ndarray
        Array of shape ``(K, 3)`` -- the reference point set.

    Returns
    -------
    float
        The minimal RMSD (same units as the inputs) after optimal superposition.

    Example
    -------
    >>> float(kabsch_rmsd(ideal_helix_ca(5), ideal_helix_ca(5)))
    0.0
    """

    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    h = Pc.T @ Qc
    u, _, vt = numpy.linalg.svd(h)
    d = numpy.sign(numpy.linalg.det(vt.T @ u.T))
    rot = vt.T @ numpy.diag([1.0, 1.0, d]) @ u.T
    Pr = Pc @ rot.T
    return float(numpy.sqrt(((Pr - Qc) ** 2).sum() / len(P)))
