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


"""
ssmutualinformation contains all the functions associated with computing things related to mutual
information. This file is not a class, but represents a set of stand alone functions, rather than
including this function in any one other class.
"""

import numpy as np
from .ssexceptions import SSException
from .ssutils import validate_weights


# ........................................................................
#
def calc_MI(X, Y, bins, weights=False, normalize=False, etol=0.0000001):
    """Mutual information :math:`I(X; Y)` between two observables.

    Computed from histogram-based estimates of :math:`p(X)`, :math:`p(Y)`,
    and the joint :math:`p(X, Y)` over the supplied bin edges:

    .. math::
        I(X; Y) = H(X) + H(Y) - H(X, Y)

    where :math:`H` is the Shannon entropy of the binned distribution.
    Optionally returns the normalised mutual information
    :math:`NMI = I / H(X, Y)` instead.

    All three histograms use the same ``bins`` edges, which must
    completely cover the range of both ``X`` and ``Y``.

    Parameters
    ----------
    X : array_like
        First observable. Must be 1D and the same length as ``Y``.
    Y : array_like
        Second observable. Must be 1D and the same length as ``X``.
    bins : array_like
        Bin edges, monotonically increasing, that cover both ``X`` and
        ``Y``. The same edges are reused for the 1D and 2D histograms.
    weights : array_like or False, optional
        Per-sample weights following the package-wide weights contract
        (one per sample, each in ``[0, 1]``, finite, summing to 1 within
        ``etol``), validated by :func:`soursop.ssutils.validate_weights`
        and forwarded to the underlying ``np.histogram*`` calls. Default
        ``False`` (uniform).
    normalize : bool, optional
        If True, divide ``I`` by the joint entropy ``H(X, Y)`` to return
        NMI in ``[0, 1]``. When ``H(X, Y) == 0`` (both variables constant)
        NMI is defined as 0. Default False.
    etol : float, optional
        Tolerance on ``|sum(weights) - 1|``. Default ``1e-7``.

    Returns
    -------
    float
        Mutual information (in nats), or NMI if ``normalize=True``. The
        absolute scale of MI depends on the chosen bin size.

    Raises
    ------
    SSException
        If ``X`` and ``Y`` differ in length or contain non-finite values,
        if ``bins`` does not cover the full data range of both, or if
        ``weights`` fails validation (see above). Previously NaN data slipped
        past the range check and gave a spurious MI, and invalid weights
        gave ``nan`` (silently turned into 0 by ``normalize=True``).

    Example
    -------
    >>> import numpy as np
    >>> from soursop.ssmutualinformation import calc_MI
    >>> rng = np.random.default_rng(0)
    >>> X = rng.uniform(-1, 1, 1000)
    >>> Y = X + 0.05 * rng.standard_normal(1000)
    >>> # bins must span both X and Y (the noise pushes Y just past +-1)
    >>> bins = np.linspace(-1.5, 1.5, 31)
    >>> round(calc_MI(X, Y, bins), 2)              # strong dependence
    2.1
    >>> round(calc_MI(X, rng.uniform(-1, 1, 1000), bins), 2)  # ~independent
    0.18
    """

    X = np.asarray(X, dtype=np.float64).ravel()
    Y = np.asarray(Y, dtype=np.float64).ravel()
    bins = np.asarray(bins, dtype=np.float64)

    if len(X) != len(Y):
        raise SSException("Error: X and Y vectors must be the same length")

    # np.min/np.max of NaN data are NaN, which made both range checks below
    # False; np.histogram then dropped the NaN samples from some histograms
    # but not others, so the marginals no longer matched the joint
    if not (np.all(np.isfinite(X)) and np.all(np.isfinite(Y))):
        raise SSException("Error: X and Y must contain only finite values")

    # the same weights contract as everywhere else in SOURSOP (previously
    # invalid weights gave nan, which normalize=True silently turned into 0)
    weights = validate_weights(weights, len(X), stride=1, etol=etol)

    if np.min(bins) > np.min(X) or np.min(bins) > np.min(Y):
        raise SSException(
            f"Error: Bins passed to calc_MI in ssmutualinformation() do not straddle the full data range. Bin minimum {np.min(bins)} is bigger than one/both of data minima: X={np.min(X)}, Y={np.min(Y)}"
        )

    if np.max(bins) < np.max(X) or np.max(bins) < np.max(Y):
        raise SSException(
            f"Error: Bins passed to calc_MI in ssmutualinformation() do not straddle the full data range. Bin max {np.max(bins)} is smaller than one/both of data maxima: X={np.max(X)}, Y={np.max(Y)}"
        )

    # the edges are given once per axis to histogram2d: a bare length-2 edge
    # array would otherwise be read as [nx, ny] bin counts
    if weights is not False and weights is not None:
        c_XY = np.histogram2d(X, Y, [bins, bins], weights=weights)[0]
        c_X = np.histogram(X, bins, weights=weights)[0]
        c_Y = np.histogram(Y, bins, weights=weights)[0]
    else:
        c_XY = np.histogram2d(X, Y, [bins, bins])[0]
        c_X = np.histogram(X, bins)[0]
        c_Y = np.histogram(Y, bins)[0]

    H_X = shan_entropy(c_X)
    H_Y = shan_entropy(c_Y)
    H_XY = shan_entropy(c_XY)

    MI = H_X + H_Y - H_XY

    if normalize:
        # NMI is defined as 0 when the joint entropy is zero (both variables
        # constant); only that case is coerced, so any other NaN still shows
        if H_XY == 0:
            return 0.0
        return MI / H_XY
    else:
        return MI


# ........................................................................
#
def shan_entropy(c):
    """Shannon entropy (in nats) of a histogram-style array.

    Treats the input as unnormalised counts: normalises by the total,
    drops zero entries (to avoid ``log(0)``), and returns
    :math:`H = -\\sum_i p_i \\ln p_i`. Works for both 1D vectors (entropy
    of a single distribution) and 2D matrices (joint entropy of a 2D
    histogram).

    Parameters
    ----------
    c : array_like
        1D vector or 2D matrix of non-negative counts (or weights).
        The values are normalised internally so they need not sum to 1.

    Returns
    -------
    float
        Shannon entropy in nats. 0 means a perfectly peaked distribution;
        the maximum is :math:`\\ln(N)` for a uniform distribution over
        ``N`` non-zero bins.

    Example
    -------
    >>> import numpy as np
    >>> from soursop.ssmutualinformation import shan_entropy
    >>> round(float(shan_entropy(np.array([1, 1, 1, 1]))), 3)   # uniform 4-bin
    1.386
    >>> float(shan_entropy(np.array([1, 0, 0, 0])))             # peaked
    -0.0
    """
    # normalize such that all elements sum up to 1 (accepting any array_like,
    # e.g. a plain list, as documented)
    c = np.asarray(c, dtype=np.float64)
    c_normalized = c / float(np.sum(c))

    # now convert into a single vector of non-zero elements
    c_normalized = c_normalized[np.nonzero(c_normalized)]

    # compute the entropy associated with this vector. The more
    # evenly distributed the greater the entropy
    H = -np.sum(c_normalized * np.log(c_normalized))
    return H
