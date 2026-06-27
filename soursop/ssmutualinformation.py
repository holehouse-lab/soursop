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


# ........................................................................
#
def calc_MI(X, Y, bins, weights=False, normalize=False):
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
        Per-sample weights. If provided, ``len(weights) == len(X)`` is
        required and the weights are forwarded to the underlying
        ``np.histogram*`` calls. Default ``False`` (uniform).
    normalize : bool, optional
        If True, divide ``I`` by the joint entropy ``H(X, Y)`` to return
        NMI in ``[0, 1]``. NaN values (e.g. from ``H(X,Y) == 0``) are
        coerced to 0. Default False.

    Returns
    -------
    float
        Mutual information (in nats), or NMI if ``normalize=True``. The
        absolute scale of MI depends on the chosen bin size.

    Raises
    ------
    SSException
        If ``X`` and ``Y`` differ in length, or ``bins`` does not cover
        the full data range of both.

    Example
    -------
    >>> import numpy as np
    >>> from soursop.ssmutualinformation import calc_MI
    >>> rng = np.random.default_rng(0)
    >>> X = rng.uniform(-1, 1, 1000)
    >>> Y = X + 0.05 * rng.standard_normal(1000)
    >>> bins = np.arange(-1.01, 1.02, 0.1)
    >>> calc_MI(X, Y, bins)              # strong dependence
    2.18
    >>> calc_MI(X, rng.uniform(-1, 1, 1000), bins)  # ~independent
    0.05
    """

    if len(X) != len(Y):
        raise SSException("Error: X and Y vectors must be the same length")

    if np.min(bins) > np.min(X) or np.min(bins) > np.min(Y):
        raise SSException(
            f"Error: Bins passed to calc_MI in ssmutualinformation() do not straddle the full data range. Bin minimum {np.min(bins)} is bigger than one/both of data minima: X={np.min(X)}, Y={np.min(Y)}"
        )

    if np.max(bins) < np.max(X) or np.max(bins) < np.max(Y):
        raise SSException(
            f"Error: Bins passed to calc_MI in ssmutualinformation() do not straddle the full data range. Bin max {np.max(bins)} is smaller than one/both of data maxima: X={np.max(X)}, Y={np.max(Y)}"
        )

    if weights is not False and weights is not None:
        c_XY = np.histogram2d(X, Y, bins, weights=weights)[0]
        c_X = np.histogram(X, bins, weights=weights)[0]
        c_Y = np.histogram(Y, bins, weights=weights)[0]
    else:
        c_XY = np.histogram2d(X, Y, bins)[0]
        c_X = np.histogram(X, bins)[0]
        c_Y = np.histogram(Y, bins)[0]

    H_X = shan_entropy(c_X)
    H_Y = shan_entropy(c_Y)
    H_XY = shan_entropy(c_XY)

    MI = H_X + H_Y - H_XY

    if normalize:
        # Uncomment to *locally* suppress numpy divide by zero runtime warning.
        # This will be cleaned up before returning nmi.
        # with np.errstate(invalid='ignore'):
        nmi = MI / H_XY
        # aforementioned cleanup
        if np.isnan(nmi):
            # print(f"Normalization led to {np.nan} - replacing with 0.")
            nmi = 0
        return nmi
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
    >>> shan_entropy(np.array([1, 1, 1, 1]))     # uniform 4-bin
    1.386
    >>> shan_entropy(np.array([1, 0, 0, 0]))     # peaked
    0.0
    """
    # normalize such that all elements sum up to 1
    c_normalized = c / float(np.sum(c))

    # now convert into a single vector of non-zero elements
    c_normalized = c_normalized[np.nonzero(c_normalized)]

    # compute the entropy associated with this vector. The more
    # evenly distributed the greater the entropy
    H = -sum(c_normalized * np.log(c_normalized))
    return H
