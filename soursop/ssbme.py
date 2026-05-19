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
ssbme - Bayesian Maximum Entropy (BME) and iterative BME (iBME) reweighting.

This module reweights molecular ensembles against experimental observables
using the Bayesian Maximum Entropy framework. BME finds the minimally biased
set of frame weights (smallest relative entropy from the prior) that brings
the ensemble averages into agreement with experimental data.

Two reweighters are provided:

``BME``
    Standard maximum-entropy reweighting. Use when the calculated observables
    are already on the same scale as the experimental values.

``iBME``
    Iterative BME. At every iteration a (weighted) linear regression of the
    ensemble-averaged calculated values against the experimental values is
    used to refit a global scale and offset before the maxent step. This is
    the method of choice when there is an unknown global scale/offset between
    the calculated and experimental data (e.g. SAXS intensities). Ported from
    the reference implementation of Bottaro et al.

``ExperimentalObservable``
    Container for a single experimental data point (value, uncertainty and
    constraint type).

The optimal regularisation strength ``theta`` can be selected automatically
with an L-curve scan via :func:`theta_scan` (also exposed as
``BME.scan_theta`` / ``iBME.scan_theta``).

Example
-------
>>> import numpy as np
>>> from soursop.ssbme import BME, iBME, ExperimentalObservable
>>> rng = np.random.default_rng(0)
>>> calc = rng.normal([24, 65], [3, 6], size=(2000, 2))
>>> obs = [ExperimentalObservable(23.0, 1.0, name="Rg"),
...        ExperimentalObservable(60.0, 2.0, name="Ree")]
>>> result = BME(obs, calc).fit(theta=2.0, auto_theta=False)
>>> reweighted_means = result.predict(calc)

References
----------
Primary reference for the BME procedure and software this module is
adapted from:

* Bottaro, S., Bengtsen, T. & Lindorff-Larsen, K. *Integrating Molecular
  Simulation and Experimental Data: A Bayesian/Maximum Entropy Reweighting
  Approach.* Methods Mol. Biol. 2112, 219-240 (2020).
  https://doi.org/10.1007/978-1-0716-0270-6_15
  (preprint: bioRxiv 2018, https://doi.org/10.1101/457952). Reference
  implementation: https://github.com/KULL-Centre/BME.

Foundational and closely related work:

* Hummer, G. & Koefinger, J. *Bayesian ensemble refinement by replica
  simulations and reweighting.* J. Chem. Phys. 143, 243150 (2015) - the
  BioEn reweighting functional that BME is mathematically equivalent to.
* Cesari, A., Gil-Ley, A. & Bussi, G. *Combining simulations and solution
  experiments as a paradigm for RNA force field refinement.* J. Chem.
  Theory Comput. 12, 6192-6200 (2016) - the maximum-entropy-with-error
  objective Gamma(lambda) minimised here.
* Cesari, A., Reisser, S. & Bussi, G. *Using the maximum entropy
  principle to combine simulations and solution experiments.*
  Computation 6, 15 (2018) - review of the method and its pitfalls.
* Rozycki, B., Kim, Y. C. & Hummer, G. *SAXS ensemble refinement of
  ESCRT-III CHMP3 conformational transitions.* Structure 19, 109-116
  (2011) - EROS; origin of the global scaling parameter theta and the
  SAXS scale/offset problem addressed by iBME.
* Bottaro, S. & Lindorff-Larsen, K. *Biophysical experiments and
  biomolecular simulations: A perfect match?* Science 361, 355-360
  (2018) - perspective on integrating simulation and experiment.

**Author(s):** Alex Holehouse
"""

from dataclasses import dataclass, field
from datetime import datetime
from typing import List, Optional, Tuple, Union

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp, rel_entr

from .ssexceptions import SSException

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Module constants
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
DEFAULT_THETA = 0.5
DEFAULT_MAX_ITERATIONS = 50000
DEFAULT_OPTIMIZER = "L-BFGS-B"
LAMBDA_INIT_SCALE = 1e-3
MIN_WEIGHT_THRESHOLD = 1e-50

#: Valid experimental constraint types.
VALID_CONSTRAINTS = {"equality", "upper", "lower"}


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
def _srel(w0: np.ndarray, w1: np.ndarray) -> float:
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
    idxs = np.where(w1 > MIN_WEIGHT_THRESHOLD)
    return float(np.sum(w1[idxs] * np.log(w1[idxs] / w0[idxs])))


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
def _weighted_linear_regression(
    x: np.ndarray,
    y: np.ndarray,
    sample_weight: np.ndarray,
    fit_intercept: bool = True,
) -> Tuple[float, float]:
    """Closed-form weighted least-squares regression of ``y`` on ``x``.

    This is a small numpy replacement for ``sklearn.linear_model``'s
    ``LinearRegression`` (SOURSOP does not depend on scikit-learn). It solves
    ``min_{a,b} sum_i s_i (y_i - (a x_i + b))^2``.

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
    x = np.asarray(x, dtype=np.float64).ravel()
    y = np.asarray(y, dtype=np.float64).ravel()
    s = np.asarray(sample_weight, dtype=np.float64).ravel()

    if fit_intercept:
        sw = np.sum(s)
        x_mean = np.sum(s * x) / sw
        y_mean = np.sum(s * y) / sw
        cov_xy = np.sum(s * (x - x_mean) * (y - y_mean))
        var_x = np.sum(s * (x - x_mean) ** 2)
        slope = cov_xy / var_x
        intercept = y_mean - slope * x_mean
    else:
        slope = np.sum(s * x * y) / np.sum(s * x * x)
        intercept = 0.0

    return float(slope), float(intercept)


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
@dataclass
class ExperimentalObservable:
    """Container for a single experimental observable.

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

    Raises
    ------
    SSException
        If ``uncertainty`` is not positive or ``constraint`` is invalid.
    """

    value: float
    uncertainty: float
    constraint: str = "equality"
    name: Optional[str] = None

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


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
@dataclass
class BMEResult:
    """Container for a BME / iBME reweighting result.

    Attributes
    ----------
    weights : numpy.ndarray
        Optimized (posterior) frame weights, summing to 1.
    initial_weights : numpy.ndarray
        Prior frame weights.
    lambdas : numpy.ndarray
        Optimized Lagrange multipliers (one per observable).
    chi_squared_initial : float
        Reduced chi-squared before reweighting.
    chi_squared_final : float
        Reduced chi-squared after reweighting.
    phi : float
        Fraction of effective frames, ``exp(-D_KL(w || w0))``.
    n_iterations : int
        Number of optimizer iterations (maxent) or iBME iterations.
    success : bool
        Whether the optimization succeeded.
    message : str
        Optimizer status message.
    theta : float
        Regularisation strength used.
    observables : list of ExperimentalObservable
        Observables used in the fit.
    calculated_values : numpy.ndarray
        Calculated values used in the fit (for iBME, the final scaled array).
    metadata : dict
        Free-form metadata (optimizer, timestamp, ...).
    scale : float, optional
        Final iBME scale factor (alpha) relative to the original input.
        ``None`` for plain BME.
    offset : float, optional
        Final iBME offset (beta) relative to the original input. ``None``
        for plain BME.
    ibme_iterations : list of dict
        Per-iteration iBME log (empty for plain BME). Each entry has
        ``iteration``, ``scale``, ``offset``, ``chi_squared`` and ``diff``.
    """

    weights: np.ndarray
    initial_weights: np.ndarray
    lambdas: np.ndarray
    chi_squared_initial: float
    chi_squared_final: float
    phi: float
    n_iterations: int
    success: bool
    message: str
    theta: float
    observables: List[ExperimentalObservable]
    calculated_values: np.ndarray
    metadata: dict
    scale: Optional[float] = None
    offset: Optional[float] = None
    ibme_iterations: List[dict] = field(default_factory=list)

    def __str__(self):
        status = "SUCCESS" if self.success else "FAILED"
        lines = [
            f"BME Result [{status}]",
            f"  Chi-squared initial: {self.chi_squared_initial:.4f}",
            f"  Chi-squared final:   {self.chi_squared_final:.4f}",
            f"  phi (effective fraction): {self.phi:.4f}",
            f"  Iterations: {self.n_iterations}",
            f"  Theta: {self.theta}",
        ]
        if self.scale is not None:
            lines.append(f"  iBME scale/offset: {self.scale:.4g} / {self.offset:.4g}")
        return "\n".join(lines)

    def __repr__(self):
        return self.__str__()

    @property
    def kl_divergence(self) -> float:
        """Relative entropy (KL divergence) of ``weights`` from prior."""
        if self.phi > 0:
            return -np.log(self.phi)
        return np.inf

    def predict(self, calculated_values: np.ndarray) -> np.ndarray:
        """Weighted average of arbitrary observables using these weights.

        Parameters
        ----------
        calculated_values : numpy.ndarray
            Per-frame observable values, shape ``(n_frames, n_observables)``.
            ``n_frames`` must match the fitted ensemble.

        Returns
        -------
        numpy.ndarray
            Weighted means, shape ``(n_observables,)``.

        Raises
        ------
        SSException
            If the number of frames does not match the fitted ensemble.
        """
        calculated_values = np.asarray(calculated_values)
        if calculated_values.shape[0] != self.weights.shape[0]:
            raise SSException(
                f"Number of frames in calculated_values "
                f"({calculated_values.shape[0]}) must match the fitted "
                f"ensemble ({self.weights.shape[0]})"
            )
        return np.sum(self.weights[:, np.newaxis] * calculated_values, axis=0)

    def diagnostics(self, warn_threshold: float = 0.5) -> dict:
        """Diagnose the reweighting result and flag potential issues.

        Parameters
        ----------
        warn_threshold : float, optional
            Minimum acceptable ``phi`` before a low-diversity warning is
            raised. Default 0.5.

        Returns
        -------
        dict
            Diagnostic metrics and a list of human-readable ``warnings``.

        Notes
        -----
        Two effective sample sizes are reported: the entropy-based
        ``neff_entropy`` (``N * phi``, the usual BME measure) and the
        Renyi-2 / participation ratio ``neff_renyi2`` (``1 / sum w^2``),
        which is more sensitive to a few dominant weights.
        """
        diag: dict = {}
        warnings: List[str] = []
        n = len(self.weights)

        neff_entropy = n * float(self.phi)
        diag["neff_entropy"] = neff_entropy
        diag["neff_entropy_fraction"] = float(self.phi)

        neff_renyi2 = 1.0 / np.sum(self.weights**2)
        diag["neff_renyi2"] = float(neff_renyi2)
        diag["neff_renyi2_fraction"] = float(neff_renyi2 / n)

        diag["weight_min"] = float(self.weights.min())
        diag["weight_max"] = float(self.weights.max())
        diag["weight_std"] = float(self.weights.std())
        if diag["weight_min"] > 0:
            diag["weight_range_orders"] = float(
                np.log10(diag["weight_max"] / diag["weight_min"])
            )
        else:
            diag["weight_range_orders"] = np.inf

        diag["chi2_improvement"] = self.chi_squared_initial - self.chi_squared_final
        diag["chi2_improvement_pct"] = (
            (diag["chi2_improvement"] / self.chi_squared_initial) * 100.0
            if self.chi_squared_initial > 0
            else np.nan
        )

        if self.phi < warn_threshold:
            warnings.append(
                f"Low Phi ({self.phi:.3f} < {warn_threshold}): significant "
                "loss of ensemble diversity. Consider increasing theta or "
                "loosening observable uncertainties."
            )
        if neff_renyi2 < 0.1 * n:
            warnings.append(
                f"Low effective sample size (1/sum w^2) ({neff_renyi2:.1f} / "
                f"{n}): only ~{diag['neff_renyi2_fraction'] * 100:.1f}% of "
                "frames are effectively used."
            )
        if diag["weight_range_orders"] > 3:
            warnings.append(
                f"Large weight range ({diag['weight_range_orders']:.1f} "
                "orders of magnitude): a few frames dominate the ensemble."
            )
        if self.chi_squared_final > 2 * len(self.observables):
            warnings.append(
                f"High final Chi-squared ({self.chi_squared_final:.2f}): "
                "poor fit; observables may be incompatible with the ensemble."
            )

        diag["warnings"] = warnings
        diag["status"] = "OK" if len(warnings) == 0 else "WARNING"
        return diag

    def print_diagnostics(self, warn_threshold: float = 0.5):
        """Print a formatted diagnostic report for this result.

        Parameters
        ----------
        warn_threshold : float, optional
            Passed through to :meth:`diagnostics`. Default 0.5.
        """
        diag = self.diagnostics(warn_threshold)
        n = len(self.weights)
        key_w, val_w = 32, 14

        print("\n" + "=" * 60)
        print("BME DIAGNOSTIC REPORT")
        print("=" * 60)
        print(f"\nOptimization Status: {self.message}")
        print(f"Success: {self.success}")
        print(f"Iterations: {self.n_iterations}")
        print("\nChi-squared:")
        print(f"  {'Initial':<{key_w}} {self.chi_squared_initial:>{val_w}.4f}")
        print(f"  {'Final':<{key_w}} {self.chi_squared_final:>{val_w}.4f}")
        print(
            f"  {'Improvement':<{key_w}} {diag['chi2_improvement']:>{val_w}.4f}"
            f" ({diag['chi2_improvement_pct']:.1f}%)"
        )
        print("\nEnsemble Diversity:")
        print(f"  {'Phi (entropy fraction)':<{key_w}} {self.phi:>{val_w}.4f}")
        print(
            f"  {'N_eff (entropy-based)':<{key_w}} "
            f"{diag['neff_entropy']:>{val_w}.1f}  / {n}"
        )
        print(
            f"  {'N_eff (1/sum w^2, Renyi-2)':<{key_w}} "
            f"{diag['neff_renyi2']:>{val_w}.1f}  / {n}"
        )
        print(f"  {'Theta':<{key_w}} {self.theta:>{val_w}.4f}")
        if self.scale is not None:
            print(f"  {'iBME scale':<{key_w}} {self.scale:>{val_w}.4g}")
            print(f"  {'iBME offset':<{key_w}} {self.offset:>{val_w}.4g}")
        print("\nWeight Distribution:")
        print(f"  {'Min':<{key_w}} {diag['weight_min']:>{val_w}.2e}")
        print(f"  {'Max':<{key_w}} {diag['weight_max']:>{val_w}.2e}")
        print(f"  {'Std Dev':<{key_w}} {diag['weight_std']:>{val_w}.2e}")
        print(
            f"  {'Range (orders of magnitude)':<{key_w}} "
            f"{diag['weight_range_orders']:>{val_w}.1f}"
        )
        if len(diag["warnings"]) > 0:
            print(f"\nWARNINGS ({len(diag['warnings'])}):")
            for i, warning in enumerate(diag["warnings"], 1):
                print(f"  {i}. {warning}")
        else:
            print(f"\nStatus: {diag['status']} - No issues detected")
        print("=" * 60)


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
@dataclass
class ThetaScanResult:
    """Container for an L-curve theta scan.

    Attributes
    ----------
    theta_values : numpy.ndarray
        Scanned theta grid.
    chi_squared_values : numpy.ndarray
        Final chi-squared at each theta.
    phi_values : numpy.ndarray
        Effective-fraction (phi) at each theta.
    kl_divergence_values : numpy.ndarray
        Relative entropy at each theta.
    results : list of BMEResult
        Per-theta fit results.
    optimal_theta : float
        Selected theta (L-curve knee).
    optimal_idx : int
        Index of ``optimal_theta`` in ``theta_values``.
    method : str
        Human-readable knee-selection method used.
    """

    theta_values: np.ndarray
    chi_squared_values: np.ndarray
    phi_values: np.ndarray
    kl_divergence_values: np.ndarray
    results: List[BMEResult]
    optimal_theta: float
    optimal_idx: int
    method: str

    def print_summary(self):
        """Print a formatted summary of the theta scan."""
        print("\n" + "=" * 60)
        print("THETA SCAN SUMMARY")
        print("=" * 60)
        print(
            f"\nScan range: {self.theta_values[0]:.4f} to {self.theta_values[-1]:.4f}"
        )
        print(f"Number of points: {len(self.theta_values)}")
        print(f"Method: {self.method}\n")
        opt = self.optimal_idx
        print(f"RECOMMENDED THETA: {self.optimal_theta:.4f}")
        print(f"Chi squared:      {self.chi_squared_values[opt]:.4f}")
        print(f"Phi (N_eff):      {self.phi_values[opt]:.4f}")
        print(f"Relative entropy: {self.kl_divergence_values[opt]:.4f}")
        print("=" * 60)


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
def _find_knee_perpendicular(x: np.ndarray, y: np.ndarray) -> int:
    """Knee index by maximum perpendicular distance to the endpoint chord."""
    x_n = (x - x.min()) / (x.max() - x.min() + 1e-10)
    y_n = (y - y.min()) / (y.max() - y.min() + 1e-10)

    p1 = np.array([x_n[0], y_n[0]])
    p2 = np.array([x_n[-1], y_n[-1]])
    line_vec = p2 - p1
    line_len = np.linalg.norm(line_vec)
    if line_len < 1e-10:
        return len(x) // 2
    line_unit = line_vec / line_len

    distances = []
    for i in range(len(x_n)):
        point = np.array([x_n[i], y_n[i]])
        vec = point - p1
        proj = p1 + np.dot(vec, line_unit) * line_unit
        distances.append(np.linalg.norm(point - proj))
    return int(np.argmax(distances))


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
def _find_knee_curvature(x: np.ndarray, y: np.ndarray) -> int:
    """Knee index by maximum Menger curvature (3-point estimate)."""
    x_n = (x - x.min()) / (x.max() - x.min() + 1e-10)
    y_n = (y - y.min()) / (y.max() - y.min() + 1e-10)
    n = len(x_n)
    curvature = np.zeros(n)
    for i in range(1, n - 1):
        p0 = np.array([x_n[i - 1], y_n[i - 1]])
        p1 = np.array([x_n[i], y_n[i]])
        p2 = np.array([x_n[i + 1], y_n[i + 1]])
        v1 = p1 - p0
        v2 = p2 - p1
        area = abs(v1[0] * v2[1] - v1[1] * v2[0]) / 2.0
        a = np.linalg.norm(p2 - p1)
        b = np.linalg.norm(p0 - p2)
        c = np.linalg.norm(p1 - p0)
        if a * b * c > 1e-10:
            curvature[i] = 4 * area / (a * b * c)
    if n > 2:
        curvature[0] = curvature[1]
        curvature[-1] = curvature[-2]
    return int(np.argmax(curvature))


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
def find_optimal_theta(
    chi_squared_values: np.ndarray,
    kl_divergence_values: np.ndarray,
    method: str = "perpendicular",
) -> Tuple[int, str]:
    """Select the L-curve knee from chi-squared vs. relative-entropy.

    Parameters
    ----------
    chi_squared_values : numpy.ndarray
        Final chi-squared per theta.
    kl_divergence_values : numpy.ndarray
        Relative entropy per theta.
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
        return _find_knee_curvature(
            chi_squared_values, kl_divergence_values
        ), "Menger curvature"
    elif method == "perpendicular":
        return _find_knee_perpendicular(
            chi_squared_values, kl_divergence_values
        ), "Perpendicular distance"
    raise SSException(
        f"Unknown method: {method}, must be 'curvature' or 'perpendicular'"
    )


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
def theta_scan(
    observables: List[ExperimentalObservable],
    calculated_values: np.ndarray,
    reweighter: str = "bme",
    theta_range: Union[Tuple[float, float], np.ndarray] = (0.01, 10.0),
    n_points: int = 15,
    log_scale: bool = True,
    initial_weights: Optional[np.ndarray] = None,
    method: str = "perpendicular",
    verbose: bool = False,
    fit_kwargs: Optional[dict] = None,
) -> ThetaScanResult:
    """Scan a grid of theta values for BME or iBME and pick the L-curve knee.

    Parameters
    ----------
    observables : list of ExperimentalObservable
        Experimental observables.
    calculated_values : numpy.ndarray
        Per-frame calculated values, shape ``(n_frames, n_observables)``.
    reweighter : str, optional
        ``"bme"`` (default) or ``"ibme"``.
    theta_range : tuple or numpy.ndarray, optional
        ``(min, max)`` for a generated grid of ``n_points``, or an explicit
        1D array of theta values. Default ``(0.01, 10.0)``.
    n_points : int, optional
        Number of grid points when ``theta_range`` is a tuple. Default 15.
    log_scale : bool, optional
        Sample logarithmically when ``theta_range`` is a tuple. Default True.
    initial_weights : numpy.ndarray, optional
        Prior weights. Uniform if None.
    method : str, optional
        Knee-selection method (see :func:`find_optimal_theta`).
    verbose : bool, optional
        Print per-theta progress. Default False.
    fit_kwargs : dict, optional
        Extra keyword arguments forwarded to the reweighter's ``fit``
        (e.g. ``ftol``/``max_ibme_iterations`` for iBME).

    Returns
    -------
    ThetaScanResult

    Raises
    ------
    SSException
        If ``reweighter`` is unknown, ``n_points`` < 1, or a log-scale grid
        has non-positive endpoints.
    """
    if reweighter not in ("bme", "ibme"):
        raise SSException(f"Unknown reweighter: {reweighter}, must be 'bme' or 'ibme'")

    if isinstance(theta_range, (tuple, list)):
        if n_points < 1:
            raise SSException(
                f"n_points must be >= 1 for a tuple theta_range, got {n_points}"
            )
        if log_scale and (theta_range[0] <= 0 or theta_range[1] <= 0):
            raise SSException(
                "theta_range endpoints must be positive when log_scale=True, "
                f"got {theta_range}"
            )
        if log_scale:
            theta_values = np.logspace(
                np.log10(theta_range[0]), np.log10(theta_range[1]), n_points
            )
        else:
            theta_values = np.linspace(theta_range[0], theta_range[1], n_points)
    else:
        theta_values = np.asarray(theta_range, dtype=np.float64)

    fit_kwargs = dict(fit_kwargs) if fit_kwargs else {}

    chi2_vals: List[float] = []
    phi_vals: List[float] = []
    kl_vals: List[float] = []
    results: List[BMEResult] = []

    for i, theta in enumerate(theta_values):
        if verbose:
            print(f"Processing theta {i + 1}/{len(theta_values)}: {theta:.4f}")
        if reweighter == "bme":
            rw = BME(observables, calculated_values, initial_weights)
            res = rw.fit(
                theta=float(theta),
                auto_theta=False,
                verbose=False,
                **fit_kwargs,
            )
        else:
            rw = iBME(observables, calculated_values, initial_weights)
            res = rw.fit(theta=float(theta), verbose=False, **fit_kwargs)

        results.append(res)
        chi2_vals.append(res.chi_squared_final)
        phi_vals.append(res.phi)
        kl_vals.append(res.kl_divergence)

    chi2_arr = np.array(chi2_vals)
    phi_arr = np.array(phi_vals)
    kl_arr = np.array(kl_vals)

    optimal_idx, method_name = find_optimal_theta(chi2_arr, kl_arr, method=method)

    return ThetaScanResult(
        theta_values=theta_values,
        chi_squared_values=chi2_arr,
        phi_values=phi_arr,
        kl_divergence_values=kl_arr,
        results=results,
        optimal_theta=float(theta_values[optimal_idx]),
        optimal_idx=optimal_idx,
        method=method_name,
    )


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
def _validate_inputs(observables, calculated_values, initial_weights):
    """Validate constructor arguments shared by :class:`BME` and :class:`iBME`.

    Raises
    ------
    SSException
        If any argument is malformed or dimensions are inconsistent.
    """
    if not isinstance(observables, (list, tuple)) or len(observables) == 0:
        raise SSException("observables must be a non-empty list")
    if not all(isinstance(obs, ExperimentalObservable) for obs in observables):
        raise SSException("All observables must be ExperimentalObservable instances")
    if not isinstance(calculated_values, np.ndarray):
        raise SSException("calculated_values must be a numpy array")
    if calculated_values.ndim != 2:
        raise SSException("calculated_values must be 2D (n_frames, n_observables)")
    if calculated_values.shape[1] != len(observables):
        raise SSException(
            f"Number of observables ({len(observables)}) must match "
            f"calculated_values columns ({calculated_values.shape[1]})"
        )
    if initial_weights is not None:
        if not isinstance(initial_weights, np.ndarray):
            raise SSException("initial_weights must be a numpy array")
        if len(initial_weights) != calculated_values.shape[0]:
            raise SSException("initial_weights length must match number of frames")


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
class BME:
    """Bayesian Maximum Entropy reweighting.

    Finds the minimally biased frame weights that bring ensemble averages
    into agreement with experimental observables, by maximising entropy
    relative to the prior subject to a Gaussian-error data restraint.

    Parameters
    ----------
    observables : list of ExperimentalObservable
        Experimental observables to fit against.
    calculated_values : numpy.ndarray
        Per-frame calculated observable values, shape
        ``(n_frames, n_observables)``.
    initial_weights : numpy.ndarray, optional
        Prior frame weights. If None, uniform weights are used. Internally
        normalized to sum to 1.

    Raises
    ------
    SSException
        If inputs are malformed (see :func:`_validate_inputs`).
    """

    def __init__(
        self,
        observables: List[ExperimentalObservable],
        calculated_values: np.ndarray,
        initial_weights: Optional[np.ndarray] = None,
    ):
        _validate_inputs(observables, calculated_values, initial_weights)

        self.observables = list(observables)
        self.calculated_values = calculated_values
        self.n_frames = calculated_values.shape[0]
        self.n_observables = len(observables)

        if initial_weights is None:
            self.initial_weights = np.ones(self.n_frames) / self.n_frames
        else:
            self.initial_weights = initial_weights / np.sum(initial_weights)

        self._lambdas = np.random.normal(
            loc=0.0, scale=LAMBDA_INIT_SCALE, size=self.n_observables
        ).astype(np.float64)
        self._bounds = [obs.get_bounds() for obs in self.observables]

        self._exp_values = np.array(
            [obs.value for obs in self.observables], dtype=np.float64
        )
        self._exp_sigma2 = np.array(
            [obs.uncertainty**2 for obs in self.observables],
            dtype=np.float64,
        )

        self._result: Optional[BMEResult] = None
        self._theta_scan_result: Optional[ThetaScanResult] = None
        self._theta: Optional[float] = None

    # ------------------------------------------------------------------
    def _compute_chi_squared(self, weights: np.ndarray) -> float:
        """Constraint-aware reduced chi-squared for a weight vector.

        ``equality`` observables always penalize deviations; ``upper`` /
        ``lower`` only penalize the disallowed side.

        Parameters
        ----------
        weights : numpy.ndarray
            Frame weights, shape ``(n_frames,)``.

        Returns
        -------
        float
            Mean of ``(diff / sigma)^2`` over observables.
        """
        chi_squared = 0.0
        for idx, obs in enumerate(self.observables):
            calc_avg = np.sum(self.calculated_values[:, idx] * weights)
            diff = calc_avg - obs.value
            if obs.constraint == "equality":
                penalize = True
            elif obs.constraint == "upper":
                penalize = diff > 0
            else:  # "lower"
                penalize = diff < 0
            if penalize:
                chi_squared += (diff / obs.uncertainty) ** 2
        return chi_squared / self.n_observables

    # ------------------------------------------------------------------
    def _objective_and_gradient(self, lambdas: np.ndarray) -> Tuple[float, np.ndarray]:
        """Maxent objective and gradient with a Gaussian-error prior.

        Implements ``L = lambda^T O_exp + (theta/2) lambda^T Sigma^2 lambda
        + log Z`` and its gradient, scaled by ``1/theta`` for numerical
        stability.

        Parameters
        ----------
        lambdas : numpy.ndarray
            Current Lagrange multipliers.

        Returns
        -------
        tuple
            ``(objective/theta, gradient/theta)``.
        """
        log_w = -np.sum(lambdas * self.calculated_values, axis=1) + np.log(
            self.initial_weights
        )
        log_z = logsumexp(log_w)
        probs = np.exp(log_w - log_z)
        avg_calc = np.sum(probs[:, np.newaxis] * self.calculated_values, axis=0)

        regularization = self._theta / 2 * np.sum(lambdas**2 * self._exp_sigma2)
        constraint = np.dot(lambdas, self._exp_values)
        objective = log_z + constraint + regularization
        gradient = (
            self._exp_values + self._theta * lambdas * self._exp_sigma2 - avg_calc
        )
        return objective / self._theta, gradient / self._theta

    # ------------------------------------------------------------------
    def _run_single_theta_optimization(
        self, max_iterations: int, optimizer: str, verbose: bool
    ) -> BMEResult:
        """Run the maxent optimization for the currently set ``self._theta``.

        Parameters
        ----------
        max_iterations : int
            Maximum optimizer iterations.
        optimizer : str
            scipy.optimize.minimize method (must support ``jac`` + bounds).
        verbose : bool
            Print progress.

        Returns
        -------
        BMEResult
        """
        chi2_initial = self._compute_chi_squared(self.initial_weights)
        if verbose:
            print("BME Optimization")
            print(f"  Theta: {self._theta}")
            print(f"  Observables: {self.n_observables}")
            print(f"  Frames: {self.n_frames}")
            print(f"  Chi-squared initial: {chi2_initial:.4f}")

        opt = minimize(
            self._objective_and_gradient,
            self._lambdas,
            options={"maxiter": max_iterations},
            method=optimizer,
            jac=True,
            bounds=self._bounds,
        )

        metadata = {
            "optimizer": optimizer,
            "max_iterations": max_iterations,
            "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }

        if opt.success:
            log_w = -np.sum(
                opt.x[np.newaxis, :] * self.calculated_values, axis=1
            ) + np.log(self.initial_weights)
            weights = np.exp(log_w - logsumexp(log_w))
            chi2_final = self._compute_chi_squared(weights)
            rel = float(np.sum(rel_entr(weights, self.initial_weights)))
            phi = float(np.exp(-rel))
            if verbose:
                print(f"  Optimization successful (iters: {opt.nit})")
                print(f"  Chi-squared final: {chi2_final:.4f}")
                print(f"  Effective fraction (phi): {phi:.4f}")
            return BMEResult(
                weights=weights,
                initial_weights=self.initial_weights.copy(),
                lambdas=opt.x.copy(),
                chi_squared_initial=chi2_initial,
                chi_squared_final=chi2_final,
                phi=phi,
                n_iterations=opt.nit,
                success=True,
                message=str(opt.message),
                theta=self._theta,
                observables=self.observables,
                calculated_values=self.calculated_values,
                metadata=metadata,
            )

        if verbose:
            print("  Optimization failed")
            print(f"  Message: {opt.message}")
        return BMEResult(
            weights=self.initial_weights.copy(),
            initial_weights=self.initial_weights.copy(),
            lambdas=opt.x.copy(),
            chi_squared_initial=chi2_initial,
            chi_squared_final=np.nan,
            phi=np.nan,
            n_iterations=getattr(opt, "nit", -1),
            success=False,
            message=str(opt.message),
            theta=self._theta,
            observables=self.observables,
            calculated_values=self.calculated_values,
            metadata=metadata,
        )

    # ------------------------------------------------------------------
    def fit(
        self,
        max_iterations: int = DEFAULT_MAX_ITERATIONS,
        optimizer: str = DEFAULT_OPTIMIZER,
        verbose: bool = True,
        *,
        theta: Optional[float] = None,
        auto_theta: bool = True,
        theta_scan_kwargs: Optional[dict] = None,
    ) -> BMEResult:
        """Fit BME weights at a fixed or automatically selected theta.

        Parameters
        ----------
        max_iterations : int, optional
            Maximum optimizer iterations. Default ``DEFAULT_MAX_ITERATIONS``.
        optimizer : str, optional
            scipy.optimize.minimize method. Default ``"L-BFGS-B"``.
        verbose : bool, optional
            Print progress. Default True.
        theta : float, optional
            If given, use this theta (no scan). Must be positive.
        auto_theta : bool, optional
            If True and ``theta`` is None, run an L-curve scan to pick theta.
            Default True.
        theta_scan_kwargs : dict, optional
            Extra keyword arguments forwarded to :meth:`scan_theta`.

        Returns
        -------
        BMEResult

        Raises
        ------
        SSException
            If ``theta`` is not positive.
        """
        if theta is not None:
            if theta <= 0:
                raise SSException(f"theta must be positive, got {theta}")
            self._theta = float(theta)
            self._theta_scan_result = None
            if verbose:
                print(f"[BME] Using manual theta = {self._theta:.4g}.")
        elif auto_theta:
            if verbose:
                print("[BME] Auto theta: running L-curve scan...")
            scan_kwargs = dict(theta_range=(0.01, 10.0), n_points=15)
            if theta_scan_kwargs:
                scan_kwargs.update(theta_scan_kwargs)
            scan = self.scan_theta(**scan_kwargs)
            opt_idx = scan.optimal_idx
            self._theta = float(scan.optimal_theta)
            self._result = scan.results[opt_idx]
            self._theta_scan_result = scan
            self._lambdas = self._result.lambdas.copy()
            if verbose:
                print(
                    f"[BME] Selected theta={self._theta:.4g} via "
                    f"{scan.method} (chi2_final="
                    f"{self._result.chi_squared_final:.3f}, "
                    f"phi={self._result.phi:.3f})."
                )
            return self._result
        else:
            if self._theta is None:
                self._theta = DEFAULT_THETA
            if verbose:
                print(
                    f"[BME] auto_theta=False, no manual theta; using "
                    f"theta={self._theta:.4g}."
                )

        self._result = self._run_single_theta_optimization(
            max_iterations=max_iterations,
            optimizer=optimizer,
            verbose=verbose,
        )
        return self._result

    # ------------------------------------------------------------------
    def scan_theta(
        self,
        theta_range: Union[Tuple[float, float], np.ndarray] = (0.01, 10.0),
        n_points: int = 15,
        log_scale: bool = True,
        method: str = "perpendicular",
        verbose: bool = False,
    ) -> ThetaScanResult:
        """Run an L-curve theta scan using this instance's data.

        Parameters
        ----------
        theta_range, n_points, log_scale, method, verbose
            See :func:`theta_scan`.

        Returns
        -------
        ThetaScanResult
        """
        scan = theta_scan(
            observables=self.observables,
            calculated_values=self.calculated_values,
            reweighter="bme",
            theta_range=theta_range,
            n_points=n_points,
            log_scale=log_scale,
            initial_weights=self.initial_weights,
            method=method,
            verbose=verbose,
        )
        self._theta_scan_result = scan
        return scan

    # ------------------------------------------------------------------
    def predict(self, calculated_values: np.ndarray) -> np.ndarray:
        """Weighted means of arbitrary observables using the fitted weights.

        Parameters
        ----------
        calculated_values : numpy.ndarray
            Per-frame values, shape ``(n_frames, n_observables)``.

        Returns
        -------
        numpy.ndarray
            Weighted means, shape ``(n_observables,)``.

        Raises
        ------
        SSException
            If :meth:`fit` has not succeeded, or frame count mismatches.
        """
        if self._result is None or not self._result.success:
            raise SSException(
                "Model has not been successfully fitted. Call fit() first."
            )
        return self._result.predict(calculated_values)

    @property
    def result(self) -> Optional[BMEResult]:
        """The most recent :class:`BMEResult`, or None."""
        return self._result

    @property
    def theta(self) -> Optional[float]:
        """Theta used in the most recent fit, or None."""
        return self._theta

    @property
    def theta_scan_result(self) -> Optional[ThetaScanResult]:
        """The most recent :class:`ThetaScanResult`, or None."""
        return self._theta_scan_result


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
class iBME:
    """Iterative Bayesian Maximum Entropy reweighting.

    At each iteration a (optionally error-weighted) linear regression of the
    ensemble-averaged calculated values against the experimental values is
    used to refit a global scale ``alpha`` and offset ``beta``; the
    calculated values are rescaled (``calc -> alpha*calc + beta``) and a
    standard BME maxent step is run. The procedure repeats until the change
    in chi-squared drops below ``ftol``. iBME is appropriate when there is an
    unknown global scale/offset between calculated and experimental data
    (e.g. SAXS). Ported from the reference implementation of Bottaro et al.

    Parameters
    ----------
    observables : list of ExperimentalObservable
        Experimental observables to fit against.
    calculated_values : numpy.ndarray
        Per-frame calculated values, shape ``(n_frames, n_observables)``.
    initial_weights : numpy.ndarray, optional
        Prior frame weights. Uniform if None.

    Raises
    ------
    SSException
        If inputs are malformed (see :func:`_validate_inputs`).
    """

    def __init__(
        self,
        observables: List[ExperimentalObservable],
        calculated_values: np.ndarray,
        initial_weights: Optional[np.ndarray] = None,
    ):
        _validate_inputs(observables, calculated_values, initial_weights)

        self.observables = list(observables)
        self.calculated_values = calculated_values
        self.n_frames = calculated_values.shape[0]
        self.n_observables = len(observables)

        if initial_weights is None:
            self.initial_weights = np.ones(self.n_frames) / self.n_frames
        else:
            self.initial_weights = initial_weights / np.sum(initial_weights)

        self._exp_values = np.array(
            [obs.value for obs in self.observables], dtype=np.float64
        )
        self._exp_sigma = np.array(
            [obs.uncertainty for obs in self.observables], dtype=np.float64
        )

        self._result: Optional[BMEResult] = None
        self._theta: Optional[float] = None
        self._theta_scan_result: Optional[ThetaScanResult] = None

    # ------------------------------------------------------------------
    def fit(
        self,
        theta: float,
        ftol: float = 0.01,
        max_ibme_iterations: int = 50,
        fit_offset: bool = True,
        lr_weights: bool = True,
        max_iterations: int = DEFAULT_MAX_ITERATIONS,
        optimizer: str = DEFAULT_OPTIMIZER,
        verbose: bool = True,
    ) -> BMEResult:
        """Run iterative BME at a fixed theta.

        Parameters
        ----------
        theta : float
            Regularisation strength. Must be positive.
        ftol : float, optional
            Convergence tolerance on ``|delta chi-squared|`` between
            iterations. Default 0.01.
        max_ibme_iterations : int, optional
            Maximum number of scale/offset+BME iterations. Default 50.
        fit_offset : bool, optional
            If True fit a scale and offset; if False fit scale only.
            Default True.
        lr_weights : bool, optional
            If True weight the linear regression by ``1/sigma^2``; otherwise
            use uniform weights. Default True.
        max_iterations : int, optional
            Maximum optimizer iterations for each inner BME step.
        optimizer : str, optional
            scipy.optimize.minimize method for each inner BME step.
        verbose : bool, optional
            Print per-iteration progress. Default True.

        Returns
        -------
        BMEResult
            Result with the final weights, ``scale``/``offset`` (net,
            relative to the original input), the ``ibme_iterations`` log,
            ``chi_squared_initial`` from the first iteration's pre-fit
            chi-squared, and ``phi`` computed against the original prior.

        Raises
        ------
        SSException
            If ``theta`` is not positive.
        """
        if theta <= 0:
            raise SSException(f"theta must be positive, got {theta}")
        self._theta = float(theta)

        w0 = self.initial_weights.copy()
        current_weights = w0.copy()
        # Working copy of the calculated values; rescaled in place each
        # iteration so the scale/offset accumulates (matches reference iBME).
        calc = self.calculated_values.astype(np.float64).copy()

        if lr_weights:
            lr_w = 1.0 / self._exp_sigma**2
        else:
            lr_w = np.ones(self.n_observables)

        # Accumulated (net) scale/offset relative to the original input.
        net_scale = 1.0
        net_offset = 0.0

        iterations: List[dict] = []
        chi2_initial = np.nan
        chi2_old = np.nan
        last_result: Optional[BMEResult] = None
        it = 0

        for it in range(max_ibme_iterations):
            calc_avg = np.sum(calc * current_weights[:, np.newaxis], axis=0)
            alpha, beta = _weighted_linear_regression(
                calc_avg, self._exp_values, lr_w, fit_intercept=fit_offset
            )
            calc = alpha * calc + beta
            net_scale *= alpha
            net_offset = alpha * net_offset + beta

            bme = BME(self.observables, calc, w0)
            res = bme.fit(
                theta=self._theta,
                auto_theta=False,
                verbose=False,
                max_iterations=max_iterations,
                optimizer=optimizer,
            )

            if it == 0:
                chi2_initial = res.chi_squared_initial

            current_weights = res.weights.copy()
            chi2_new = res.chi_squared_final
            diff = abs(chi2_old - chi2_new)
            chi2_old = chi2_new
            last_result = res

            iterations.append(
                {
                    "iteration": it,
                    "scale": alpha,
                    "offset": beta,
                    "chi_squared": chi2_new,
                    "diff": diff,
                }
            )
            if verbose:
                print(
                    f"Iteration {it:3d}  scale={alpha:8.4f}  "
                    f"offset={beta:8.4f}  chi2={chi2_new:8.4f}  "
                    f"diff={diff:8.2e}"
                )

            if np.isfinite(diff) and diff < ftol:
                if verbose:
                    print(
                        f"iBME converged below tolerance {ftol:.2e} after "
                        f"{it + 1} iterations"
                    )
                break

        phi = float(np.exp(-_srel(w0, current_weights)))
        success = last_result is not None and last_result.success
        metadata = {
            "optimizer": optimizer,
            "max_iterations": max_iterations,
            "max_ibme_iterations": max_ibme_iterations,
            "ftol": ftol,
            "fit_offset": fit_offset,
            "lr_weights": lr_weights,
            "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }

        self._result = BMEResult(
            weights=current_weights,
            initial_weights=w0,
            lambdas=(
                last_result.lambdas.copy()
                if last_result is not None
                else np.zeros(self.n_observables)
            ),
            chi_squared_initial=chi2_initial,
            chi_squared_final=chi2_old,
            phi=phi,
            n_iterations=it + 1,
            success=success,
            message=(
                last_result.message if last_result is not None else "no iterations run"
            ),
            theta=self._theta,
            observables=self.observables,
            calculated_values=calc,
            metadata=metadata,
            scale=net_scale,
            offset=net_offset,
            ibme_iterations=iterations,
        )
        return self._result

    # ------------------------------------------------------------------
    def scan_theta(
        self,
        theta_range: Union[Tuple[float, float], np.ndarray] = (0.01, 10.0),
        n_points: int = 15,
        log_scale: bool = True,
        method: str = "perpendicular",
        verbose: bool = False,
        fit_kwargs: Optional[dict] = None,
    ) -> ThetaScanResult:
        """Run an L-curve theta scan using iterative BME.

        Parameters
        ----------
        theta_range, n_points, log_scale, method, verbose
            See :func:`theta_scan`.
        fit_kwargs : dict, optional
            Extra keyword arguments forwarded to :meth:`fit` for each theta
            (e.g. ``ftol``, ``max_ibme_iterations``, ``fit_offset``).

        Returns
        -------
        ThetaScanResult
        """
        scan = theta_scan(
            observables=self.observables,
            calculated_values=self.calculated_values,
            reweighter="ibme",
            theta_range=theta_range,
            n_points=n_points,
            log_scale=log_scale,
            initial_weights=self.initial_weights,
            method=method,
            verbose=verbose,
            fit_kwargs=fit_kwargs,
        )
        self._theta_scan_result = scan
        return scan

    # ------------------------------------------------------------------
    def predict(self, calculated_values: np.ndarray) -> np.ndarray:
        """Weighted means of arbitrary observables using the fitted weights.

        Parameters
        ----------
        calculated_values : numpy.ndarray
            Per-frame values, shape ``(n_frames, n_observables)``.

        Returns
        -------
        numpy.ndarray
            Weighted means, shape ``(n_observables,)``.

        Raises
        ------
        SSException
            If :meth:`fit` has not succeeded, or frame count mismatches.
        """
        if self._result is None or not self._result.success:
            raise SSException(
                "Model has not been successfully fitted. Call fit() first."
            )
        return self._result.predict(calculated_values)

    def get_ibme_weights(self) -> np.ndarray:
        """Final iBME weights.

        Returns
        -------
        numpy.ndarray

        Raises
        ------
        SSException
            If :meth:`fit` has not been called.
        """
        if self._result is None:
            raise SSException("iBME weights not available. Call fit() first.")
        return self._result.weights

    def get_ibme_stats(self) -> List[dict]:
        """Per-iteration iBME log (scale, offset, chi-squared, diff).

        Returns
        -------
        list of dict

        Raises
        ------
        SSException
            If :meth:`fit` has not been called.
        """
        if self._result is None:
            raise SSException("iBME stats not available. Call fit() first.")
        return self._result.ibme_iterations

    @property
    def result(self) -> Optional[BMEResult]:
        """The most recent :class:`BMEResult`, or None."""
        return self._result

    @property
    def theta(self) -> Optional[float]:
        """Theta used in the most recent fit, or None."""
        return self._theta

    @property
    def theta_scan_result(self) -> Optional[ThetaScanResult]:
        """The most recent :class:`ThetaScanResult`, or None."""
        return self._theta_scan_result
