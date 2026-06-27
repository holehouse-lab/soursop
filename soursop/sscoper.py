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
sscoper - Convex Optimization for Ensemble Reweighting (COPER) and iCOPER.

This module reweights molecular ensembles against experimental observables
using the **COPER** method of Leung *et al.* (Leung, Bignucolo, Aregger,
Dames, Mazur, Bernèche & Grzesiek, *J. Chem. Theory Comput.* **12**,
383-394, 2016). COPER solves the *primal constrained convex* problem
directly over the **N frame weights** in two steps:

1. **Chi-squared minimisation / feasibility**: minimise chi-squared(w) over
   the simplex (``0 <= w_i <= 1``, ``sum w_i = 1``). If the minimum is below
   the chosen limit, the minimiser is a feasible interior point; otherwise
   the data cannot be satisfied by any reweighting and the problem **has
   no solution** (a finding of the paper).
2. **Entropy maximisation**: maximise the Shannon entropy
   ``S = -sum w_i ln w_i`` (generalised here to the relative entropy
   ``-KL(w || w0)`` for an arbitrary prior ``w0``), subject to a *hard*
   constraint ``chi-squared <= limit`` plus the simplex, starting from
   the step-1 point.

There is no regularisation parameter ``theta`` (contrast :mod:`soursop.ssbme`,
which solves the dual penalty ``L = (m/2) chi^2 - theta * S``); the knob
is the chi-squared limit, default ``1`` (i.e. agreement with experiment
within error). The entropy reduction ``delta_S = S(w) - S(w0) <= 0`` is a
well-defined measure of the information content of the data, equal in
magnitude to the mean free energy change ``<delta_G>/kT`` (paper eq 10).
Each frame's reweighting factor is ``r_i = w_i / w0_i`` (paper eq 9a).

Two reweighters are provided:

``COPER``
    Standard COPER reweighting. Use when the calculated observables are
    already on the same scale as the experimental values. Per-data-type
    chi-squared constraints (``chi^2_alpha <= limit`` for each group, as
    in the paper's RDC + J-coupling fits) are supported via the
    :attr:`~soursop.ssutils.ExperimentalObservable.group` field.

``iCOPER``
    Iterative COPER. At every iteration a (weighted) linear regression of
    the ensemble-averaged calculated values against the experimental
    values is used to refit a global scale and offset before the
    constrained-maxent step. The analogue of :class:`soursop.ssbme.iBME`
    for data with an unknown global scale/offset (e.g. SAXS intensities).

``ExperimentalObservable``
    Container for a single experimental data point. Re-exported from
    :mod:`soursop.ssutils` so the syntax is identical to BME / iBME.

The chi-squared limit can be selected automatically with
:func:`chi2_limit_scan` (also exposed as ``COPER.scan_chi2_limit`` /
``iCOPER.scan_chi2_limit``), the COPER analogue of BME's L-curve
``theta_scan``.

Example
-------
>>> import numpy as np
>>> from soursop.sscoper import COPER, ExperimentalObservable
>>> rng = np.random.default_rng(0)
>>> calc = rng.normal([24, 65], [3, 6], size=(600, 2))
>>> obs = [ExperimentalObservable(23.0, 1.0, name="Rg"),
...        ExperimentalObservable(60.0, 2.0, name="Ree")]
>>> result = COPER(obs, calc).fit(chi2_limit=1.0, verbose=False)
>>> reweighted_means = result.predict(calc)

References
----------
Primary reference for the COPER method this module is adapted from:

* Leung, H. T. A., Bignucolo, O., Aregger, R., Dames, S. A., Mazur, A.,
  Bernèche, S. & Grzesiek, S. *A Rigorous and Efficient Method To
  Reweight Very Large Conformational Ensembles Using Average
  Experimental Data and To Determine Their Relative Information
  Content.* J. Chem. Theory Comput. **12**, 383-394 (2016).
  https://doi.org/10.1021/acs.jctc.5b00759

Closely related and foundational work:

* Kullback, S. & Leibler, R. A. *On Information and Sufficiency.* Ann.
  Math. Stat. **22**, 79-86 (1951).
* Bottaro, S., Bengtsen, T. & Lindorff-Larsen, K. *Integrating Molecular
  Simulation and Experimental Data: A Bayesian/Maximum Entropy
  Reweighting Approach.* Methods Mol. Biol. **2112**, 219-240 (2020) -
  the BME / iBME method implemented in :mod:`soursop.ssbme`, which
  differs from COPER in formulating the problem as a penalty
  (``L = (m/2) chi^2 - theta * S``) rather than a hard chi-squared
  constraint.
* Hummer, G. & Köfinger, J. *Bayesian ensemble refinement by replica
  simulations and reweighting.* J. Chem. Phys. **143**, 243150 (2015).
* Różycki, B., Kim, Y. C. & Hummer, G. *SAXS ensemble refinement of
  ESCRT-III CHMP3 conformational transitions.* Structure **19**,
  109-116 (2011) - EROS; SAXS scale/offset problem that motivates
  iCOPER.

**Author(s):** Alex Holehouse
"""

from dataclasses import dataclass, field
from datetime import datetime
from typing import List, Optional, Tuple, Union

import numpy as np
from scipy.optimize import (
    NonlinearConstraint,
    minimize,
)
from scipy.sparse.linalg import LinearOperator

from .ssexceptions import SSException
from .ssutils import (
    MIN_WEIGHT_THRESHOLD,
    ExperimentalObservable,
    constraint_chi_squared,
    find_optimal_theta,
    relative_entropy,
    validate_reweighting_inputs,
    weighted_linear_regression,
)

# Internal aliases for symbols whose canonical home is soursop.ssutils but
# which sscoper uses pervasively under shorter / iBME-style names.
_srel = relative_entropy
_weighted_linear_regression = weighted_linear_regression
_validate_inputs = validate_reweighting_inputs

__all__ = [
    "COPER",
    "iCOPER",
    "COPERResult",
    "COPERScanResult",
    "chi2_limit_scan",
    "ExperimentalObservable",
]


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Module constants
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
DEFAULT_CHI2_LIMIT = 1.0
DEFAULT_MAX_ITERATIONS = 2000
DEFAULT_OPTIMIZER = "trust-constr"
#: Numerical slack applied to the feasibility check (``chi2_min`` may exceed
#: ``chi2_limit`` by this much and still be considered feasible).
FEASIBILITY_TOLERANCE = 1e-6
#: Default convergence tolerances forwarded to scipy's trust-constr; the
#: paper's IPOPT settings were ``1e-3`` (chi^2-min) and ``1e-5`` (entropy
#: max). Slightly looser than scipy's ``1e-8`` defaults, which keeps run
#: times reasonable on large ensembles without affecting the qualitative
#: solution.
_TRUST_CONSTR_OPTIONS = {"xtol": 1e-7, "gtol": 1e-6, "verbose": 0}


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
def _shannon_entropy(weights: np.ndarray) -> float:
    """Shannon entropy ``-sum w_i ln w_i`` with a floor for zero weights.

    Parameters
    ----------
    weights : numpy.ndarray
        Frame weights summing to 1.

    Returns
    -------
    float
        ``-sum_i w_i ln(w_i)``, treating weights at or below
        :data:`MIN_WEIGHT_THRESHOLD` as contributing zero
        (``0 * log 0 = 0``).
    """
    mask = weights > MIN_WEIGHT_THRESHOLD
    if not np.any(mask):
        return 0.0
    w = weights[mask]
    return float(-np.sum(w * np.log(w)))


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
@dataclass
class COPERResult:
    """Container for a COPER / iCOPER reweighting result.

    Attributes
    ----------
    weights : numpy.ndarray
        Optimized (posterior) frame weights, summing to 1.
    initial_weights : numpy.ndarray
        Prior frame weights.
    chi_squared_initial : float
        Reduced chi-squared at the prior (max over groups when per-data-type
        constraints are used).
    chi_squared_min : float
        Reduced chi-squared at the end of the step-1 minimisation (max over
        groups). Determines feasibility: ``feasible`` is True iff
        ``chi_squared_min <= chi2_limit`` (up to ``FEASIBILITY_TOLERANCE``).
    chi_squared_final : float
        Reduced chi-squared after the step-2 entropy maximisation (max over
        groups). Equal to ``chi_squared_min`` when ``feasible`` is False.
    chi2_limit : float
        The chi-squared limit used (``chi^2_alpha <= chi2_limit`` per group).
    feasible : bool
        Whether the data can be satisfied by any reweighting at this limit.
    entropy : float
        Shannon entropy ``-sum w * ln w`` at the optimized weights.
    entropy_initial : float
        Shannon entropy at the prior.
    delta_S : float
        Entropy change ``S(w) - S(w0) <= 0`` (equal to ``-KL(w || w0)``).
        Measures the information content of the experimental data
        relative to the prior ensemble.
    mean_delta_G_kT : float
        Mean free-energy change in ``kT`` units, ``<delta_G>/kT = -delta_S``
        (Leung et al. eq 10).
    phi : float
        Fraction of effective frames, ``exp(-KL(w || w0)) in (0, 1]``.
    n_iterations : int
        Optimizer iterations of the step-2 (entropy) optimisation, or
        step-1 when infeasible.
    success : bool
        Whether the optimisation succeeded.
    message : str
        Optimizer / feasibility status message.
    observables : list of ExperimentalObservable
        Observables used in the fit.
    calculated_values : numpy.ndarray
        Calculated values used in the fit (for iCOPER, the final scaled
        array).
    metadata : dict
        Free-form metadata (optimizer, timestamp, group structure, ...).
    scale : float, optional
        Final iCOPER scale factor relative to the original input. ``None``
        for plain COPER.
    offset : float, optional
        Final iCOPER offset relative to the original input. ``None`` for
        plain COPER.
    icoper_iterations : list of dict
        Per-iteration iCOPER log (empty for plain COPER). Each entry has
        ``iteration``, ``scale``, ``offset``, ``chi_squared``, ``feasible``
        and ``diff``.
    """

    weights: np.ndarray
    initial_weights: np.ndarray
    chi_squared_initial: float
    chi_squared_min: float
    chi_squared_final: float
    chi2_limit: float
    feasible: bool
    entropy: float
    entropy_initial: float
    delta_S: float
    mean_delta_G_kT: float
    phi: float
    n_iterations: int
    success: bool
    message: str
    observables: List[ExperimentalObservable]
    calculated_values: np.ndarray
    metadata: dict
    scale: Optional[float] = None
    offset: Optional[float] = None
    icoper_iterations: List[dict] = field(default_factory=list)

    def __str__(self):
        status = "SUCCESS" if self.success else "FAILED"
        feas = "feasible" if self.feasible else "INFEASIBLE"
        lines = [
            f"COPER Result [{status}, {feas}]",
            f"  Chi-squared initial: {self.chi_squared_initial:.4f}",
            f"  Chi-squared min:     {self.chi_squared_min:.4f}",
            f"  Chi-squared final:   {self.chi_squared_final:.4f}",
            f"  Chi-squared limit:   {self.chi2_limit:.4f}",
            f"  phi (effective fraction): {self.phi:.4f}",
            f"  delta_S:             {self.delta_S:.4f}",
            f"  <delta_G>/kT:        {self.mean_delta_G_kT:.4f}",
            f"  Iterations: {self.n_iterations}",
        ]
        if self.scale is not None:
            lines.append(f"  iCOPER scale/offset: {self.scale:.4g} / {self.offset:.4g}")
        return "\n".join(lines)

    def __repr__(self):
        return self.__str__()

    @property
    def kl_divergence(self) -> float:
        """Relative entropy (KL divergence) of ``weights`` from prior."""
        if self.phi > 0:
            return -float(np.log(self.phi))
        return float(np.inf)

    @property
    def reweighting_factors(self) -> np.ndarray:
        """Per-frame reweighting factors ``r_i = w_i / w0_i``
        (Leung et al. eq 9a).
        """
        return self.weights / self.initial_weights

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
        ``neff_entropy`` (``N * phi``, the COPER ``<delta_G>/kT`` measure)
        and the Renyi-2 / participation ratio ``neff_renyi2``
        (``1 / sum w^2``), which is more sensitive to a few dominant
        weights.
        """
        diag: dict = {}
        warnings: List[str] = []
        n = len(self.weights)

        diag["neff_entropy"] = n * float(self.phi)
        diag["neff_entropy_fraction"] = float(self.phi)
        diag["neff_renyi2"] = float(1.0 / np.sum(self.weights**2))
        diag["neff_renyi2_fraction"] = float(diag["neff_renyi2"] / n)

        diag["weight_min"] = float(self.weights.min())
        diag["weight_max"] = float(self.weights.max())
        diag["weight_std"] = float(self.weights.std())
        if diag["weight_min"] > 0:
            diag["weight_range_orders"] = float(
                np.log10(diag["weight_max"] / diag["weight_min"])
            )
        else:
            diag["weight_range_orders"] = float(np.inf)

        rfac = self.reweighting_factors
        diag["rfac_min"] = float(rfac.min())
        diag["rfac_max"] = float(rfac.max())

        diag["chi2_improvement"] = self.chi_squared_initial - self.chi_squared_final
        diag["chi2_improvement_pct"] = (
            (diag["chi2_improvement"] / self.chi_squared_initial) * 100.0
            if self.chi_squared_initial > 0
            else float(np.nan)
        )

        if not self.feasible:
            warnings.append(
                f"Infeasible: minimum achievable chi-squared "
                f"({self.chi_squared_min:.3f}) exceeds the limit "
                f"({self.chi2_limit:.3f}). The data cannot be reproduced by "
                "any reweighting of the prior ensemble; check sampling, the "
                "force field, or the experimental uncertainties."
            )
        if self.phi < warn_threshold:
            warnings.append(
                f"Low Phi ({self.phi:.3f} < {warn_threshold}): significant "
                "loss of ensemble diversity. Consider loosening the chi^2 "
                "limit or the observable uncertainties, or improving "
                "sampling."
            )
        if diag["neff_renyi2"] < 0.1 * n:
            warnings.append(
                f"Low effective sample size (1/sum w^2) "
                f"({diag['neff_renyi2']:.1f} / {n}): only "
                f"~{diag['neff_renyi2_fraction'] * 100:.1f}% of frames "
                "are effectively used."
            )
        if diag["weight_range_orders"] > 3:
            warnings.append(
                f"Large weight range ({diag['weight_range_orders']:.1f} "
                "orders of magnitude): a few frames dominate the ensemble."
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
        print("COPER DIAGNOSTIC REPORT")
        print("=" * 60)
        print(f"\nOptimization Status: {self.message}")
        print(f"Success: {self.success}    Feasible: {self.feasible}")
        print(f"Iterations: {self.n_iterations}")
        print("\nChi-squared (max over groups):")
        print(f"  {'Initial':<{key_w}} {self.chi_squared_initial:>{val_w}.4f}")
        print(f"  {'Step-1 minimum':<{key_w}} {self.chi_squared_min:>{val_w}.4f}")
        print(f"  {'Final':<{key_w}} {self.chi_squared_final:>{val_w}.4f}")
        print(f"  {'Limit':<{key_w}} {self.chi2_limit:>{val_w}.4f}")
        print(
            f"  {'Improvement':<{key_w}} {diag['chi2_improvement']:>{val_w}.4f}"
            f" ({diag['chi2_improvement_pct']:.1f}%)"
        )
        print("\nInformation content / diversity:")
        print(f"  {'Phi (entropy fraction)':<{key_w}} {self.phi:>{val_w}.4f}")
        print(
            f"  {'N_eff (entropy-based)':<{key_w}} "
            f"{diag['neff_entropy']:>{val_w}.1f}  / {n}"
        )
        print(
            f"  {'N_eff (1/sum w^2, Renyi-2)':<{key_w}} "
            f"{diag['neff_renyi2']:>{val_w}.1f}  / {n}"
        )
        print(f"  {'delta_S':<{key_w}} {self.delta_S:>{val_w}.4f}")
        print(f"  {'<delta_G>/kT':<{key_w}} {self.mean_delta_G_kT:>{val_w}.4f}")
        if self.scale is not None:
            print(f"  {'iCOPER scale':<{key_w}} {self.scale:>{val_w}.4g}")
            print(f"  {'iCOPER offset':<{key_w}} {self.offset:>{val_w}.4g}")
        print("\nWeight distribution:")
        print(f"  {'Min':<{key_w}} {diag['weight_min']:>{val_w}.2e}")
        print(f"  {'Max':<{key_w}} {diag['weight_max']:>{val_w}.2e}")
        print(f"  {'Std Dev':<{key_w}} {diag['weight_std']:>{val_w}.2e}")
        print(
            f"  {'Range (orders of magnitude)':<{key_w}} "
            f"{diag['weight_range_orders']:>{val_w}.1f}"
        )
        print(f"  {'Reweighting factor min':<{key_w}} {diag['rfac_min']:>{val_w}.4g}")
        print(f"  {'Reweighting factor max':<{key_w}} {diag['rfac_max']:>{val_w}.4g}")

        if len(diag["warnings"]) > 0:
            print(f"\nWARNINGS ({len(diag['warnings'])}):")
            for i, warning in enumerate(diag["warnings"], 1):
                print(f"  {i}. {warning}")
        else:
            print(f"\nStatus: {diag['status']} - No issues detected")
        print("=" * 60)


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
@dataclass
class COPERScanResult:
    """Container for a chi-squared-limit scan.

    Attributes
    ----------
    chi2_limits : numpy.ndarray
        Scanned chi-squared-limit grid.
    chi_squared_values : numpy.ndarray
        Final chi-squared (``chi_squared_final``) per limit.
    phi_values : numpy.ndarray
        Effective-fraction (phi) per limit.
    kl_divergence_values : numpy.ndarray
        Relative entropy per limit.
    feasible_mask : numpy.ndarray of bool
        ``feasible`` flag per limit.
    results : list of COPERResult
        Per-limit fit results.
    optimal_chi2_limit : float
        Selected chi-squared limit (L-curve knee).
    optimal_idx : int
        Index of ``optimal_chi2_limit`` in ``chi2_limits``.
    method : str
        Human-readable knee-selection method used.
    """

    chi2_limits: np.ndarray
    chi_squared_values: np.ndarray
    phi_values: np.ndarray
    kl_divergence_values: np.ndarray
    feasible_mask: np.ndarray
    results: List[COPERResult]
    optimal_chi2_limit: float
    optimal_idx: int
    method: str

    def print_summary(self):
        """Print a formatted summary of the chi-squared-limit scan."""
        print("\n" + "=" * 60)
        print("COPER CHI^2-LIMIT SCAN SUMMARY")
        print("=" * 60)
        print(f"\nScan range: {self.chi2_limits[0]:.4f} to {self.chi2_limits[-1]:.4f}")
        print(f"Number of points: {len(self.chi2_limits)}")
        print(f"Feasible at: {int(self.feasible_mask.sum())} / {len(self.chi2_limits)}")
        print(f"Method: {self.method}\n")
        opt = self.optimal_idx
        print(f"RECOMMENDED LIMIT: {self.optimal_chi2_limit:.4f}")
        print(f"Chi squared:      {self.chi_squared_values[opt]:.4f}")
        print(f"Phi (N_eff):      {self.phi_values[opt]:.4f}")
        print(f"Relative entropy: {self.kl_divergence_values[opt]:.4f}")
        print("=" * 60)


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
def chi2_limit_scan(
    observables: List[ExperimentalObservable],
    calculated_values: np.ndarray,
    reweighter: str = "coper",
    chi2_limits: Union[Tuple[float, float], np.ndarray] = (0.25, 4.0),
    n_points: int = 8,
    log_scale: bool = True,
    initial_weights: Optional[np.ndarray] = None,
    method: str = "perpendicular",
    verbose: bool = False,
    fit_kwargs: Optional[dict] = None,
) -> COPERScanResult:
    """Scan a grid of chi-squared limits for COPER or iCOPER.

    The L-curve analogue of :func:`soursop.ssbme.theta_scan`: tightening
    the chi-squared limit (the paper's error-scaling range of roughly
    0.25-4) trades fit quality for ensemble diversity. The knee of the
    chi-squared / relative-entropy curve gives a principled choice.

    Parameters
    ----------
    observables : list of ExperimentalObservable
        Experimental observables.
    calculated_values : numpy.ndarray
        Per-frame calculated values, shape ``(n_frames, n_observables)``.
    reweighter : str, optional
        ``"coper"`` (default) or ``"icoper"``.
    chi2_limits : tuple or numpy.ndarray, optional
        ``(min, max)`` for a generated grid of ``n_points``, or an explicit
        1D array of limits. Default ``(0.25, 4.0)`` (the paper's error-
        scaling range).
    n_points : int, optional
        Number of grid points when ``chi2_limits`` is a tuple. Default 8.
    log_scale : bool, optional
        Sample logarithmically when ``chi2_limits`` is a tuple. Default True.
    initial_weights : numpy.ndarray, optional
        Prior weights. Uniform if None.
    method : str, optional
        Knee-selection method (see
        :func:`soursop.ssutils.find_optimal_theta`).
    verbose : bool, optional
        Print per-limit progress. Default False.
    fit_kwargs : dict, optional
        Extra keyword arguments forwarded to the reweighter's ``fit``
        (e.g. ``ftol`` / ``max_icoper_iterations`` for iCOPER).

    Returns
    -------
    COPERScanResult

    Raises
    ------
    SSException
        If ``reweighter`` is unknown, ``n_points`` < 1, or a log-scale grid
        has non-positive endpoints.
    """
    if reweighter not in ("coper", "icoper"):
        raise SSException(
            f"Unknown reweighter: {reweighter}, must be 'coper' or 'icoper'"
        )

    if isinstance(chi2_limits, (tuple, list)):
        if n_points < 1:
            raise SSException(
                f"n_points must be >= 1 for a tuple chi2_limits, got {n_points}"
            )
        if log_scale and (chi2_limits[0] <= 0 or chi2_limits[1] <= 0):
            raise SSException(
                "chi2_limits endpoints must be positive when log_scale=True, "
                f"got {chi2_limits}"
            )
        if log_scale:
            limits = np.logspace(
                np.log10(chi2_limits[0]), np.log10(chi2_limits[1]), n_points
            )
        else:
            limits = np.linspace(chi2_limits[0], chi2_limits[1], n_points)
    else:
        limits = np.asarray(chi2_limits, dtype=np.float64)

    fit_kwargs = dict(fit_kwargs) if fit_kwargs else {}

    chi2_vals: List[float] = []
    phi_vals: List[float] = []
    kl_vals: List[float] = []
    feas: List[bool] = []
    results: List[COPERResult] = []

    for i, lim in enumerate(limits):
        if verbose:
            print(f"Processing chi2_limit {i + 1}/{len(limits)}: {lim:.4f}")
        if reweighter == "coper":
            rw = COPER(observables, calculated_values, initial_weights)
        else:
            rw = iCOPER(observables, calculated_values, initial_weights)
        res = rw.fit(chi2_limit=float(lim), verbose=False, **fit_kwargs)
        results.append(res)
        chi2_vals.append(res.chi_squared_final)
        phi_vals.append(res.phi)
        kl_vals.append(res.kl_divergence)
        feas.append(bool(res.feasible))

    chi2_arr = np.array(chi2_vals)
    phi_arr = np.array(phi_vals)
    kl_arr = np.array(kl_vals)
    feas_arr = np.array(feas, dtype=bool)

    optimal_idx, method_name = find_optimal_theta(chi2_arr, kl_arr, method=method)

    return COPERScanResult(
        chi2_limits=limits,
        chi_squared_values=chi2_arr,
        phi_values=phi_arr,
        kl_divergence_values=kl_arr,
        feasible_mask=feas_arr,
        results=results,
        optimal_chi2_limit=float(limits[optimal_idx]),
        optimal_idx=optimal_idx,
        method=method_name,
    )


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
class COPER:
    """Convex Optimization for Ensemble Reweighting (COPER).

    Reweights an ensemble by *maximising entropy subject to a hard
    chi-squared constraint*, solved as a primal convex optimisation over
    the N frame weights via SciPy's ``trust-constr`` interior-point
    method. The two-step procedure of Leung et al. (2016):

    1. minimise the per-group chi-squared over the simplex to test
       feasibility, then
    2. maximise the (relative) entropy subject to
       ``chi^2_alpha <= chi2_limit`` for every observable group, starting
       from the step-1 point.

    Parameters
    ----------
    observables : list of ExperimentalObservable
        Experimental observables to fit against. Observables sharing the
        same :attr:`~soursop.ssutils.ExperimentalObservable.group` label
        contribute to a single chi-squared term ``chi^2_alpha``;
        observables without a group are pooled into one default group.
    calculated_values : numpy.ndarray
        Per-frame calculated observable values, shape
        ``(n_frames, n_observables)``.
    initial_weights : numpy.ndarray, optional
        Prior frame weights. If None, uniform weights are used. Internally
        normalised to sum to 1.

    Raises
    ------
    SSException
        If inputs are malformed
        (see :func:`soursop.ssutils.validate_reweighting_inputs`).

    Notes
    -----
    The simplex (``w_i >= 0``, ``sum_i w_i = 1``) is enforced implicitly via
    a softmax reparameterisation (``w = softmax(z)``), so the optimiser only
    has to handle the few chi-squared constraints rather than ``N`` bound
    constraints plus a sum-to-one equality. Both the entropy objective and the
    chi-squared constraints supply **exact analytic Hessians as matrix-free**
    ``LinearOperator`` objects (the entropy Hessian is effectively diagonal;
    each chi-squared Hessian is rank ``<= m`` via the calculated-value matrix), so
    no dense ``N x N`` matrix is ever formed and memory stays ``O(N)`` rather
    than ``O(N^2)``. ``trust-constr`` then converges in tens of iterations and
    scales to ~10^5 frames (a binding 10^5-conformer fit takes seconds-to-tens-
    of-seconds; the equivalent dense Hessian would need tens of GB). The
    problem is convex with a unique global solution, so neither the
    reparameterisation nor the exact Hessians change the optimum.
    """

    def __init__(
        self,
        observables: List[ExperimentalObservable],
        calculated_values: np.ndarray,
        initial_weights: Optional[np.ndarray] = None,
    ):
        _validate_inputs(observables, calculated_values, initial_weights)

        self.observables = list(observables)
        self.calculated_values = np.asarray(calculated_values, dtype=np.float64)
        self.n_frames = self.calculated_values.shape[0]
        self.n_observables = len(observables)

        if initial_weights is None:
            self.initial_weights = np.ones(self.n_frames) / self.n_frames
        else:
            self.initial_weights = np.asarray(
                initial_weights, dtype=np.float64
            ) / float(np.sum(initial_weights))

        # Build {group_name: [observable indices]}; ungrouped observables
        # are pooled into the default group "all".
        groups: dict = {}
        for i, obs in enumerate(self.observables):
            key = obs.group if obs.group is not None else "all"
            groups.setdefault(key, []).append(i)
        self._group_names = list(groups.keys())
        self._group_indices = [np.asarray(v, dtype=int) for v in groups.values()]

        self._result: Optional[COPERResult] = None
        self._chi2_limit: Optional[float] = None
        self._scan_result: Optional[COPERScanResult] = None

    # ------------------------------------------------------------------
    def _chi2_per_group_vec(self, weights: np.ndarray) -> np.ndarray:
        """Per-group constraint-aware reduced chi-squared at ``weights``."""
        return np.array(
            [
                constraint_chi_squared(
                    weights, self.calculated_values, self.observables, indices
                )
                for indices in self._group_indices
            ],
            dtype=np.float64,
        )

    def _chi2_per_group_jac(self, weights: np.ndarray) -> np.ndarray:
        """Jacobian of the per-group chi-squared vector, shape ``(G, N)``.

        Row ``alpha`` is the analytic gradient of ``chi^2_alpha(w)``. For an
        ``upper`` / ``lower`` constraint the gradient is zero whenever the
        constraint is not active (i.e. on the unpenalised side).
        """
        jac = np.zeros((len(self._group_indices), self.n_frames))
        for a, indices in enumerate(self._group_indices):
            m = len(indices)
            for idx in indices:
                obs = self.observables[idx]
                calc_col = self.calculated_values[:, idx]
                calc_avg = float(np.sum(calc_col * weights))
                diff = calc_avg - obs.value
                if obs.constraint == "equality":
                    penalize = True
                elif obs.constraint == "upper":
                    penalize = diff > 0
                else:  # "lower"
                    penalize = diff < 0
                if penalize:
                    jac[a] += (2.0 / m) * (diff / obs.uncertainty**2) * calc_col
        return jac

    def _chi2_sum_obj(self, weights: np.ndarray) -> float:
        """Sum of per-group chi-squared (objective of the step-1 problem)."""
        return float(np.sum(self._chi2_per_group_vec(weights)))

    def _chi2_sum_grad(self, weights: np.ndarray) -> np.ndarray:
        """Analytic gradient of :meth:`_chi2_sum_obj`."""
        return self._chi2_per_group_jac(weights).sum(axis=0)

    def _entropy_obj(self, weights: np.ndarray) -> float:
        """Relative entropy ``sum w * ln(w/w0)`` (minimise = maximise -KL)."""
        safe = np.maximum(weights, MIN_WEIGHT_THRESHOLD)
        return float(np.sum(safe * np.log(safe / self.initial_weights)))

    def _entropy_grad(self, weights: np.ndarray) -> np.ndarray:
        """Analytic gradient of :meth:`_entropy_obj`."""
        safe = np.maximum(weights, MIN_WEIGHT_THRESHOLD)
        return np.log(safe / self.initial_weights) + 1.0

    # ------------------------------------------------------------------
    # Softmax reparameterisation helpers.
    #
    # The simplex (w_i >= 0, sum_i w_i = 1) is enforced implicitly by
    # optimising over an unconstrained vector ``z`` with ``w = softmax(z)``.
    # This removes the N bound constraints and the sum-to-one equality from
    # the problem handed to the optimiser, leaving only the (few) chi-squared
    # constraints. trust-constr then converges in tens of iterations instead
    # of the hundreds-to-thousands needed when it must also juggle N active
    # bound constraints - a several-orders-of-magnitude speed-up for typical
    # ensembles, with an identical optimum (same convex problem).
    @staticmethod
    def _softmax(z: np.ndarray) -> np.ndarray:
        """Numerically stable softmax mapping ``z`` -> simplex weights."""
        z = z - np.max(z)
        ez = np.exp(z)
        return ez / np.sum(ez)

    @staticmethod
    def _grad_w_to_z(w: np.ndarray, grad_w: np.ndarray) -> np.ndarray:
        """Map a w-space gradient of a scalar ``g`` to z-space (``w = softmax(z)``).

        ``dg/dz_j = w_j (dg/dw_j - sum_k w_k dg/dw_k)``. Any additive constant
        in ``grad_w`` cancels, so the w-space gradients of the entropy and
        chi-squared objectives can be reused verbatim.
        """
        return w * (grad_w - np.dot(grad_w, w))

    def _chi2_per_group_jac_z(self, w: np.ndarray) -> np.ndarray:
        """Per-group chi-squared Jacobian in z-space, shape ``(G, N)``."""
        jac_w = self._chi2_per_group_jac(w)
        dots = jac_w @ w
        return w[None, :] * (jac_w - dots[:, None])

    # ------------------------------------------------------------------
    # Exact, matrix-free Hessians in z-space.
    #
    # Both the entropy objective and each chi-squared constraint have cheap
    # analytic Hessians, supplied to trust-constr as ``LinearOperator``s so the
    # dense N x N matrix is never formed. This keeps memory at O(N) (rather than
    # the O(N^2) of a BFGS approximation) and lets the solver use its iterative
    # trust-region subproblem solver, so COPER scales to very large ensembles.
    #
    # For g(w) with w = softmax(z), the Hessian-vector product is
    #   H_z v = u o (gw - a) + w o (Hw u) - (b + c) w,
    # where u = w o v - w (w . v) is the softmax-Jacobian applied to v,
    # gw = grad_w g, Hw = hess_w g, a = w . gw, b = u . gw, c = w . (Hw u).
    # (Verified against finite differences.)
    @staticmethod
    def _hvp_z(w, v, grad_w, hw_matvec, a=None):
        """Hessian-vector product in z-space for one scalar function of w."""
        if a is None:
            a = float(np.dot(w, grad_w))
        u = w * v - w * float(np.dot(w, v))
        hu = hw_matvec(u)
        b = float(np.dot(u, grad_w))
        c = float(np.dot(w, hu))
        return u * (grad_w - a) + w * hu - (b + c) * w

    def _entropy_hess_z(self, z: np.ndarray) -> LinearOperator:
        """LinearOperator for the z-space Hessian of the entropy objective."""
        w = self._softmax(z)
        safe = np.maximum(w, MIN_WEIGHT_THRESHOLD)
        grad_w = np.log(safe / self.initial_weights) + 1.0
        a = float(np.dot(w, grad_w))

        # H_w = diag(1/w); H_w u = u / w.
        def hw(u):
            return u / safe

        n = w.size

        def mv(v):
            return self._hvp_z(w, np.asarray(v).ravel(), grad_w, hw, a)

        return LinearOperator((n, n), matvec=mv, rmatvec=mv, dtype=np.float64)

    def _group_chi2_grad_and_hw(self, w: np.ndarray):
        """Per-group ``(grad_w, Hw-matvec)`` for the chi-squared constraints.

        Only the *active* observables (those actually penalised at ``w`` - all
        of them for ``equality``, the violated side for ``upper`` / ``lower``)
        contribute, matching :meth:`_chi2_per_group_jac`. The Hessian
        ``H_w = (2/m) sum_k active  calc_k outer calc_k / sigma_k^2`` is applied
        as a rank-(<=k) matrix-free product, ``O(N * k)`` per call.
        """
        per_group = []
        for indices in self._group_indices:
            m = len(indices)
            cols, inv_sig2, diffs = [], [], []
            for idx in indices:
                obs = self.observables[idx]
                calc_col = self.calculated_values[:, idx]
                diff = float(np.sum(calc_col * w)) - obs.value
                if obs.constraint == "equality":
                    penalize = True
                elif obs.constraint == "upper":
                    penalize = diff > 0
                else:  # "lower"
                    penalize = diff < 0
                if penalize:
                    cols.append(calc_col)
                    inv_sig2.append(1.0 / obs.uncertainty**2)
                    diffs.append(diff)
            if cols:
                A = np.column_stack(cols)  # (N, k) active columns
                inv_sig2 = np.asarray(inv_sig2)
                grad_w = (2.0 / m) * (A @ (np.asarray(diffs) * inv_sig2))

                def hw(u, A=A, inv_sig2=inv_sig2, m=m):
                    return (2.0 / m) * (A @ ((A.T @ u) * inv_sig2))
            else:
                grad_w = np.zeros(self.n_frames)

                def hw(u):
                    return np.zeros_like(u)

            per_group.append((grad_w, hw))
        return per_group

    def _constraint_hess_z(self, z: np.ndarray, lagrange: np.ndarray) -> LinearOperator:
        """LinearOperator for the constraint Hessian ``sum_a v_a H_z[chi2_a]``."""
        w = self._softmax(z)
        per_group = self._group_chi2_grad_and_hw(w)
        terms = [
            (float(lagrange[a]), grad_w, float(np.dot(w, grad_w)), hw)
            for a, (grad_w, hw) in enumerate(per_group)
            if lagrange[a] != 0.0
        ]
        n = w.size

        def mv(v):
            v = np.asarray(v).ravel()
            out = np.zeros(n)
            for coef, grad_w, a, hw in terms:
                out += coef * self._hvp_z(w, v, grad_w, hw, a)
            return out

        return LinearOperator((n, n), matvec=mv, rmatvec=mv, dtype=np.float64)

    # ------------------------------------------------------------------
    def _make_result(
        self,
        *,
        weights: np.ndarray,
        chi_squared_initial: float,
        chi_squared_min: float,
        chi_squared_final: float,
        chi2_limit: float,
        feasible: bool,
        success: bool,
        message: str,
        n_iterations: int,
        metadata: dict,
    ) -> COPERResult:
        kl = _srel(self.initial_weights, weights)
        entropy = _shannon_entropy(weights)
        entropy_initial = _shannon_entropy(self.initial_weights)
        delta_S = entropy - entropy_initial
        # <delta_G>/kT = KL(w || w0); for a uniform w0 this equals -delta_S
        # (matches Leung et al. eq 10).
        mean_delta_G_kT = float(kl)
        phi = float(np.exp(-kl))
        return COPERResult(
            weights=weights,
            initial_weights=self.initial_weights.copy(),
            chi_squared_initial=float(chi_squared_initial),
            chi_squared_min=float(chi_squared_min),
            chi_squared_final=float(chi_squared_final),
            chi2_limit=float(chi2_limit),
            feasible=bool(feasible),
            entropy=float(entropy),
            entropy_initial=float(entropy_initial),
            delta_S=float(delta_S),
            mean_delta_G_kT=mean_delta_G_kT,
            phi=phi,
            n_iterations=int(n_iterations),
            success=bool(success),
            message=str(message),
            observables=self.observables,
            calculated_values=self.calculated_values,
            metadata=metadata,
        )

    # ------------------------------------------------------------------
    def fit(
        self,
        chi2_limit: float = DEFAULT_CHI2_LIMIT,
        max_iterations: int = DEFAULT_MAX_ITERATIONS,
        optimizer: str = DEFAULT_OPTIMIZER,
        verbose: bool = True,
    ) -> COPERResult:
        """Fit COPER weights at a given chi-squared limit.

        Parameters
        ----------
        chi2_limit : float, optional
            Hard upper bound on each per-group chi-squared. Default ``1.0``
            (the paper's ``chi^2 <= 1`` convention).
        max_iterations : int, optional
            Maximum optimizer iterations (each of the two steps). Default
            from the module.
        optimizer : str, optional
            scipy.optimize.minimize method used for the constrained
            entropy-maximisation step (must support nonlinear constraints).
            Default ``"trust-constr"``. The feasibility step is always solved
            with the unconstrained ``"L-BFGS-B"``.
        verbose : bool, optional
            Print progress. Default True.

        Returns
        -------
        COPERResult

        Raises
        ------
        SSException
            If ``chi2_limit`` is not positive.

        Notes
        -----
        The result is also stored on the instance (``self.result``,
        ``self.chi2_limit``). When the problem is infeasible the returned
        result has ``feasible=False`` and ``success=False``; its
        ``weights`` are the step-1 chi-squared minimiser (the closest the
        prior ensemble can come to satisfying the data).
        """
        if chi2_limit <= 0:
            raise SSException(f"chi2_limit must be positive, got {chi2_limit}")
        self._chi2_limit = float(chi2_limit)
        n = self.n_frames

        # Optimise over z with w = softmax(z); z0 = log(w0) gives softmax(z0) = w0.
        w0 = self.initial_weights.copy()
        z0 = np.log(np.maximum(w0, MIN_WEIGHT_THRESHOLD))
        chi2_initial = float(np.max(self._chi2_per_group_vec(w0)))

        metadata = {
            "optimizer": optimizer,
            "max_iterations": max_iterations,
            "groups": dict(
                zip(
                    self._group_names,
                    [list(map(int, g)) for g in self._group_indices],
                )
            ),
            "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }

        if verbose:
            print("COPER Optimization")
            print(
                f"  Frames: {n}, Observables: {self.n_observables}, "
                f"Groups: {len(self._group_indices)}"
            )
            print(f"  Chi-squared limit:   {chi2_limit}")
            print(f"  Chi-squared initial (max over groups): {chi2_initial:.4f}")

        # ---------- Step 1: chi-squared minimisation / feasibility ----------
        # Unconstrained in z (the softmax handles the simplex), so a plain
        # quasi-Newton solve suffices and is fast.
        opt1 = minimize(
            fun=lambda z: self._chi2_sum_obj(self._softmax(z)),
            x0=z0,
            method="L-BFGS-B",
            jac=lambda z: self._grad_w_to_z(
                self._softmax(z), self._chi2_sum_grad(self._softmax(z))
            ),
            options={"maxiter": max_iterations},
        )
        w_min = self._softmax(opt1.x)
        chi2_min = float(np.max(self._chi2_per_group_vec(w_min)))
        feasible = chi2_min <= chi2_limit + FEASIBILITY_TOLERANCE

        if verbose:
            print(f"  Step 1 chi^2 (max over groups): {chi2_min:.4f}")
            print(f"  Feasible (<= limit + tol): {feasible}")

        if not feasible:
            self._result = self._make_result(
                weights=w_min,
                chi_squared_initial=chi2_initial,
                chi_squared_min=chi2_min,
                chi_squared_final=chi2_min,
                chi2_limit=chi2_limit,
                feasible=False,
                success=False,
                message=(
                    f"Data infeasible at chi2_limit={chi2_limit:.4g}: "
                    f"minimum max-group chi^2 ({chi2_min:.4g}) "
                    "exceeds the limit."
                ),
                n_iterations=int(getattr(opt1, "nit", -1)),
                metadata=metadata,
            )
            return self._result

        # ---------- Step 2: entropy maximisation (KL minimisation) ----------
        # Only the per-group chi-squared upper bounds remain as explicit
        # constraints (the simplex is implicit in the softmax), so trust-constr
        # converges quickly.
        nl_constraint = NonlinearConstraint(
            fun=lambda z: self._chi2_per_group_vec(self._softmax(z)),
            lb=-np.inf,
            ub=chi2_limit,
            jac=lambda z: self._chi2_per_group_jac_z(self._softmax(z)),
            hess=self._constraint_hess_z,
        )
        # Start the entropy maximisation from the uniform prior (z0), the
        # unconstrained max-entropy point, rather than from the step-1
        # chi-squared minimum (a heavily-reweighted, low-entropy point). The
        # problem is convex so the optimum is unique, but starting at the
        # max-entropy end and tightening toward the constraint is far better
        # conditioned and converges reliably to the true (highest-entropy)
        # solution.
        opt2 = minimize(
            fun=lambda z: self._entropy_obj(self._softmax(z)),
            x0=z0,
            method=optimizer,
            jac=lambda z: self._grad_w_to_z(
                self._softmax(z), self._entropy_grad(self._softmax(z))
            ),
            hess=self._entropy_hess_z,
            constraints=[nl_constraint],
            options={"maxiter": max_iterations, **_TRUST_CONSTR_OPTIONS},
        )
        w_opt = self._softmax(opt2.x)
        chi2_final = float(np.max(self._chi2_per_group_vec(w_opt)))

        if verbose:
            print(f"  Step 2 chi^2 (max over groups): {chi2_final:.4f}")
            print(f"  Optimization successful: {bool(opt2.success)}")

        self._result = self._make_result(
            weights=w_opt,
            chi_squared_initial=chi2_initial,
            chi_squared_min=chi2_min,
            chi_squared_final=chi2_final,
            chi2_limit=chi2_limit,
            feasible=True,
            success=bool(opt2.success),
            message=str(opt2.message),
            n_iterations=int(opt2.nit),
            metadata=metadata,
        )
        return self._result

    # ------------------------------------------------------------------
    def scan_chi2_limit(
        self,
        chi2_limits: Union[Tuple[float, float], np.ndarray] = (0.25, 4.0),
        n_points: int = 8,
        log_scale: bool = True,
        method: str = "perpendicular",
        verbose: bool = False,
    ) -> COPERScanResult:
        """Run a chi-squared-limit scan using this instance's data."""
        scan = chi2_limit_scan(
            observables=self.observables,
            calculated_values=self.calculated_values,
            reweighter="coper",
            chi2_limits=chi2_limits,
            n_points=n_points,
            log_scale=log_scale,
            initial_weights=self.initial_weights,
            method=method,
            verbose=verbose,
        )
        self._scan_result = scan
        return scan

    def predict(self, calculated_values: np.ndarray) -> np.ndarray:
        """Weighted means of arbitrary observables using the fitted weights.

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
    def result(self) -> Optional[COPERResult]:
        """The most recent :class:`COPERResult`, or None."""
        return self._result

    @property
    def chi2_limit(self) -> Optional[float]:
        """The chi-squared limit used in the most recent fit, or None."""
        return self._chi2_limit

    @property
    def scan_result(self) -> Optional[COPERScanResult]:
        """The most recent :class:`COPERScanResult`, or None."""
        return self._scan_result


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
class iCOPER:
    """Iterative COPER for data with an unknown global scale and offset.

    At each iteration a (optionally inverse-variance-weighted) linear
    regression of the ensemble-averaged calculated values against the
    experimental values is used to refit a global scale ``alpha`` and
    offset ``beta``; the calculated values are rescaled
    (``calc -> alpha * calc + beta``) and a standard :class:`COPER` step
    is run. The procedure repeats until the change in the final
    chi-squared drops below ``ftol``. This is the COPER analogue of
    :class:`soursop.ssbme.iBME`.

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
        If inputs are malformed.
    """

    def __init__(
        self,
        observables: List[ExperimentalObservable],
        calculated_values: np.ndarray,
        initial_weights: Optional[np.ndarray] = None,
    ):
        _validate_inputs(observables, calculated_values, initial_weights)

        self.observables = list(observables)
        self.calculated_values = np.asarray(calculated_values, dtype=np.float64)
        self.n_frames = self.calculated_values.shape[0]
        self.n_observables = len(observables)

        if initial_weights is None:
            self.initial_weights = np.ones(self.n_frames) / self.n_frames
        else:
            self.initial_weights = np.asarray(
                initial_weights, dtype=np.float64
            ) / float(np.sum(initial_weights))

        self._exp_values = np.array(
            [obs.value for obs in self.observables], dtype=np.float64
        )
        self._exp_sigma = np.array(
            [obs.uncertainty for obs in self.observables], dtype=np.float64
        )

        self._result: Optional[COPERResult] = None
        self._chi2_limit: Optional[float] = None
        self._scan_result: Optional[COPERScanResult] = None

    # ------------------------------------------------------------------
    def fit(
        self,
        chi2_limit: float = DEFAULT_CHI2_LIMIT,
        ftol: float = 0.01,
        max_icoper_iterations: int = 50,
        fit_offset: bool = True,
        lr_weights: bool = True,
        max_iterations: int = DEFAULT_MAX_ITERATIONS,
        optimizer: str = DEFAULT_OPTIMIZER,
        verbose: bool = True,
    ) -> COPERResult:
        """Run iterative COPER at a fixed chi-squared limit.

        Parameters
        ----------
        chi2_limit : float, optional
            Hard upper bound on each per-group chi-squared (passed to the
            inner :class:`COPER` step). Default 1.0.
        ftol : float, optional
            Convergence tolerance on ``|delta chi-squared|`` between
            iCOPER iterations. Default 0.01.
        max_icoper_iterations : int, optional
            Maximum number of scale/offset + COPER iterations. Default 50.
        fit_offset : bool, optional
            If True fit a scale and offset; if False fit scale only.
            Default True.
        lr_weights : bool, optional
            If True weight the linear regression by ``1/sigma^2``; otherwise
            use uniform weights. Default True.
        max_iterations : int, optional
            Maximum inner-optimizer iterations per COPER step.
        optimizer : str, optional
            scipy.optimize.minimize method for each inner COPER step.
        verbose : bool, optional
            Print per-iteration progress. Default True.

        Returns
        -------
        COPERResult
            Result with the final weights, ``scale``/``offset`` (net,
            relative to the original input), the ``icoper_iterations``
            log, ``chi_squared_initial`` from the first iteration's
            pre-fit chi-squared, and ``phi`` computed against the original
            prior.

        Raises
        ------
        SSException
            If ``chi2_limit`` is not positive.
        """
        if chi2_limit <= 0:
            raise SSException(f"chi2_limit must be positive, got {chi2_limit}")
        self._chi2_limit = float(chi2_limit)

        w0 = self.initial_weights.copy()
        current_weights = w0.copy()
        calc = self.calculated_values.copy()
        if lr_weights:
            lr_w = 1.0 / self._exp_sigma**2
        else:
            lr_w = np.ones(self.n_observables)

        net_scale = 1.0
        net_offset = 0.0

        iterations: List[dict] = []
        chi2_initial = float("nan")
        chi2_old = float("nan")
        last_result: Optional[COPERResult] = None
        it = 0

        for it in range(max_icoper_iterations):
            calc_avg = np.sum(calc * current_weights[:, np.newaxis], axis=0)
            alpha, beta = _weighted_linear_regression(
                calc_avg, self._exp_values, lr_w, fit_intercept=fit_offset
            )
            calc = alpha * calc + beta
            net_scale *= alpha
            net_offset = alpha * net_offset + beta

            coper = COPER(self.observables, calc, w0)
            res = coper.fit(
                chi2_limit=self._chi2_limit,
                max_iterations=max_iterations,
                optimizer=optimizer,
                verbose=False,
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
                    "feasible": res.feasible,
                    "diff": diff,
                }
            )
            if verbose:
                print(
                    f"Iteration {it:3d}  scale={alpha:8.4f}  "
                    f"offset={beta:8.4f}  chi2={chi2_new:8.4f}  "
                    f"feasible={res.feasible}  diff={diff:8.2e}"
                )

            if np.isfinite(diff) and diff < ftol:
                if verbose:
                    print(
                        f"iCOPER converged below tolerance {ftol:.2e} after "
                        f"{it + 1} iterations"
                    )
                break

        if last_result is None:
            raise SSException("iCOPER: no iterations were run")

        metadata = {
            "optimizer": optimizer,
            "max_iterations": max_iterations,
            "max_icoper_iterations": max_icoper_iterations,
            "ftol": ftol,
            "fit_offset": fit_offset,
            "lr_weights": lr_weights,
            "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }

        # Recompute entropy / phi against the original prior w0.
        kl = _srel(w0, current_weights)
        entropy = _shannon_entropy(current_weights)
        entropy_initial = _shannon_entropy(w0)
        delta_S = entropy - entropy_initial
        phi = float(np.exp(-kl))

        self._result = COPERResult(
            weights=current_weights,
            initial_weights=w0,
            chi_squared_initial=float(chi2_initial),
            chi_squared_min=float(last_result.chi_squared_min),
            chi_squared_final=float(chi2_old),
            chi2_limit=float(chi2_limit),
            feasible=bool(last_result.feasible),
            entropy=float(entropy),
            entropy_initial=float(entropy_initial),
            delta_S=float(delta_S),
            mean_delta_G_kT=float(kl),
            phi=phi,
            n_iterations=int(it + 1),
            success=bool(last_result.success),
            message=str(last_result.message),
            observables=self.observables,
            calculated_values=calc,
            metadata=metadata,
            scale=float(net_scale),
            offset=float(net_offset),
            icoper_iterations=iterations,
        )
        return self._result

    # ------------------------------------------------------------------
    def scan_chi2_limit(
        self,
        chi2_limits: Union[Tuple[float, float], np.ndarray] = (0.25, 4.0),
        n_points: int = 8,
        log_scale: bool = True,
        method: str = "perpendicular",
        verbose: bool = False,
        fit_kwargs: Optional[dict] = None,
    ) -> COPERScanResult:
        """Run a chi-squared-limit scan using iCOPER."""
        scan = chi2_limit_scan(
            observables=self.observables,
            calculated_values=self.calculated_values,
            reweighter="icoper",
            chi2_limits=chi2_limits,
            n_points=n_points,
            log_scale=log_scale,
            initial_weights=self.initial_weights,
            method=method,
            verbose=verbose,
            fit_kwargs=fit_kwargs,
        )
        self._scan_result = scan
        return scan

    def predict(self, calculated_values: np.ndarray) -> np.ndarray:
        """Weighted means of arbitrary observables using the fitted weights.

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

    def get_icoper_weights(self) -> np.ndarray:
        """Final iCOPER weights.

        Raises
        ------
        SSException
            If :meth:`fit` has not been called.
        """
        if self._result is None:
            raise SSException("iCOPER weights not available. Call fit() first.")
        return self._result.weights

    def get_icoper_stats(self) -> List[dict]:
        """Per-iteration iCOPER log
        (``iteration``, ``scale``, ``offset``, ``chi_squared``,
        ``feasible``, ``diff``).

        Raises
        ------
        SSException
            If :meth:`fit` has not been called.
        """
        if self._result is None:
            raise SSException("iCOPER stats not available. Call fit() first.")
        return self._result.icoper_iterations

    @property
    def result(self) -> Optional[COPERResult]:
        """The most recent :class:`COPERResult`, or None."""
        return self._result

    @property
    def chi2_limit(self) -> Optional[float]:
        """The chi-squared limit used in the most recent fit, or None."""
        return self._chi2_limit

    @property
    def scan_result(self) -> Optional[COPERScanResult]:
        """The most recent :class:`COPERScanResult`, or None."""
        return self._scan_result
