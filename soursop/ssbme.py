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
from .ssutils import (
    MIN_WEIGHT_THRESHOLD,
    VALID_CONSTRAINTS,
    ExperimentalObservable,
    constraint_chi_squared,
    find_optimal_theta,
    relative_entropy,
    validate_reweighting_inputs,
    weighted_linear_regression,
)

# These reweighting primitives were factored out to soursop.ssutils so
# they can be shared with soursop.sscoper (COPER / iCOPER). The aliases
# below preserve this module's public API
# (``ExperimentalObservable``, ``find_optimal_theta``,
# ``MIN_WEIGHT_THRESHOLD``, ``VALID_CONSTRAINTS``) and its internal call
# sites (``_srel``, ``_weighted_linear_regression``, ``_validate_inputs``)
# unchanged.
_srel = relative_entropy
_weighted_linear_regression = weighted_linear_regression
_validate_inputs = validate_reweighting_inputs

__all__ = [
    "BME",
    "iBME",
    "BMECustom",
    "BMEResult",
    "BMECustomResult",
    "ThetaScanResult",
    "theta_scan",
    "ExperimentalObservable",
    "find_optimal_theta",
    "MIN_WEIGHT_THRESHOLD",
    "VALID_CONSTRAINTS",
]

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Module constants
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
DEFAULT_THETA = 0.5
DEFAULT_MAX_ITERATIONS = 50000
DEFAULT_OPTIMIZER = "L-BFGS-B"
LAMBDA_INIT_SCALE = 1e-3


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
        return constraint_chi_squared(weights, self.calculated_values, self.observables)

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


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
@dataclass
class BMECustomResult:
    """Container for a :class:`BMECustom` reweighting result.

    Mirrors :class:`BMEResult` but for the raw vector/matrix + custom-cost
    variant: the scalar fit metric is a general ``cost`` (the reduced
    chi-squared by default, or whatever the user's ``cost_function``
    returns) rather than a constraint-aware chi-squared, and the
    experimental data are held as plain arrays rather than a list of
    :class:`~soursop.ssutils.ExperimentalObservable`.

    Attributes
    ----------
    weights : numpy.ndarray
        Optimized (posterior) frame weights, summing to 1.
    initial_weights : numpy.ndarray
        Prior frame weights.
    cost_initial : float
        Cost (goodness-of-fit) at the prior weights.
    cost_final : float
        Cost at the optimized weights.
    phi : float
        Fraction of effective frames, ``exp(-D_KL(w || w0)) in (0, 1]``.
    n_iterations : int
        Optimizer iterations.
    success : bool
        Whether the optimization succeeded.
    message : str
        Optimizer status message.
    theta : float
        Entropy-penalty strength used.
    experiment : numpy.ndarray
        Experimental vector, shape ``(m,)``.
    calculated_values : numpy.ndarray
        Per-conformer calculated matrix used in the fit, shape ``(n, m)``.
    metadata : dict
        Free-form metadata (optimizer, whether a custom cost was used, ...).
    """

    weights: np.ndarray
    initial_weights: np.ndarray
    cost_initial: float
    cost_final: float
    phi: float
    n_iterations: int
    success: bool
    message: str
    theta: float
    experiment: np.ndarray
    calculated_values: np.ndarray
    metadata: dict

    def __str__(self):
        status = "SUCCESS" if self.success else "FAILED"
        return "\n".join(
            [
                f"BMECustom Result [{status}]",
                f"  Cost initial: {self.cost_initial:.4f}",
                f"  Cost final:   {self.cost_final:.4f}",
                f"  phi (effective fraction): {self.phi:.4f}",
                f"  Iterations: {self.n_iterations}",
                f"  Theta: {self.theta}",
            ]
        )

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
        """Per-frame reweighting factors ``r_i = w_i / w0_i``."""
        return self.weights / self.initial_weights

    def predict(self, calculated_values: np.ndarray) -> np.ndarray:
        """Weighted average of arbitrary observables using these weights.

        Parameters
        ----------
        calculated_values : numpy.ndarray
            Per-frame values, shape ``(n_frames, ...)``; ``n_frames`` must
            match the fitted ensemble.

        Returns
        -------
        numpy.ndarray
            Weighted average over frames (axis 0).

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
        weights = self.weights.reshape((-1,) + (1,) * (calculated_values.ndim - 1))
        return np.sum(weights * calculated_values, axis=0)

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
        """
        diag: dict = {}
        warnings: List[str] = []
        n = len(self.weights)

        diag["neff_entropy"] = n * float(self.phi)
        diag["neff_entropy_fraction"] = float(self.phi)
        diag["neff_renyi2"] = float(1.0 / np.sum(self.weights**2))
        diag["neff_renyi2_fraction"] = float(diag["neff_renyi2"] / n)
        diag["cost_improvement"] = self.cost_initial - self.cost_final

        if self.phi < warn_threshold:
            warnings.append(
                f"Low Phi ({self.phi:.3f} < {warn_threshold}): significant "
                "loss of ensemble diversity. Consider increasing theta or "
                "loosening the experimental uncertainties."
            )
        if diag["neff_renyi2"] < 0.1 * n:
            warnings.append(
                f"Low effective sample size (1/sum w^2) "
                f"({diag['neff_renyi2']:.1f} / {n}): only "
                f"~{diag['neff_renyi2_fraction'] * 100:.1f}% of frames "
                "are effectively used."
            )

        diag["warnings"] = warnings
        diag["status"] = "OK" if len(warnings) == 0 else "WARNING"
        return diag

    def print_diagnostics(self, warn_threshold: float = 0.5):
        """Print a formatted diagnostic report for this result."""
        diag = self.diagnostics(warn_threshold)
        n = len(self.weights)
        key_w, val_w = 32, 14
        print("\n" + "=" * 60)
        print("BME (custom) DIAGNOSTIC REPORT")
        print("=" * 60)
        print(f"\nOptimization Status: {self.message}")
        print(f"Success: {self.success}")
        print(f"Iterations: {self.n_iterations}")
        print("\nCost:")
        print(f"  {'Initial':<{key_w}} {self.cost_initial:>{val_w}.4f}")
        print(f"  {'Final':<{key_w}} {self.cost_final:>{val_w}.4f}")
        print(f"  {'Improvement':<{key_w}} {diag['cost_improvement']:>{val_w}.4f}")
        print("\nEnsemble diversity:")
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
        if len(diag["warnings"]) > 0:
            print(f"\nWARNINGS ({len(diag['warnings'])}):")
            for i, warning in enumerate(diag["warnings"], 1):
                print(f"  {i}. {warning}")
        else:
            print(f"\nStatus: {diag['status']} - No issues detected")
        print("=" * 60)


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
class BMECustom:
    """BME reweighting against an experimental vector with a custom cost.

    A "BME variant" that reweights an ensemble against a whole experimental
    *profile* (e.g. a SAXS or PRE curve), with an **optional user-supplied
    goodness-of-fit function**. It is the penalty form of Bayesian Ensemble
    Refinement: it minimises

    .. math::

        \\mathcal{L}(w) = \\text{cost}(w) + \\theta\\, D_{\\mathrm{KL}}(w\\,\\|\\,w^0)

    directly over the ``n`` frame weights (subject to the simplex
    ``w_i \\ge 0``, ``\\sum_i w_i = 1``) using SciPy's ``trust-constr``. The
    entropy penalty ``theta`` plays the same role as in :class:`BME`: large
    ``theta`` keeps the ensemble close to the prior, small ``theta`` fits
    the data harder.

    Unlike :class:`BME` (which exploits a closed-form exponential solution
    special to the Gaussian chi-squared), an arbitrary cost has no such
    closed form, so the weights are optimised directly. The default cost
    reproduces :class:`BME`'s reduced chi-squared exactly.

    Parameters
    ----------
    experiment : numpy.ndarray
        Experimental vector, shape ``(m,)`` (e.g. one value per SAXS
        ``q``-point or per PRE residue).
    calculated_values : numpy.ndarray
        Per-conformer calculated values for the same observable, shape
        ``(n_frames, m)``.
    uncertainty : float or numpy.ndarray, optional
        Experimental uncertainty used by the **default** chi-squared cost:
        a scalar (broadcast to all ``m`` points) or a length-``m`` vector.
        Ignored when ``cost_function`` is supplied. Defaults to ``1.0``.
    cost_function : callable, optional
        Custom goodness-of-fit, called as
        ``cost_function(experiment, calculated_values, weights) -> float``,
        where ``experiment`` is ``(m,)``, ``calculated_values`` is
        ``(n_frames, m)`` and ``weights`` is the current ``(n_frames,)``
        weight vector. Lower is better. If ``None`` (default), the reduced
        chi-squared ``mean_k ((<calc>_k - V_k)/sigma_k)^2`` is used, where
        ``<calc>_k = sum_i w_i calc[i, k]``.
    initial_weights : numpy.ndarray, optional
        Prior frame weights. If None, uniform weights are used. Internally
        normalised to sum to 1.

    Raises
    ------
    SSException
        If inputs are malformed.

    Notes
    -----
    With a custom ``cost_function`` the optimiser uses a finite-difference
    gradient, so very large ensembles or expensive cost functions can be
    slow; the default chi-squared path uses an analytic gradient. The
    problem is convex (unique optimum) when the cost is convex in ``w``
    (the default chi-squared is); an arbitrary cost is optimised locally.

    Examples
    --------
    >>> import numpy as np
    >>> from soursop.ssbme import BMECustom
    >>> # SAXS-like profile: m points, n conformers
    >>> bme = BMECustom(I_exp, calc_I, uncertainty=sigma_exp)
    >>> result = bme.fit(theta=1.0)
    >>> weights = result.weights
    >>>
    >>> # custom cost: chi-squared on a log scale
    >>> def log_chi2(exp, calc, w):
    ...     avg = w @ calc
    ...     return float(np.mean((np.log(avg) - np.log(exp)) ** 2))
    >>> result = BMECustom(I_exp, calc_I, cost_function=log_chi2).fit(theta=1.0)
    """

    def __init__(
        self,
        experiment: np.ndarray,
        calculated_values: np.ndarray,
        uncertainty=None,
        cost_function=None,
        initial_weights: Optional[np.ndarray] = None,
    ):
        experiment = np.asarray(experiment, dtype=np.float64)
        calculated_values = np.asarray(calculated_values, dtype=np.float64)

        if experiment.ndim != 1:
            raise SSException(
                f"experiment must be a 1D vector, got shape {experiment.shape}"
            )
        if calculated_values.ndim != 2:
            raise SSException(
                "calculated_values must be 2D (n_frames, m), got shape "
                f"{calculated_values.shape}"
            )
        self.m = experiment.shape[0]
        self.n_frames = calculated_values.shape[0]
        if calculated_values.shape[1] != self.m:
            raise SSException(
                f"calculated_values columns ({calculated_values.shape[1]}) "
                f"must match the experiment length ({self.m})"
            )
        if cost_function is not None and not callable(cost_function):
            raise SSException("cost_function must be callable or None")

        # Uncertainty (used only by the default chi-squared cost).
        if uncertainty is None:
            sigma = np.ones(self.m, dtype=np.float64)
        else:
            sigma = np.asarray(uncertainty, dtype=np.float64)
            if sigma.ndim == 0:
                sigma = np.full(self.m, float(sigma))
            elif sigma.shape != (self.m,):
                raise SSException(
                    "uncertainty must be a scalar or a length-m vector "
                    f"(m={self.m}), got shape {sigma.shape}"
                )
            if np.any(sigma <= 0):
                raise SSException("uncertainty values must be positive")

        self.experiment = experiment
        self.calculated_values = calculated_values
        self.uncertainty = sigma
        self.cost_function = cost_function

        if initial_weights is None:
            self.initial_weights = np.ones(self.n_frames) / self.n_frames
        else:
            initial_weights = np.asarray(initial_weights, dtype=np.float64)
            if len(initial_weights) != self.n_frames:
                raise SSException("initial_weights length must match number of frames")
            self.initial_weights = initial_weights / np.sum(initial_weights)

        self._result: Optional[BMECustomResult] = None
        self._theta: Optional[float] = None
        self._theta_scan_result: Optional[ThetaScanResult] = None

    # ------------------------------------------------------------------
    def _cost(self, weights: np.ndarray) -> float:
        """Evaluate the (user or default) cost at ``weights``."""
        if self.cost_function is None:
            avg = weights @ self.calculated_values
            diff = (avg - self.experiment) / self.uncertainty
            return float(np.mean(diff**2))
        return float(
            self.cost_function(self.experiment, self.calculated_values, weights)
        )

    def _default_cost_grad(self, weights: np.ndarray) -> np.ndarray:
        """Analytic gradient of the default chi-squared cost."""
        avg = weights @ self.calculated_values
        coef = (2.0 / self.m) * (avg - self.experiment) / self.uncertainty**2
        return self.calculated_values @ coef

    # ------------------------------------------------------------------
    def fit(
        self,
        theta: float = 1.0,
        max_iterations: int = 2000,
        optimizer: str = "L-BFGS-B",
        verbose: bool = True,
    ) -> BMECustomResult:
        """Reweight by minimising ``cost(w) + theta * D_KL(w || w0)``.

        Internally the simplex ``w_i >= 0``, ``sum_i w_i = 1`` is enforced
        via a softmax reparameterisation (``w = softmax(z)``), reducing the
        problem to an unconstrained minimisation over ``z`` solved by
        L-BFGS-B. This avoids boundary kinks in ``log w`` and gives clean
        convergence on both the analytic-cost and custom-cost paths.

        Parameters
        ----------
        theta : float, optional
            Entropy-penalty strength (must be positive). Default 1.0.
        max_iterations : int, optional
            Maximum optimizer iterations. Default 2000.
        optimizer : str, optional
            scipy.optimize.minimize method for the unconstrained problem
            in the softmax-reparameterised variable ``z``. Default
            ``"L-BFGS-B"``.
        verbose : bool, optional
            Print progress. Default True.

        Returns
        -------
        BMECustomResult

        Raises
        ------
        SSException
            If ``theta`` is not positive.
        """
        if theta <= 0:
            raise SSException(f"theta must be positive, got {theta}")
        self._theta = float(theta)
        n = self.n_frames
        w0 = self.initial_weights.copy()
        log_w0 = np.log(np.maximum(w0, MIN_WEIGHT_THRESHOLD))

        def softmax(z):
            # Numerically stable softmax.
            z = z - np.max(z)
            ez = np.exp(z)
            return ez / np.sum(ez)

        def objective(z):
            w = softmax(z)
            # KL(w || w0) = sum_i w_i (log w_i - log w0_i); with softmax
            # log w_i = z_i - logsumexp(z).
            log_w = z - logsumexp(z)
            kl = float(np.sum(w * (log_w - log_w0)))
            return self._cost(w) + theta * kl

        cost_initial = self._cost(w0)
        if verbose:
            print("BMECustom Optimization")
            print(f"  Frames: {n}, Observables (vector length): {self.m}")
            custom = self.cost_function is not None
            print(f"  Cost function: {'custom' if custom else 'default chi^2'}")
            print(f"  Theta: {theta}")
            print(f"  Cost initial: {cost_initial:.4f}")

        # Initial parameter: z = log w0 -> softmax(z) = w0.
        z0 = log_w0.copy()

        opt = minimize(
            objective,
            z0,
            method=optimizer,
            jac="2-point",
            options={"maxiter": max_iterations, "ftol": 1e-9, "gtol": 1e-7},
        )

        w_opt = softmax(opt.x)
        cost_final = self._cost(w_opt)
        phi = float(np.exp(-_srel(w0, w_opt)))

        # L-BFGS-B reports success when it converges by ``ftol`` / ``gtol``.
        # Fall back to "did we end up no worse than the prior?" so a stalled
        # but improving fit isn't reported as a failure.
        success = bool(opt.success) or cost_final <= cost_initial + 1e-9

        if verbose:
            print(f"  Cost final: {cost_final:.4f}")
            print(f"  Effective fraction (phi): {phi:.4f}")
            print(f"  Optimization successful: {success}")
            print(f"  Optimizer message: {opt.message}")

        self._result = BMECustomResult(
            weights=w_opt,
            initial_weights=w0,
            cost_initial=float(cost_initial),
            cost_final=float(cost_final),
            phi=phi,
            n_iterations=int(getattr(opt, "nit", -1)),
            success=success,
            message=str(opt.message),
            theta=self._theta,
            experiment=self.experiment,
            calculated_values=self.calculated_values,
            metadata={
                "optimizer": optimizer,
                "max_iterations": max_iterations,
                "custom_cost": self.cost_function is not None,
                "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            },
        )
        return self._result

    # ------------------------------------------------------------------
    def scan_theta(
        self,
        theta_range: Union[Tuple[float, float], np.ndarray] = (0.01, 100.0),
        n_points: int = 12,
        log_scale: bool = True,
        method: str = "perpendicular",
        verbose: bool = False,
    ) -> ThetaScanResult:
        """Scan ``theta`` and pick the cost vs. relative-entropy knee.

        The L-curve analogue of :meth:`BME.scan_theta`. The returned
        :class:`ThetaScanResult` stores the cost in its
        ``chi_squared_values`` field (the generic fit metric for this
        variant).

        Parameters
        ----------
        theta_range, n_points, log_scale, method, verbose
            See :func:`theta_scan`. Default range ``(0.01, 100.0)``.

        Returns
        -------
        ThetaScanResult
        """
        if isinstance(theta_range, (tuple, list)):
            if n_points < 1:
                raise SSException(
                    f"n_points must be >= 1 for a tuple theta_range, got {n_points}"
                )
            if log_scale and (theta_range[0] <= 0 or theta_range[1] <= 0):
                raise SSException(
                    "theta_range endpoints must be positive when "
                    f"log_scale=True, got {theta_range}"
                )
            if log_scale:
                thetas = np.logspace(
                    np.log10(theta_range[0]), np.log10(theta_range[1]), n_points
                )
            else:
                thetas = np.linspace(theta_range[0], theta_range[1], n_points)
        else:
            thetas = np.asarray(theta_range, dtype=np.float64)

        costs: List[float] = []
        phis: List[float] = []
        kls: List[float] = []
        results: List[BMECustomResult] = []
        for i, t in enumerate(thetas):
            if verbose:
                print(f"Processing theta {i + 1}/{len(thetas)}: {t:.4f}")
            res = self.fit(theta=float(t), verbose=False)
            results.append(res)
            costs.append(res.cost_final)
            phis.append(res.phi)
            kls.append(res.kl_divergence)

        cost_arr = np.array(costs)
        kl_arr = np.array(kls)
        optimal_idx, method_name = find_optimal_theta(cost_arr, kl_arr, method=method)
        scan = ThetaScanResult(
            theta_values=thetas,
            chi_squared_values=cost_arr,
            phi_values=np.array(phis),
            kl_divergence_values=kl_arr,
            results=results,
            optimal_theta=float(thetas[optimal_idx]),
            optimal_idx=optimal_idx,
            method=method_name,
        )
        self._theta_scan_result = scan
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
    def result(self) -> Optional[BMECustomResult]:
        """The most recent :class:`BMECustomResult`, or None."""
        return self._result

    @property
    def theta(self) -> Optional[float]:
        """Theta used in the most recent fit, or None."""
        return self._theta

    @property
    def theta_scan_result(self) -> Optional[ThetaScanResult]:
        """The most recent :class:`ThetaScanResult`, or None."""
        return self._theta_scan_result
