"""
Tests for soursop.sscoper - Convex Optimization for Ensemble Reweighting
(COPER) and its iterative variant (iCOPER).

Covers:
  * ExperimentalObservable (shared via ssutils) with the new ``group`` field;
    BME ignores ``group``.
  * COPER feasibility: valid simplex, chi^2_final <= chi2_limit, phi in
    (0, 1], delta_S <= 0, reweighting_factors shape.
  * COPER infeasibility: targets outside ensemble range -> feasible=False,
    success=False, weights equal the chi^2-min point, informative message.
  * Constraint handling: an already-satisfied upper bound barely perturbs
    weights (phi ~ 1).
  * Per-data-type chi^2 (``group=`` field) constrains each group.
  * predict() shape + equals sum(w * calc).
  * chi2_limit_scan returns arrays of correct length, |KL| decreases as the
    chi^2 limit is loosened, optimal_idx valid.
  * iCOPER recovers an unknown global scale on calc = alpha * true + beta,
    converges within max_icoper_iterations, and logs per-iteration stats.
  * Sanity parity with BME on the same satisfiable problem (both produce
    valid weights and reduce chi^2; the methods need not agree numerically).
"""

import numpy as np
import pytest

from soursop.ssbme import BME
from soursop.sscoper import (
    COPER,
    COPERResult,
    COPERScanResult,
    ExperimentalObservable,
    chi2_limit_scan,
    iCOPER,
)
from soursop.ssexceptions import SSException


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
def _make_ensemble(n_frames=200, seed=42):
    """Two-observable Gaussian ensemble with known means."""
    rng = np.random.default_rng(seed)
    return rng.normal(loc=[24.0, 65.0], scale=[3.0, 6.0], size=(n_frames, 2))


# --------------------------------------------------------------------------
# ExperimentalObservable (shared with ssbme via ssutils)
# --------------------------------------------------------------------------
class TestExperimentalObservable:
    def test_group_field_default_none(self):
        obs = ExperimentalObservable(23.0, 1.0, name="Rg")
        assert obs.group is None

    def test_group_field_kept(self):
        obs = ExperimentalObservable(23.0, 1.0, name="Rg", group="size")
        assert obs.group == "size"

    def test_bad_uncertainty_raises(self):
        with pytest.raises(SSException):
            ExperimentalObservable(10.0, 0.0)

    def test_bme_ignores_group(self):
        """BME must still work on observables carrying a ``group=`` tag."""
        np.random.seed(0)
        calc = _make_ensemble()
        obs = [
            ExperimentalObservable(21.0, 1.0, name="Rg", group="size"),
            ExperimentalObservable(58.0, 2.0, name="Ree", group="size"),
        ]
        res = BME(obs, calc).fit(theta=2.0, auto_theta=False, verbose=False)
        assert res.success
        assert np.isclose(np.sum(res.weights), 1.0)


# --------------------------------------------------------------------------
# COPER core
# --------------------------------------------------------------------------
class TestCOPER:
    def test_feasible_fit(self):
        np.random.seed(0)
        calc = _make_ensemble()
        # Targets a little off the ensemble means: feasible at chi^2 <= 1.
        obs = [
            ExperimentalObservable(23.0, 1.0, name="Rg"),
            ExperimentalObservable(63.0, 2.0, name="Ree"),
        ]
        res = COPER(obs, calc).fit(chi2_limit=1.0, verbose=False)

        assert isinstance(res, COPERResult)
        assert res.success
        assert res.feasible
        assert np.isclose(np.sum(res.weights), 1.0)
        assert np.all(res.weights >= -1e-12)
        assert res.chi_squared_final <= 1.0 + 1e-3
        assert 0.0 < res.phi <= 1.0
        assert res.delta_S <= 1e-9
        assert res.mean_delta_G_kT >= -1e-9
        # reweighting factors shape
        assert res.reweighting_factors.shape == (calc.shape[0],)

    def test_predict_matches_weighted_mean(self):
        np.random.seed(0)
        calc = _make_ensemble()
        obs = [
            ExperimentalObservable(23.0, 1.0, name="Rg"),
            ExperimentalObservable(63.0, 2.0, name="Ree"),
        ]
        coper = COPER(obs, calc)
        res = coper.fit(chi2_limit=1.0, verbose=False)
        pred = coper.predict(calc)
        expected = np.sum(res.weights[:, None] * calc, axis=0)
        assert pred.shape == (2,)
        assert np.allclose(pred, expected)

    def test_predict_before_fit_raises(self):
        obs = [
            ExperimentalObservable(23.0, 1.0),
            ExperimentalObservable(63.0, 2.0),
        ]
        with pytest.raises(SSException):
            COPER(obs, _make_ensemble()).predict(_make_ensemble())

    def test_infeasible_returns_chi2_min(self):
        """Targets far outside the ensemble range -> infeasible."""
        np.random.seed(0)
        calc = _make_ensemble()
        obs = [
            ExperimentalObservable(100.0, 1.0, name="Rg-impossible"),
            ExperimentalObservable(63.0, 2.0, name="Ree"),
        ]
        res = COPER(obs, calc).fit(chi2_limit=1.0, verbose=False)
        assert not res.feasible
        assert not res.success
        assert "infeasible" in res.message.lower()
        # chi^2_final equals chi^2_min for infeasible results.
        assert res.chi_squared_final == res.chi_squared_min
        # No reweighting can satisfy the constraint, but the chi^2-min point
        # is still a valid probability vector.
        assert np.isclose(np.sum(res.weights), 1.0)
        assert np.all(res.weights >= -1e-12)

    def test_upper_constraint_already_satisfied(self):
        np.random.seed(0)
        calc = _make_ensemble()
        # Equality observable already very close to the ensemble mean, plus
        # a slack upper bound -> phi should stay close to 1.
        obs = [
            ExperimentalObservable(24.0, 1.0, name="Rg"),
            ExperimentalObservable(200.0, 5.0, constraint="upper", name="Ree<=200"),
        ]
        res = COPER(obs, calc).fit(chi2_limit=5.0, verbose=False)
        assert res.feasible
        assert res.phi > 0.95

    def test_per_data_type_chi2(self):
        """With two groups, each chi^2_alpha is constrained <= limit."""
        np.random.seed(0)
        calc = _make_ensemble()
        obs = [
            ExperimentalObservable(23.0, 1.0, name="Rg", group="size"),
            ExperimentalObservable(63.0, 2.0, name="Ree", group="size"),
            ExperimentalObservable(22.0, 1.5, name="Rg_alt", group="extra"),
        ]
        # Three-column "calc": reuse Rg column for the extra group.
        calc3 = np.column_stack([calc[:, 0], calc[:, 1], calc[:, 0]])
        coper = COPER(obs, calc3)
        res = coper.fit(chi2_limit=1.0, verbose=False)
        assert res.feasible
        # Every per-group chi^2 must respect the limit.
        per_group = coper._chi2_per_group_vec(res.weights)
        assert per_group.shape == (2,)
        assert np.all(per_group <= 1.0 + 1e-3)

    def test_bad_chi2_limit_raises(self):
        obs = [
            ExperimentalObservable(23.0, 1.0),
            ExperimentalObservable(63.0, 2.0),
        ]
        with pytest.raises(SSException):
            COPER(obs, _make_ensemble()).fit(chi2_limit=0.0, verbose=False)


# --------------------------------------------------------------------------
# iCOPER
# --------------------------------------------------------------------------
class TestICOPER:
    def _scaled_problem(self, seed=1):
        """calc = a * true + b with known (a, b); experimental ~ mean(true)."""
        rng = np.random.default_rng(seed)
        true = rng.normal(loc=[30.0, 80.0], scale=[4.0, 7.0], size=(200, 2))
        a, b = 3.0, -10.0
        calc = a * true + b
        obs = [
            ExperimentalObservable(29.0, 0.8, name="o1"),
            ExperimentalObservable(78.0, 1.5, name="o2"),
        ]
        return true, calc, obs, a, b

    def test_recovers_scale(self):
        np.random.seed(0)
        _, calc, obs, a, _b = self._scaled_problem()
        ic = iCOPER(obs, calc)
        res = ic.fit(
            chi2_limit=1.0,
            ftol=1e-3,
            max_icoper_iterations=30,
            verbose=False,
        )
        assert res.success
        assert res.scale == pytest.approx(1.0 / a, rel=0.1)
        # After applying the recovered net mapping the reweighted ensemble
        # should reproduce the experimental targets within ~uncertainty.
        scaled = res.scale * calc + res.offset
        pred = np.sum(res.weights[:, None] * scaled, axis=0)
        exp_vals = np.array([o.value for o in obs])
        assert np.allclose(pred, exp_vals, atol=1.0)
        assert np.isclose(np.sum(res.weights), 1.0)
        assert 0.0 < res.phi <= 1.0

    def test_converges_and_logs_stats(self):
        np.random.seed(0)
        _, calc, obs, _, _ = self._scaled_problem()
        ic = iCOPER(obs, calc)
        res = ic.fit(
            chi2_limit=1.0,
            ftol=1e-3,
            max_icoper_iterations=30,
            verbose=False,
        )
        assert res.n_iterations <= 30
        stats = ic.get_icoper_stats()
        assert len(stats) == res.n_iterations
        # Each iteration entry must carry the expected keys.
        for entry in stats:
            for key in (
                "iteration",
                "scale",
                "offset",
                "chi_squared",
                "feasible",
                "diff",
            ):
                assert key in entry

    def test_stats_before_fit_raises(self):
        _, calc, obs, _, _ = self._scaled_problem()
        with pytest.raises(SSException):
            iCOPER(obs, calc).get_icoper_weights()


# --------------------------------------------------------------------------
# chi2_limit_scan
# --------------------------------------------------------------------------
class TestChi2LimitScan:
    def test_scan_shapes_and_kl_trend(self):
        np.random.seed(0)
        calc = _make_ensemble(n_frames=150)
        obs = [
            ExperimentalObservable(22.5, 1.0, name="Rg"),
            ExperimentalObservable(62.0, 2.0, name="Ree"),
        ]
        scan = chi2_limit_scan(
            obs,
            calc,
            reweighter="coper",
            chi2_limits=(0.5, 4.0),
            n_points=5,
        )
        assert isinstance(scan, COPERScanResult)
        assert scan.chi2_limits.shape == (5,)
        assert scan.chi_squared_values.shape == (5,)
        assert scan.phi_values.shape == (5,)
        assert scan.feasible_mask.shape == (5,)
        assert 0 <= scan.optimal_idx < 5
        assert scan.optimal_chi2_limit == scan.chi2_limits[scan.optimal_idx]
        # A looser chi^2 limit allows a more uniform reweighted ensemble,
        # so the relative entropy at the largest limit is no larger than
        # at the smallest limit.
        feas_idx = np.where(scan.feasible_mask)[0]
        if len(feas_idx) >= 2:
            assert scan.kl_divergence_values[feas_idx[-1]] <= (
                scan.kl_divergence_values[feas_idx[0]] + 1e-6
            )

    def test_unknown_reweighter_raises(self):
        obs = [
            ExperimentalObservable(23.0, 1.0),
            ExperimentalObservable(63.0, 2.0),
        ]
        with pytest.raises(SSException):
            chi2_limit_scan(obs, _make_ensemble(), reweighter="bogus")


# --------------------------------------------------------------------------
# Cross-method sanity vs BME
# --------------------------------------------------------------------------
class TestCOPERvsBMESanity:
    def test_both_reduce_chi2(self):
        """COPER and BME should both reduce chi^2 on the same satisfiable data."""
        np.random.seed(0)
        calc = _make_ensemble()
        obs = [
            ExperimentalObservable(22.5, 1.0, name="Rg"),
            ExperimentalObservable(62.0, 2.0, name="Ree"),
        ]
        # COPER
        coper_res = COPER(obs, calc).fit(chi2_limit=1.0, verbose=False)
        # BME
        bme_res = BME(obs, calc).fit(theta=2.0, auto_theta=False, verbose=False)
        # Both succeeded and produced valid weight distributions.
        assert coper_res.feasible and coper_res.success
        assert bme_res.success
        assert np.isclose(np.sum(coper_res.weights), 1.0)
        assert np.isclose(np.sum(bme_res.weights), 1.0)
        # Both reduced chi-squared.
        assert coper_res.chi_squared_final <= coper_res.chi_squared_initial + 1e-6
        assert bme_res.chi_squared_final <= bme_res.chi_squared_initial + 1e-6
