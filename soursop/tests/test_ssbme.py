"""
Tests for soursop.ssbme - Bayesian Maximum Entropy (BME) and iterative
BME (iBME) ensemble reweighting.

Covers:
  * ExperimentalObservable validation -> SSException
  * BME maxent fit reduces chi-squared; weights are a valid distribution
  * constraint handling (upper/lower already satisfied -> no perturbation)
  * iBME recovers an unknown global scale/offset and beats plain BME on
    scale-mismatched data; converges within max_ibme_iterations
  * theta_scan for both reweighters (shapes, phi-vs-theta trend, knee)
  * BMEResult.predict / reweighter.predict shape + value
"""

import numpy as np
import pytest

from soursop.ssbme import (
    BME,
    BMECustom,
    BMECustomResult,
    BMEResult,
    ExperimentalObservable,
    ThetaScanResult,
    iBME,
    theta_scan,
)
from soursop.ssexceptions import SSException


# --------------------------------------------------------------------------
# Synthetic ensemble helpers
# --------------------------------------------------------------------------
def _make_ensemble(n_frames=2000, seed=42):
    """Two-observable Gaussian ensemble with known means."""
    rng = np.random.default_rng(seed)
    calc = rng.normal(loc=[24.0, 65.0], scale=[3.0, 6.0], size=(n_frames, 2))
    return calc


# --------------------------------------------------------------------------
# ExperimentalObservable
# --------------------------------------------------------------------------
class TestExperimentalObservable:
    def test_valid(self):
        obs = ExperimentalObservable(23.0, 1.0, name="Rg")
        assert obs.constraint == "equality"
        assert obs.get_bounds() == (None, None)

    def test_constraint_normalization_and_bounds(self):
        assert ExperimentalObservable(1, 1, "UPPER").get_bounds() == (0.0, None)
        assert ExperimentalObservable(1, 1, " Lower ").get_bounds() == (None, 0.0)

    def test_bad_uncertainty_raises(self):
        with pytest.raises(SSException):
            ExperimentalObservable(10.0, 0.0)
        with pytest.raises(SSException):
            ExperimentalObservable(10.0, -1.0)

    def test_bad_constraint_raises(self):
        with pytest.raises(SSException):
            ExperimentalObservable(10.0, 1.0, constraint="nonsense")
        with pytest.raises(SSException):
            ExperimentalObservable(10.0, 1.0, constraint=5)


# --------------------------------------------------------------------------
# Input validation
# --------------------------------------------------------------------------
class TestInputValidation:
    def test_empty_observables(self):
        with pytest.raises(SSException):
            BME([], _make_ensemble())

    def test_non_observable_entries(self):
        with pytest.raises(SSException):
            BME([1, 2], _make_ensemble())

    def test_shape_mismatch(self):
        obs = [ExperimentalObservable(23.0, 1.0)]
        with pytest.raises(SSException):
            BME(obs, _make_ensemble())  # 2 columns, 1 observable

    def test_initial_weights_length(self):
        obs = [ExperimentalObservable(23.0, 1.0), ExperimentalObservable(60.0, 2.0)]
        with pytest.raises(SSException):
            BME(obs, _make_ensemble(), initial_weights=np.ones(10))


# --------------------------------------------------------------------------
# BME core
# --------------------------------------------------------------------------
class TestBME:
    def test_fit_reduces_chi2_and_valid_weights(self):
        np.random.seed(0)
        calc = _make_ensemble()
        # Pull both observables away from the ensemble means.
        obs = [
            ExperimentalObservable(21.0, 1.0, name="Rg"),
            ExperimentalObservable(58.0, 2.0, name="Ree"),
        ]
        res = BME(obs, calc).fit(theta=2.0, auto_theta=False, verbose=False)

        assert isinstance(res, BMEResult)
        assert res.success
        assert np.isclose(np.sum(res.weights), 1.0)
        assert np.all(res.weights >= 0)
        assert res.chi_squared_final < res.chi_squared_initial
        assert 0.0 < res.phi <= 1.0

    def test_predict_matches_weighted_mean(self):
        np.random.seed(0)
        calc = _make_ensemble()
        obs = [ExperimentalObservable(21.0, 1.0), ExperimentalObservable(58.0, 2.0)]
        bme = BME(obs, calc)
        res = bme.fit(theta=2.0, auto_theta=False, verbose=False)

        pred = bme.predict(calc)
        expected = np.sum(res.weights[:, None] * calc, axis=0)
        assert pred.shape == (2,)
        assert np.allclose(pred, expected)
        # Reweighting should pull the Rg mean toward the experimental value.
        assert abs(pred[0] - 21.0) < abs(np.mean(calc[:, 0]) - 21.0)

    def test_predict_before_fit_raises(self):
        obs = [ExperimentalObservable(21.0, 1.0), ExperimentalObservable(58.0, 2.0)]
        with pytest.raises(SSException):
            BME(obs, _make_ensemble()).predict(_make_ensemble())

    def test_upper_constraint_already_satisfied(self):
        """An upper bound far above the ensemble mean barely perturbs it."""
        np.random.seed(0)
        calc = _make_ensemble()
        obs = [
            ExperimentalObservable(24.0, 1.0, name="Rg"),
            ExperimentalObservable(200.0, 5.0, constraint="upper", name="Ree<=200"),
        ]
        res = BME(obs, calc).fit(theta=5.0, auto_theta=False, verbose=False)
        # Ree mean (~65) is well below 200 -> contributes ~0 to chi2;
        # weights stay close to uniform and phi stays high.
        assert res.phi > 0.95
        assert np.allclose(res.weights, res.initial_weights, atol=1e-3)


# --------------------------------------------------------------------------
# iBME
# --------------------------------------------------------------------------
class TestIBME:
    def _scaled_problem(self, seed=1):
        """calc = a*true + b with known (a, b); exp ~ mean(true)."""
        rng = np.random.default_rng(seed)
        true = rng.normal(loc=[30.0, 80.0], scale=[4.0, 7.0], size=(1500, 2))
        a, b = 3.0, -10.0
        calc = a * true + b
        # Targets slightly off the unweighted mean so reweighting matters.
        obs = [
            ExperimentalObservable(28.5, 0.8, name="o1"),
            ExperimentalObservable(78.0, 1.5, name="o2"),
        ]
        return true, calc, obs, a, b

    def test_recovers_scale_and_offset(self):
        np.random.seed(0)
        true, calc, obs, a, b = self._scaled_problem()
        ib = iBME(obs, calc)
        res = ib.fit(theta=5.0, ftol=1e-3, max_ibme_iterations=50, verbose=False)

        assert res.success
        # Net mapping should invert the global scale: alpha ~ 1/a.
        assert res.scale == pytest.approx(1.0 / a, rel=0.05)
        # Applying the recovered scale/offset to the calculated values
        # and reweighting should reproduce the experimental targets.
        scaled = res.scale * calc + res.offset
        pred = np.sum(res.weights[:, None] * scaled, axis=0)
        exp_vals = np.array([o.value for o in obs])
        assert np.allclose(pred, exp_vals, atol=0.5)
        assert np.isclose(np.sum(res.weights), 1.0)
        assert 0.0 < res.phi <= 1.0

    def test_converges_before_max_iterations(self):
        np.random.seed(0)
        _, calc, obs, _, _ = self._scaled_problem()
        ib = iBME(obs, calc)
        res = ib.fit(theta=5.0, ftol=1e-3, max_ibme_iterations=50, verbose=False)
        assert res.n_iterations < 50
        assert len(ib.get_ibme_stats()) == res.n_iterations

    def test_beats_plain_bme_on_scaled_data(self):
        np.random.seed(0)
        _, calc, obs, _, _ = self._scaled_problem()

        bme_res = BME(obs, calc).fit(theta=5.0, auto_theta=False, verbose=False)
        ibme_res = iBME(obs, calc).fit(theta=5.0, ftol=1e-3, verbose=False)

        # With an unknown global scale/offset, iBME should fit the data
        # markedly better than plain BME.
        assert ibme_res.chi_squared_final < bme_res.chi_squared_final

    def test_stats_before_fit_raises(self):
        _, calc, obs, _, _ = self._scaled_problem()
        with pytest.raises(SSException):
            iBME(obs, calc).get_ibme_weights()


# --------------------------------------------------------------------------
# theta_scan
# --------------------------------------------------------------------------
class TestThetaScan:
    def test_bme_scan_shapes_and_knee(self):
        np.random.seed(0)
        calc = _make_ensemble(n_frames=800)
        obs = [ExperimentalObservable(21.0, 1.0), ExperimentalObservable(58.0, 2.0)]
        scan = theta_scan(
            obs, calc, reweighter="bme", theta_range=(0.05, 20.0), n_points=8
        )

        assert scan.theta_values.shape == (8,)
        assert scan.chi_squared_values.shape == (8,)
        assert scan.phi_values.shape == (8,)
        assert 0 <= scan.optimal_idx < 8
        assert scan.optimal_theta == scan.theta_values[scan.optimal_idx]
        # Larger theta = more regularisation = higher effective fraction.
        assert scan.phi_values[-1] > scan.phi_values[0]

    def test_ibme_scan_runs(self):
        np.random.seed(0)
        rng = np.random.default_rng(2)
        true = rng.normal([30.0, 80.0], [4.0, 7.0], size=(600, 2))
        calc = 2.0 * true + 5.0
        obs = [ExperimentalObservable(28.5, 1.0), ExperimentalObservable(78.0, 2.0)]
        scan = theta_scan(
            obs,
            calc,
            reweighter="ibme",
            theta_range=(0.1, 10.0),
            n_points=5,
            fit_kwargs={"ftol": 1e-2, "max_ibme_iterations": 30},
        )
        assert scan.theta_values.shape == (5,)
        assert all(r.scale is not None for r in scan.results)

    def test_unknown_reweighter_raises(self):
        obs = [ExperimentalObservable(21.0, 1.0), ExperimentalObservable(58.0, 2.0)]
        with pytest.raises(SSException):
            theta_scan(obs, _make_ensemble(), reweighter="bogus")


# --------------------------------------------------------------------------
# BMECustom (raw vector + matrix, optional custom cost)
# --------------------------------------------------------------------------
def _make_profile_problem(n_frames=200, m=8, seed=3):
    """Synthetic SAXS-like 'profile': m experimental points, n conformers."""
    rng = np.random.default_rng(seed)
    q = np.linspace(0.05, 0.3, m)
    # True intensity (Debye-like decay) + per-conformer log-normal noise.
    true = 100.0 * np.exp(-30.0 * q**2)
    noise = rng.lognormal(mean=0.0, sigma=0.15, size=(n_frames, m))
    calc = true[None, :] * noise  # (n, m)
    # Experimental targets: a known shifted shape that *is* reachable.
    exp = 100.0 * np.exp(-25.0 * q**2)
    sigma = 0.05 * exp
    return exp, calc, sigma


class TestBMECustom:
    def test_default_cost_reduces_chi2(self):
        np.random.seed(0)
        exp, calc, sigma = _make_profile_problem()
        res = BMECustom(exp, calc, uncertainty=sigma).fit(theta=0.5, verbose=False)
        assert isinstance(res, BMECustomResult)
        assert res.success
        assert np.isclose(np.sum(res.weights), 1.0)
        assert np.all(res.weights >= -1e-12)
        assert res.cost_final < res.cost_initial
        assert 0.0 < res.phi <= 1.0
        assert res.reweighting_factors.shape == (calc.shape[0],)

    def test_custom_cost_matches_default(self):
        """A user cost identical to the default chi^2 must give the same fit."""
        np.random.seed(0)
        exp, calc, sigma = _make_profile_problem()

        def reduced_chi2(experiment, calc_matrix, weights):
            avg = weights @ calc_matrix
            return float(np.mean(((avg - experiment) / sigma) ** 2))

        res_def = BMECustom(exp, calc, uncertainty=sigma).fit(theta=0.5, verbose=False)
        res_cus = BMECustom(exp, calc, cost_function=reduced_chi2).fit(
            theta=0.5, verbose=False
        )
        # Same global optimum (convex problem); both paths converge to the
        # same cost within optimiser tolerance. The weight vectors are
        # close but not identical because the default uses an analytic
        # gradient while the custom path uses a finite-difference gradient.
        assert res_def.cost_final == pytest.approx(res_cus.cost_final, rel=5e-3)
        assert res_def.phi == pytest.approx(res_cus.phi, rel=5e-2)

    def test_custom_cost_runs_and_reduces(self):
        """A non-default cost (log-scale chi^2) should still drive the fit."""
        np.random.seed(0)
        exp, calc, _ = _make_profile_problem()

        def log_chi2(experiment, calc_matrix, weights):
            avg = weights @ calc_matrix
            return float(np.mean((np.log(avg) - np.log(experiment)) ** 2))

        res = BMECustom(exp, calc, cost_function=log_chi2).fit(theta=0.5, verbose=False)
        assert res.success
        assert res.cost_final < res.cost_initial
        assert np.isclose(np.sum(res.weights), 1.0)

    def test_scalar_uncertainty(self):
        np.random.seed(0)
        exp, calc, _ = _make_profile_problem()
        # Scalar uncertainty broadcasts to all m points.
        res = BMECustom(exp, calc, uncertainty=5.0).fit(theta=0.5, verbose=False)
        assert res.success
        assert np.isclose(np.sum(res.weights), 1.0)

    def test_predict_matches_weighted_mean(self):
        np.random.seed(0)
        exp, calc, sigma = _make_profile_problem()
        bme = BMECustom(exp, calc, uncertainty=sigma)
        res = bme.fit(theta=0.5, verbose=False)
        pred = bme.predict(calc)
        expected = res.weights @ calc
        assert pred.shape == (calc.shape[1],)
        assert np.allclose(pred, expected)

    def test_predict_before_fit_raises(self):
        exp, calc, sigma = _make_profile_problem()
        with pytest.raises(SSException):
            BMECustom(exp, calc, uncertainty=sigma).predict(calc)

    def test_invalid_inputs_raise(self):
        exp, calc, sigma = _make_profile_problem()
        # Wrong vector dim
        with pytest.raises(SSException):
            BMECustom(calc, calc)
        # Wrong matrix dim
        with pytest.raises(SSException):
            BMECustom(exp, exp)
        # Mismatched m
        with pytest.raises(SSException):
            BMECustom(exp[:-1], calc)
        # Non-callable cost
        with pytest.raises(SSException):
            BMECustom(exp, calc, cost_function=42)
        # Bad uncertainty shape
        with pytest.raises(SSException):
            BMECustom(exp, calc, uncertainty=np.ones(exp.shape[0] + 1))
        # Non-positive uncertainty
        with pytest.raises(SSException):
            BMECustom(exp, calc, uncertainty=0.0)
        # Bad theta
        with pytest.raises(SSException):
            BMECustom(exp, calc, uncertainty=sigma).fit(theta=0.0, verbose=False)

    def test_scan_theta(self):
        np.random.seed(0)
        exp, calc, sigma = _make_profile_problem(n_frames=120)
        bme = BMECustom(exp, calc, uncertainty=sigma)
        scan = bme.scan_theta(theta_range=(0.05, 5.0), n_points=4)
        assert isinstance(scan, ThetaScanResult)
        assert scan.theta_values.shape == (4,)
        # chi_squared_values is overloaded to store the cost for BMECustom.
        assert scan.chi_squared_values.shape == (4,)
        assert 0 <= scan.optimal_idx < 4
