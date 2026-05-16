"""
Tests for the consistent deterministic frame-`weights` capability and the
shared soursop.ssutils weight-vector validation / reduction helpers.

Covers the extremes:
  * uniform weights  == unweighted ensemble mean
  * one-hot weights  == that single frame's value
  * invalid weights  -> SSException (wrong length / <0 / >1 / sum!=1 /
    non-finite / non-numeric)
  * stride + weights  -> subsample then renormalise (sum == 1)
  * deterministic-everywhere edge cases that intentionally raise
"""

import numpy as np
import pytest

from soursop import ssutils
from soursop.ssexceptions import SSException


# --------------------------------------------------------------------------
# ssutils.validate_weights
# --------------------------------------------------------------------------
class TestValidateWeights:

    def test_false_and_none_are_noops(self):
        assert ssutils.validate_weights(False, 10) is False
        assert ssutils.validate_weights(None, 10) is False

    def test_uniform_ok_and_sums_to_one(self):
        w = ssutils.validate_weights(np.full(10, 0.1), 10)
        assert isinstance(w, np.ndarray)
        assert w.dtype == np.float64
        assert abs(w.sum() - 1.0) < 1e-12

    def test_wrong_length_raises(self):
        with pytest.raises(SSException):
            ssutils.validate_weights(np.full(9, 1.0 / 9), 10)

    def test_negative_element_raises(self):
        bad = np.full(10, 1.0 / 9)
        bad[0] = -1.0 / 9
        with pytest.raises(SSException):
            ssutils.validate_weights(bad, 10)

    def test_element_greater_than_one_raises(self):
        bad = np.zeros(10)
        bad[0] = 1.5
        with pytest.raises(SSException):
            ssutils.validate_weights(bad, 10)

    def test_sum_not_one_raises(self):
        with pytest.raises(SSException):
            ssutils.validate_weights(np.full(10, 0.05), 10)   # sums to 0.5

    def test_nonfinite_raises(self):
        bad = np.full(10, 0.1)
        bad[3] = np.nan
        with pytest.raises(SSException):
            ssutils.validate_weights(bad, 10)
        bad[3] = np.inf
        with pytest.raises(SSException):
            ssutils.validate_weights(bad, 10)

    def test_non_numeric_raises(self):
        with pytest.raises(SSException):
            ssutils.validate_weights("not a vector", 10)

    def test_etol_band(self):
        w = np.full(10, 0.1)
        w[0] += 1e-9          # sum within default etol (1e-7)
        assert ssutils.validate_weights(w, 10) is not False
        w2 = np.full(10, 0.1)
        w2[0] += 1e-3         # sum outside default etol
        with pytest.raises(SSException):
            ssutils.validate_weights(w2, 10)

    def test_stride_subsample_then_renormalise(self):
        n, stride = 12, 3
        raw = np.arange(1, n + 1, dtype=float)
        raw = raw / raw.sum()                       # valid full-length vector
        out = ssutils.validate_weights(raw, n, stride=stride)
        expected = raw[::stride] / raw[::stride].sum()
        assert len(out) == len(raw[::stride])
        assert abs(out.sum() - 1.0) < 1e-12
        assert np.allclose(out, expected, rtol=0, atol=1e-12)


# --------------------------------------------------------------------------
# deterministic reduction helpers
# --------------------------------------------------------------------------
class TestWeightedHelpers:

    def test_weighted_mean_matches_numpy_average(self):
        x = np.array([1.0, 3.0, 5.0, 7.0])
        w = np.array([0.1, 0.2, 0.3, 0.4])
        assert np.isclose(ssutils.weighted_mean(x, w), np.average(x, weights=w))

    def test_uniform_weighted_mean_equals_plain_mean(self):
        x = np.random.RandomState(0).rand(20)
        w = np.full(20, 1.0 / 20)
        assert np.isclose(ssutils.weighted_mean(x, w), x.mean(), rtol=1e-12)

    def test_weighted_rms(self):
        x = np.array([3.0, 4.0])
        w = np.array([0.5, 0.5])
        assert np.isclose(ssutils.weighted_rms(x, w), np.sqrt((9 + 16) / 2))

    def test_weighted_std_zero_for_constant(self):
        assert ssutils.weighted_std(np.full(6, 2.0), np.full(6, 1.0 / 6)) == 0.0

    def test_uniform_weighted_std_equals_population_std(self):
        x = np.random.RandomState(1).rand(15)
        w = np.full(15, 1.0 / 15)
        assert np.isclose(ssutils.weighted_std(x, w), x.std(), rtol=1e-12)

    def test_weighted_corr_perfect_correlation(self):
        a = np.array([1.0, 2.0, 3.0, 4.0])
        b = 2.0 * a + 1.0
        w = np.full(4, 0.25)
        assert np.isclose(ssutils.weighted_corr(a, b, w), 1.0, rtol=1e-9)


# --------------------------------------------------------------------------
# per-frame getters: uniform == unweighted mean ; one-hot == that frame
# --------------------------------------------------------------------------
class TestPerFrameGetters:

    def _scalar_getters(self, P):
        return {
            'rg':   lambda **k: P.get_radius_of_gyration(**k),
            'asph': lambda **k: P.get_asphericity(verbose=False, **k),
            'e2e':  lambda **k: P.get_end_to_end_distance(**k),
            'Rh':   lambda **k: P.get_hydrodynamic_radius(**k),
        }

    def test_uniform_equals_unweighted_mean(self, NTL9_CP):
        P = NTL9_CP
        n = P.n_frames
        uni = np.full(n, 1.0 / n)
        for name, fn in self._scalar_getters(P).items():
            base = np.asarray(fn())
            wmean = fn(weights=uni)
            assert np.allclose(wmean, base.mean(), rtol=1e-5, atol=1e-6), name

    def test_one_hot_equals_single_frame(self, GS6_CP):
        P = GS6_CP
        n = P.n_frames
        for name, fn in self._scalar_getters(P).items():
            base = np.asarray(fn())
            for k in range(n):
                oneh = np.zeros(n)
                oneh[k] = 1.0
                assert np.isclose(fn(weights=oneh), base[k], rtol=1e-9, atol=1e-9), (name, k)

    def test_gyration_tensor_uniform_and_onehot(self, GS6_CP):
        P = GS6_CP
        n = P.n_frames
        gt = np.asarray(P.get_gyration_tensor(verbose=False))
        uni = np.full(n, 1.0 / n)
        assert np.allclose(P.get_gyration_tensor(verbose=False, weights=uni),
                           gt.mean(axis=0), rtol=1e-5, atol=1e-6)
        oneh = np.zeros(n); oneh[1] = 1.0
        assert np.allclose(P.get_gyration_tensor(verbose=False, weights=oneh),
                           gt[1], rtol=1e-9, atol=1e-9)

    def test_invalid_weights_through_public_method_raises(self, GS6_CP):
        P = GS6_CP
        with pytest.raises(SSException):
            P.get_radius_of_gyration(weights=np.full(P.n_frames + 1, 1.0 / (P.n_frames + 1)))
        with pytest.raises(SSException):
            P.get_asphericity(verbose=False, weights=np.full(P.n_frames, 0.5 / P.n_frames))


# --------------------------------------------------------------------------
# collapsing functions: uniform weights == unweighted result
# --------------------------------------------------------------------------
class TestCollapsingWeighted:

    def test_internal_scaling_mean_vals_uniform(self, NTL9_CP):
        P = NTL9_CP
        n = P.n_frames
        uni = np.full(n, 1.0 / n)
        s0, m0 = P.get_internal_scaling(mean_vals=True, verbose=False)
        s1, m1 = P.get_internal_scaling(mean_vals=True, weights=uni, verbose=False)
        assert np.allclose(np.asarray(m0, float), np.asarray(m1, float), rtol=1e-5, atol=1e-5)

    def test_internal_scaling_RMS_uniform(self, NTL9_CP):
        P = NTL9_CP
        n = P.n_frames
        uni = np.full(n, 1.0 / n)
        _, r0 = P.get_internal_scaling_RMS(verbose=False)
        _, r1 = P.get_internal_scaling_RMS(weights=uni, verbose=False)
        assert np.allclose(np.asarray(r0, float), np.asarray(r1, float), rtol=1e-5, atol=1e-5)

    def test_scaling_exponent_uniform(self, NTL9_CP):
        P = NTL9_CP
        n = P.n_frames
        uni = np.full(n, 1.0 / n)
        np.random.seed(42)
        a = P.get_scaling_exponent(verbose=False)
        np.random.seed(42)
        b = P.get_scaling_exponent(weights=uni, verbose=False)
        # nu and A0 point estimates should agree under uniform weights
        assert np.allclose([a[0], a[1]], [b[0], b[1]], rtol=1e-4, atol=1e-4)

    def test_end_to_end_vs_rg_correlation_uniform(self, NTL9_CP):
        P = NTL9_CP
        n = P.n_frames
        uni = np.full(n, 1.0 / n)
        assert np.isclose(P.get_end_to_end_vs_rg_correlation(weights=uni),
                          P.get_end_to_end_vs_rg_correlation(), rtol=1e-5, atol=1e-6)

    def test_local_collapse_uniform(self, NTL9_CP):
        P = NTL9_CP
        n = P.n_frames
        uni = np.full(n, 1.0 / n)
        m0, s0, _, _ = P.get_local_collapse(window_size=5, verbose=False)
        m1, s1, _, _ = P.get_local_collapse(window_size=5, verbose=False, weights=uni)
        assert np.allclose(np.asarray(m0), np.asarray(m1), rtol=1e-5, atol=1e-5)
        assert np.allclose(np.asarray(s0), np.asarray(s1), rtol=1e-5, atol=1e-5)

    def test_angle_decay_uniform(self, NTL9_CP):
        P = NTL9_CP
        n = P.n_frames
        uni = np.full(n, 1.0 / n)
        d0 = P.get_angle_decay()
        d1 = P.get_angle_decay(weights=uni)
        assert np.allclose(d0, d1, rtol=1e-5, atol=1e-5)

    def test_sasa_uniform(self, NTL9_CP):
        P = NTL9_CP
        n = P.n_frames
        uni = np.full(n, 1.0 / n)
        base = np.asarray(P.get_all_SASA(mode='residue', stride=1))   # (n_frames, n_res)
        w = P.get_all_SASA(mode='residue', stride=1, weights=uni)     # (n_res,)
        assert np.allclose(w, base.mean(axis=0), rtol=1e-5, atol=1e-5)

    def test_interchain_distance_map_uniform(self, GMX_2CHAINS):
        T = GMX_2CHAINS
        n = T.proteinTrajectoryList[0].n_frames
        uni = np.full(n, 1.0 / n)
        d0, s0 = T.get_interchain_distance_map(0, 1, mode='CA')
        d1, s1 = T.get_interchain_distance_map(0, 1, mode='CA', weights=uni)
        assert np.allclose(d0, d1, rtol=1e-5, atol=1e-5)
        assert np.allclose(s0, s1, rtol=1e-5, atol=1e-5)


# --------------------------------------------------------------------------
# deterministic-everywhere: intentionally-unsupported combinations raise
# --------------------------------------------------------------------------
class TestDeterministicEdgeCases:

    def test_internal_scaling_mean_vals_false_with_weights_raises(self, GS6_CP):
        P = GS6_CP
        uni = np.full(P.n_frames, 1.0 / P.n_frames)
        with pytest.raises(SSException):
            P.get_internal_scaling(mean_vals=False, weights=uni, verbose=False)

    def test_local_heterogeneity_with_weights_raises(self, GS6_CP):
        P = GS6_CP
        uni = np.full(P.n_frames, 1.0 / P.n_frames)
        with pytest.raises(SSException):
            P.get_local_heterogeneity(fragment_size=3, stride=1, verbose=False, weights=uni)


def _uniform(n):
    return np.full(n, 1.0 / n)


def _one_hot(n, k):
    w = np.zeros(n)
    w[k] = 1.0
    return w


# --------------------------------------------------------------------------
# get_local_to_global_correlation: this is the function whose latent
# "weights silently discarded" bug was fixed - so we explicitly assert
# (a) uniform weights reproduce the unweighted result and (b) non-uniform
# weights actually CHANGE the result (i.e. the weights are now used).
# --------------------------------------------------------------------------
class TestLocalToGlobalCorrelationWeights:

    @staticmethod
    def _l2g(P, **kw):
        # seed so the stochastic pair selection is identical across calls
        np.random.seed(0)
        return P.get_local_to_global_correlation(
            n_cycles=15, max_num_pairs=4, stride=20, verbose=False, **kw
        )

    def test_uniform_matches_unweighted(self, CTL9_CP):
        P = CTL9_CP
        uni = _uniform(P.n_frames)
        base = self._l2g(P)
        w = self._l2g(P, weights=uni)
        assert np.allclose(base[2], w[2], rtol=1e-5, atol=1e-6)   # mean_corr
        assert np.allclose(base[3], w[3], rtol=1e-5, atol=1e-6)   # std_corr

    def test_nonuniform_changes_result(self, CTL9_CP):
        # regression guard for the fixed bug: a non-uniform weight MUST
        # change the correlation (previously the weights were discarded
        # and overwritten with a uniform vector).
        P = CTL9_CP
        n = P.n_frames
        w = np.arange(1, n + 1, dtype=float)
        w = w / w.sum()
        base = self._l2g(P)
        wt = self._l2g(P, weights=w)
        assert not np.allclose(base[2], wt[2], rtol=1e-4, atol=1e-4)


# --------------------------------------------------------------------------
# SASA summaries + the method-level stride+weights renormalisation path
# --------------------------------------------------------------------------
class TestSASASummariesWeights:

    def test_site_accessibility_uniform(self, NTL9_CP):
        P = NTL9_CP
        uni = _uniform(P.n_frames)
        base = P.get_site_accessibility([1, 2, 3], mode='resid', stride=1)
        wtd = P.get_site_accessibility([1, 2, 3], mode='resid', stride=1, weights=uni)
        for key in base:
            assert np.allclose(base[key], wtd[key], rtol=1e-5, atol=1e-6), key

    def test_site_accessibility_one_hot(self, GS6_CP):
        P = GS6_CP
        n = P.n_frames
        per_res = np.transpose(P.get_all_SASA(mode='residue', stride=1))  # (n_res, n_frames)
        lookup = P.get_amino_acid_sequence()
        for k in range(n):
            wtd = P.get_site_accessibility([1, 2], mode='resid', stride=1,
                                           weights=_one_hot(n, k))
            for i in (1, 2):
                mean_k, std_k = wtd[lookup[i]]
                assert np.isclose(mean_k, per_res[i][k], rtol=1e-9, atol=1e-9)
                assert np.isclose(std_k, 0.0, atol=1e-9)

    def test_regional_SASA_uniform_and_one_hot(self, GS6_CP):
        P = GS6_CP
        n = P.n_frames
        total = np.transpose(P.get_all_SASA(stride=1))      # (n_res, n_frames)
        base = P.get_regional_SASA(1, 5, stride=1)
        assert np.isclose(base, P.get_regional_SASA(1, 5, stride=1, weights=_uniform(n)),
                          rtol=1e-5, atol=1e-6)
        for k in range(n):
            expected = sum(float(total[i][k]) for i in range(1, 5))
            got = P.get_regional_SASA(1, 5, stride=1, weights=_one_hot(n, k))
            # SASA arrays are float32 (mdtraj); allow float32 round-off
            assert np.isclose(got, expected, rtol=1e-5, atol=1e-3)

    def test_all_SASA_mode_all_uniform_and_one_hot(self, GS6_CP):
        P = GS6_CP
        n = P.n_frames
        base = P.get_all_SASA(mode='all', stride=1)          # 3-tuple, each (n_frames, n_res)
        wtd_uni = P.get_all_SASA(mode='all', stride=1, weights=_uniform(n))
        assert isinstance(wtd_uni, tuple) and len(wtd_uni) == 3
        for b, w in zip(base, wtd_uni):
            assert np.allclose(w, np.asarray(b).mean(axis=0), rtol=1e-5, atol=1e-6)
        for k in range(n):
            wtd_k = P.get_all_SASA(mode='all', stride=1, weights=_one_hot(n, k))
            for b, w in zip(base, wtd_k):
                assert np.allclose(w, np.asarray(b)[k], rtol=1e-9, atol=1e-9)

    def test_all_SASA_stride_plus_weights_renormalises(self, NTL9_CP):
        # end-to-end check of the stride -> subsample -> renormalise path
        # through a real public method (not just the validator).
        P = NTL9_CP
        n = P.n_frames
        stride = 2
        full = np.arange(1, n + 1, dtype=float)
        full = full / full.sum()                       # valid full-length vector
        per_frame = np.asarray(P.get_all_SASA(mode='residue', stride=stride))  # (n_str, n_res)
        sub = full[::stride]
        sub = sub / sub.sum()
        expected = np.average(per_frame, axis=0, weights=sub)
        got = P.get_all_SASA(mode='residue', stride=stride, weights=full)
        assert np.allclose(got, expected, rtol=1e-6, atol=1e-8)


# --------------------------------------------------------------------------
# SSTrajectory.get_overall_* delegating getters
# --------------------------------------------------------------------------
class TestOverallTrajectoryWeights:

    def _getters(self, T):
        return {
            'rg':   T.get_overall_radius_of_gyration,
            'asph': T.get_overall_asphericity,
            'Rh':   T.get_overall_hydrodynamic_radius,
        }

    def test_uniform_equals_unweighted_mean(self, NTL9_CO):
        T = NTL9_CO
        n = T.proteinTrajectoryList[0].n_frames
        uni = _uniform(n)
        for name, fn in self._getters(T).items():
            base = np.asarray(fn())
            assert np.allclose(fn(weights=uni), base.mean(), rtol=1e-5, atol=1e-6), name

    def test_one_hot_equals_single_frame(self, GS6_CO):
        T = GS6_CO
        n = T.proteinTrajectoryList[0].n_frames
        for name, fn in self._getters(T).items():
            base = np.asarray(fn())
            for k in range(n):
                assert np.isclose(fn(weights=_one_hot(n, k)), base[k],
                                  rtol=1e-9, atol=1e-9), (name, k)


# --------------------------------------------------------------------------
# Pre-existing deterministically-weighted methods: confirm the routed
# validation did not change behaviour and uniform == unweighted.
# --------------------------------------------------------------------------
class TestPreexistingDeterministicUniformInvariance:

    def test_distance_map_uniform(self, NTL9_CP):
        P = NTL9_CP
        uni = _uniform(P.n_frames)
        for mode in ('CA', 'COM'):
            base = P.get_distance_map(mode=mode, verbose=False)[0]
            wtd = P.get_distance_map(mode=mode, weights=uni, verbose=False)[0]
            assert np.allclose(base, wtd, rtol=1e-5, atol=1e-5), mode

    def test_contact_map_uniform(self, NTL9_CP):
        P = NTL9_CP
        uni = _uniform(P.n_frames)
        base = P.get_contact_map(mode='ca')[0]
        wtd = P.get_contact_map(mode='ca', weights=uni)[0]
        assert np.allclose(base, wtd, rtol=1e-5, atol=1e-5)

    def test_Q_protein_average_with_weights_raises(self, NTL9_CP):
        # get_Q(protein_average=True) intentionally refuses weights:
        # frame-averaged Q must be reweighted outside SOURSOP. This is a
        # documented guard - lock it in.
        P = NTL9_CP
        uni = _uniform(P.n_frames)
        with pytest.raises(SSException):
            P.get_Q(protein_average=True, weights=uni)

    def test_dihedral_mutual_information_uniform(self, NTL9_CP):
        P = NTL9_CP
        uni = _uniform(P.n_frames)
        base = P.get_dihedral_mutual_information()
        wtd = P.get_dihedral_mutual_information(weights=uni)
        assert np.allclose(base, wtd, rtol=1e-4, atol=1e-4, equal_nan=True)
