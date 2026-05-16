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


# --------------------------------------------------------------------------
# get_Q: the weighted (protein_average=False) breakdown must accept a
# full-length (len == n_frames) weight vector AND work together with a
# frame stride (previously stride != 1 + weights was hard-blocked).
# --------------------------------------------------------------------------
class TestQWeightsAndStride:

    def _shape0(self, P, **kw):
        out = P.get_Q(protein_average=False, **kw)
        return np.asarray(out[0], dtype=float), np.asarray(out[1])

    def test_full_length_weights_stride1(self, NTL9_CP):
        P = NTL9_CP
        w = _uniform(P.n_frames)                       # one weight per frame
        frac, pairs = self._shape0(P, weights=w)
        base_frac, base_pairs = self._shape0(P)        # unweighted breakdown
        assert frac.shape == base_frac.shape == (len(base_pairs),)
        assert np.all(np.isfinite(frac))
        assert frac.min() >= -1e-9 and frac.max() <= 1.0 + 1e-9

    def test_stride_plus_weights_now_supported(self, NTL9_CP):
        # the headline fix: stride != 1 together with a full-length
        # weights vector used to raise; it must now run and be well-formed.
        P = NTL9_CP
        w = _uniform(P.n_frames)
        frac, pairs = self._shape0(P, stride=2, weights=w)
        assert frac.shape == (len(pairs),)
        assert np.all(np.isfinite(frac))
        assert frac.min() >= -1e-9 and frac.max() <= 1.0 + 1e-9

    def test_weights_actually_used(self, NTL9_CP):
        # a non-uniform weight must change the per-contact breakdown
        # (guards against the weights being silently ignored), and the
        # result must be deterministic.
        P = NTL9_CP
        n = P.n_frames
        uni = _uniform(n)
        nonuni = np.arange(1, n + 1, dtype=float)
        nonuni = nonuni / nonuni.sum()
        f_uni, _ = self._shape0(P, weights=uni)
        f_non, _ = self._shape0(P, weights=nonuni)
        f_non2, _ = self._shape0(P, weights=nonuni)
        # get_Q is not bit-reproducible (mdtraj superpose alignment is
        # stable only to ~1e-6 relative); the key assertion is that a
        # non-uniform weight visibly changes the per-contact breakdown.
        assert np.allclose(f_non, f_non2, rtol=1e-4, atol=1e-6)
        assert not np.allclose(f_uni, f_non, rtol=1e-3, atol=1e-3)

    def test_invalid_weights_raise(self, NTL9_CP):
        P = NTL9_CP
        n = P.n_frames
        with pytest.raises(SSException):                # wrong length
            P.get_Q(protein_average=False, weights=_uniform(n - 1))
        with pytest.raises(SSException):                # wrong length + stride
            P.get_Q(protein_average=False, stride=2, weights=_uniform(n - 1))
        with pytest.raises(SSException):                # sum != 1
            P.get_Q(protein_average=False, weights=np.full(n, 0.5 / n))


# --------------------------------------------------------------------------
# get_contact_map: weights + stride != 1 must work together (the explicit
# "stride must be 1 if weights given" hard-block was removed).
# --------------------------------------------------------------------------
class TestContactMapWeightsAndStride:

    def test_full_length_weights_stride1(self, NTL9_CP):
        P = NTL9_CP
        uni = _uniform(P.n_frames)
        base = P.get_contact_map(mode='ca')[0]
        wtd = P.get_contact_map(mode='ca', weights=uni)[0]
        assert wtd.shape == base.shape
        assert np.allclose(base, wtd, rtol=1e-5, atol=1e-5)
        assert wtd.min() >= -1e-9 and wtd.max() <= 1.0 + 1e-9

    def test_stride_plus_weights_matches_unweighted_stride(self, NTL9_CP):
        # the headline fix: stride != 1 + a full-length weight vector used
        # to raise. With uniform weights it must now run AND reproduce the
        # unweighted strided contact map (uniform weights over the strided
        # frames == the plain strided mean).
        P = NTL9_CP
        uni = _uniform(P.n_frames)
        base = P.get_contact_map(mode='ca', stride=2)[0]
        wtd = P.get_contact_map(mode='ca', stride=2, weights=uni)[0]
        assert wtd.shape == base.shape
        assert np.allclose(base, wtd, rtol=1e-5, atol=1e-5)

    def test_weights_actually_used(self, NTL9_CP):
        P = NTL9_CP
        n = P.n_frames
        uni = _uniform(n)
        nonuni = np.arange(1, n + 1, dtype=float)
        nonuni = nonuni / nonuni.sum()
        c_uni = P.get_contact_map(mode='ca', weights=uni)[0]
        c_non = P.get_contact_map(mode='ca', weights=nonuni)[0]
        c_non2 = P.get_contact_map(mode='ca', weights=nonuni)[0]
        assert np.array_equal(c_non, c_non2)               # deterministic
        assert not np.allclose(c_uni, c_non, rtol=1e-6, atol=1e-6)

    def test_invalid_weights_raise(self, NTL9_CP):
        P = NTL9_CP
        n = P.n_frames
        with pytest.raises(SSException):                   # wrong length
            P.get_contact_map(mode='ca', weights=_uniform(n - 1))
        with pytest.raises(SSException):                   # wrong length + stride
            P.get_contact_map(mode='ca', stride=2, weights=_uniform(n - 1))
        with pytest.raises(SSException):                   # sum != 1
            P.get_contact_map(mode='ca', weights=np.full(n, 0.5 / n))


# --------------------------------------------------------------------------
# Additional opt-in reweighted getters: the two inter-residue distance
# functions, get_t, and DSSP / BBSEG (return_per_frame=False only).
# --------------------------------------------------------------------------
class TestExtraReweightedGetters:

    def test_inter_residue_COM_distance(self, NTL9_CP):
        P = NTL9_CP
        n = P.n_frames
        base = np.asarray(P.get_inter_residue_COM_distance(1, 20))
        assert np.isclose(P.get_inter_residue_COM_distance(1, 20, weights=_uniform(n)),
                          base.mean(), rtol=1e-6, atol=1e-9)

    def test_inter_residue_COM_distance_one_hot(self, GS6_CP):
        P = GS6_CP
        n = P.n_frames
        base = np.asarray(P.get_inter_residue_COM_distance(1, 5))
        for k in range(n):
            assert np.isclose(P.get_inter_residue_COM_distance(1, 5, weights=_one_hot(n, k)),
                              base[k], rtol=1e-9, atol=1e-9)

    def test_inter_residue_atomic_distance(self, NTL9_CP):
        P = NTL9_CP
        n = P.n_frames
        for mode in ('atom', 'ca', 'closest', 'closest-heavy'):
            base = np.asarray(P.get_inter_residue_atomic_distance(1, 20, mode=mode))
            assert np.isclose(P.get_inter_residue_atomic_distance(1, 20, mode=mode, weights=_uniform(n)),
                              base.mean(), rtol=1e-5, atol=1e-6), mode

    def test_inter_residue_atomic_distance_one_hot(self, GS6_CP):
        P = GS6_CP
        n = P.n_frames
        base = np.asarray(P.get_inter_residue_atomic_distance(1, 5))
        for k in range(n):
            assert np.isclose(P.get_inter_residue_atomic_distance(1, 5, weights=_one_hot(n, k)),
                              base[k], rtol=1e-9, atol=1e-9)

    def test_inter_residue_distances_stride_plus_weights(self, NTL9_CP):
        # full-length weights + stride must work and equal the unweighted
        # strided mean (uniform over strided frames == plain strided mean)
        P = NTL9_CP
        n = P.n_frames
        uni = _uniform(n)
        for fn in (lambda **k: P.get_inter_residue_COM_distance(1, 20, stride=2, **k),
                   lambda **k: P.get_inter_residue_atomic_distance(1, 20, stride=2, **k)):
            base = np.asarray(fn()).mean()
            assert np.isclose(fn(weights=uni), base, rtol=1e-5, atol=1e-6)

    def test_get_t(self, NTL9_CP):
        P = NTL9_CP
        n = P.n_frames
        base = np.asarray(P.get_t())
        assert np.isclose(P.get_t(weights=_uniform(n)), base.mean(), rtol=1e-6, atol=1e-9)

    def test_get_t_one_hot(self, GS6_CP):
        P = GS6_CP
        n = P.n_frames
        base = np.asarray(P.get_t())
        for k in range(n):
            assert np.isclose(P.get_t(weights=_one_hot(n, k)), base[k], rtol=1e-9, atol=1e-9)

    def test_dssp_uniform_and_onehot(self, NTL9_CP):
        P = NTL9_CP
        n = P.n_frames
        r0, h0, e0, c0 = P.get_secondary_structure_DSSP()
        r1, h1, e1, c1 = P.get_secondary_structure_DSSP(weights=_uniform(n))
        assert r0 == r1
        assert np.allclose(h0, h1, rtol=1e-5, atol=1e-6)
        assert np.allclose(e0, e1, rtol=1e-5, atol=1e-6)
        assert np.allclose(c0, c1, rtol=1e-5, atol=1e-6)
        # one-hot -> the fractions are exactly that frame's 0/1 masks
        _, hk, ek, ck = P.get_secondary_structure_DSSP(weights=_one_hot(n, 0))
        for arr in (hk, ek, ck):
            assert set(np.unique(arr)).issubset({0.0, 1.0})

    def test_dssp_per_frame_with_weights_raises(self, GS6_CP):
        P = GS6_CP
        with pytest.raises(SSException):
            P.get_secondary_structure_DSSP(return_per_frame=True,
                                           weights=_uniform(P.n_frames))

    def test_bbseg_uniform(self, NTL9_CP):
        P = NTL9_CP
        n = P.n_frames
        r0, cls0 = P.get_secondary_structure_BBSEG()
        r1, cls1 = P.get_secondary_structure_BBSEG(weights=_uniform(n))
        assert r0 == r1
        for c in range(9):
            assert np.allclose(cls0[c], cls1[c], rtol=1e-5, atol=1e-6), c

    def test_bbseg_per_frame_with_weights_raises(self, GS6_CP):
        P = GS6_CP
        with pytest.raises(SSException):
            P.get_secondary_structure_BBSEG(return_per_frame=True,
                                            weights=_uniform(P.n_frames))

    def test_invalid_weights_raise(self, GS6_CP):
        P = GS6_CP
        n = P.n_frames
        with pytest.raises(SSException):
            P.get_inter_residue_COM_distance(1, 5, weights=_uniform(n - 1))
        with pytest.raises(SSException):
            P.get_inter_residue_atomic_distance(1, 5, weights=np.full(n, 0.5 / n))
        with pytest.raises(SSException):
            P.get_t(weights=np.full(n, 2.0 / n))


# --------------------------------------------------------------------------
# Remaining weights-bearing entry points: get_polymer_scaled_distance_map
# and the low-level ssmutualinformation.calc_MI (the function whose
# weighted path was bug-fixed).
# --------------------------------------------------------------------------
class TestPolymerScaledDistanceMapWeights:

    def test_uniform_equals_unweighted_explicit_model(self, NTL9_CP):
        # with nu/A0 supplied the only weighted dependence is the internal
        # get_distance_map, so uniform weights must reproduce the
        # unweighted deviation matrix.
        P = NTL9_CP
        uni = _uniform(P.n_frames)
        m0, nu0, a0, _ = P.get_polymer_scaled_distance_map(nu=0.5, A0=5.5,
                                                           min_separation=5,
                                                           verbose=False)
        m1, nu1, a1, _ = P.get_polymer_scaled_distance_map(nu=0.5, A0=5.5,
                                                           min_separation=5,
                                                           weights=uni,
                                                           verbose=False)
        assert (nu0, a0) == (nu1, a1) == (0.5, 5.5)
        assert np.allclose(m0, m1, rtol=1e-5, atol=1e-5)

    def test_default_model_with_weights_runs(self, NTL9_CP):
        # default nu/A0 -> internally fits via get_scaling_exponent; with
        # weights it must still run and return a finite matrix.
        P = NTL9_CP
        uni = _uniform(P.n_frames)
        np.random.seed(42)
        m, nu, a0, _ = P.get_polymer_scaled_distance_map(min_separation=5,
                                                         weights=uni,
                                                         verbose=False)
        assert np.all(np.isfinite(np.asarray(m, dtype=float)))
        assert np.isfinite(nu) and np.isfinite(a0)

    def test_invalid_weights_raise(self, NTL9_CP):
        P = NTL9_CP
        n = P.n_frames
        with pytest.raises(SSException):
            P.get_polymer_scaled_distance_map(nu=0.5, A0=5.5,
                                              weights=_uniform(n - 1),
                                              verbose=False)
        with pytest.raises(SSException):
            P.get_polymer_scaled_distance_map(nu=0.5, A0=5.5,
                                              weights=np.full(n, 0.5 / n),
                                              verbose=False)


class TestCalcMIWeights:

    def test_uniform_weights_equal_unweighted(self):
        from soursop.ssmutualinformation import calc_MI
        rng = np.random.RandomState(0)
        X = rng.uniform(-1, 1, 200)
        Y = 0.7 * X + 0.3 * rng.uniform(-1, 1, 200)
        bins = np.linspace(-1.05, 1.05, 12)
        base = calc_MI(X, Y, bins)                       # unweighted
        w = np.full(len(X), 1.0 / len(X))
        wtd = calc_MI(X, Y, bins, weights=w)             # uniform weights
        assert np.isclose(base, wtd, rtol=1e-9, atol=1e-9)

    def test_array_weights_do_not_raise(self):
        # regression guard: `if weights:` used to raise ValueError on an
        # array ("truth value of an array ... is ambiguous").
        from soursop.ssmutualinformation import calc_MI
        rng = np.random.RandomState(1)
        X = rng.uniform(-1, 1, 100)
        Y = rng.uniform(-1, 1, 100)
        bins = np.linspace(-1.05, 1.05, 10)
        w = np.full(len(X), 1.0 / len(X))
        val = calc_MI(X, Y, bins, weights=w)
        assert np.isfinite(val)
        assert calc_MI(X, Y, bins, weights=False) is not None
