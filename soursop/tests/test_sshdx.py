"""
Tests for soursop.sshdx - HDX protection factors (Best-Vendruscolo).

Covers:
  * shape, range and integer dtype of the per-residue per-frame ``N_c``
    and ``N_h`` arrays;
  * exclusion of immediate sequence neighbours (|i - j| <= 2 by default);
  * consistency between ``compute_Nc`` / ``compute_Nh`` residue lists;
  * the Best-Vendruscolo combination ``lnP = beta_c*Nc + beta_h*Nh +
    beta_0`` (and ``beta_0`` shifts ``lnP`` exactly);
  * the package-wide ``weights`` / ``stride`` contract (uniform = mean,
    one-hot = single frame, invalid weights raise ``SSException``).

The folded fixture (CTL9) exercises the H-bond path; GS6 (disordered)
is used for the input-validation and stride tests because it is fast.
"""

import numpy as np
import pytest

from soursop import sshdx
from soursop.ssexceptions import SSException


# --------------------------------------------------------------------------
# compute_Nc
# --------------------------------------------------------------------------
class TestComputeNc:
    def test_shape_and_range(self, GS6_CP):
        res, Nc = sshdx.compute_Nc(GS6_CP)
        assert Nc.dtype.kind == "i"
        assert Nc.shape == (GS6_CP.n_frames, len(res))
        assert Nc.min() >= 0
        # Defaults: 6.5 A cutoff, exclude |i-j|<=2. On GS6 (a hexapeptide
        # of residues 0..7) this is a small number per residue.
        assert Nc.max() < 1000  # absurdly loose, just a finiteness check

    def test_returns_residue_indices(self, GS6_CP):
        res, _ = sshdx.compute_Nc(GS6_CP)
        # Should be a strictly increasing sequence of integer residue indices.
        assert np.all(np.diff(res) > 0)

    def test_stride(self, GS6_CP):
        _, Nc = sshdx.compute_Nc(GS6_CP)
        _, Nc2 = sshdx.compute_Nc(GS6_CP, stride=2)
        assert Nc2.shape[1] == Nc.shape[1]
        assert np.array_equal(Nc2, Nc[::2])

    def test_exclude_neighbours_monotonic(self, GS6_CP):
        """Larger exclusion radius -> fewer or equal contacts."""
        _, Nc_default = sshdx.compute_Nc(GS6_CP, exclude_neighbours=2)
        _, Nc_strict = sshdx.compute_Nc(GS6_CP, exclude_neighbours=4)
        # Per-residue per-frame: stricter exclusion can only remove
        # contacts, not add them.
        assert np.all(Nc_strict <= Nc_default)

    def test_cutoff_monotonic(self, GS6_CP):
        """Smaller contact cutoff -> fewer or equal contacts."""
        _, Nc_default = sshdx.compute_Nc(GS6_CP, contact_cutoff=0.65)
        _, Nc_tight = sshdx.compute_Nc(GS6_CP, contact_cutoff=0.45)
        assert np.all(Nc_tight <= Nc_default)


# --------------------------------------------------------------------------
# compute_Nh
# --------------------------------------------------------------------------
class TestComputeNh:
    def test_shape_and_range_on_gs6(self, GS6_CP):
        """GS6 is disordered so backbone H-bonds should be sparse."""
        res, Nh = sshdx.compute_Nh(GS6_CP)
        assert Nh.dtype.kind == "i"
        assert Nh.shape == (GS6_CP.n_frames, len(res))
        # Backbone amide H can donate at most one H-bond at a time.
        assert Nh.min() >= 0
        assert Nh.max() <= 1

    def test_nh_detects_folded_hbonds(self, CTL9_CP):
        """A folded protein (CTL9) MUST contain backbone H-bonds."""
        res, Nh = sshdx.compute_Nh(CTL9_CP)
        assert Nh.shape == (CTL9_CP.n_frames, len(res))
        assert int(Nh.sum()) > 0, "no backbone H-bonds in folded CTL9 - sanity check"


# --------------------------------------------------------------------------
# compute_protection_factors
# --------------------------------------------------------------------------
class TestProtectionFactors:
    def test_shape_and_formula(self, GS6_CP):
        """lnP must equal beta_c*Nc + beta_h*Nh + beta_0 exactly."""
        res, lnP = sshdx.compute_protection_factors(GS6_CP)
        res_nc, Nc = sshdx.compute_Nc(GS6_CP)
        res_nh, Nh = sshdx.compute_Nh(GS6_CP)
        assert np.array_equal(res, res_nc)
        assert np.array_equal(res, res_nh)
        assert lnP.shape == (GS6_CP.n_frames, len(res))
        expected = (
            sshdx.DEFAULT_BETA_C * Nc + sshdx.DEFAULT_BETA_H * Nh + sshdx.DEFAULT_BETA_0
        )
        assert np.allclose(lnP, expected)

    def test_beta_0_shifts_lnP(self, GS6_CP):
        _, lnP = sshdx.compute_protection_factors(GS6_CP, beta_0=0.0)
        _, lnP_shifted = sshdx.compute_protection_factors(GS6_CP, beta_0=1.5)
        assert np.allclose(lnP_shifted, lnP + 1.5)

    def test_uniform_weights_match_unweighted_mean(self, GS6_CP):
        _, lnP = sshdx.compute_protection_factors(GS6_CP)
        n = GS6_CP.n_frames
        _, lnP_mean = sshdx.compute_protection_factors(
            GS6_CP, weights=np.full(n, 1.0 / n)
        )
        assert lnP_mean.shape == (lnP.shape[1],)
        assert np.allclose(lnP_mean, lnP.mean(axis=0))

    def test_one_hot_weights_pick_a_single_frame(self, GS6_CP):
        _, lnP = sshdx.compute_protection_factors(GS6_CP)
        n = GS6_CP.n_frames
        w = np.zeros(n)
        w[2] = 1.0
        _, lnP_frame = sshdx.compute_protection_factors(GS6_CP, weights=w)
        assert np.allclose(lnP_frame, lnP[2])

    def test_invalid_weights_raise(self, GS6_CP):
        n = GS6_CP.n_frames
        bad = np.full(n, 2.0 / n)  # doesn't sum to 1
        with pytest.raises(SSException):
            sshdx.compute_protection_factors(GS6_CP, weights=bad)

    def test_ctl9_lnP_range_is_physical(self, CTL9_CP):
        """Folded CTL9 must show a wide range of lnP between exposed and buried residues."""
        _, lnP = sshdx.compute_protection_factors(CTL9_CP)
        assert lnP.min() == pytest.approx(0.0, abs=1e-12)  # fully exposed
        assert lnP.max() > 5.0  # something is well-protected
