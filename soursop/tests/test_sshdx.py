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
def _bruteforce_distance_nh(protein, cutoff=0.24, exclude_neighbours=None):
    """Best-Vendruscolo N_h straight from the coordinates: protein oxygens
    within ``cutoff`` (nm) of each amide H."""
    top = protein.traj.topology
    res_NH, _, h_idx = sshdx._backbone_nh_map(top)
    oxygens = [a for a in top.atoms if a.element.symbol == "O"]
    o_idx = np.array([a.index for a in oxygens])
    o_res = np.array([a.residue.index for a in oxygens])
    xyz = protein.traj.xyz
    Nh = np.zeros((protein.n_frames, len(res_NH)), dtype=int)
    for k, (r, h) in enumerate(zip(res_NH, h_idx)):
        keep = np.ones(len(o_idx), dtype=bool)
        if exclude_neighbours is not None:
            keep = np.abs(o_res - r) > exclude_neighbours
        d = np.linalg.norm(xyz[:, o_idx[keep], :] - xyz[:, [h], :], axis=2)
        Nh[:, k] = (d < cutoff).sum(axis=1)
    return res_NH, Nh


class TestComputeNh:
    def test_shape_and_range_on_gs6(self, GS6_CP):
        """GS6 is disordered so H-bonds should be sparse."""
        res, Nh = sshdx.compute_Nh(GS6_CP)
        assert Nh.dtype.kind == "i"
        assert Nh.shape == (GS6_CP.n_frames, len(res))
        assert Nh.min() >= 0
        # a single amide H rarely sits within 2.4 A of more than two oxygens
        assert Nh.max() <= 3

    def test_default_is_best_vendruscolo_distance_count(self, NTL9_CP, GS6_CP):
        """The default N_h is the HDXer/Best-Vendruscolo definition: every
        protein oxygen within 2.4 A of the amide H, no sequence exclusion."""
        for protein in (NTL9_CP, GS6_CP):
            res, Nh = sshdx.compute_Nh(protein)
            res_bf, Nh_bf = _bruteforce_distance_nh(protein)
            assert np.array_equal(res, res_bf)
            assert np.array_equal(Nh, Nh_bf)
        # on the folded NTL9 fixture this definition finds many more H-bonds
        # than the backbone-only geometric one (119 vs 21 over 10 frames)
        _, Nh_dist = sshdx.compute_Nh(NTL9_CP)
        _, Nh_wn = sshdx.compute_Nh(NTL9_CP, hbond_method="wernet-nilsson")
        assert Nh_dist.sum() == 119
        assert Nh_wn.sum() == 21

    def test_distance_method_honours_cutoff_and_exclusion(self, NTL9_CP):
        _, Nh = sshdx.compute_Nh(NTL9_CP)
        _, Nh_tight = sshdx.compute_Nh(NTL9_CP, hbond_cutoff=0.20)
        assert np.all(Nh_tight <= Nh)
        _, Nh_excl = sshdx.compute_Nh(NTL9_CP, exclude_neighbours=2)
        assert np.all(Nh_excl <= Nh)
        assert np.array_equal(
            Nh_excl, _bruteforce_distance_nh(NTL9_CP, exclude_neighbours=2)[1]
        )
        with pytest.raises(SSException):
            sshdx.compute_Nh(NTL9_CP, hbond_cutoff=0.0)
        with pytest.raises(SSException):
            sshdx.compute_Nh(NTL9_CP, hbond_method="baker-hubbard")

    def test_wernet_nilsson_method_is_backbone_only(self, CTL9_CP):
        """The legacy definition: backbone donor H to backbone carbonyl O,
        |i - j| > 2, one H-bond per amide at most."""
        res, Nh = sshdx.compute_Nh(CTL9_CP, hbond_method="wernet-nilsson")
        assert Nh.shape == (CTL9_CP.n_frames, len(res))
        assert Nh.max() <= 1
        assert int(Nh.sum()) > 0, "no backbone H-bonds in folded CTL9 - sanity check"

        # independent reconstruction from mdtraj's H-bond list
        import mdtraj as md

        top = CTL9_CP.traj.topology
        res_NH, _, h_idx = sshdx._backbone_nh_map(top)
        h_to_res = dict(zip(h_idx.tolist(), res_NH.tolist()))
        o_to_res = {a.index: a.residue.index for a in top.atoms if a.name == "O"}
        expected = np.zeros_like(Nh)
        for f, hb in enumerate(md.wernet_nilsson(CTL9_CP.traj, periodic=False)):
            for _, h, a in hb:
                i, j = h_to_res.get(int(h)), o_to_res.get(int(a))
                if i is not None and j is not None and abs(i - j) > 2:
                    expected[f, list(res_NH).index(i)] += 1
        assert np.array_equal(Nh, expected)

    def test_hydrogen_free_topology_raises(self, SIGA_CG_CO):
        cg = SIGA_CG_CO.proteinTrajectoryList[0]
        with pytest.raises(SSException, match="no residues with a backbone amide"):
            sshdx.compute_Nh(cg)
        with pytest.raises(SSException, match="no residues with a backbone amide"):
            sshdx.compute_Nc(cg)


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

    def test_hbond_kwargs_forwarded(self, NTL9_CP):
        _, lnP_default = sshdx.compute_protection_factors(NTL9_CP)
        _, lnP_wn = sshdx.compute_protection_factors(
            NTL9_CP, hbond_method="wernet-nilsson"
        )
        _, Nc = sshdx.compute_Nc(NTL9_CP)
        _, Nh_wn = sshdx.compute_Nh(NTL9_CP, hbond_method="wernet-nilsson")
        assert np.allclose(lnP_wn, 0.35 * Nc + 2.0 * Nh_wn)
        assert not np.allclose(lnP_default, lnP_wn)
        # hbond_exclude_neighbours reaches compute_Nh
        _, lnP_excl = sshdx.compute_protection_factors(
            NTL9_CP, hbond_exclude_neighbours=2
        )
        _, Nh_excl = sshdx.compute_Nh(NTL9_CP, exclude_neighbours=2)
        assert np.allclose(lnP_excl, 0.35 * Nc + 2.0 * Nh_excl)

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
