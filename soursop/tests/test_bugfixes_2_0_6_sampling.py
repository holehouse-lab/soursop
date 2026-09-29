"""
Regression tests for the sssampling fixes that went into 2.0.6.

Covers: finite relative entropy via pseudocount smoothing, compute_pdf
normalising with the bins actually passed, correct inverse-CDF resampling of
the EV reference, seeded (reproducible) EV references, strict bwidth
validation, quality_plot on uncapped chains, frame-count mismatch handling,
per-trajectory topology lists on the sequential loader, SSException for bad
selectors, and the non-standard residue error on the EV path.
"""

import os

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

import soursop
from soursop.ssexceptions import SSException
from soursop.sssampling import (
    PrecomputedDihedralInterface,
    SamplingQuality,
    rel_entropy,
    smooth_pdf,
)
from soursop.sstools import find_trajectory_files
from soursop.sstrajectory import SSTrajectory


def _sampling_paths():
    data_dir = soursop.get_data("test_data")
    wt_traj, wt_top = find_trajectory_files(
        os.path.join(data_dir, "sampling_quality/WT"), 3
    )
    ev_traj, ev_top = find_trajectory_files(
        os.path.join(data_dir, "sampling_quality/EV"), 3
    )
    return wt_traj, wt_top, ev_traj, ev_top


def _ntl9_paths():
    data_dir = soursop.get_data("test_data")
    return (
        os.path.join(data_dir, "ntl9_AA.xtc"),
        os.path.join(data_dir, "ntl9_AA.pdb"),
    )


@pytest.fixture(scope="module")
def wt_vs_ev_1d():
    wt_traj, wt_top, ev_traj, ev_top = _sampling_paths()
    return SamplingQuality(
        wt_traj,
        ev_traj,
        top_file=wt_top,
        ref_top=ev_top,
        method="1D angle distributions",
        force_sequential=True,
    )


@pytest.fixture(scope="module")
def wt_vs_ev_2d():
    wt_traj, wt_top, ev_traj, ev_top = _sampling_paths()
    return SamplingQuality(
        wt_traj,
        ev_traj,
        top_file=wt_top,
        ref_top=ev_top,
        method="2D angle distributions",
        force_sequential=True,
    )


@pytest.fixture(scope="module")
def wt_vs_precomputed():
    wt_traj, wt_top, _, _ = _sampling_paths()
    return SamplingQuality(
        wt_traj,
        top_file=wt_top,
        method="1D angle distributions",
        force_sequential=True,
        seed=0,
    )


# --------------------------------------------------------------------------
# 1. relative entropy is finite by default
# --------------------------------------------------------------------------
class TestRelativeEntropy:
    def test_smooth_pdf_rows_sum_to_one(self):
        pdf = np.array([[1.0, 0.0, 0.0], [0.5, 0.5, 0.0]])
        out = smooth_pdf(pdf, 1e-3)
        assert np.allclose(out.sum(axis=-1), 1.0)
        assert np.all(out > 0)
        assert np.array_equal(smooth_pdf(pdf, 0), pdf)
        with pytest.raises(SSException):
            smooth_pdf(pdf, -1.0)

    def test_rel_entropy_is_inf_on_disjoint_support(self):
        assert np.isinf(rel_entropy(np.array([1.0, 0.0]), np.array([0.0, 1.0])))

    def test_compute_dihedral_rel_entropy_finite(self, wt_vs_ev_1d):
        re_smoothed = wt_vs_ev_1d.compute_dihedral_rel_entropy()
        assert re_smoothed.shape == (2, 3, 2)
        assert np.all(np.isfinite(re_smoothed))
        assert np.all(re_smoothed >= 0)

        re_raw = wt_vs_ev_1d.compute_dihedral_rel_entropy(pseudocount=0)
        assert np.any(np.isinf(re_raw))

    def test_compute_dihedral_rel_entropy_finite_ev_path(self, wt_vs_precomputed):
        re_smoothed = wt_vs_precomputed.compute_dihedral_rel_entropy()
        assert np.all(np.isfinite(re_smoothed))
        assert np.all(re_smoothed >= 0)

    def test_all_to_all_relative_entropy_finite(self, wt_vs_ev_1d):
        phi_df, psi_df = wt_vs_ev_1d.get_all_to_all_trj_comparisons(
            metric="relative entropy"
        )
        assert np.all(np.isfinite(phi_df.values))
        assert np.all(np.isfinite(psi_df.values))
        assert np.all(phi_df.values >= 0)
        assert np.all(psi_df.values >= 0)

        phi_raw, psi_raw = wt_vs_ev_1d.get_all_to_all_trj_comparisons(
            metric="relative entropy", pseudocount=0
        )
        assert np.any(np.isinf(phi_raw.values)) or np.any(np.isinf(psi_raw.values))


# --------------------------------------------------------------------------
# 2. compute_pdf uses the bins it is given
# --------------------------------------------------------------------------
class TestComputePdfNormalisation:
    @pytest.mark.parametrize("width", [30.0, 15.0, 7.5])
    def test_rows_sum_to_one_for_custom_bins(self, wt_vs_ev_1d, width):
        bins = np.arange(-180.0, 180.0 + width, width)
        pdf_3d = wt_vs_ev_1d.compute_pdf(wt_vs_ev_1d.phi_angles, bins=bins)
        pdf_2d = wt_vs_ev_1d.compute_pdf(wt_vs_ev_1d.phi_angles[0], bins=bins)
        assert pdf_3d.shape[-1] == len(bins) - 1
        assert np.allclose(pdf_3d.sum(axis=-1), 1.0)
        assert np.allclose(pdf_2d.sum(axis=-1), 1.0)

    def test_default_bins_unchanged(self, wt_vs_ev_1d):
        # the default 15 degree grid must give the same numbers as before
        bins = wt_vs_ev_1d.get_degree_bins()
        expected = (
            np.apply_along_axis(
                lambda col: np.histogram(col, bins=bins, density=True)[0],
                axis=1,
                arr=wt_vs_ev_1d.psi_angles[0],
            )
            * 15.0
        )
        assert np.allclose(
            expected, wt_vs_ev_1d.compute_pdf(wt_vs_ev_1d.psi_angles[0], bins=bins)
        )


# --------------------------------------------------------------------------
# 3. inverse-CDF resampling reproduces the reference histogram
# --------------------------------------------------------------------------
class TestEVResampling:
    def test_resampled_histogram_matches_reference(self):
        bins = np.arange(-180.0, 181.0, 15.0)
        ev = PrecomputedDihedralInterface(
            "GGGGGGGG", bins=bins, num_trajs=1, nsamples=200000, seed=1
        )
        reference = ev.gather_phi_reference_dihedrals("GGGGGGGG")
        sampled = ev.ref_phi_angles[0]

        assert sampled.shape == (8, 200000)
        assert sampled.min() >= -180.0 and sampled.max() <= 180.0

        for i in range(reference.shape[0]):
            ref_hist = np.histogram(reference[i], bins=bins)[0] / reference.shape[1]
            smp_hist = np.histogram(sampled[i], bins=bins)[0] / sampled.shape[1]
            assert np.max(np.abs(ref_hist - smp_hist)) < 0.005
            # first and last bins explicitly, the old interp never reached them
            assert abs(ref_hist[0] - smp_hist[0]) < 0.005
            assert abs(ref_hist[-1] - smp_hist[-1]) < 0.005


# --------------------------------------------------------------------------
# 4. seeded EV reference is reproducible
# --------------------------------------------------------------------------
class TestSeed:
    def test_interface_seed(self):
        bins = np.arange(-180.0, 181.0, 15.0)
        a = PrecomputedDihedralInterface("AAAA", bins, 2, 100, seed=7)
        b = PrecomputedDihedralInterface("AAAA", bins, 2, 100, seed=7)
        c = PrecomputedDihedralInterface("AAAA", bins, 2, 100, seed=8)
        assert np.array_equal(a.ref_phi_angles, b.ref_phi_angles)
        assert np.array_equal(a.ref_psi_angles, b.ref_psi_angles)
        assert not np.array_equal(a.ref_phi_angles, c.ref_phi_angles)

    def test_sampling_quality_seed(self):
        wt_traj, wt_top, _, _ = _sampling_paths()
        kw = dict(top_file=wt_top, force_sequential=True)
        a = SamplingQuality(wt_traj, seed=3, **kw)
        b = SamplingQuality(wt_traj, seed=3, **kw)
        c = SamplingQuality(wt_traj, seed=4, **kw)
        assert np.array_equal(a.ref_phi_angles, b.ref_phi_angles)
        assert np.array_equal(
            a.compute_dihedral_hellingers(), b.compute_dihedral_hellingers()
        )
        assert not np.array_equal(a.ref_phi_angles, c.ref_phi_angles)


# --------------------------------------------------------------------------
# 5. bwidth validation and exact bin grid
# --------------------------------------------------------------------------
class TestBwidth:
    def test_seven_and_a_half_degrees(self):
        wt_traj, wt_top, _, _ = _sampling_paths()
        sq = SamplingQuality(
            wt_traj,
            top_file=wt_top,
            bwidth=np.deg2rad(7.5),
            force_sequential=True,
            seed=0,
        )
        assert len(sq.bins) == 49
        assert sq.bins[0] == -180.0
        assert sq.bins[-1] == 180.0
        assert np.allclose(np.diff(sq.bins), 7.5)
        h = sq.compute_dihedral_hellingers()
        assert np.all(np.isfinite(h))

    def test_small_divisor_no_longer_crashes(self):
        # 0.4 degrees divides 360 exactly (900 bins); the old code rounded the
        # width to 0 and blew up in np.arange.
        wt_traj, wt_top, _, _ = _sampling_paths()
        sq = SamplingQuality(
            wt_traj,
            top_file=wt_top,
            bwidth=np.deg2rad(0.4),
            method="1D angle distributions",
            force_sequential=True,
            seed=0,
        )
        assert len(sq.bins) == 901
        assert sq.bins[-1] == 180.0

    @pytest.mark.parametrize("deg", [7.0, 100.0, 0.7])
    def test_non_divisor_raises(self, deg):
        wt_traj, wt_top, _, _ = _sampling_paths()
        with pytest.raises(SSException, match="whole number of bins"):
            SamplingQuality(
                wt_traj,
                top_file=wt_top,
                bwidth=np.deg2rad(deg),
                force_sequential=True,
            )

    @pytest.mark.parametrize("bwidth", [0.0, -0.1, 7.0])
    def test_out_of_range_raises(self, bwidth):
        wt_traj, wt_top, _, _ = _sampling_paths()
        with pytest.raises(SSException, match="between 0 and 2\\*pi"):
            SamplingQuality(
                wt_traj, top_file=wt_top, bwidth=bwidth, force_sequential=True
            )


# --------------------------------------------------------------------------
# 6. quality_plot on uncapped chains and figname handling
# --------------------------------------------------------------------------
class TestQualityPlot:
    def test_uncapped_reference_plot(self, tmp_path):
        xtc, pdb = _ntl9_paths()
        sq = SamplingQuality(
            [xtc, xtc],
            [xtc, xtc],
            top_file=pdb,
            ref_top=pdb,
            force_sequential=True,
        )
        fig, axd = sq.quality_plot(save_dir=str(tmp_path))
        assert set(axd) == {"A", "B", "C", "D"}
        assert axd["A"].get_title() == "Comparison to reference ensemble"
        assert os.path.isfile(os.path.join(tmp_path, "2D_hellingers.pdf"))
        plt.close(fig)

    def test_uncapped_ev_plot_title(self):
        xtc, pdb = _ntl9_paths()
        sq = SamplingQuality([xtc, xtc], top_file=pdb, force_sequential=True, seed=0)
        fig, axd = sq.quality_plot()
        assert axd["A"].get_title() == "Comparison to the Excluded Volume Limit"
        plt.close(fig)


# --------------------------------------------------------------------------
# 7. unequal frame counts
# --------------------------------------------------------------------------
class TestFrameCountMismatch:
    def test_mismatch_raises_without_truncate(self, tmp_path):
        xtc, pdb = _ntl9_paths()
        short = str(tmp_path / "short.xtc")
        SSTrajectory(xtc, pdb).traj[0:5].save_xtc(short)

        with pytest.raises(SSException, match="truncate=True"):
            SamplingQuality([xtc, short], top_file=pdb, force_sequential=True)

        with pytest.raises(SSException, match="truncate=True"):
            SamplingQuality(
                [xtc],
                [short],
                top_file=pdb,
                ref_top=pdb,
                force_sequential=True,
            )

    def test_mismatch_works_with_truncate(self, tmp_path):
        xtc, pdb = _ntl9_paths()
        short = str(tmp_path / "short.xtc")
        SSTrajectory(xtc, pdb).traj[0:5].save_xtc(short)

        sq = SamplingQuality(
            [xtc, short],
            top_file=pdb,
            force_sequential=True,
            truncate=True,
            verbose=False,
            seed=0,
        )
        assert sq.phi_angles.shape[-1] == 5
        assert sq.ref_phi_angles.shape[-1] == 5


# --------------------------------------------------------------------------
# 9. bad metrics / selectors raise SSException; recompute is honoured
# --------------------------------------------------------------------------
class TestSelectors:
    def test_2d_comparison_bad_metric(self, wt_vs_ev_2d):
        with pytest.raises(SSException):
            wt_vs_ev_2d.get_all_to_all_2d_trj_comparison(metric="relative entropy")

    def test_2d_comparison_recompute(self, wt_vs_ev_2d):
        a = wt_vs_ev_2d.get_all_to_all_2d_trj_comparison()
        b = wt_vs_ev_2d.get_all_to_all_2d_trj_comparison(recompute=True)
        assert isinstance(a, np.ndarray)
        assert a.shape == (3, 2)
        assert np.allclose(a, b)

    def test_bad_selectors_raise_ssexception(self, wt_vs_ev_1d):
        with pytest.raises(SSException):
            wt_vs_ev_1d.trj_pdfs(dihedral="nope")
        with pytest.raises(SSException):
            wt_vs_ev_1d.ref_pdfs(dihedral="nope")
        with pytest.raises(SSException):
            wt_vs_ev_1d.get_all_to_all_trj_comparisons(metric="nope")


# --------------------------------------------------------------------------
# 10. per-trajectory topology lists on the sequential path
# --------------------------------------------------------------------------
class TestTopologyLists:
    def test_list_topologies_sequential(self):
        wt_traj, wt_top, ev_traj, ev_top = _sampling_paths()
        sq = SamplingQuality(
            wt_traj,
            ev_traj,
            top_file=list(wt_top),
            ref_top=list(ev_top),
            method="1D angle distributions",
            force_sequential=True,
        )
        assert sq.phi_angles.shape == (3, 2, 500)

    def test_single_trajectory_list_topology(self):
        wt_traj, wt_top, ev_traj, ev_top = _sampling_paths()
        sq = SamplingQuality(
            wt_traj[:1],
            ev_traj[:1],
            top_file=wt_top[:1],
            ref_top=ev_top[:1],
            force_sequential=True,
        )
        assert sq.phi_angles.shape == (1, 2, 500)

    def test_wrong_length_topology_list_raises(self):
        wt_traj, wt_top, _, _ = _sampling_paths()
        with pytest.raises(SSException, match="top_file"):
            SamplingQuality(wt_traj, top_file=wt_top[:2], force_sequential=True)


# --------------------------------------------------------------------------
# 11. non-standard residue on the EV path
# --------------------------------------------------------------------------
class TestNonStandardResidue:
    def test_named_in_error(self):
        bins = np.arange(-180.0, 181.0, 15.0)
        with pytest.raises(SSException, match="'X' at position 1"):
            PrecomputedDihedralInterface("AXA", bins, 1, 10)


# ---------------------------------------------------------------------------
# uncapped chains: phi and psi columns must refer to the same residue
# ---------------------------------------------------------------------------
def test_uncapped_chain_phi_psi_columns_are_residue_aligned():
    import soursop
    from soursop.sssampling import SamplingQuality

    data = soursop.get_data("test_data")
    xtc = os.path.join(data, "ntl9_AA.xtc")
    pdb = os.path.join(data, "ntl9_AA.pdb")
    sq = SamplingQuality(
        [xtc],
        top_file=pdb,
        method="2D angle distributions",
        seed=0,
        force_sequential=True,
    )
    protein = sq.trajs[0].proteinTrajectoryList[0]
    n = protein.n_residues
    # ntl9 is uncapped: residue 0 has no phi and residue n-1 has no psi, so
    # the aligned set is the interior residues 1..n-2
    assert sq.residue_indices == list(range(1, n - 1))
    assert sq.phi_angles.shape[1] == sq.psi_angles.shape[1] == n - 2
    assert sq.ref_phi_angles.shape[1] == sq.ref_psi_angles.shape[1] == n - 2

    # column k of phi and psi must be the dihedrals of residue k+1
    phi_atoms, phi = protein.get_angles("phi")
    psi_atoms, psi = protein.get_angles("psi")
    phi_by_res = {q[1].residue.index: row for q, row in zip(phi_atoms, phi)}
    psi_by_res = {q[1].residue.index: row for q, row in zip(psi_atoms, psi)}
    for k, r in enumerate(sq.residue_indices):
        assert np.allclose(sq.phi_angles[0, k], phi_by_res[r])
        assert np.allclose(sq.psi_angles[0, k], psi_by_res[r])

    # the 2D Hellinger distance must be computable and finite
    H = sq.compute_dihedral_hellingers()
    assert H.shape == (1, n - 2)
    assert np.all(np.isfinite(H))


def test_capped_chain_residue_indices_unchanged():
    import soursop
    from soursop.sssampling import SamplingQuality
    from soursop.sstools import find_trajectory_files

    data = soursop.get_data("test_data")
    wt, wp = find_trajectory_files(os.path.join(data, "sampling_quality/WT"), 3)
    sq = SamplingQuality(
        wt[:1],
        top_file=wp[0],
        method="2D angle distributions",
        seed=0,
        force_sequential=True,
    )
    protein = sq.trajs[0].proteinTrajectoryList[0]
    # capped <AA>: both non-cap residues carry phi and psi
    assert sq.residue_indices == [1, 2]
    assert sq.phi_angles.shape[1] == 2
    assert protein.ncap and protein.ccap
