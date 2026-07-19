"""Regression tests for the correctness fixes shipped in SOURSOP 2.0.3.

Each test pins a specific bug that was fixed; the test name references the
behaviour that previously regressed. See CHANGELOG.md (2.0.3) for details.
"""

import os

import numpy as np
import pytest

import soursop
from soursop.sstrajectory import SSTrajectory
from soursop.ssexceptions import SSException

DATA = soursop.get_data("test_data")


def _prot(name):
    return SSTrajectory(
        os.path.join(DATA, f"{name}.xtc"), os.path.join(DATA, f"{name}.pdb")
    ).proteinTrajectoryList[0]


# ---------------------------------------------------------------------------
# ssprotein.get_local_to_global_correlation - terminal residue was dropped
# ---------------------------------------------------------------------------
@pytest.mark.filterwarnings("ignore::RuntimeWarning")
def test_l2g_includes_terminal_residue():
    # n_pairs (counts per separation) and array lengths are independent of the
    # random pair sampling, so no seed is needed. The RuntimeWarning filter
    # covers the (separate, pre-existing) undefined-correlation edge when a
    # sampled distance series has zero variance.
    P = _prot("ntl9_AA")
    raw, n_pairs, mean_corr, std_corr = P.get_local_to_global_correlation(verbose=False)
    # the largest sequence separation must contribute at least one pair
    # (the end-to-end pair); previously it contributed zero because the
    # terminal residue was excluded from all pairs.
    assert n_pairs[-1] >= 1
    assert len(mean_corr) == len(n_pairs) == len(std_corr)


# ---------------------------------------------------------------------------
# ssprotein COM-mode distance calculations - stride was silently ignored
# ---------------------------------------------------------------------------
def test_com_distance_map_honours_stride():
    P = _prot("ctl9_AA")
    m1, _ = P.get_distance_map(mode="COM", stride=1, verbose=False)
    m3, _ = P.get_distance_map(mode="COM", stride=3, verbose=False)
    assert not np.allclose(m1, m3)


def test_calculate_all_ca_distances_com_stride_shape():
    P = _prot("ctl9_AA")
    full = P.calculate_all_CA_distances(5, mode="COM", stride=1)
    strided = P.calculate_all_CA_distances(5, mode="COM", stride=3)
    assert strided.shape[0] < full.shape[0]


def test_com_distance_map_weights_and_stride_no_crash():
    P = _prot("ctl9_AA")
    w = np.full(P.n_frames, 1.0 / P.n_frames)
    # previously raised ValueError (weights strided to n<n_frames vs unstrided data)
    dmap, std = P.get_distance_map(mode="COM", stride=3, weights=w, verbose=False)
    assert dmap.shape[0] == len(P.resid_with_CA)
    assert std is None


# ---------------------------------------------------------------------------
# ssprotein.get_clusters - degenerate frame count now errors cleanly
# ---------------------------------------------------------------------------
def test_get_clusters_too_few_frames_raises():
    P = _prot("ntl9_AA")
    with pytest.raises(SSException):
        P.get_clusters(stride=100000)


# ---------------------------------------------------------------------------
# ssprotein.get_angle_decay - single-frame trajectories crashed
# ---------------------------------------------------------------------------
def test_angle_decay_single_frame():
    P = _prot("ntl9_AA")
    single = SSTrajectory(TRJ=P.traj[0:1]).proteinTrajectoryList[0]
    out = single.get_angle_decay()  # previously raised AxisError on 1 frame
    assert out.ndim == 2 and out.shape[1] >= 2
    assert np.all(np.isfinite(out))


# ---------------------------------------------------------------------------
# ssprotein.get_secondary_structure_DSSP - fractions must sum to 1
# ---------------------------------------------------------------------------
def test_dssp_fractions_sum_to_one_uncapped():
    # ntl9 is uncapped, so DSSP emits 'NA' at the termini; those must fold into
    # coil so the H/E/C partition sums to 1.
    P = _prot("ntl9_AA")
    _, H, E, C = P.get_secondary_structure_DSSP()
    assert np.allclose(np.array(H) + np.array(E) + np.array(C), 1.0)

    _, H, E, C = P.get_secondary_structure_DSSP(return_per_frame=True)
    assert np.allclose(np.array(H) + np.array(E) + np.array(C), 1.0)


# ---------------------------------------------------------------------------
# ssprotein.get_distance_map - weighted std map is None (not a NaN array)
# ---------------------------------------------------------------------------
def test_distance_map_weighted_std_is_none():
    P = _prot("gs6_AA")
    w = np.full(P.n_frames, 1.0 / P.n_frames)
    _, std = P.get_distance_map(weights=w, verbose=False)
    assert std is None


# ---------------------------------------------------------------------------
# ssprotein.get_RMSD - NumPy integer frame indices were rejected
# ---------------------------------------------------------------------------
def test_get_rmsd_accepts_numpy_int():
    P = _prot("ntl9_AA")
    r_py = P.get_RMSD(0, frame2=3)
    r_np = P.get_RMSD(0, frame2=np.int64(3))
    assert len(r_py) == 1 and len(r_np) == 1
    assert np.allclose(np.asarray(r_py), np.asarray(r_np))


# ---------------------------------------------------------------------------
# ssprotein.get_asphericity - finite (guarded) output
# ---------------------------------------------------------------------------
def test_asphericity_finite():
    P = _prot("ntl9_AA")
    assert np.all(np.isfinite(P.get_asphericity(verbose=False)))


# ---------------------------------------------------------------------------
# ssdata - KM2 / KM3 are distinct valid residue names
# ---------------------------------------------------------------------------
def test_km2_km3_are_valid_residue_names():
    from soursop.ssdata import ALL_VALID_RESIDUE_NAMES

    assert "KM2" in ALL_VALID_RESIDUE_NAMES
    assert "KM3" in ALL_VALID_RESIDUE_NAMES
    assert "KM2KM3" not in ALL_VALID_RESIDUE_NAMES


# ---------------------------------------------------------------------------
# sstools.find_trajectory_files - replicate-index range matching
# ---------------------------------------------------------------------------
def test_find_trajectory_files_replicate_range(tmp_path):
    from soursop.sstools import find_trajectory_files

    for rep in ["rep0", "rep3", "rep7", "rep15"]:
        d = tmp_path / "setA" / rep
        d.mkdir(parents=True)
        (d / "__traj.xtc").write_text("x")
        (d / "__START.pdb").write_text("x")

    trajs, tops = find_trajectory_files(str(tmp_path), num_replicates=5)
    got = sorted(os.path.basename(os.path.dirname(t)) for t in trajs)
    # rep7 (7) and rep15 (15) are out of range for num_replicates=5
    assert got == ["rep0", "rep3"]
    assert len(trajs) == len(tops) == 2


# ---------------------------------------------------------------------------
# ssnmr.compute_random_coil_chemical_shifts - phospho + perdeuteration guard
# and out-of-range character handling
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("seq", ["AA(SEP)AA", "AA(TPO)AA", "AA(PTR)AA"])
def test_perdeuteration_phospho_raises_ssexception(seq):
    from soursop.ssnmr import compute_random_coil_chemical_shifts

    # all three phospho residues must raise the documented SSException rather
    # than an IndexError (pSer/pThr previously slipped past the guard).
    with pytest.raises(SSException):
        compute_random_coil_chemical_shifts(seq, use_perdeuteration=True)


@pytest.mark.parametrize("seq", ["AAcAA", "AA]AA"])
def test_out_of_range_character_does_not_crash(seq):
    from soursop.ssnmr import compute_random_coil_chemical_shifts

    # characters that map to codes 26-35 previously raised IndexError; they are
    # now skipped, and the skipped residue does not desync the output labels.
    out = compute_random_coil_chemical_shifts(seq)
    assert isinstance(out, list)
    # 'AAcAA' -> the 'c' is skipped, leaving 4 valid residues (AAAA)
    assert all(d["Res"] in ("A",) for d in out)


# ---------------------------------------------------------------------------
# sssampling - EV reference path works for uncapped chains
# ---------------------------------------------------------------------------
def test_sampling_ev_reference_uncapped():
    from soursop.sssampling import SamplingQuality

    sq = SamplingQuality(
        [os.path.join(DATA, "ntl9_AA.xtc")],
        reference_list=None,
        top_file=os.path.join(DATA, "ntl9_AA.pdb"),
        method="1D angle distributions",
    )
    h = sq.compute_dihedral_hellingers()
    assert np.all(np.isfinite(h))


# ---------------------------------------------------------------------------
# ssbme.BMECustom - default chi^2 path (analytic gradient) fits cleanly
# ---------------------------------------------------------------------------
def test_bmecustom_default_path_fits():
    from soursop.ssbme import BMECustom

    rng = np.random.default_rng(1)
    calc = rng.normal(size=(60, 5))
    exp = calc.mean(axis=0) + 0.1
    res = BMECustom(exp, calc, uncertainty=0.5).fit(theta=1.0, verbose=False)
    w = res.weights
    assert abs(w.sum() - 1.0) < 1e-8
    assert w.min() >= -1e-12
    # reweighting should not increase the cost relative to the uniform prior
    prior = np.full(len(w), 1.0 / len(w))
    bc = BMECustom(exp, calc, uncertainty=0.5)
    assert bc._cost(w) <= bc._cost(prior) + 1e-9
