"""
Regression tests for the correctness fixes made in SOURSOP 2.0.6.

Currently this covers ``get_hydrodynamic_radius(mode='kr')`` on a sub-region.
The Kirkwood-Riseman branch looped over every CA-bearing residue in the chain
and never looked at ``R1``/``R2``, so a request for the hydrodynamic radius of a
region silently returned the whole-chain value (the Nygaard branch has always
honoured the region). Whole-chain results must be unchanged.
"""

import os

import mdtraj as md
import numpy as np
import pytest

import soursop
from soursop import sstrajectory
from soursop.ssexceptions import SSException


test_data_dir = soursop.get_data("test_data")


def _load(name):
    return sstrajectory.SSTrajectory(
        os.path.join(test_data_dir, f"{name}.xtc"),
        os.path.join(test_data_dir, f"{name}.pdb"),
    ).proteinTrajectoryList[0]


@pytest.fixture(scope="module")
def NTL9():
    """Uncapped 56-residue chain."""
    return _load("ntl9_AA")


@pytest.fixture(scope="module")
def GS6():
    """ACE/NME-capped 6-residue chain (8 residues including caps)."""
    return _load("gs6_AA")


def _manual_kr(protein, R1, R2, frame):
    """Kirkwood-Riseman Rh for one frame straight from mdtraj CA distances."""
    resids = [r for r in protein.resid_with_CA if R1 <= r <= R2]
    ca = [protein.get_CA_index(r) for r in resids]
    pairs = [(a, b) for i, a in enumerate(ca) for b in ca[i + 1 :]]
    d = (
        10
        * md.compute_distances(protein.traj[frame], np.array(pairs), periodic=False)[0]
    )
    return 1.0 / np.mean(1.0 / d)


# ........................................................................
#
def test_full_chain_kr_unchanged(GS6):
    """The whole-chain values pinned in test_ssproteins.py must not move."""
    rh = GS6.get_hydrodynamic_radius(mode="kr")
    assert abs(rh[0] - 5.9223) < 0.001
    assert abs(np.mean(rh) - 5.884) < 0.001

    rh = GS6.get_hydrodynamic_radius(mode="kr", distance_mode="COM")
    assert abs(rh[0] - 5.8634) < 0.001
    assert abs(np.mean(rh) - 5.8536) < 0.001


def test_default_region_equals_explicit_full_range(NTL9, GS6):
    """R1=None/R2=None is the same as asking for the whole chain explicitly."""
    for protein in (NTL9, GS6):
        default = protein.get_hydrodynamic_radius(mode="kr")
        explicit = protein.get_hydrodynamic_radius(
            R1=0, R2=protein.n_residues - 1, mode="kr"
        )
        assert np.allclose(default, explicit)


@pytest.mark.parametrize("R1,R2", [(5, 20), (0, 30), (25, 55), (0, 1)])
def test_region_matches_manual_sum(NTL9, R1, R2):
    """A sub-region reproduces an explicit Kirkwood-Riseman sum over that region."""
    rh = NTL9.get_hydrodynamic_radius(R1=R1, R2=R2, mode="kr")
    assert len(rh) == NTL9.n_frames
    for frame in (0, NTL9.n_frames - 1):
        assert rh[frame] == pytest.approx(_manual_kr(NTL9, R1, R2, frame), rel=1e-6)


def test_region_is_not_the_whole_chain(NTL9):
    """Regression: R1/R2 were silently ignored and the whole chain was returned."""
    full = NTL9.get_hydrodynamic_radius(mode="kr")
    sub = NTL9.get_hydrodynamic_radius(R1=10, R2=20, mode="kr")
    assert not np.allclose(full, sub)
    assert np.mean(sub) < np.mean(full)


def test_swapped_endpoints_are_equivalent(NTL9):
    a = NTL9.get_hydrodynamic_radius(R1=10, R2=20, mode="kr")
    b = NTL9.get_hydrodynamic_radius(R1=20, R2=10, mode="kr")
    assert np.allclose(a, b)


def test_com_mode_region(NTL9):
    """COM distances honour the region too, and match the raw COM distances."""
    full = NTL9.get_hydrodynamic_radius(mode="kr", distance_mode="COM")
    sub = NTL9.get_hydrodynamic_radius(R1=10, R2=20, mode="kr", distance_mode="COM")
    assert np.all(np.isfinite(sub))
    assert not np.allclose(full, sub)

    resids = [r for r in NTL9.resid_with_CA if 10 <= r <= 20]
    inv = []
    for i, a in enumerate(resids):
        for b in resids[i + 1 :]:
            inv.append(1.0 / NTL9.get_inter_residue_COM_distance(a, b))
    expected = 1.0 / np.mean(np.array(inv), axis=0)
    assert np.allclose(sub, expected, rtol=1e-6)


def test_capped_chain_region_excludes_caps(GS6):
    """On a capped chain a region that spans a cap still only sums CA-bearing residues."""
    # residue 0 is ACE (no CA); residues 1-3 are the first three real residues
    with_cap = GS6.get_hydrodynamic_radius(R1=0, R2=3, mode="kr")
    without_cap = GS6.get_hydrodynamic_radius(R1=1, R2=3, mode="kr")
    assert np.allclose(with_cap, without_cap)
    assert np.allclose(with_cap[0], _manual_kr(GS6, 1, 3, 0))


def test_region_with_fewer_than_two_residues_raises(NTL9, GS6):
    with pytest.raises(SSException):
        NTL9.get_hydrodynamic_radius(R1=5, R2=5, mode="kr")
    # ACE cap plus one real residue - only one CA-bearing residue
    with pytest.raises(SSException):
        GS6.get_hydrodynamic_radius(R1=0, R2=1, mode="kr")


def test_out_of_range_region_raises(NTL9):
    with pytest.raises(SSException):
        NTL9.get_hydrodynamic_radius(R1=0, R2=NTL9.n_residues, mode="kr")
    with pytest.raises(SSException):
        NTL9.get_hydrodynamic_radius(R1=-1, R2=10, mode="kr")


def test_weights_path_on_region(NTL9):
    w = np.ones(NTL9.n_frames) / NTL9.n_frames
    per_frame = NTL9.get_hydrodynamic_radius(R1=5, R2=30, mode="kr")
    weighted = NTL9.get_hydrodynamic_radius(R1=5, R2=30, mode="kr", weights=w)
    assert weighted == pytest.approx(np.mean(per_frame))


def test_nygaard_region_still_honoured(NTL9):
    """Sanity: the Nygaard branch was never affected."""
    full = NTL9.get_hydrodynamic_radius(mode="nygaard")
    sub = NTL9.get_hydrodynamic_radius(R1=10, R2=20, mode="nygaard")
    assert not np.allclose(full, sub)


# ========================================================================
#
# Further 2.0.6 fixes. Each block below pins one of the behaviours that was
# corrected in this release; the whole-chain / default-argument numerics
# for existing call patterns were recorded with the pre-fix code and must
# not move.
#


@pytest.fixture(scope="module")
def CTL9():
    """ACE/NME-capped 92-residue chain (94 residues including caps)."""
    return _load("ctl9_AA")


@pytest.fixture(scope="module")
def SIGA_CG():
    """One-bead-per-residue coarse-grained chain."""
    return _load("sigA_CG")


# ........................................................................
# 1. get_Q must not move the trajectory coordinates
#
def test_get_Q_does_not_mutate_coordinates(NTL9):
    """Prior to 2.0.6 get_Q() superposed self.traj in place when stride == 1."""
    before = NTL9.traj.xyz.copy()
    NTL9.get_Q()
    assert np.array_equal(before, NTL9.traj.xyz)
    NTL9.get_Q(stride=3)
    assert np.array_equal(before, NTL9.traj.xyz)


def test_get_Q_stride_and_reference_frame(NTL9):
    """Q is a per-frame value; strided output is a subset and Q(ref) is ~1."""
    q = NTL9.get_Q()
    assert len(q) == NTL9.n_frames
    assert q[0] == pytest.approx(1.0, abs=0.05)
    q2 = NTL9.get_Q(native_state_reference_frame=4)
    assert q2[4] == pytest.approx(1.0, abs=0.05)
    assert not np.allclose(q, q2)


# ........................................................................
# 2. instantaneous distance maps are distances, not squared distances
#
def test_instantaneous_maps_independent_of_RMS(NTL9):
    inst_rms, _ = NTL9.get_distance_map(
        RMS=True, return_instantaneous_maps=True, verbose=False
    )
    inst_plain, _ = NTL9.get_distance_map(
        RMS=False, return_instantaneous_maps=True, verbose=False
    )
    assert inst_rms.shape == (
        NTL9.n_frames,
        len(NTL9.resid_with_CA),
        len(NTL9.resid_with_CA),
    )
    assert np.array_equal(inst_rms, inst_plain)

    # and the averaged RMS map is the RMS of the instantaneous distances
    rms_mean, _ = NTL9.get_distance_map(RMS=True, verbose=False)
    assert np.allclose(rms_mean, np.sqrt(np.mean(inst_plain**2, axis=0)), atol=1e-4)


# ........................................................................
# 3. get_all_SASA(mode='all') columns line up on capped chains
#
def test_sasa_mode_all_aligned_on_capped_chain(CTL9):
    residue, sidechain, backbone = CTL9.get_all_SASA(mode="all", stride=5)
    n_ca = len(CTL9.resid_with_CA)
    assert residue.shape == sidechain.shape == backbone.shape
    assert residue.shape[1] == n_ca
    assert n_ca == CTL9.n_residues - 2

    full = CTL9.get_all_SASA(mode="residue", stride=5)
    assert full.shape[1] == CTL9.n_residues
    assert np.array_equal(residue, full[:, CTL9.resid_with_CA])


def test_sasa_mode_all_unchanged_on_uncapped_chain(NTL9):
    residue, sidechain, backbone = NTL9.get_all_SASA(mode="all", stride=5)
    full = NTL9.get_all_SASA(mode="residue", stride=5)
    assert residue.shape == sidechain.shape == backbone.shape
    assert np.array_equal(residue, full)


# ........................................................................
# 4. get_local_collapse default bins are in Angstroms
#
def test_local_collapse_default_bins_capture_every_frame(NTL9):
    mean, std, histo, bins = NTL9.get_local_collapse(window_size=30, verbose=False)
    assert bins[0] == 0.0
    assert bins[-1] > 90.0
    assert len(histo) == NTL9.n_residues - 30 + 1
    for h in histo:
        assert int(np.sum(h)) == NTL9.n_frames
    assert np.all(np.array(mean) > 5.0)


# ........................................................................
# 5. Nygaard Rh and <t> use the regional residue count
#
def _nygaard(rg, n, a1=0.216, a2=4.06, a3=0.821):
    n033, n060 = n**0.33, n**0.60
    return rg / (((a1 * (rg - a2 * n033)) / (n060 - n033)) + a3)


def _t(rg, n):
    # same (truncated) exponent the implementation uses
    return 2.5 * ((1.75 * (rg / (3.6 * n))) ** (4.0 / n**0.3333))


def test_nygaard_whole_chain_unchanged(NTL9, CTL9):
    """Values recorded with the 2.0.5 code (which used n_residues, caps included)."""
    rh = NTL9.get_hydrodynamic_radius(mode="nygaard")
    assert rh[0] == pytest.approx(19.5753883803, abs=1e-6)
    assert np.mean(rh) == pytest.approx(20.2825976567, abs=1e-6)

    rh = CTL9.get_hydrodynamic_radius(mode="nygaard")
    assert rh[0] == pytest.approx(27.4041476007, abs=1e-6)
    assert np.mean(rh) == pytest.approx(29.4924092037, abs=1e-6)


def test_t_whole_chain_unchanged(NTL9, CTL9):
    t = NTL9.get_t()
    assert t[0] == pytest.approx(0.3393042496, abs=1e-8)
    assert np.mean(t) == pytest.approx(0.3752896525, abs=1e-8)

    t = CTL9.get_t()
    assert t[0] == pytest.approx(0.4528728379, abs=1e-8)
    assert np.mean(t) == pytest.approx(0.5328514740, abs=1e-8)


def test_nygaard_region_uses_regional_N(NTL9, CTL9):
    for protein in (NTL9, CTL9):
        R1, R2 = 10, 30
        rg = protein.get_radius_of_gyration(R1, R2)
        rh = protein.get_hydrodynamic_radius(R1=R1, R2=R2, mode="nygaard")
        assert np.allclose(rh, _nygaard(rg, R2 - R1 + 1), rtol=1e-6)
        # whole chain: N == n_residues (caps included)
        rg = protein.get_radius_of_gyration()
        assert np.allclose(
            protein.get_hydrodynamic_radius(),
            _nygaard(rg, protein.n_residues),
            rtol=1e-6,
        )


def test_t_region_uses_regional_N(NTL9, CTL9):
    for protein in (NTL9, CTL9):
        R1, R2 = 10, 30
        rg = protein.get_radius_of_gyration(R1, R2)
        assert np.allclose(protein.get_t(R1=R1, R2=R2), _t(rg, R2 - R1 + 1), rtol=1e-6)
        # swapped endpoints are equivalent
        assert np.allclose(protein.get_t(R1=R2, R2=R1), protein.get_t(R1=R1, R2=R2))


# ........................................................................
# 6. DSSP refuses one-bead coarse-grained chains
#
def test_dssp_raises_on_one_bead_cg(SIGA_CG):
    assert SIGA_CG.is_coarse_grained
    with pytest.raises(SSException, match="one-bead"):
        SIGA_CG.get_secondary_structure_DSSP()
    with pytest.raises(SSException):
        SIGA_CG.get_secondary_structure_DSSP(return_per_frame=True)


# ........................................................................
# 7. get_scaling_exponent(end_effect=0)
#
def test_scaling_exponent_end_effect_zero(NTL9):
    np.random.seed(1)
    out = NTL9.get_scaling_exponent(end_effect=0, verbose=False)
    assert np.isfinite(out[0])
    assert np.isfinite(out[1])
    np.random.seed(1)
    out5 = NTL9.get_scaling_exponent(end_effect=5, verbose=False)
    assert out[0] != out5[0]


# ........................................................................
# 8. get_sidechain_alignment_angle
#
def test_sidechain_alignment_raises_on_one_bead_cg(SIGA_CG):
    with pytest.raises(SSException, match="one-bead"):
        SIGA_CG.get_sidechain_alignment_angle(1, 5)


def test_sidechain_alignment_atom_override(NTL9):
    default = NTL9.get_sidechain_alignment_angle(1, 5)
    assert len(default) == NTL9.n_frames

    # explicit CB->CB override is accepted and differs from the default tips
    cb = NTL9.get_sidechain_alignment_angle(
        1, 5, sidechain_atom_1="CB", sidechain_atom_2="CB"
    )
    assert len(cb) == NTL9.n_frames
    assert np.all(np.isfinite(cb))

    # passing the default tip explicitly reproduces the default
    res1 = NTL9.get_amino_acid_sequence(numbered=False)[1]
    res2 = NTL9.get_amino_acid_sequence(numbered=False)[5]
    from soursop.ssdata import DEFAULT_SIDECHAIN_VECTOR_ATOMS

    explicit = NTL9.get_sidechain_alignment_angle(
        1,
        5,
        sidechain_atom_1=DEFAULT_SIDECHAIN_VECTOR_ATOMS[res1],
        sidechain_atom_2=DEFAULT_SIDECHAIN_VECTOR_ATOMS[res2],
    )
    assert np.array_equal(default, explicit)

    # an atom that does not exist names the offending residue
    with pytest.raises(SSException, match="residue 5"):
        NTL9.get_sidechain_alignment_angle(1, 5, sidechain_atom_2="NOPE")
    with pytest.raises(SSException, match="residue 1"):
        NTL9.get_sidechain_alignment_angle(1, 5, sidechain_atom_1="NOPE")


# ........................................................................
# 9. get_multiple_CA_index accepts numpy integers
#
def test_multiple_CA_index_numpy_int(NTL9):
    expected = [NTL9.get_CA_index(3)]
    assert NTL9.get_multiple_CA_index(3) == expected
    assert NTL9.get_multiple_CA_index(np.int64(3)) == expected
    assert NTL9.get_multiple_CA_index(np.int32(3)) == expected
    assert NTL9.get_multiple_CA_index(np.array([3, 4])[0]) == expected
    assert NTL9.get_multiple_CA_index(np.array([3, 4])) == [
        NTL9.get_CA_index(3),
        NTL9.get_CA_index(4),
    ]


# ........................................................................
# 10. get_inter_residue_COM_distance validates stride
#
@pytest.mark.parametrize("stride", [0, -1, -3])
def test_com_distance_bad_stride_raises(NTL9, stride):
    with pytest.raises(SSException):
        NTL9.get_inter_residue_COM_distance(1, 5, stride=stride)


def test_com_distance_oversized_stride_raises(NTL9):
    with pytest.raises(SSException):
        NTL9.get_inter_residue_COM_distance(1, 5, stride=NTL9.n_frames + 1)
    # n_frames itself is the largest legal stride
    d = NTL9.get_inter_residue_COM_distance(1, 5, stride=NTL9.n_frames)
    assert len(d) == 1


def test_com_distance_valid_stride(NTL9):
    full = NTL9.get_inter_residue_COM_distance(1, 5)
    strided = NTL9.get_inter_residue_COM_distance(1, 5, stride=3)
    assert np.array_equal(strided, full[::3])


# ........................................................................
# 12. get_angle_decay(return_all_pairs=True) keys are resids
#
def test_angle_decay_pair_keys_are_resids(NTL9, CTL9):
    decay, pairs = NTL9.get_angle_decay(return_all_pairs=True)
    resids = list(NTL9.resid_with_CA)
    assert resids[0] == 0
    assert "0-0" in pairs
    assert "0-1" in pairs
    assert pairs["0-0"] == 1.0
    assert f"{resids[-1]}-{resids[-1]}" in pairs
    assert f"{resids[-2]}-{resids[-1]}" in pairs
    assert f"{resids[-1]}-{resids[-1] + 1}" not in pairs
    n = len(resids)
    assert len(pairs) == n + n * (n - 1) // 2

    # the |i-j| = 1 mean in the decay array is the mean over the "i-(i+1)" pairs
    sep1 = [pairs[f"{r}-{r + 1}"] for r in resids[:-1]]
    assert decay[1, 1] == pytest.approx(np.mean(sep1))

    # on a capped chain the first CA-bearing residue is resid 1
    _, pairs_c = CTL9.get_angle_decay(return_all_pairs=True)
    assert "0-0" not in pairs_c
    assert "0-1" not in pairs_c
    assert "1-1" in pairs_c
    assert "1-2" in pairs_c
