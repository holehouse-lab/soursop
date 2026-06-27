"""Tests for SWAN 2-bead (CA backbone / CB sidechain) coarse-grained support.

SWAN trajectories carry a single ``CA`` bead per residue plus a single ``CB``
bead per non-glycine residue (and nothing else). SOURSOP auto-detects this
model on load (``SSTrajectory.swan_trajectory``) and switches sidechain-vector
and secondary-structure analyses to CA/CB definitions. These tests verify
detection, the SWAN-specific behaviour, and that the broad geometric API runs
end-to-end on the three example SWAN trajectories.
"""

import os

import numpy as np
import pytest

import soursop
from soursop import sstrajectory
from soursop.ssexceptions import SSException


test_data_dir = soursop.get_data("test_data")


# The example helix trajectory imposes two helical subregions (SWAN HelixSpec,
# 1-based inclusive residue numbering): residues 7-21 are 100% helical and
# residues 32-46 are ~50% helical. In SOURSOP's 0-based residue indexing these
# are slices [6:21] and [31:46].
HELIX_SEG_FULL = slice(6, 21)  # 1-based residues 7-21, expected ~1.0
HELIX_SEG_HALF = slice(31, 46)  # 1-based residues 32-46, expected ~0.5
HELIX_TOL = 0.05


# ----------------------------------------------------------------------------
# Detection
# ----------------------------------------------------------------------------
def test_swan_detected_on_load(SWAN_HELIX_CO, SWAN_ASH1_CO, SWAN_ASYN_CO):
    assert SWAN_HELIX_CO.swan_trajectory is True
    assert SWAN_ASH1_CO.swan_trajectory is True
    assert SWAN_ASYN_CO.swan_trajectory is True


def test_swan_flag_propagates_to_protein(SWAN_HELIX_CP, SWAN_ASH1_CP, SWAN_ASYN_CP):
    assert SWAN_HELIX_CP.is_swan is True
    assert SWAN_ASH1_CP.is_swan is True
    assert SWAN_ASYN_CP.is_swan is True


def test_non_swan_not_detected(CTL9_CO, SIGA_CG_CO):
    # all-atom and one-bead CA-only coarse-grained trajectories are NOT SWAN
    assert CTL9_CO.swan_trajectory is False
    assert SIGA_CG_CO.swan_trajectory is False
    assert CTL9_CO.proteinTrajectoryList[0].is_swan is False
    assert SIGA_CG_CO.proteinTrajectoryList[0].is_swan is False


def test_swan_flag_can_be_forced():
    # forcing swan_trajectory=False on a real SWAN trajectory disables detection
    traj = sstrajectory.SSTrajectory(
        os.path.join(test_data_dir, "swan_trajectories/helix.xtc"),
        os.path.join(test_data_dir, "swan_trajectories/helix.pdb"),
        swan_trajectory=False,
    )
    assert traj.swan_trajectory is False
    assert traj.proteinTrajectoryList[0].is_swan is False


# ----------------------------------------------------------------------------
# Secondary structure (helicity / beta)
# ----------------------------------------------------------------------------
def test_helix_helicity_reproduced(SWAN_HELIX_CP):
    resids, H, E, C = SWAN_HELIX_CP.get_secondary_structure_DSSP()

    # the fully-helical subregion should be ~100% helical
    assert H[HELIX_SEG_FULL].mean() == pytest.approx(1.0, abs=HELIX_TOL)

    # the ~50% subregion should be ~50% helical
    assert H[HELIX_SEG_HALF].mean() == pytest.approx(0.5, abs=HELIX_TOL)

    # H / E / C are a partition at every residue
    assert np.allclose(H + E + C, 1.0)


def test_helix_coil_regions_have_low_helicity(SWAN_HELIX_CP):
    resids, H, E, C = SWAN_HELIX_CP.get_secondary_structure_DSSP()
    coil = np.concatenate([H[0:6], H[21:31], H[46:]])
    assert coil.mean() < HELIX_TOL


def test_idrs_have_no_helicity_or_beta(SWAN_ASH1_CP, SWAN_ASYN_CP):
    for protein in (SWAN_ASH1_CP, SWAN_ASYN_CP):
        resids, H, E, C = protein.get_secondary_structure_DSSP()
        assert H.mean() < HELIX_TOL
        assert E.mean() < HELIX_TOL
        assert np.allclose(H + E + C, 1.0)


def test_dssp_per_frame_shape(SWAN_HELIX_CP):
    resids, H, E, C = SWAN_HELIX_CP.get_secondary_structure_DSSP(return_per_frame=True)
    n_frames = SWAN_HELIX_CP.n_frames
    assert H.shape == (n_frames, len(resids))
    assert E.shape == (n_frames, len(resids))
    assert C.shape == (n_frames, len(resids))
    # binary masks that partition every (frame, residue)
    assert np.array_equal(H + E + C, np.ones_like(H))


# ----------------------------------------------------------------------------
# Sidechain vectors (CA -> CB for every residue)
# ----------------------------------------------------------------------------
def test_sidechain_alignment_angle(SWAN_HELIX_CP):
    seq = SWAN_HELIX_CP.get_amino_acid_sequence(numbered=False)
    non_gly = [i for i, name in enumerate(seq) if name != "GLY"][:2]
    angle = SWAN_HELIX_CP.get_sidechain_alignment_angle(non_gly[0], non_gly[1])
    assert angle.shape == (SWAN_HELIX_CP.n_frames,)
    assert np.isfinite(angle).all()
    # angle between two unit vectors is in [0, 180] degrees
    assert angle.min() >= -1e-6
    assert angle.max() <= 180.0 + 1e-6


def test_sidechain_alignment_angle_glycine_raises(SWAN_HELIX_CP):
    seq = SWAN_HELIX_CP.get_amino_acid_sequence(numbered=False)
    gly = [i for i, name in enumerate(seq) if name == "GLY"][0]
    non_gly = [i for i, name in enumerate(seq) if name != "GLY"][0]
    with pytest.raises(SSException):
        SWAN_HELIX_CP.get_sidechain_alignment_angle(gly, non_gly)


# ----------------------------------------------------------------------------
# Dihedral / BBSEG functions are undefined for SWAN and must raise clearly
# ----------------------------------------------------------------------------
def test_get_angles_raises(SWAN_ASH1_CP):
    for angle in ("phi", "psi", "omega", "chi1"):
        with pytest.raises(SSException):
            SWAN_ASH1_CP.get_angles(angle)


def test_bbseg_raises(SWAN_ASH1_CP):
    with pytest.raises(SSException):
        SWAN_ASH1_CP.get_secondary_structure_BBSEG()


def test_dihedral_mutual_information_raises(SWAN_ASH1_CP):
    with pytest.raises(SSException):
        SWAN_ASH1_CP.get_dihedral_mutual_information()


# ----------------------------------------------------------------------------
# Coverage smoke test: the broad geometric API runs end-to-end on SWAN models
# ----------------------------------------------------------------------------
@pytest.fixture(params=["helix", "ash1", "asyn"])
def swan_protein(request, SWAN_HELIX_CP, SWAN_ASH1_CP, SWAN_ASYN_CP):
    return {"helix": SWAN_HELIX_CP, "ash1": SWAN_ASH1_CP, "asyn": SWAN_ASYN_CP}[
        request.param
    ]


def test_swan_protein_geometric_api(swan_protein):
    p = swan_protein
    n = p.n_residues

    # scalar / per-frame observables
    assert np.all(np.isfinite(p.get_radius_of_gyration()))
    assert np.all(np.isfinite(p.get_end_to_end_distance()))
    assert np.all(np.isfinite(p.get_asphericity()))
    assert np.all(np.isfinite(p.get_hydrodynamic_radius()))
    assert np.all(np.isfinite(p.get_center_of_mass()))
    assert np.all(np.isfinite(p.get_molecular_volume()))

    # gyration tensor
    gt = p.get_gyration_tensor()
    assert np.all(np.isfinite(gt))

    # CA-based maps and distances
    dmap, _ = p.get_distance_map()
    assert dmap.shape == (n, n)
    cmap = p.get_contact_map(mode="ca")[0]
    assert cmap.shape == (n, n)
    assert np.all(np.isfinite(p.get_inter_residue_atomic_distance(1, n - 2)))
    assert np.all(np.isfinite(p.get_inter_residue_COM_distance(1, n - 2)))

    # per-residue COM (CA and full residue)
    assert p.get_residue_COM(1).shape == (p.n_frames, 3)
    assert p.get_residue_COM(1, "CA").shape == (p.n_frames, 3)

    # internal scaling and scaling exponent
    seps, dists = p.get_internal_scaling()
    assert len(seps) == len(dists)


def test_swan_trajectory_overall_api(SWAN_HELIX_CO, SWAN_ASH1_CO, SWAN_ASYN_CO):
    for traj in (SWAN_HELIX_CO, SWAN_ASH1_CO, SWAN_ASYN_CO):
        assert np.isfinite(traj.get_overall_radius_of_gyration()).all()
        assert np.isfinite(traj.get_overall_asphericity()).all()
        assert np.isfinite(traj.get_overall_hydrodynamic_radius()).all()
