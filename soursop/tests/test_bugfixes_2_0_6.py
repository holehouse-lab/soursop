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
    d = 10 * md.compute_distances(protein.traj[frame], np.array(pairs), periodic=False)[0]
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
