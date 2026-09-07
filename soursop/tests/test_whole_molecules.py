"""
Tests for the whole-molecule input contract.

SOURSOP applies no periodic-boundary corrections, so it relies on the chains
it is given being whole. ``SSTrajectory`` checks for the signature of a
chain wrapped across the box (two consecutive residues further apart than
half the shortest box vector) on load and exposes the same test as
``check_molecules_whole``.
"""

import os
import warnings

import numpy as np
import pytest

import soursop
from soursop.ssexceptions import SSException
from soursop.sstrajectory import SSTrajectory

test_data_dir = soursop.get_data("test_data")

BOXED_FIXTURES = [
    ("gs6_AA.pdb", "gs6_AA.xtc"),
    ("ntl9_AA.pdb", "ntl9_AA.xtc"),
    ("gromacs1chain/top.pdb", "gromacs1chain/traj.xtc"),
    ("gromacs2chains/top.pdb", "gromacs2chains/traj.xtc"),
    ("swan_trajectories/helix.pdb", "swan_trajectories/helix.xtc"),
]


def _load(pdb, xtc, **kwargs):
    return SSTrajectory(
        os.path.join(test_data_dir, xtc), os.path.join(test_data_dir, pdb), **kwargs
    )


def _wrapped_copy(T, first_shifted_resid=30, axis=0):
    """A copy of ``T.traj`` with the C-terminal part of the chain moved by one
    box vector, i.e. what a simulation engine's wrapping does to a chain that
    straddles the boundary."""
    t = T.traj[:]
    sel = t.topology.select(f"resid {first_shifted_resid} to {t.n_residues - 1}")
    t.xyz[:, sel, axis] += t.unitcell_lengths[:, [axis]]
    return t


@pytest.mark.parametrize("pdb,xtc", BOXED_FIXTURES)
def test_whole_fixtures_are_not_flagged(pdb, xtc):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        T = _load(pdb, xtc)
    assert not any("split across" in str(w.message) for w in caught)
    for entry in T.check_molecules_whole():
        assert entry["tested"]
        assert not entry["split"]
        assert entry["n_frames_split"] == 0
        assert entry["max_ca_distance"] < entry["half_box"]


def test_chain_larger_than_half_the_box_is_still_whole():
    # gromacs1chain spans more than half of its 7.4 nm box but is a whole
    # molecule: consecutive residues stay ~3.8 A apart, so it must pass
    T = _load("gromacs1chain/top.pdb", "gromacs1chain/traj.xtc")
    P = T.proteinTrajectoryList[0]
    assert np.max(P.get_end_to_end_distance()) > 0.5 * np.min(P.unitcell)
    entry = T.check_molecules_whole()[0]
    assert entry["tested"] and not entry["split"]
    assert entry["max_ca_distance"] < 5.0


def test_no_box_cannot_be_tested():
    T = _load("sigA_CG.pdb", "sigA_CG.xtc")
    assert T.traj.unitcell_lengths is None
    entry = T.check_molecules_whole()[0]
    assert not entry["tested"] and not entry["split"]
    assert entry["half_box"] is None and entry["max_ca_distance"] is None


def test_zero_length_box_cannot_be_tested():
    T = _load("ntl9_AA.pdb", "ntl9_AA.xtc")
    t = T.traj[:]
    t.unitcell_lengths = np.zeros_like(t.unitcell_lengths)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        T0 = SSTrajectory(TRJ=t)
    assert not any("split across" in str(w.message) for w in caught)
    assert not T0.check_molecules_whole()[0]["tested"]


class TestWrappedChain:
    def test_warns_on_load_with_details(self):
        T = _load("ntl9_AA.pdb", "ntl9_AA.xtc")
        wrapped = _wrapped_copy(T, first_shifted_resid=30)
        with pytest.warns(
            UserWarning, match="split across the periodic boundary"
        ) as rec:
            T2 = SSTrajectory(TRJ=wrapped)
        message = str(rec[0].message)
        assert "residues 29-30" in message
        assert f"{T.n_frames} of {T.n_frames} frames" in message
        assert "make_molecules_whole" in message

        entry = T2.check_molecules_whole()[0]
        assert entry["tested"] and entry["split"]
        assert entry["n_frames_split"] == T.n_frames
        assert entry["worst_pair"] == (29, 30)
        # the broken pair is roughly one box vector apart
        box = 10 * T.traj.unitcell_lengths[0, 0]
        assert abs(entry["max_ca_distance"] - box) < 10.0
        assert entry["max_ca_distance"] > entry["half_box"]

    def test_partial_wrap_counts_only_affected_frames(self):
        T = _load("ntl9_AA.pdb", "ntl9_AA.xtc")
        t = T.traj[:]
        sel = t.topology.select("resid 40 to 55")
        # wrap the tail in three frames only
        for f in (1, 4, 7):
            t.xyz[f, sel, 2] += t.unitcell_lengths[f, 2]
        with pytest.warns(UserWarning, match="3 of 10 frames"):
            T2 = SSTrajectory(TRJ=t)
        entry = T2.check_molecules_whole()[0]
        assert entry["n_frames_split"] == 3
        assert entry["worst_pair"] == (39, 40)
        assert entry["worst_frame"] in (1, 4, 7)

    def test_second_chain_is_reported_by_index(self):
        T = _load("gromacs2chains/top.pdb", "gromacs2chains/traj.xtc")
        t = T.traj[:]
        chain1 = t.topology.chain(1)
        n_res = chain1.n_residues
        sel = [a.index for r in list(chain1.residues)[n_res // 2 :] for a in r.atoms]
        t.xyz[:, sel, 1] += t.unitcell_lengths[:, [1]]
        with pytest.warns(UserWarning, match="protein 1:"):
            T2 = SSTrajectory(TRJ=t)
        report = T2.check_molecules_whole()
        assert not report[0]["split"]
        assert report[1]["split"]
        # resids are reported in the chain's own (0-based) numbering
        assert report[1]["worst_pair"] == (n_res // 2 - 1, n_res // 2)

    def test_raise_mode_and_off_switch(self):
        T = _load("ntl9_AA.pdb", "ntl9_AA.xtc")
        wrapped = _wrapped_copy(T)
        with pytest.raises(SSException, match="split across the periodic boundary"):
            SSTrajectory(TRJ=wrapped, check_whole_molecules="raise")
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            T2 = SSTrajectory(TRJ=wrapped, check_whole_molecules=False)
        assert not any("split across" in str(w.message) for w in caught)
        # the method itself is always available
        assert T2.check_molecules_whole()[0]["split"]
        with pytest.raises(SSException, match="check_whole_molecules"):
            SSTrajectory(TRJ=wrapped, check_whole_molecules="maybe")

    def test_check_passes_through_the_file_loader(self, tmp_path):
        # the same check runs when reading from disk, e.g. via parallel loaders
        T = _load("ntl9_AA.pdb", "ntl9_AA.xtc")
        wrapped = _wrapped_copy(T)
        xtc = str(tmp_path / "wrapped.xtc")
        pdb = str(tmp_path / "wrapped.pdb")
        wrapped.save_xtc(xtc)
        wrapped[0].save_pdb(pdb)
        with pytest.warns(UserWarning, match="split across the periodic boundary"):
            SSTrajectory(xtc, pdb)

    def test_chunking_gives_the_same_report(self):
        T = _load("ntl9_AA.pdb", "ntl9_AA.xtc")
        T2 = SSTrajectory(TRJ=_wrapped_copy(T), check_whole_molecules=False)
        a = T2.check_molecules_whole(chunk_size=3)
        b = T2.check_molecules_whole(chunk_size=5000)
        assert a == b
