"""
Reading an XTC file whose final frame was only partly written.

A trajectory that is still being written (or whose run stopped mid-write)
ends part-way through a frame, which makes ``mdtraj.load`` fail outright.
``SSTrajectory`` reads the complete frames instead, drops the incomplete final
frame and warns. Damage anywhere else in the file still raises.

The truncated files are built from the bundled GS6 trajectory in pytest's
temporary directory.
"""

import os
import warnings

import mdtraj as md
import numpy as np
import pytest

import soursop
from soursop.ssexceptions import SoursopWarning
from soursop.sstrajectory import SSTrajectory

test_data_dir = soursop.get_data("test_data")
GS6_XTC = os.path.join(test_data_dir, "gs6_AA.xtc")
GS6_PDB = os.path.join(test_data_dir, "gs6_AA.pdb")


def _frame_offsets(xtc):
    with md.formats.XTCTrajectoryFile(xtc) as f:
        return [int(o) for o in f.offsets]


@pytest.fixture
def truncated_xtc(tmp_path):
    """GS6 with its final frame cut off half-way (header intact)."""
    raw = open(GS6_XTC, "rb").read()
    last = _frame_offsets(GS6_XTC)[-1]
    path = tmp_path / "unfinished.xtc"
    path.write_bytes(raw[: last + (len(raw) - last) // 2])
    return str(path)


def test_md_load_fails_on_the_truncated_file(truncated_xtc):
    # the premise: mdtraj on its own cannot read any of the file
    with pytest.raises(RuntimeError):
        md.load(truncated_xtc, top=GS6_PDB)


def test_complete_frames_are_read_and_final_frame_dropped(truncated_xtc):
    full = md.load(GS6_XTC, top=GS6_PDB)
    with pytest.warns(SoursopWarning, match="final frame is incomplete"):
        T = SSTrajectory(truncated_xtc, GS6_PDB)
    assert T.n_frames == full.n_frames - 1
    np.testing.assert_array_equal(T.traj.xyz, full.xyz[:-1])
    np.testing.assert_array_equal(T.traj.time, full.time[:-1])
    np.testing.assert_allclose(T.traj.unitcell_vectors, full.unitcell_vectors[:-1])
    # analyses run as normal on the recovered frames
    np.testing.assert_allclose(
        T.proteinTrajectoryList[0].get_radius_of_gyration(),
        SSTrajectory(GS6_XTC, GS6_PDB)
        .proteinTrajectoryList[0]
        .get_radius_of_gyration()[:-1],
    )


def test_list_of_files_is_handled(truncated_xtc):
    # sssampling passes trajectories to SSTrajectory as a list of paths
    with pytest.warns(SoursopWarning, match="final frame is incomplete"):
        T = SSTrajectory([truncated_xtc], GS6_PDB)
    assert T.n_frames == md.load(GS6_XTC, top=GS6_PDB).n_frames - 1


def test_damage_before_the_final_frame_still_raises(tmp_path):
    raw = open(GS6_XTC, "rb").read()
    offsets = _frame_offsets(GS6_XTC)
    # frame 1 cut in half, every other frame intact
    half = (offsets[2] - offsets[1]) // 2
    path = tmp_path / "corrupt.xtc"
    path.write_bytes(raw[: offsets[1] + half] + raw[offsets[2] :])
    with pytest.raises(RuntimeError):
        SSTrajectory(str(path), GS6_PDB)


def test_complete_file_loads_without_a_warning():
    with warnings.catch_warnings():
        warnings.simplefilter("error", SoursopWarning)
        T = SSTrajectory(GS6_XTC, GS6_PDB)
    assert T.n_frames == md.load(GS6_XTC, top=GS6_PDB).n_frames
