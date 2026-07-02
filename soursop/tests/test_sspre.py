# Import package, test suite, and other packages as needed
import os

import numpy as np
import soursop
from soursop.sspre import (
    SSPRE,
    _fibonacci_cap_directions,
    _rotation_from_z,
)


test_data_dir = soursop.get_data("test_data")


def test_pre_profiles_ctl9(CTL9_CP):

    # check we can correctly compute PRE profiles

    # PRE settings for Stenzowski et al 2020 (note in original calulations
    # we used a tau_c of 5 ns (ultimately used 4 ns for the paper but the saved
    # profiles here used 5)
    PRE = SSPRE(CTL9_CP, 5, 12, 14, 850000000)

    # for each label position
    for P in [61, 74, 96, 109, 119, 149]:
        corrected_val = P - 58

        # calculate PRE profile using the classic point-at-CB model
        # (use_label=False), which is what the 2019 reference data used;
        # use_label now defaults to True (the calibrated label-cloud model).
        vals = PRE.generate_PRE_profile(corrected_val, use_label=False)

        # get previously-calculated profile
        file_path = os.path.join(test_data_dir, f"pre_data/PRE_{P}_unweighted.csv")
        pre_test_data = np.loadtxt(file_path)

        # ensure the largest differences is less than 0.02. There is some numerical inprecision
        # and some small changes in how the calculations are done between 2019 and the current
        # version of SOURSOP but these do not materially influence the profiles
        assert max(np.abs(pre_test_data - vals[0])) < 0.02


# ----------------------------------------------------------------------
# Coarse-grained label-cloud model (use_label=True)
# ----------------------------------------------------------------------


def test_fibonacci_cap_directions():
    # full sphere: unit vectors covering both hemispheres
    dirs = _fibonacci_cap_directions(200, 180.0)
    assert dirs.shape == (200, 3)
    np.testing.assert_allclose(np.linalg.norm(dirs, axis=1), 1.0, atol=1e-8)
    assert dirs[:, 2].min() < -0.9 and dirs[:, 2].max() > 0.9

    # narrow cap: every direction within the requested half-angle of +z
    half = 40.0
    cap = _fibonacci_cap_directions(200, half)
    ang = np.degrees(np.arccos(np.clip(cap[:, 2], -1, 1)))
    assert ang.max() <= half + 1e-6

    # determinism
    np.testing.assert_array_equal(
        _fibonacci_cap_directions(50, 130.0), _fibonacci_cap_directions(50, 130.0)
    )


def test_rotation_from_z():
    # maps +z onto the target axis, and is a proper rotation
    for axis in [
        np.array([0.0, 0.0, 1.0]),
        np.array([0.0, 0.0, -1.0]),
        np.array([1.0, 2.0, -3.0]),
    ]:
        u = axis / np.linalg.norm(axis)
        R = _rotation_from_z(u)
        np.testing.assert_allclose(R @ np.array([0.0, 0.0, 1.0]), u, atol=1e-8)
        np.testing.assert_allclose(R @ R.T, np.eye(3), atol=1e-8)
        assert abs(np.linalg.det(R) - 1.0) < 1e-8


def test_use_label_is_default(CTL9_CP):
    # as of SOURSOP 2.0.2 use_label defaults to True: the no-argument call must
    # match the explicit label-cloud call and differ from the point model
    PRE = SSPRE(CTL9_CP, 5, 12, 14, 850000000)
    default, _ = PRE.generate_PRE_profile(38)
    cloud, _ = PRE.generate_PRE_profile(38, use_label=True)
    point, _ = PRE.generate_PRE_profile(38, use_label=False)
    np.testing.assert_array_equal(np.asarray(default), np.asarray(cloud))
    assert not np.allclose(np.asarray(default), np.asarray(point))


def test_label_cloud_basic(CTL9_CP):
    # the label-cloud model runs, returns finite profiles of the right length,
    # is bounded in [0, 1], and is reproducible
    PRE = SSPRE(CTL9_CP, 5, 12, 14, 850000000)
    point, _ = PRE.generate_PRE_profile(38, use_label=False)
    # use_label defaults to True, but we pass it explicitly for clarity
    cloud, gamma = PRE.generate_PRE_profile(38, use_label=True)

    assert len(cloud) == len(point)
    cloud = np.asarray(cloud)
    assert np.all(np.isfinite(cloud))
    assert cloud.min() >= -1e-9 and cloud.max() <= 1.0 + 1e-9

    # deterministic (Fibonacci cloud, no RNG)
    cloud2, _ = PRE.generate_PRE_profile(38, use_label=True)
    np.testing.assert_array_equal(cloud, np.asarray(cloud2))

    # the cloud broadens quenching relative to the point model, but the two
    # remain strongly correlated
    r = np.corrcoef(cloud, np.asarray(point))[0, 1]
    assert r > 0.8
    assert not np.allclose(cloud, np.asarray(point))


def test_label_cloud_bead_radius(CTL9_CP):
    # steric exclusion is always applied; changing the bead radius changes the
    # profile (a larger bead prunes more clashing conformers), and every radius
    # produces a valid, finite profile
    PRE = SSPRE(CTL9_CP, 5, 12, 14, 850000000)
    small, _ = PRE.generate_PRE_profile(38, use_label=True, label_bead_radius=4.5)
    large, _ = PRE.generate_PRE_profile(38, use_label=True, label_bead_radius=7.0)
    small = np.asarray(small)
    large = np.asarray(large)
    assert np.all(np.isfinite(small)) and np.all(np.isfinite(large))
    assert small.min() >= -1e-9 and large.max() <= 1.0 + 1e-9
    # a different bead radius yields a different (sterically re-pruned) profile
    assert not np.allclose(small, large)


def test_label_cloud_soft_steric(CTL9_CP):
    # the soft steric wall runs, is deterministic and bounded, differs from the
    # hard cutoff, and converges to the hard cutoff as the wall stiffens
    PRE = SSPRE(CTL9_CP, 5, 12, 14, 850000000)
    hard, _ = PRE.generate_PRE_profile(
        38, use_label=True, label_steric="hard", label_bead_radius=6.0
    )
    soft, _ = PRE.generate_PRE_profile(
        38,
        use_label=True,
        label_steric="soft",
        label_bead_radius=6.0,
        label_wall_stiffness=6.0,
    )
    hard = np.asarray(hard)
    soft = np.asarray(soft)
    assert np.all(np.isfinite(soft))
    assert soft.min() >= -1e-9 and soft.max() <= 1.0 + 1e-9
    # deterministic
    soft2, _ = PRE.generate_PRE_profile(
        38,
        use_label=True,
        label_steric="soft",
        label_bead_radius=6.0,
        label_wall_stiffness=6.0,
    )
    np.testing.assert_array_equal(soft, np.asarray(soft2))
    # soft differs from hard, but a very stiff wall converges to the hard cutoff
    assert not np.allclose(soft, hard)
    stiff, _ = PRE.generate_PRE_profile(
        38,
        use_label=True,
        label_steric="soft",
        label_bead_radius=6.0,
        label_wall_stiffness=120.0,
    )
    assert np.max(np.abs(np.asarray(stiff) - hard)) < 0.02


def test_label_cloud_bad_steric_raises(CTL9_CP):
    # an invalid steric mode is rejected
    import pytest
    from soursop.ssexceptions import SSException

    PRE = SSPRE(CTL9_CP, 5, 12, 14, 850000000)
    with pytest.raises(SSException):
        PRE.generate_PRE_profile(38, use_label=True, label_steric="medium")
