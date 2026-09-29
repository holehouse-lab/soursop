"""
Regression tests for the 2.0.6 bug fixes in sssampling, ssbme and sscoper.

Each test is tied to a specific fix:

  * EV reference lookup uses the residue's own table (sssampling).
  * trj_pdfs / ref_pdfs no longer share the 'joint' cache key.
  * fractional_helicity honours the instance's proteinID.
  * BMEResult / COPERResult.predict handle 1D input.
  * COPER / iCOPER tolerate zero prior weights.
  * SamplingQuality raises on mismatched trajectory / reference counts.
  * BMECustom's theta is documented relative to BME's (factor m/2).
  * iBME / iCOPER refuse an underdetermined scale + offset fit.
  * quality_plot resolves the method / dihedral branch lazily.
"""

import os

import matplotlib
import numpy as np
import pytest

import soursop
from soursop.ssbme import BME, BMECustom, ExperimentalObservable, iBME
from soursop.sscoper import COPER, iCOPER
from soursop.ssexceptions import SSException
from soursop.sssampling import (
    EV_RESIDUE_MAPPER,
    PHI_EV_ANGLES_DICT,
    PSI_EV_ANGLES_DICT,
    PrecomputedDihedralInterface,
    SamplingQuality,
)
from soursop.sstools import find_trajectory_files
from soursop.sstrajectory import SSTrajectory

matplotlib.use("Agg")

DATA_DIR = soursop.get_data("test_data")


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
def _sampling_paths():
    wt_data = os.path.join(DATA_DIR, "sampling_quality/WT")
    ev_data = os.path.join(DATA_DIR, "sampling_quality/EV")
    wt_traj, wt_top = find_trajectory_files(wt_data, 3)
    ev_traj, ev_top = find_trajectory_files(ev_data, 3)
    return wt_traj, wt_top, ev_traj, ev_top


def _make_ensemble(n_frames=200, seed=42):
    rng = np.random.default_rng(seed)
    return rng.normal(loc=[24.0, 65.0], scale=[3.0, 6.0], size=(n_frames, 2))


def _two_observables():
    return [
        ExperimentalObservable(25.0, 0.5, name="a"),
        ExperimentalObservable(64.0, 1.0, name="b"),
    ]


# --------------------------------------------------------------------------
# 1. EV reference dihedral lookup
# --------------------------------------------------------------------------
class TestEVReferenceLookup:
    def test_gpg_phi_and_psi_rows(self):
        bins = np.arange(-180, 181, 15)
        ev = PrecomputedDihedralInterface("GPG", bins=bins, num_trajs=1, nsamples=10)

        phi = ev.gather_phi_reference_dihedrals("GPG")
        psi = ev.gather_psi_reference_dihedrals("GPG")
        assert phi.shape[0] == 3 and psi.shape[0] == 3

        # position 0: glycine with the default alanine preceding context
        assert np.array_equal(
            phi[0], PHI_EV_ANGLES_DICT["GLY"][EV_RESIDUE_MAPPER["ALA"]]
        )
        # position 1: proline preceded by glycine -> the PRO table
        assert np.array_equal(
            phi[1], PHI_EV_ANGLES_DICT["PRO"][EV_RESIDUE_MAPPER["GLY"]]
        )
        # position 2: glycine preceded by proline -> the GLY table, PRO context
        assert np.array_equal(
            phi[2], PHI_EV_ANGLES_DICT["GLY"][EV_RESIDUE_MAPPER["PRO"]]
        )

        # psi uses the following residue as context
        assert np.array_equal(
            psi[0], PSI_EV_ANGLES_DICT["GLY"][EV_RESIDUE_MAPPER["PRO"]]
        )
        assert np.array_equal(
            psi[1], PSI_EV_ANGLES_DICT["PRO"][EV_RESIDUE_MAPPER["GLY"]]
        )
        assert np.array_equal(
            psi[2], PSI_EV_ANGLES_DICT["GLY"][EV_RESIDUE_MAPPER["ALA"]]
        )

    def test_residue_identity_matters(self):
        """The old lookup keyed everything on the neighbour, so 'AP' and 'AA'
        gave the same phi row at position 1. They must differ now."""
        bins = np.arange(-180, 181, 15)
        ev = PrecomputedDihedralInterface("AA", bins=bins, num_trajs=1, nsamples=10)
        phi_ap = ev.gather_phi_reference_dihedrals("AP")[1]
        phi_aa = ev.gather_phi_reference_dihedrals("AA")[1]
        assert not np.array_equal(phi_ap, phi_aa)


# --------------------------------------------------------------------------
# 2. trj_pdfs / ref_pdfs cache collision
# --------------------------------------------------------------------------
class TestJointPdfCache:
    def test_trj_and_ref_joint_pdfs_are_distinct(self):
        wt_traj, wt_top, ev_traj, ev_top = _sampling_paths()
        sq = SamplingQuality(
            wt_traj,
            ev_traj,
            top_file=wt_top,
            ref_top=ev_top,
            method="2D angle distributions",
        )
        trj_joint = sq.trj_pdfs("joint")
        ref_joint = sq.ref_pdfs("joint")

        assert trj_joint.shape == ref_joint.shape
        assert not np.allclose(trj_joint, ref_joint)

        # the cached reference must match a fresh recompute
        ref_fresh = sq.ref_pdfs("joint", recompute=True)
        assert np.allclose(ref_joint, ref_fresh)
        # and the trajectory cache must be untouched by the reference calls
        assert np.allclose(sq.trj_pdfs("joint"), trj_joint)

        # order independence: reference first, then trajectory
        sq2 = SamplingQuality(
            wt_traj,
            ev_traj,
            top_file=wt_top,
            ref_top=ev_top,
            method="2D angle distributions",
        )
        assert np.allclose(sq2.ref_pdfs("joint"), ref_joint)
        assert np.allclose(sq2.trj_pdfs("joint"), trj_joint)


# --------------------------------------------------------------------------
# 3. fractional_helicity honours proteinID
# --------------------------------------------------------------------------
class TestHelicityProteinID:
    def test_helicity_uses_selected_chain(self):
        top = os.path.join(DATA_DIR, "gromacs2chains/top.pdb")
        trj = os.path.join(DATA_DIR, "gromacs2chains/traj.xtc")

        sq = SamplingQuality([trj], top_file=top, proteinID=1)
        trj_h, ref_h = sq.fractional_helicity()

        ref = SSTrajectory(trj, pdb_filename=top)
        expected_1 = ref.proteinTrajectoryList[1].get_secondary_structure_DSSP()[1]
        expected_0 = ref.proteinTrajectoryList[0].get_secondary_structure_DSSP()[1]

        assert trj_h.shape == (1, len(expected_1))
        assert np.allclose(trj_h[0], expected_1)
        # sanity: the two chains do not have identical DSSP helicity, so the
        # test would have caught the old chain-0 default
        assert len(expected_0) != len(expected_1) or not np.allclose(
            expected_0, expected_1
        )
        assert np.all(ref_h == 0)


# --------------------------------------------------------------------------
# 4. 1D input to predict
# --------------------------------------------------------------------------
class TestPredict1D:
    def test_bme_predict_1d(self):
        calc = _make_ensemble()
        res = BME(_two_observables(), calc).fit(
            theta=2.0, auto_theta=False, verbose=False
        )
        col = calc[:, 0]
        out = res.predict(col)
        assert out.shape == (1,)
        assert np.isclose(out[0], np.sum(res.weights * col))
        # consistent with the 2D path
        assert np.isclose(out[0], res.predict(calc)[0])

    def test_coper_predict_1d(self):
        calc = _make_ensemble()
        res = COPER(_two_observables(), calc).fit(chi2_limit=1.0, verbose=False)
        col = calc[:, 1]
        out = res.predict(col)
        assert out.shape == (1,)
        assert np.isclose(out[0], np.sum(res.weights * col))
        assert np.isclose(out[0], res.predict(calc)[1])

    def test_predict_frame_mismatch_raises(self):
        calc = _make_ensemble()
        res = BME(_two_observables(), calc).fit(
            theta=2.0, auto_theta=False, verbose=False
        )
        with pytest.raises(SSException):
            res.predict(calc[:-1, 0])


# --------------------------------------------------------------------------
# 5. zero prior weights in COPER / iCOPER
# --------------------------------------------------------------------------
class TestZeroPriorWeights:
    @staticmethod
    def _prior(n):
        w0 = np.ones(n)
        w0[:5] = 0.0
        return w0

    def test_coper_zero_prior_frames(self):
        calc = _make_ensemble()
        obs = _two_observables()
        w0 = self._prior(calc.shape[0])

        res = COPER(obs, calc, initial_weights=w0).fit(chi2_limit=1.0, verbose=False)
        assert res.success and res.feasible
        assert np.all(np.isfinite(res.weights))
        assert np.isclose(np.sum(res.weights), 1.0)
        assert np.all(res.weights[:5] == 0.0)
        assert np.all(res.weights[5:] > 0.0)
        assert res.chi_squared_final <= 1.0 + 1e-4
        assert 0.0 < res.phi <= 1.0
        assert np.isfinite(res.mean_delta_G_kT)
        assert res.metadata["n_zero_prior_frames"] == 5
        assert np.all(res.reweighting_factors[:5] == 0.0)

        # identical to fitting the sub-ensemble directly
        sub = COPER(obs, calc[5:]).fit(chi2_limit=1.0, verbose=False)
        assert np.allclose(res.weights[5:], sub.weights, atol=1e-6)

        # qualitatively the same as BME with the same prior
        bme = BME(obs, calc, initial_weights=w0).fit(
            theta=1.0, auto_theta=False, verbose=False
        )
        assert np.all(bme.weights[:5] == 0.0)
        assert np.corrcoef(bme.weights[5:], res.weights[5:])[0, 1] > 0.5

    def test_icoper_zero_prior_frames(self):
        calc = _make_ensemble()
        obs = _two_observables()
        w0 = self._prior(calc.shape[0])

        res = iCOPER(obs, 1.1 * calc + 0.5, initial_weights=w0).fit(
            chi2_limit=1.0, verbose=False, max_icoper_iterations=10
        )
        assert np.all(np.isfinite(res.weights))
        assert np.all(res.weights[:5] == 0.0)
        assert np.isclose(np.sum(res.weights), 1.0)
        assert np.isfinite(res.phi)

    def test_negative_prior_raises(self):
        calc = _make_ensemble()
        w0 = np.ones(calc.shape[0])
        w0[0] = -1.0
        with pytest.raises(SSException):
            COPER(_two_observables(), calc, initial_weights=w0)


# --------------------------------------------------------------------------
# 6. mismatched trajectory / reference counts
# --------------------------------------------------------------------------
class TestReferenceCountMismatch:
    def test_one_trajectory_many_references_raises(self):
        wt_traj, wt_top, ev_traj, ev_top = _sampling_paths()
        with pytest.raises(SSException, match="one-to-one"):
            SamplingQuality(
                wt_traj[:1],
                ev_traj,
                top_file=wt_top,
                ref_top=ev_top,
                method="1D angle distributions",
            )

    def test_two_trajectories_three_references_raises(self):
        wt_traj, wt_top, ev_traj, ev_top = _sampling_paths()
        with pytest.raises(SSException, match="one-to-one"):
            SamplingQuality(
                wt_traj[:2],
                ev_traj,
                top_file=wt_top,
                ref_top=ev_top,
                method="1D angle distributions",
            )

    def test_one_to_one_still_works(self):
        wt_traj, wt_top, ev_traj, ev_top = _sampling_paths()
        # the single-trajectory branch takes one topology, not a list
        sq = SamplingQuality(
            wt_traj[:1],
            ev_traj[:1],
            top_file=wt_top[0],
            ref_top=ev_top[0],
            method="1D angle distributions",
        )
        assert len(sq.trajs) == 1 and len(sq.ref_trajs) == 1


# --------------------------------------------------------------------------
# 7. BMECustom theta vs BME theta
# --------------------------------------------------------------------------
class TestBMECustomThetaScale:
    def test_theta_conversion_reproduces_bme(self):
        rng = np.random.default_rng(1)
        n, m = 300, 6
        calc = rng.normal(loc=np.linspace(10, 20, m), scale=2.0, size=(n, m))
        exp = np.linspace(10, 20, m) + 1.0
        sig = np.full(m, 0.5)
        obs = [ExperimentalObservable(v, s) for v, s in zip(exp, sig)]

        for theta in (0.5, 2.0):
            rb = BME(obs, calc).fit(theta=theta, auto_theta=False, verbose=False)
            rc = BMECustom(exp, calc, uncertainty=sig).fit(
                theta=2.0 * theta / m, verbose=False
            )
            assert np.allclose(rb.weights, rc.weights, atol=1e-4)
            assert np.isclose(rb.phi, rc.phi, rtol=1e-3)
            assert np.isclose(rb.chi_squared_final, rc.cost_final, rtol=1e-2)

            # and the same nominal theta is NOT the same problem (documented)
            rc_same = BMECustom(exp, calc, uncertainty=sig).fit(
                theta=theta, verbose=False
            )
            assert not np.isclose(rb.phi, rc_same.phi, rtol=1e-3)


# --------------------------------------------------------------------------
# 8. underdetermined scale + offset fit
# --------------------------------------------------------------------------
class TestUnderdeterminedOffsetFit:
    def test_ibme_single_observable_offset_raises(self):
        calc = _make_ensemble()[:, :1]
        obs = [ExperimentalObservable(25.0, 0.5)]
        with pytest.raises(SSException, match="two observables"):
            iBME(obs, calc).fit(theta=1.0, fit_offset=True, verbose=False)
        # scale-only is fine
        res = iBME(obs, calc).fit(theta=1.0, fit_offset=False, verbose=False)
        assert np.isfinite(res.scale)

    def test_icoper_single_observable_offset_raises(self):
        calc = _make_ensemble()[:, :1]
        obs = [ExperimentalObservable(25.0, 0.5)]
        with pytest.raises(SSException, match="two observables"):
            iCOPER(obs, calc).fit(chi2_limit=1.0, fit_offset=True, verbose=False)
        res = iCOPER(obs, calc).fit(chi2_limit=1.0, fit_offset=False, verbose=False)
        assert np.isfinite(res.scale)


# --------------------------------------------------------------------------
# 9. quality_plot branch selection
# --------------------------------------------------------------------------
class TestQualityPlotBranches:
    def test_2d_method_rejects_phi_psi_cleanly(self):
        wt_traj, wt_top, ev_traj, ev_top = _sampling_paths()
        sq = SamplingQuality(
            wt_traj,
            ev_traj,
            top_file=wt_top,
            ref_top=ev_top,
            method="2D angle distributions",
        )
        for dihedral in ("phi", "psi"):
            with pytest.raises(SSException, match="not available"):
                sq.quality_plot(dihedral=dihedral)
        fig, axd = sq.quality_plot(dihedral="2D")
        assert set(axd.keys()) == {"A", "B", "C", "D"}
        matplotlib.pyplot.close(fig)

    def test_1d_method_phi_psi_work_and_2d_rejected(self):
        wt_traj, wt_top, ev_traj, ev_top = _sampling_paths()
        sq = SamplingQuality(
            wt_traj,
            ev_traj,
            top_file=wt_top,
            ref_top=ev_top,
            method="1D angle distributions",
        )
        with pytest.raises(SSException, match="not available"):
            sq.quality_plot(dihedral="2D")
        for dihedral in ("phi", "psi"):
            fig, axd = sq.quality_plot(dihedral=dihedral)
            assert set(axd.keys()) == {"A", "B", "C", "D"}
            matplotlib.pyplot.close(fig)
