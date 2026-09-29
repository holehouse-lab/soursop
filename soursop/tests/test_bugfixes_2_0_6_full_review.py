"""
Regression tests for the full-codebase review of SOURSOP 2.0.6.

Each test pins one fix from the "Bug fixes - full-codebase review" section of
the CHANGELOG. They are grouped by module and use the session fixtures from
conftest.py; nothing here mutates a shared fixture (copies of trajectories are
taken wherever coordinates or topologies are edited).
"""

import os
import pathlib
import warnings

import mdtraj as md
import numpy as np
import pytest

import soursop
from soursop import ssnmr, ssutils
from soursop.ssbme import BME, BMECustom, ExperimentalObservable, iBME
from soursop.sscoper import COPER, chi2_limit_scan, iCOPER
from soursop.ssexceptions import SoursopWarning, SSException
from soursop.sshdx import compute_Nc, compute_protection_factors
from soursop.ssmutualinformation import calc_MI, shan_entropy
from soursop.ssprotein import SSProtein
from soursop.sssampling import SamplingQuality
from soursop.sstools import get_distance_periodic
from soursop.sstrajectory import SSTrajectory

test_data_dir = soursop.get_data("test_data")


def _path(*parts):
    return os.path.join(test_data_dir, *parts)


def _uniform(n):
    return np.full(n, 1.0 / n)


# =============================================================================
# ssutils / sstools
# =============================================================================
class TestSharedHelpers:
    def test_validate_weights_float_stride_raises_ssexception(self):
        with pytest.raises(SSException):
            ssutils.validate_weights(_uniform(10), 10, stride=2.0)

    def test_validate_weights_checks_sum_before_stride(self):
        w = np.full(10, 0.5)  # sums to 5
        with pytest.raises(SSException, match="do not sum to 1"):
            ssutils.validate_weights(w, 10, stride=2)

    def test_stride_with_weights_emits_soursop_warning(self):
        with pytest.warns(SoursopWarning, match="stride with weights"):
            out = ssutils.validate_weights(_uniform(10), 10, stride=2)
        assert np.isclose(out.sum(), 1.0) and len(out) == 5

    def test_weighted_corr_degenerate_weights_raise(self):
        one_hot = np.zeros(5)
        one_hot[2] = 1.0
        with pytest.raises(SSException, match="at least two frames"):
            ssutils.weighted_corr(np.arange(5.0), np.arange(5.0) ** 2, one_hot)

    def test_build_scan_grid(self):
        np.testing.assert_allclose(
            ssutils.build_scan_grid([0.5, 1.0, 2.0, 4.0], 8, True, "x"),
            [0.5, 1.0, 2.0, 4.0],
        )
        assert len(ssutils.build_scan_grid((0.5, 4.0), 8, True, "x")) == 8
        for bad in [(0.1, 1.0, 2.0), np.array([]), [1.0, -1.0], [[1.0, 2.0]]]:
            with pytest.raises(SSException):
                ssutils.build_scan_grid(bad, 8, True, "x")

    def test_curvature_knee_is_interior(self):
        x = np.array([0.0, 0.1, 1.0])
        y = np.array([1.0, 0.1, 0.0])
        idx, _ = ssutils.find_optimal_theta(x, y, method="curvature")
        assert idx == 1

    def test_distance_periodic_beyond_one_and_a_half_boxes(self):
        d = get_distance_periodic([[0.0, 0.0, 0.0]], [[21.0, 0.0, 0.0]], 10.0)
        np.testing.assert_allclose(d, [1.0])
        # per-axis (slab) box
        d = get_distance_periodic([[0.0, 0.0, 0.0]], [[0.0, 0.0, 29.0]], [10, 10, 30])
        np.testing.assert_allclose(d, [1.0])

    def test_degenerate_regression_raises(self):
        with pytest.raises(SSException):
            ssutils.weighted_linear_regression(
                np.ones(4), np.arange(4.0), np.ones(4), fit_intercept=True
            )

    def test_observable_rejects_nan(self):
        with pytest.raises(SSException):
            ExperimentalObservable(1.0, np.nan)
        with pytest.raises(SSException):
            ExperimentalObservable(np.nan, 1.0)


# =============================================================================
# ssprotein
# =============================================================================
class TestSSProtein:
    def test_site_accessibility_keys_use_soursop_resid(self, GS6_CP):
        out = GS6_CP.get_site_accessibility([1, 2], mode="resid", stride=1)
        names = GS6_CP.get_amino_acid_sequence(numbered=False)
        assert sorted(out) == sorted([f"{names[1]}-1", f"{names[2]}-2"])

    def test_internal_scaling_stride_with_weights(self, NTL9_CP):
        n = NTL9_CP.n_frames
        base = NTL9_CP.get_internal_scaling(mean_vals=True, stride=2, verbose=False)
        with pytest.warns(SoursopWarning):
            wtd = NTL9_CP.get_internal_scaling(
                mean_vals=True, stride=2, weights=_uniform(n), verbose=False
            )
        np.testing.assert_allclose(base[1], wtd[1], rtol=1e-10)
        none = NTL9_CP.get_internal_scaling(
            mean_vals=True, stride=2, weights=None, verbose=False
        )
        np.testing.assert_allclose(base[1], none[1])

    def test_internal_scaling_rejects_cap_endpoint(self, GS6_CP):
        assert GS6_CP.ncap
        with pytest.raises(SSException, match="not a CA-bearing residue"):
            GS6_CP.get_internal_scaling(R1=0, mode="COM", verbose=False)

    def test_local_heterogeneity_weights_none_is_unweighted(self, GS6_CP):
        GS6_CP.get_local_heterogeneity(
            fragment_size=4, stride=1, weights=None, verbose=False
        )

    def test_local_to_global_validation(self, NTL9_CP):
        n = NTL9_CP.n_frames
        one_hot = np.zeros(n)
        one_hot[0] = 1.0
        with pytest.raises(SSException):
            NTL9_CP.get_local_to_global_correlation(
                stride=1, weights=one_hot, verbose=False
            )
        with pytest.raises(SSException):
            NTL9_CP.get_local_to_global_correlation(stride=1.5, verbose=False)
        with pytest.raises(SSException, match="distinct pairs"):
            NTL9_CP.get_local_to_global_correlation(
                stride=1, max_num_pairs=10**6, n_cycles=1, verbose=False
            )

    def test_local_to_global_pairs_without_replacement(self, GS6_CP):
        # GS6 has 6 CA residues -> 15 pairs; asking for 15 distinct pairs must
        # work (and every selection is then the full set, so the correlation is
        # identical across cycles)
        np.random.seed(0)
        raw, ks, mean_c, std_c = GS6_CP.get_local_to_global_correlation(
            stride=1, max_num_pairs=16, n_cycles=3, verbose=False
        )
        assert np.isclose(std_c[-1], 0.0, atol=1e-12)

    def test_local_collapse_weighted_histogram_counts(self, NTL9_CP):
        n = NTL9_CP.n_frames
        _, _, h0, _ = NTL9_CP.get_local_collapse(window_size=10, verbose=False)
        _, _, h1, _ = NTL9_CP.get_local_collapse(
            window_size=10, weights=_uniform(n), verbose=False
        )
        np.testing.assert_allclose(h0[0], h1[0])
        with pytest.raises(SSException, match="uneven"):
            NTL9_CP.get_local_collapse(bins=[0, 2, 3, 6], verbose=False)
        with pytest.raises(SSException):
            NTL9_CP.get_local_collapse(window_size=2.5, verbose=False)

    def test_sidechain_alignment_self_angle_and_degenerate_tip(self, NTL9_CP):
        theta = NTL9_CP.get_sidechain_alignment_angle(5, 5)
        assert np.all(theta < 1e-5)
        with pytest.raises(SSException, match="zero length"):
            NTL9_CP.get_sidechain_alignment_angle(5, 5, sidechain_atom_1="CA")

    def test_scaling_exponent_fit_reaches_largest_separation(self, NTL9_CP):
        out = NTL9_CP.get_scaling_exponent(end_effect=0, verbose=False, n_bootstrap=5)
        assert isinstance(out, list) and len(out) == 10
        fit_region, all_points = out[8], out[9]
        assert fit_region[0].max() == all_points[0].max()
        with pytest.raises(SSException):
            NTL9_CP.get_scaling_exponent(num_fitting_points=2, verbose=False)
        with pytest.raises(SSException):
            NTL9_CP.get_scaling_exponent(
                fraction_override=True, fraction_of_points=1.5, verbose=False
            )

    def test_regional_sasa_rejects_float(self, GS6_CP):
        with pytest.raises(SSException):
            GS6_CP.get_regional_SASA(1.5, 4)

    def test_molecular_volume_too_few_atoms(self, GS6_CP):
        with pytest.warns(SoursopWarning, match="no atoms named CA"):
            P = SSProtein(GS6_CP.traj.atom_slice([0, 1, 2]))
        with pytest.raises(SSException, match="four atoms"):
            P.get_molecular_volume()

    def test_residue_mass_raises_on_cg(self, SIGA_CG_CO):
        with pytest.raises(SSException, match="coarse-grained"):
            SIGA_CG_CO.proteinTrajectoryList[0].get_residue_mass(3)

    def test_unrecognised_two_bead_chain_dssp_and_angles_raise(self, SWAN_HELIX_CO):
        P = SSProtein(SWAN_HELIX_CO.proteinTrajectoryList[0].traj, swan=False)
        with pytest.raises(SSException, match="'NA'"):
            P.get_secondary_structure_DSSP()
        with pytest.raises(SSException, match="backbone"):
            P.get_angles("phi")

    def test_two_bead_window_longer_than_chain_warns(self, SWAN_HELIX_CP):
        n = len(SWAN_HELIX_CP.resid_with_CA)
        with pytest.warns(SoursopWarning, match="helix_window"):
            SWAN_HELIX_CP.get_secondary_structure_DSSP(helix_window=n + 1)
        with pytest.raises(SSException):
            SWAN_HELIX_CP.get_secondary_structure_DSSP(helix_window=2)


# =============================================================================
# sstrajectory
# =============================================================================
class TestSSTrajectory:
    def test_swan_detection_is_per_protein(self, SWAN_HELIX_CO):
        top = md.Topology()
        chain = top.add_chain()
        res = top.add_residue("NA", chain)
        top.add_atom("NA", md.element.sodium, res)
        t = SWAN_HELIX_CO.traj
        ion = md.Trajectory(
            np.zeros((t.n_frames, 1, 3)),
            top,
            unitcell_lengths=t.unitcell_lengths,
            unitcell_angles=t.unitcell_angles,
        )
        T = SSTrajectory(TRJ=t.stack(ion), check_whole_molecules=False)
        assert T.swan_trajectory
        assert T.proteinTrajectoryList[0].is_swan

    def test_whole_molecule_check_catches_split_sidechain(self, NTL9_CO):
        t = NTL9_CO.traj[:]
        res = t.topology.residue(10)
        side = [
            a.index for a in res.atoms if a.name not in ("N", "CA", "C", "O", "H", "HA")
        ]
        t.xyz[:, side, 0] += t.unitcell_lengths[:, [0]]
        with pytest.raises(SSException, match="split across the periodic boundary"):
            SSTrajectory(TRJ=t, check_whole_molecules="raise")

    def test_whole_molecule_check_ignores_chain_boundary_in_group(self, GMX_2CHAINS):
        n_res = GMX_2CHAINS.traj.n_residues
        with warnings.catch_warnings():
            warnings.simplefilter("error", SoursopWarning)
            T = SSTrajectory(
                TRJ=GMX_2CHAINS.traj, protein_grouping=[list(range(n_res))]
            )
        assert not T.check_molecules_whole()[0]["split"]
        with pytest.raises(SSException):
            T.check_molecules_whole(chunk_size=0)

    def test_overall_rg_respects_protein_grouping(self, NTL9_CO):
        T = SSTrajectory(TRJ=NTL9_CO.traj, protein_grouping=[list(range(30))])
        np.testing.assert_allclose(
            T.get_overall_radius_of_gyration(),
            T.proteinTrajectoryList[0].get_radius_of_gyration(),
        )

    def test_protein_grouping_out_of_range_raises(self, NTL9_CO):
        with pytest.raises(SSException, match="valid indices"):
            SSTrajectory(TRJ=NTL9_CO.traj, protein_grouping=[[0, 1, 2, 999]])
        with pytest.raises(SSException):
            SSTrajectory(TRJ=NTL9_CO.traj, protein_grouping=[[True, 2]])

    def test_periodic_interchain_invariant_to_lattice_shift(self, GMX_2CHAINS):
        t = GMX_2CHAINS.traj[:]
        chain1 = [a.index for a in t.topology.chain(1).atoms]
        t.xyz[:, chain1, 0] += 2 * t.unitcell_lengths[:, [0]]
        T2 = SSTrajectory(TRJ=t, check_whole_molecules=False)
        d1 = GMX_2CHAINS.get_interchain_distance(0, 1, 10, 20, periodic=True)
        d2 = T2.get_interchain_distance(0, 1, 10, 20, periodic=True)
        np.testing.assert_allclose(d1, d2, atol=1e-3)
        m1, _ = GMX_2CHAINS.get_interchain_distance_map(0, 1, periodic=True)
        m2, _ = T2.get_interchain_distance_map(0, 1, periodic=True)
        np.testing.assert_allclose(m1, m2, atol=1e-3)

    def test_interchain_input_validation(self, GMX_2CHAINS):
        with pytest.raises(SSException, match="out of range"):
            GMX_2CHAINS.get_interchain_distance(-1, 0, 1, 1)
        with pytest.raises(SSException, match="integer"):
            GMX_2CHAINS.get_interchain_distance(0, 1, 2.7, 1)
        with pytest.raises(SSException):
            GMX_2CHAINS.get_interchain_contact_map(0, 1, stride=True)
        with pytest.raises(SSException):
            GMX_2CHAINS.get_interchain_distance_map(5, 0)

    def test_pdblead_prepends_only_the_first_model(self, GS6_CO, tmp_path):
        pdb = str(tmp_path / "three_models.pdb")
        GS6_CO.traj[:3].save_pdb(pdb)
        T = SSTrajectory(pdb, pdb, pdblead=True, check_whole_molecules=False)
        assert T.n_frames == 4


# =============================================================================
# ssbme
# =============================================================================
def _bme_problem(scale=1.0):
    rng = np.random.default_rng(1)
    calc = np.column_stack([rng.normal(20, 3, 500), rng.normal(5, 1, 500)])
    obs = [
        ExperimentalObservable(18 * scale, 0.5 * scale),
        ExperimentalObservable(5.5 * scale, 0.2 * scale),
    ]
    return obs, calc * scale


class TestBME:
    def test_units_do_not_change_the_answer(self):
        results = []
        for scale in (1.0, 1e-6):
            obs, calc = _bme_problem(scale)
            results.append(
                BME(obs, calc).fit(theta=2.0, auto_theta=False, verbose=False)
            )
        assert results[1].n_iterations > 0
        np.testing.assert_allclose(results[0].weights, results[1].weights, atol=1e-8)

    def test_unbounded_optimizer_with_one_sided_observable_raises(self):
        obs, calc = _bme_problem()
        obs[0] = ExperimentalObservable(25.0, 0.5, constraint="upper")
        with pytest.raises(SSException, match="cannot handle bounds"):
            BME(obs, calc).fit(theta=2.0, optimizer="BFGS", verbose=False)

    def test_nan_theta_raises(self):
        obs, calc = _bme_problem()
        with pytest.raises(SSException):
            BME(obs, calc).fit(theta=np.nan, verbose=False)

    def test_column_vector_prior_raises(self):
        obs, calc = _bme_problem()
        with pytest.raises(SSException, match="1D"):
            BME(obs, calc, np.ones((500, 1)))

    def test_nan_calc_on_zero_prior_frames_is_allowed(self):
        obs, calc = _bme_problem()
        calc = calc.copy()
        calc[:10] = np.nan
        w0 = np.ones(500)
        w0[:10] = 0
        res = BME(obs, calc, w0).fit(theta=2.0, auto_theta=False, verbose=False)
        assert res.success and np.all(res.weights[:10] == 0)

    def test_bmecustom_zero_prior_frame(self):
        rng = np.random.default_rng(0)
        x = rng.normal(size=(200, 1))
        w0 = np.exp(-1.5 * x[:, 0])
        w0[0] = 0.0
        w0 /= w0.sum()
        res = BMECustom(np.array([x[:, 0].mean()]), x, initial_weights=w0).fit(
            theta=0.1, verbose=False
        )
        assert res.weights[0] == 0.0
        assert np.isfinite(res.kl_divergence) and 0 < res.phi < 1

    def test_bmecustom_success_and_failed_diagnostics(self):
        rng = np.random.default_rng(0)
        calc = rng.normal(size=(400, 20))
        exp = calc.mean(axis=0) + 0.5
        res = BMECustom(exp, calc).fit(theta=0.05, max_iterations=1, verbose=False)
        assert not res.success
        assert res.diagnostics()["status"] == "FAILED"

    def test_bmecustom_scan_leaves_state_untouched(self):
        obs, calc = _bme_problem()
        bc = BMECustom(np.array([18.0, 5.5]), calc, uncertainty=np.array([0.5, 0.2]))
        bc.scan_theta(theta_range=(0.1, 10.0), n_points=3)
        assert bc.result is None and bc.theta is None

    def test_ibme_unconverged_is_not_success(self):
        rng = np.random.default_rng(3)
        calc = rng.normal([24, 65, 12], [3, 6, 2], size=(400, 3))
        obs = [
            ExperimentalObservable(v, s) for v, s in zip([23, 60, 11.5], [1, 2, 0.5])
        ]
        with pytest.warns(SoursopWarning, match="did not converge"):
            res = iBME(obs, 1.1 * calc + 0.3).fit(
                theta=1.0, max_ibme_iterations=1, verbose=False
            )
        assert not res.success and res.metadata["converged"] is False
        assert res.n_iterations == 1


# =============================================================================
# sscoper
# =============================================================================
class TestCOPER:
    def test_large_single_group_problem_is_feasible(self):
        rng = np.random.default_rng(0)
        F = rng.normal(0.0, 1.0, size=(10000, 1))
        obs = [ExperimentalObservable(F.mean() + 2.0, 1.0)]
        r = COPER(obs, F).fit(chi2_limit=1.0, verbose=False)
        assert r.feasible and r.success
        assert r.chi_squared_final <= 1.0 + 1e-5
        # exact dual solution of this problem
        assert abs(r.kl_divergence - 0.510476) < 1e-3

    def test_clearly_infeasible_problem_is_certified(self):
        rng = np.random.default_rng(0)
        F = rng.normal(0.0, 1.0, size=(500, 1))
        r = COPER([ExperimentalObservable(10.0, 0.1)], F).fit(verbose=False)
        assert not r.feasible and r.metadata["infeasibility_certified"]

    def test_feasible_prior_is_returned_unchanged(self):
        rng = np.random.default_rng(0)
        F = rng.normal(0, 1, (500, 2))
        obs = [
            ExperimentalObservable(F[:, 0].mean() + 0.9, 1.0),
            ExperimentalObservable(F[:, 1].mean(), 1.0),
        ]
        c = COPER(obs, F)
        r = c.fit(chi2_limit=1.0, verbose=False)
        np.testing.assert_array_equal(r.weights, c.initial_weights)

    def test_optimizer_and_limit_validation(self):
        rng = np.random.default_rng(0)
        F = rng.normal(0, 1, (100, 1))
        obs = [ExperimentalObservable(1.0, 1.0)]
        with pytest.raises(SSException, match="nonlinear"):
            COPER(obs, F).fit(optimizer="L-BFGS-B", verbose=False)
        with pytest.raises(SSException):
            COPER(obs, F).fit(chi2_limit=np.nan, verbose=False)

    def test_scan_list_is_explicit_values(self):
        rng = np.random.default_rng(0)
        F = rng.normal(0, 1, (300, 1))
        obs = [ExperimentalObservable(F.mean() + 0.5, 0.2)]
        scan = chi2_limit_scan(obs, F, chi2_limits=[0.5, 1.0, 2.0, 4.0])
        np.testing.assert_allclose(scan.chi2_limits, [0.5, 1.0, 2.0, 4.0])

    def test_icoper_iterates_to_self_consistency(self):
        rng = np.random.default_rng(0)
        N, M = 400, 12
        base = rng.normal(0, 1, size=(N, M)) + np.linspace(1, 3, M)
        tilt = np.exp(0.8 * base[:, 0] - 0.5 * base[:, 5])
        tilt /= tilt.sum()
        exp_vals = 2.5 * (tilt @ base) + 0.7
        obs = [ExperimentalObservable(v, 0.05) for v in exp_vals]
        r = iCOPER(obs, base).fit(chi2_limit=1.0, verbose=False)
        assert r.success and r.metadata["converged"]
        assert r.n_iterations > 2
        assert abs(r.scale - 2.5) < 0.1

    def test_group_named_all_is_rejected_when_mixed(self):
        F = np.random.default_rng(0).normal(size=(50, 2))
        obs = [
            ExperimentalObservable(0.0, 1.0, group="all"),
            ExperimentalObservable(0.0, 1.0),
        ]
        with pytest.raises(SSException, match="'all'"):
            COPER(obs, F)


# =============================================================================
# sssampling
# =============================================================================
class TestSamplingQuality:
    NTL9 = (_path("ntl9_AA.xtc"), _path("ntl9_AA.pdb"))

    def test_truncate_keeps_protein_grouping_and_protein_id(self):
        kw = dict(protein_grouping=[list(range(20)), list(range(20, 56))])
        a = SamplingQuality(
            [self.NTL9[0]], top_file=self.NTL9[1], seed=0, proteinID=1, **kw
        )
        b = SamplingQuality(
            [self.NTL9[0]],
            top_file=self.NTL9[1],
            seed=0,
            proteinID=1,
            truncate=True,
            verbose=False,
            **kw,
        )
        assert b.proteinID == 1
        np.testing.assert_allclose(a.phi_angles, b.phi_angles)

    def test_helicity_cache_is_per_protein(self):
        kw = dict(protein_grouping=[list(range(20)), list(range(20, 56))])
        sq = SamplingQuality([self.NTL9[0]], top_file=self.NTL9[1], seed=0, **kw)
        assert sq.compute_frac_helicity()[0].shape[1] == 20
        assert sq.compute_frac_helicity(proteinID=1)[0].shape[1] == 36

    @pytest.mark.parametrize("name", ["ACE_NH2", "FOR_NME", "FOR_NH2", "UCAP_NH2"])
    def test_for_and_nh2_caps_use_the_ev_reference(self, name):
        from soursop.ssdata import PHI_EV_ANGLES_DICT

        pdb = _path("cap_tests", f"{name}.pdb")
        sq = SamplingQuality([pdb], top_file=pdb, seed=0, check_whole_molecules=False)
        assert sq.phi_angles.shape[1] == sq.ref_phi_angles.shape[1]
        assert len(PHI_EV_ANGLES_DICT) == 20

    def test_pathlib_topology_and_bwidth(self):
        sq = SamplingQuality(
            [self.NTL9[0]], top_file=pathlib.Path(self.NTL9[1]), seed=0
        )
        assert sq.phi_angles.shape[1] == 54
        with pytest.raises(SSException, match="at least two bins"):
            SamplingQuality([self.NTL9[0]], top_file=self.NTL9[1], bwidth=2 * np.pi)

    def test_argument_validation(self):
        for bad in (
            dict(proteinID=3),
            dict(seed=-1),
            dict(reference_list=[self.NTL9[0]]),
        ):
            with pytest.raises(SSException):
                SamplingQuality([self.NTL9[0]], top_file=self.NTL9[1], **bad)


# =============================================================================
# sshdx
# =============================================================================
class TestHDX:
    def test_unrecognised_cap_is_not_reported(self, GS6_CP):
        t = GS6_CP.traj[:]
        for residue in t.topology.residues:
            if residue.name == "NME":
                residue.name = "NAC"
        P = SSProtein(t)
        res, _ = compute_protection_factors(P)
        names = [t.topology.residue(int(r)).name for r in res]
        assert "NAC" not in names

    def test_deuterium_is_not_a_heavy_atom(self, NTL9_CP):
        _, nc_h = compute_Nc(NTL9_CP)
        t = NTL9_CP.traj[:]
        for atom in t.topology.atoms:
            if atom.element.symbol == "H":
                atom.element = md.element.deuterium
        _, nc_d = compute_Nc(SSProtein(t))
        np.testing.assert_array_equal(nc_h, nc_d)

    def test_cutoff_and_exclusion_validation(self, NTL9_CP):
        with pytest.warns(SoursopWarning, match="nanometres"):
            compute_Nc(NTL9_CP, contact_cutoff=6.5, stride=5)
        with pytest.raises(SSException):
            compute_Nc(NTL9_CP, contact_cutoff=-1)
        with pytest.raises(SSException):
            compute_Nc(NTL9_CP, exclude_neighbours=None)


# =============================================================================
# ssnmr / sspre / ssmutualinformation
# =============================================================================
class TestNMRPREandMI:
    def test_noe_atom_index_validation(self, NTL9_CP):
        n = NTL9_CP.traj.n_atoms
        with pytest.raises(SSException, match="outside"):
            ssnmr.compute_NOE_distances(NTL9_CP, [[0, n]])
        with pytest.raises(SSException, match="integer"):
            ssnmr.compute_NOE_distances(NTL9_CP, [[0.7, 10.9]])

    def test_noe_average_stride_and_power(self, NTL9_CP):
        n = NTL9_CP.n_frames
        d = ssnmr.compute_NOE_distances(NTL9_CP, [[0, 100]], stride=2)
        with pytest.warns(SoursopWarning):
            wtd = ssnmr.noe_ensemble_average(d, weights=_uniform(n), stride=2)
        np.testing.assert_allclose(wtd, ssnmr.noe_ensemble_average(d))
        with pytest.raises(SSException):
            ssnmr.noe_ensemble_average(d, power=0)

    def test_chemical_shift_strings_have_three_decimals(self):
        out = ssnmr.compute_random_coil_chemical_shifts(
            "GP(TPO)CG", temperature=0, pH=14, asFloat=False
        )
        for row in out:
            for key in ("CA", "CB", "CO", "N", "HN", "HA"):
                value = row[key]
                if "*" not in value:
                    assert len(value.split(".")[1]) == 3, (key, value)

    def test_chemical_shift_input_checks(self):
        with pytest.raises(SSException):
            ssnmr.compute_random_coil_chemical_shifts("AGA", pH=np.nan)
        with pytest.warns(SoursopWarning, match="phospho"):
            ssnmr.compute_random_coil_chemical_shifts("ApSGA")

    def test_pre_two_bead_default_target(self, SWAN_HELIX_CP):
        from soursop.sspre import SSPRE

        pre = SSPRE(
            SWAN_HELIX_CP, tau_c=5.0, t_delay=10.0, R_2D=10.0, W_H=2 * np.pi * 600e6
        )
        label = SWAN_HELIX_CP.resid_with_CA[5]
        profile, gamma = pre.generate_PRE_profile(label, n_label_conformers=20)
        assert len(profile) == len(SWAN_HELIX_CP.resid_with_CA)
        with pytest.raises(SSException):
            pre.generate_PRE_profile(label, label_bead_radius=-5.5)

    def test_pre_linear_frequency_warns(self, NTL9_CP):
        from soursop.sspre import SSPRE

        with pytest.warns(SoursopWarning, match="W_H"):
            SSPRE(NTL9_CP, tau_c=5.0, t_delay=10.0, R_2D=10.0, W_H=1.2e9)

    def test_calc_mi_validation(self):
        rng = np.random.default_rng(0)
        X = rng.normal(size=500)
        Y = rng.normal(size=500)
        bins = np.linspace(-6, 6, 25)
        Xn = X.copy()
        Xn[0] = np.nan
        with pytest.raises(SSException):
            calc_MI(Xn, Y, bins)
        w = np.ones(500)
        w[3] = -1
        with pytest.raises(SSException):
            calc_MI(X, Y, bins, weights=w)
        # a two-edge (single-bin) grid gives zero MI rather than crashing
        assert calc_MI(X, Y, np.array([-6.0, 6.0])) == pytest.approx(0.0)
        assert shan_entropy([1, 1, 1, 1]) == pytest.approx(np.log(4))


# =============================================================================
# etol on every weights method
# =============================================================================
def test_every_public_weights_function_takes_etol():
    """Guard: a public function or method that takes ``weights`` must also
    take ``etol`` so callers can set the sum-to-one tolerance."""
    import importlib
    import inspect

    missing = []
    for name in ("ssprotein", "sstrajectory", "ssnmr", "sshdx", "ssmutualinformation"):
        mod = importlib.import_module(f"soursop.{name}")
        for obj_name, obj in inspect.getmembers(mod):
            if inspect.isclass(obj) and obj.__module__ == mod.__name__:
                funcs = [
                    (f"{obj_name}.{n}", f)
                    for n, f in inspect.getmembers(obj, inspect.isfunction)
                    if not n.startswith("_")
                ]
            elif inspect.isfunction(obj) and obj.__module__ == mod.__name__:
                funcs = [(obj_name, obj)]
            else:
                continue
            for fname, f in funcs:
                params = inspect.signature(f).parameters
                if "weights" in params and "etol" not in params:
                    missing.append(f"{name}.{fname}")
    assert missing == []


def _almost_normalised(n):
    # sums to 1 + 1e-4: rejected by the default etol (1e-7), accepted by 1e-3
    w = np.full(n, 1.0 / n)
    w[0] += 1e-4
    return w


@pytest.mark.parametrize(
    "call",
    [
        lambda P, **kw: P.get_contact_map(**kw),
        lambda P, **kw: P.get_Q(protein_average=False, **kw),
        lambda P, **kw: P.get_dihedral_mutual_information(**kw),
        lambda P, **kw: P.get_local_to_global_correlation(
            stride=1, n_cycles=2, max_num_pairs=3, verbose=False, **kw
        ),
    ],
    ids=["contact_map", "Q", "dihedral_MI", "local_to_global"],
)
def test_new_etol_arguments_are_honoured(NTL9_CP, call):
    w = _almost_normalised(NTL9_CP.n_frames)
    with pytest.raises(SSException, match="do not sum to 1"):
        call(NTL9_CP, weights=w)
    call(NTL9_CP, weights=w, etol=1e-3)


def test_calc_mi_etol():
    rng = np.random.default_rng(0)
    X = rng.normal(size=200)
    Y = rng.normal(size=200)
    bins = np.linspace(-6, 6, 13)
    w = _almost_normalised(200)
    with pytest.raises(SSException, match="do not sum to 1"):
        calc_MI(X, Y, bins, weights=w)
    assert np.isfinite(calc_MI(X, Y, bins, weights=w, etol=1e-3))


def test_distance_map_weighted_std_matches_interchain(NTL9_CO):
    P = NTL9_CO.proteinTrajectoryList[0]
    w = np.random.default_rng(1).random(P.n_frames)
    w /= w.sum()
    _, std = P.get_distance_map(weights=w, verbose=False)
    _, std_i = NTL9_CO.get_interchain_distance_map(0, 0, weights=w)
    iu = np.triu_indices_from(std, 1)
    # the inter-chain path works from float32 centres of mass
    np.testing.assert_allclose(std[iu], std_i[iu], rtol=1e-5, atol=1e-5)
