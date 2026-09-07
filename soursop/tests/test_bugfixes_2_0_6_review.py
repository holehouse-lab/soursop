"""
Regression tests for the fixes made in the full 2.0.6 code review.

Each test pins one behaviour that was wrong before the review: residue-name
aliases in the one-letter sequence, the units of the RMS distance-map standard
deviation, the DSSP cap guard, the BME diagnostics threshold and auto-theta
optimizer forwarding, the COPER multi-group feasibility check, the shared
residue axis of the PENGUIN summary figure, and the stride validation of the
stand-alone HDX / NMR functions.
"""

import os
import platform
import warnings

import mdtraj as md
import numpy as np
import pytest

import soursop
from soursop import sshdx, ssnmr, ssutils
from soursop.ssbme import BME, ExperimentalObservable
from soursop.sscoper import COPER
from soursop.ssdata import RESIDUE_NAME_ALIASES, THREE_TO_ONE, normalize_residue_name
from soursop.ssexceptions import SSException
from soursop.sspre import SSPRE
from soursop.ssprotein import SSProtein
from soursop.sstools import find_trajectory_files

test_data_dir = soursop.get_data("test_data")


def _renamed_protein(protein, renames):
    """Copy ``protein``'s trajectory with residue ``index -> name`` renamed."""
    top = protein.traj.topology.copy()
    residues = list(top.residues)
    for idx, name in renames.items():
        residues[idx].name = name
    return SSProtein(md.Trajectory(protein.traj.xyz, top))


# ---------------------------------------------------------------------------
# 1. force-field residue names in the one-letter sequence
# ---------------------------------------------------------------------------
class TestResidueNameAliases:
    def test_aliases_map_onto_canonical_names(self):
        assert normalize_residue_name("HIE") == "HIS"
        assert normalize_residue_name("HSD") == "HIS"
        assert normalize_residue_name("ASH") == "ASP"
        assert normalize_residue_name("CYX") == "CYS"
        assert normalize_residue_name("NMA") == "NME"
        assert normalize_residue_name("ALA") == "ALA"
        # every alias resolves to something with a one-letter code
        for alias, canonical in RESIDUE_NAME_ALIASES.items():
            assert canonical in THREE_TO_ONE, alias

    def test_oneletter_sequence_with_amber_names(self, GS6_CP):
        # residue 0 is ACE, 1..6 are GSGSGS, 7 is NME
        renamed = _renamed_protein(GS6_CP, {1: "HIE", 2: "CYX", 3: "ASH"})
        assert renamed.get_amino_acid_sequence(oneletter=True) == "<HCDSGS>"
        assert renamed.get_amino_acid_sequence(oneletter=True, numbered=False)[:4] == [
            "<",
            "H",
            "C",
            "D",
        ]
        # the three-letter sequence keeps the force-field names
        assert renamed.get_amino_acid_sequence(numbered=False)[1:4] == [
            "HIE",
            "CYX",
            "ASH",
        ]

    def test_unknown_residue_raises_ssexception(self, GS6_CP):
        renamed = _renamed_protein(GS6_CP, {2: "XXX"})
        with pytest.raises(SSException, match="XXX"):
            renamed.get_amino_acid_sequence(oneletter=True)

    def test_three_letter_sequence_unchanged(self, GS6_CP):
        assert GS6_CP.get_amino_acid_sequence()[:3] == ["ACE-1", "GLY-2", "SER-3"]
        assert GS6_CP.get_amino_acid_sequence(oneletter=True) == "<GSGSGS>"


# ---------------------------------------------------------------------------
# 2. RMS distance map: the std map is the spread of the distance, in Angstroms
# ---------------------------------------------------------------------------
class TestDistanceMapStd:
    def test_rms_std_is_std_of_distance(self, GS6_CP):
        _, std_plain = GS6_CP.get_distance_map(verbose=False)
        rms_map, std_rms = GS6_CP.get_distance_map(RMS=True, verbose=False)
        # the same quantity whether or not RMS is requested
        np.testing.assert_allclose(std_rms, std_plain)

        resids = GS6_CP.resid_with_CA
        for i, j in [(0, 5), (1, 3), (2, 5)]:
            d = GS6_CP.get_inter_residue_atomic_distance(resids[i], resids[j])
            assert std_rms[i, j] == pytest.approx(np.std(d), rel=1e-6)
            assert rms_map[i, j] == pytest.approx(np.sqrt(np.mean(d**2)), rel=1e-6)
            # and it is no longer the std of the squared distance
            assert std_rms[i, j] != pytest.approx(np.std(d**2), rel=1e-3)


# ---------------------------------------------------------------------------
# 3. DSSP refuses a region that includes a cap
# ---------------------------------------------------------------------------
class TestDSSPCapGuard:
    def test_explicit_cap_region_raises(self, GS6_CP):
        assert GS6_CP.ncap and GS6_CP.ccap
        with pytest.raises(SSException, match="without a CA atom"):
            GS6_CP.get_secondary_structure_DSSP(R1=0)
        with pytest.raises(SSException, match="without a CA atom"):
            GS6_CP.get_secondary_structure_DSSP(R2=GS6_CP.n_residues - 1)

    def test_defaults_and_ca_regions_unchanged(self, GS6_CP, NTL9_CP):
        resids, H, E, C = GS6_CP.get_secondary_structure_DSSP()
        assert resids == GS6_CP.resid_with_CA
        assert len(H) == 6
        # explicit CA-bearing endpoints on the capped chain are fine
        resids, H, E, C = GS6_CP.get_secondary_structure_DSSP(R1=1, R2=6)
        assert resids == list(range(1, 7))
        # an uncapped chain may ask for its full range explicitly
        n = NTL9_CP.n_residues
        resids, H, E, C = NTL9_CP.get_secondary_structure_DSSP(R1=0, R2=n - 1)
        assert resids == list(range(n))


# ---------------------------------------------------------------------------
# 4. interchain distance map indexing still matches the intra-chain map
# ---------------------------------------------------------------------------
def test_interchain_map_matches_protein_map(GS6_CO):
    P = GS6_CO.proteinTrajectoryList[0]
    dmap, smap = GS6_CO.get_interchain_distance_map(0, 0)
    intra, intra_std = P.get_distance_map(verbose=False)
    n = len(P.resid_with_CA)
    assert dmap.shape == (n, n)
    iu = np.triu_indices(n, k=1)
    # the two paths differ only in float32 round-off (COM vs compute_distances)
    np.testing.assert_allclose(dmap[iu], intra[iu], rtol=1e-5)
    np.testing.assert_allclose(smap[iu], intra_std[iu], rtol=1e-4)
    np.testing.assert_allclose(dmap, dmap.T)


# ---------------------------------------------------------------------------
# 5. BME diagnostics: the high chi-squared warning uses the reduced chi-squared
# ---------------------------------------------------------------------------
class TestBMEDiagnostics:
    def _poorly_fitting_result(self, n_obs):
        rng = np.random.default_rng(0)
        calc = rng.normal(0.0, 1.0, size=(400, n_obs))
        # every observable sits sqrt(6) sigma away -> reduced chi2 ~ 6
        obs = [ExperimentalObservable(np.sqrt(6.0), 1.0) for _ in range(n_obs)]
        # an enormous theta pins the weights to the prior, so chi2_final ~ 6
        return BME(obs, calc).fit(theta=1e6, auto_theta=False, verbose=False)

    def test_warning_fires_between_two_and_two_m(self):
        # with five observables reduced chi2 ~ 6 lies between 2 and 2*m = 10:
        # the old ``> 2 * m`` threshold stayed silent here
        res = self._poorly_fitting_result(5)
        assert 2.0 < res.chi_squared_final < 10.0
        warnings_ = res.diagnostics()["warnings"]
        assert any("reduced chi-squared" in w for w in warnings_)

    def test_no_warning_for_a_good_fit(self):
        rng = np.random.default_rng(1)
        calc = rng.normal(0.0, 1.0, size=(400, 3))
        obs = [ExperimentalObservable(0.0, 1.0) for _ in range(3)]
        res = BME(obs, calc).fit(theta=1.0, auto_theta=False, verbose=False)
        assert res.chi_squared_final < 2.0
        assert not any("chi-squared" in w for w in res.diagnostics()["warnings"])

    def test_auto_theta_forwards_optimizer_settings(self):
        rng = np.random.default_rng(2)
        calc = rng.normal([24, 65], [3, 6], size=(300, 2))
        obs = [ExperimentalObservable(23.0, 1.0), ExperimentalObservable(60.0, 2.0)]
        bme = BME(obs, calc)
        res = bme.fit(max_iterations=321, optimizer="SLSQP", verbose=False)
        assert res.metadata["optimizer"] == "SLSQP"
        assert res.metadata["max_iterations"] == 321
        assert all(
            r.metadata["optimizer"] == "SLSQP" for r in bme.theta_scan_result.results
        )


# ---------------------------------------------------------------------------
# 6. random-coil shifts: skipped characters are reported
# ---------------------------------------------------------------------------
class TestRandomCoilSkippedCharacters:
    def test_unrecognised_letter_warns_and_is_dropped(self):
        with pytest.warns(UserWarning, match="skipped 1 unrecognised"):
            out = ssnmr.compute_random_coil_chemical_shifts("AXA")
        assert [r["Res"] for r in out] == ["A", "A"]

    def test_unknown_parenthesised_code_warns(self):
        with pytest.warns(UserWarning, match=r"\(FOO\)"):
            ssnmr.compute_random_coil_chemical_shifts("A(FOO)A")

    def test_whitespace_is_silent(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            out = ssnmr.compute_random_coil_chemical_shifts(" AS GA ")
        assert [r["Res"] for r in out] == ["A", "S", "G", "A"]


# ---------------------------------------------------------------------------
# 7. COPER: several groups, feasibility must not depend on the summed chi2
# ---------------------------------------------------------------------------
class TestCOPERMultiGroupFeasibility:
    # six frames whose (x, y) values span the triangle (0,0)-(2,0)-(0,2)
    CALC = np.array(
        [[0.0, 0.0], [2.0, 0.0], [0.0, 2.0], [1.0, 0.0], [0.0, 1.0], [0.5, 0.5]]
    )

    def _fit(self, targets, sigmas):
        obs = [
            ExperimentalObservable(targets[0], sigmas[0], group="a"),
            ExperimentalObservable(targets[1], sigmas[1], group="b"),
        ]
        return COPER(obs, self.CALC).fit(chi2_limit=1.0, verbose=False), obs

    @pytest.mark.parametrize(
        "targets,sigmas",
        [((0.782, 2.160), (0.765, 0.243)), ((1.966, 1.295), (0.787, 0.537))],
    )
    def test_feasible_two_group_problem_is_not_reported_infeasible(
        self, targets, sigmas
    ):
        # for both cases a weight vector with chi2_a, chi2_b < 1 exists, but the
        # minimiser of chi2_a + chi2_b over-fits one group and violates the
        # other, so the summed check alone said "infeasible"
        res, obs = self._fit(targets, sigmas)
        assert res.feasible
        assert res.success
        assert res.metadata["feasibility_check"] == "per-group violation"
        assert res.weights.sum() == pytest.approx(1.0)
        avg = res.weights @ self.CALC
        for k, o in enumerate(obs):
            assert ((avg[k] - o.value) / o.uncertainty) ** 2 <= 1.0 + 1e-4

    def test_genuinely_infeasible_problem_still_reported(self):
        # both targets far outside the reachable triangle
        res, _ = self._fit((5.0, 5.0), (0.1, 0.1))
        assert not res.feasible
        assert not res.success
        assert res.chi_squared_min > 1.0

    def test_single_group_uses_summed_check(self):
        obs = [ExperimentalObservable(0.9, 0.3), ExperimentalObservable(0.9, 0.3)]
        res = COPER(obs, self.CALC).fit(chi2_limit=1.0, verbose=False)
        assert res.feasible
        assert res.metadata["feasibility_check"] == "summed chi-squared"


# ---------------------------------------------------------------------------
# 8. PENGUIN summary figure: one residue axis for all four panels
# ---------------------------------------------------------------------------
class TestQualityPlotResidueAxis:
    def test_uncapped_chain_panels_share_positions(self):
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from soursop.sssampling import SamplingQuality

        xtc = os.path.join(test_data_dir, "ntl9_AA.xtc")
        pdb = os.path.join(test_data_dir, "ntl9_AA.pdb")
        sq = SamplingQuality([xtc, xtc], top_file=pdb, force_sequential=True, seed=0)
        n = sq.trajs[0].proteinTrajectoryList[0].n_residues
        fig, axd = sq.quality_plot()
        try:
            # dihedral panels: interior residues 1..n-2 sit at positions 2..n-1
            np.testing.assert_array_equal(
                axd["A"].lines[0].get_xdata(), np.arange(2, n)
            )
            np.testing.assert_array_equal(
                axd["B"].lines[0].get_xdata(), np.arange(2, n)
            )
            # helicity panel: every residue, positions 1..n
            np.testing.assert_array_equal(
                axd["D"].lines[0].get_xdata(), np.arange(1, n + 1)
            )
            assert axd["A"].get_xlim() == axd["D"].get_xlim() == (0.0, n + 1.0)
        finally:
            plt.close(fig)

    def test_capped_chain_positions_unchanged(self):
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from soursop.sssampling import SamplingQuality

        wt, wp = find_trajectory_files(
            os.path.join(test_data_dir, "sampling_quality/WT"), 3
        )
        sq = SamplingQuality(wt, top_file=wp, force_sequential=True, seed=0)
        fig, axd = sq.quality_plot()
        try:
            np.testing.assert_array_equal(axd["A"].lines[0].get_xdata(), [1, 2])
            np.testing.assert_array_equal(axd["D"].lines[0].get_xdata(), [1, 2])
        finally:
            plt.close(fig)


# ---------------------------------------------------------------------------
# 9. stride validation in the stand-alone HDX / NMR functions
# ---------------------------------------------------------------------------
class TestStandaloneStrideValidation:
    @pytest.mark.parametrize("bad", [0, -1, 2.5, True])
    def test_validate_stride_rejects(self, bad):
        with pytest.raises(SSException):
            ssutils.validate_stride(bad, 10)

    def test_validate_stride_accepts(self):
        assert ssutils.validate_stride(np.int64(3), 10) == 3
        assert ssutils.validate_stride(10, 10) == 10
        with pytest.raises(SSException, match="larger than the number of frames"):
            ssutils.validate_stride(11, 10)

    def test_nmr_and_hdx_reject_zero_stride(self, NTL9_CP):
        with pytest.raises(SSException):
            ssnmr.compute_J3_HN_HA(NTL9_CP, stride=0)
        with pytest.raises(SSException):
            ssnmr.compute_NOE_distances(NTL9_CP, [[0, 5]], stride=0)
        with pytest.raises(SSException):
            sshdx.compute_Nc(NTL9_CP, stride=0)
        with pytest.raises(SSException):
            sshdx.compute_Nh(NTL9_CP, stride=-2)
        # a valid stride still works and subsamples
        _, J = ssnmr.compute_J3_HN_HA(NTL9_CP, stride=2)
        assert J.shape[0] == len(range(0, NTL9_CP.n_frames, 2))


# ---------------------------------------------------------------------------
# 10. odds and ends
# ---------------------------------------------------------------------------
def test_locate_libraries_unsupported_platform(monkeypatch):
    monkeypatch.setattr(platform, "system", lambda: "FreeBSD")
    with pytest.raises(SSException, match="not supported on this platform"):
        ssutils._locate_libraries("openblas")


def test_sspre_repr_reports_angular_frequency(GS6_CP):
    spre = SSPRE(GS6_CP, tau_c=5.0, t_delay=10.0, R_2D=10.0, W_H=2 * np.pi * 600e6)
    assert "rad/s" in repr(spre)


# ---------------------------------------------------------------------------
# 11. plain (non minimum-image) distances in the stand-alone NMR / HDX functions
# ---------------------------------------------------------------------------
class TestNoMinimumImage:
    def test_noe_distances_match_protein_distances(self):
        # gromacs1chain: the chain spans more than half of its 7.4 nm box, so
        # mdtraj's default minimum-image distances were wrong for far pairs
        import mdtraj as md
        from soursop.sstrajectory import SSTrajectory

        T = SSTrajectory(
            os.path.join(test_data_dir, "gromacs1chain/traj.xtc"),
            os.path.join(test_data_dir, "gromacs1chain/top.pdb"),
        )
        P = T.proteinTrajectoryList[0]
        first, last = P.resid_with_CA[0], P.resid_with_CA[-1]
        i, j = P.get_CA_index(first), P.get_CA_index(last)
        d = ssnmr.compute_NOE_distances(P, [[i, j]])[:, 0]
        expected = P.get_inter_residue_atomic_distance(first, last)
        np.testing.assert_allclose(d, expected, rtol=1e-5)
        wrong = 10 * md.compute_distances(P.traj, [[i, j]], periodic=True)[:, 0]
        assert not np.allclose(d, wrong)

    def test_hdx_contacts_ignore_the_box(self):
        import mdtraj as md
        from soursop.ssprotein import SSProtein

        # shrink the box so that minimum imaging would fold the chain back on
        # itself; the counts must not change
        base = SSProtein(
            md.load(
                os.path.join(test_data_dir, "ntl9_AA.xtc"),
                top=os.path.join(test_data_dir, "ntl9_AA.pdb"),
            )
        )
        small = base.traj[:]
        small.unitcell_lengths = np.full((small.n_frames, 3), 1.5, dtype=np.float32)
        small.unitcell_angles = np.full((small.n_frames, 3), 90.0, dtype=np.float32)
        P_small = SSProtein(small)
        for fn in (sshdx.compute_Nc, sshdx.compute_Nh):
            _, a = fn(base)
            _, b = fn(P_small)
            assert np.array_equal(a, b)
        _, a = sshdx.compute_Nh(base, hbond_method="wernet-nilsson")
        _, b = sshdx.compute_Nh(P_small, hbond_method="wernet-nilsson")
        assert np.array_equal(a, b)


# ---------------------------------------------------------------------------
# 12. PRE label-cloud model: coarse-grained fallback and no silent NaN profiles
# ---------------------------------------------------------------------------
class TestPRELabelCloudGuards:
    def _spre(self, protein):
        return SSPRE(protein, tau_c=5.0, t_delay=10.0, R_2D=10.0, W_H=2 * np.pi * 600e6)

    def test_one_bead_chain_uses_the_bead(self, SIGA_CG_CO):
        cg = SIGA_CG_CO.proteinTrajectoryList[0]
        spre = self._spre(cg)
        label = cg.resid_with_CA[5]
        profile, gamma = spre.generate_PRE_profile(label)
        assert len(profile) == len(gamma) == len(cg.resid_with_CA)
        assert np.all(np.isfinite(gamma))
        explicit = spre.generate_PRE_profile(
            label, spin_label_atom="CA", target_relaxation_atom="CA"
        )
        np.testing.assert_allclose(profile, explicit[0])

    def test_bad_inputs_raise_instead_of_nan(self, GS6_CP):
        spre = self._spre(GS6_CP)
        # GS6 is <GSGSGS>: resid_with_CA[1] is a serine, which carries a CB
        label = GS6_CP.resid_with_CA[1]
        assert GS6_CP.get_amino_acid_sequence(oneletter=True)[label] == "S"
        with pytest.raises(SSException, match="n_label_conformers"):
            spre.generate_PRE_profile(label, n_label_conformers=0)
        with pytest.raises(SSException, match="target relaxation atom"):
            spre.generate_PRE_profile(label, target_relaxation_atom="XX")
        with pytest.raises(SSException, match="sterically excluded"):
            spre.generate_PRE_profile(
                label, label_steric="hard", label_bead_radius=60.0
            )
        with pytest.raises(SSException):
            spre.generate_PRE_profile(GS6_CP.n_residues + 3)

    def test_missing_target_is_nan_per_residue(self, GS6_CP):
        # OG exists on serine only, so the glycines come back as nan
        spre = self._spre(GS6_CP)
        label = GS6_CP.resid_with_CA[1]
        profile, gamma = spre.generate_PRE_profile(label, target_relaxation_atom="OG")
        names = GS6_CP.get_amino_acid_sequence(oneletter=True, numbered=False)
        for k, r in enumerate(GS6_CP.resid_with_CA):
            assert np.isnan(gamma[k]) == (names[r] != "S")


# ---------------------------------------------------------------------------
# 13. previously untested public API (values, not just shapes)
# ---------------------------------------------------------------------------
class TestPreviouslyUntestedAPI:
    def test_reweighter_scan_methods_and_predict(self):
        from soursop.ssbme import BMECustom, iBME
        from soursop.sscoper import iCOPER

        rng = np.random.default_rng(3)
        calc = rng.normal([24, 65, 12], [3, 6, 2], size=(400, 3))
        truth = np.array([23.0, 60.0, 11.5])
        obs = [ExperimentalObservable(v, s) for v, s in zip(truth, [1.0, 2.0, 0.5])]

        coper = COPER(obs, calc)
        scan = coper.scan_chi2_limit(chi2_limits=(0.5, 2.0), n_points=4)
        assert len(scan.results) == 4 and 0.5 <= scan.optimal_chi2_limit <= 2.0
        assert coper.scan_result is scan
        scan.print_summary()

        ib = iBME(obs, 1.1 * calc + 0.3)
        res = ib.fit(theta=1.0, verbose=False)
        assert res.success
        pred = ib.predict(calc)
        assert pred.shape == (3,)
        np.testing.assert_allclose(pred, res.weights @ calc)
        scan = ib.scan_theta(theta_range=(0.5, 5.0), n_points=3)
        assert len(scan.results) == 3 and ib.theta_scan_result is scan

        ic = iCOPER(obs, 1.1 * calc + 0.3)
        res = ic.fit(chi2_limit=1.0, verbose=False)
        assert res.feasible
        np.testing.assert_allclose(ic.predict(calc), res.weights @ calc)
        assert res.diagnostics()["status"] in ("OK", "WARNING")
        res.print_diagnostics()
        scan = ic.scan_chi2_limit(chi2_limits=(0.5, 2.0), n_points=3)
        assert len(scan.results) == 3 and ic.scan_result is scan

        bc = BMECustom(truth, calc, uncertainty=np.array([1.0, 2.0, 0.5]))
        scan = bc.scan_theta(theta_range=(0.1, 10.0), n_points=4)
        assert len(scan.results) == 4 and bc.theta_scan_result is scan
        diag = scan.results[-1].diagnostics()
        assert diag["neff_entropy"] <= 400 + 1e-9
        scan.results[-1].print_diagnostics()
        np.testing.assert_allclose(bc.predict(calc), scan.results[-1].weights @ calc)

    def test_get_distance_periodic_minimum_image(self):
        from soursop.sstools import chunks, get_distance_periodic

        a = np.array([[0.1, 0.1, 0.1], [1.0, 1.0, 1.0]])
        b = np.array([[9.9, 9.9, 9.9], [1.0, 1.0, 2.0]])
        d = get_distance_periodic(a, b, 10.0)
        np.testing.assert_allclose(d, [np.sqrt(3 * 0.2**2), 1.0])
        with pytest.raises(SSException):
            get_distance_periodic(a, b[:1], 10.0)
        with pytest.raises(SSException):
            get_distance_periodic(a, b, 10.0, box_shape="sphere")
        assert list(chunks([1, 2, 3, 4, 5], 2)) == [[1, 2], [3, 4]]

    def test_residue_mass_and_atom_indices(self, GS6_CP):
        top = GS6_CP.topology
        for r in GS6_CP.resid_with_CA:
            expected = sum(a.element.mass for a in top.residue(r).atoms)
            assert GS6_CP.get_residue_mass(r) == pytest.approx(expected)
            assert list(GS6_CP.get_all_atomic_indices(r)) == [
                a.index for a in top.residue(r).atoms
            ]
        names = GS6_CP.get_amino_acid_sequence(oneletter=True, numbered=False)
        masses = {names[r]: GS6_CP.get_residue_mass(r) for r in GS6_CP.resid_with_CA}
        assert masses["G"] < masses["S"]
        with pytest.raises(SSException):
            GS6_CP.get_residue_mass(GS6_CP.n_residues)

    def test_get_clusters_partitions_all_frames(self, NTL9_CP):
        members, trajs, dmats, centroids, frames = NTL9_CP.get_clusters(
            n_clusters=3, stride=1
        )
        assert sum(members) == NTL9_CP.n_frames
        assert sorted(np.concatenate(frames).tolist()) == list(range(NTL9_CP.n_frames))
        for m, t, dm, c, f in zip(members, trajs, dmats, centroids, frames):
            assert t.n_frames == m == len(f) == dm.shape[0] == dm.shape[1]
            assert 0 <= c < m
            # mdtraj's aligned self-RMSD is float32 noise (~0.01 A), not exactly 0
            np.testing.assert_allclose(np.diag(dm), 0.0, atol=0.05)
            np.testing.assert_allclose(dm, dm.T, atol=0.05)
        with pytest.raises(SSException):
            NTL9_CP.get_clusters(stride=0)


# ---------------------------------------------------------------------------
# 14. ssprotein fixes from the cross-validation pass
# ---------------------------------------------------------------------------
class TestSecondaryStructureRegions:
    def test_bbseg_region_keeps_its_endpoints(self, NTL9_CP, CTL9_CP):
        # the topology used to be sliced to the region before phi/psi were
        # computed, which lost phi(R1) and psi(R2)
        resids, classes = NTL9_CP.get_secondary_structure_BBSEG(R1=6, R2=16)
        assert resids == list(range(6, 17))
        all_resids, all_classes = NTL9_CP.get_secondary_structure_BBSEG()
        assert all_resids == list(range(1, NTL9_CP.n_residues - 1))
        for c in range(9):
            np.testing.assert_allclose(classes[c], all_classes[c][5:16])
        # on a capped chain every real residue is classified
        resids, _ = CTL9_CP.get_secondary_structure_BBSEG(R1=1, R2=92)
        assert resids == list(range(1, 93))

    def test_dssp_region_is_a_slice_of_the_whole_chain(self, CTL9_CP):
        resids, H, E, C = CTL9_CP.get_secondary_structure_DSSP()
        sub_resids, sH, sE, sC = CTL9_CP.get_secondary_structure_DSSP(R1=10, R2=40)
        assert sub_resids == list(range(10, 41))
        np.testing.assert_allclose(sH, H[9:40])
        np.testing.assert_allclose(sE, E[9:40])
        np.testing.assert_allclose(sC, C[9:40])


class TestContactMapResidues:
    def test_for_nh2_caps_are_excluded(self):
        from soursop.sstrajectory import SSTrajectory

        pdb = os.path.join(test_data_dir, "cap_tests/FOR_NH2.pdb")
        P = SSTrajectory(pdb, pdb).proteinTrajectoryList[0]
        assert P.ncap and P.ccap
        cmap, corder = P.get_contact_map()
        n = len(P.resid_with_CA)
        assert cmap.shape == (n, n) and corder.shape == (n,)

    def test_one_bead_sidechain_modes_raise(self, SIGA_CG_CO):
        cg = SIGA_CG_CO.proteinTrajectoryList[0]
        with pytest.raises(SSException):
            cg.get_contact_map(mode="sidechain")
        cmap, _ = cg.get_contact_map(mode="ca", stride=10)
        assert cmap.shape == (cg.n_residues, cg.n_residues)


class TestNygaardHydrodynamicRadius:
    def test_uses_calpha_rg_and_residue_count(self, GS6_CP, NTL9_CP):
        import mdtraj as md

        for P in (GS6_CP, NTL9_CP):
            ca = [P.get_CA_index(r) for r in P.resid_with_CA]
            rg = 10 * md.compute_rg(P.traj.atom_slice(ca))
            N = len(P.resid_with_CA)
            n033, n060 = N**0.33, N**0.60
            expected = rg / (((0.216 * (rg - 4.06 * n033)) / (n060 - n033)) + 0.821)
            np.testing.assert_allclose(P.get_hydrodynamic_radius(), expected)
        # GS6: 6 CA-bearing residues (caps no longer counted), CA-only Rg
        rh = GS6_CP.get_hydrodynamic_radius()
        assert rh[0] == pytest.approx(13.7213, abs=1e-3)
        assert np.mean(rh) == pytest.approx(13.8392, abs=1e-3)
        # explicit cap-inclusive endpoints give the same answer as the default
        np.testing.assert_allclose(
            GS6_CP.get_hydrodynamic_radius(R1=0, R2=GS6_CP.n_residues - 1), rh
        )


def test_gyration_tensor_trace_is_rg_squared(NTL9_CP, GS6_CP):
    for P in (NTL9_CP, GS6_CP):
        T = P.get_gyration_tensor(verbose=False)
        rg = P.get_radius_of_gyration()
        np.testing.assert_allclose(np.trace(T, axis1=1, axis2=2), rg**2, rtol=1e-6)
        T = P.get_gyration_tensor(R1=1, R2=3, verbose=False)
        np.testing.assert_allclose(
            np.trace(T, axis1=1, axis2=2),
            P.get_radius_of_gyration(1, 3) ** 2,
            rtol=1e-6,
        )


def test_local_heterogeneity_default_bins_cover_the_data(CTL9_CP):
    mean, std, histo, bins = CTL9_CP.get_local_heterogeneity(
        fragment_size=40, stride=200, verbose=False
    )
    assert bins[-1] > 99.0
    n_ref = len(range(0, CTL9_CP.n_frames, 200))
    # every pooled RMSD value lands in the histogram
    assert all(h.sum() == n_ref * CTL9_CP.n_frames for h in histo)
    assert max(mean) > 10.0  # the old 0-10 A default would have dropped these


def test_sasa_terminal_groups_are_backbone(NTL9_CP):
    import mdtraj as md

    residue, sidechain, backbone = NTL9_CP.get_all_SASA(stride=5, mode="all")
    np.testing.assert_allclose(sidechain + backbone, residue, atol=1e-3)
    # independent split of the per-atom SASA for the N-terminal methionine
    top = NTL9_CP.topology
    areas = 100 * md.shrake_rupley(NTL9_CP.traj[::5], mode="atom")
    res0 = top.residue(NTL9_CP.resid_with_CA[0])
    bb_names = {"N", "CA", "C", "O", "H", "HA", "H1", "H2", "H3"}
    sc_atoms = [a.index for a in res0.atoms if a.name not in bb_names]
    np.testing.assert_allclose(
        sidechain[:, 0], areas[:, sc_atoms].sum(axis=1), rtol=1e-5
    )
    # and the C-terminal carboxylate oxygens are backbone, not sidechain
    last = top.residue(NTL9_CP.resid_with_CA[-1])
    oxt = [
        a.index for a in last.atoms if a.name in ("OXT", "1OXT", "2OXT", "OT1", "OT2")
    ]
    assert len(oxt) > 0
    sc_last = [
        a.index
        for a in last.atoms
        if a.name not in {"N", "CA", "C", "O", "H", "HA"} and a.index not in oxt
    ]
    np.testing.assert_allclose(
        sidechain[:, -1], areas[:, sc_last].sum(axis=1), rtol=1e-5
    )


class TestSilentGarbageNowRaises:
    def test_q_without_native_contacts(self, GS6_CP, NTL9_CP):
        with pytest.raises(SSException, match="no heavy-atom contacts"):
            GS6_CP.get_Q()
        with pytest.raises(SSException, match="native_state_reference_frame"):
            NTL9_CP.get_Q(native_state_reference_frame=NTL9_CP.n_frames)
        assert np.all(np.isfinite(NTL9_CP.get_Q()))

    def test_residue_com_with_missing_atom(self, NTL9_CP):
        with pytest.raises(SSException, match="XX"):
            NTL9_CP.get_residue_COM(3, atom_name="XX")

    def test_one_bead_dihedrals_raise(self, SIGA_CG_CO):
        cg = SIGA_CG_CO.proteinTrajectoryList[0]
        with pytest.raises(SSException):
            cg.get_angles("phi")
        with pytest.raises(SSException):
            cg.get_dihedral_mutual_information()

    def test_rmsd_frame_validation(self, NTL9_CP):
        n = NTL9_CP.n_frames
        with pytest.raises(SSException):
            NTL9_CP.get_RMSD(n + 3)
        with pytest.raises(SSException):
            NTL9_CP.get_RMSD(0, frame2=n + 3)
        with pytest.raises(SSException):
            NTL9_CP.get_RMSD(-1)
        assert NTL9_CP.get_RMSD(0, frame2=n - 1).shape == (1,)

    def test_region_and_resid_validation(self, NTL9_CP):
        with pytest.raises(SSException):
            NTL9_CP.get_regional_SASA(50, 70)
        with pytest.raises(SSException):
            NTL9_CP.get_regional_SASA(10, 3)
        with pytest.raises(SSException):
            NTL9_CP.get_site_accessibility([NTL9_CP.n_residues], mode="resid")
        with pytest.raises(SSException):
            NTL9_CP.get_inter_residue_atomic_distance(0, -1)
        with pytest.raises(SSException):
            NTL9_CP.get_inter_residue_atomic_distance(
                0, NTL9_CP.n_residues, mode="closest"
            )
        with pytest.raises(SSException):
            NTL9_CP.get_clusters(stride=0)

    def test_non_integer_stride_raises(self, NTL9_CP):
        for call in (
            lambda: NTL9_CP.get_distance_map(stride=2.5, verbose=False),
            lambda: NTL9_CP.get_RMSD(0, stride=2.5),
            lambda: NTL9_CP.get_contact_map(stride=2.5),
            lambda: NTL9_CP.get_all_SASA(stride=2.5),
            lambda: NTL9_CP.get_inter_residue_COM_distance(0, 5, stride=2.5),
        ):
            with pytest.raises(SSException):
                call()
        # numpy integers are fine
        assert NTL9_CP.get_inter_residue_COM_distance(
            0, 5, stride=np.int64(2)
        ).shape == (5,)

    def test_angle_decay_on_two_bead_chain(self, SWAN_HELIX_CP):
        with pytest.raises(SSException, match="get_angle_decay"):
            SWAN_HELIX_CP.get_angle_decay()

    def test_sidechain_alignment_self_angle_is_finite(self, CTL9_CP):
        a = CTL9_CP.get_sidechain_alignment_angle(5, 5)
        assert np.all(np.isfinite(a))
        assert np.max(a) < 0.1


# ---------------------------------------------------------------------------
# 15. reweighting robustness from the cross-validation pass
# ---------------------------------------------------------------------------
class TestReweightingInputs:
    def _data(self):
        rng = np.random.default_rng(4)
        calc = rng.normal([24, 65, 12], [3, 6, 2], size=(300, 3))
        obs = [
            ExperimentalObservable(23.0, 1.0),
            ExperimentalObservable(60.0, 2.0),
            ExperimentalObservable(11.5, 0.5),
        ]
        return obs, calc

    def test_negative_priors_rejected_everywhere(self):
        from soursop.ssbme import BMECustom, iBME
        from soursop.sscoper import iCOPER

        obs, calc = self._data()
        w0 = np.ones(300)
        w0[:5] = -0.5
        for cls in (BME, iBME, COPER, iCOPER):
            with pytest.raises(SSException, match="non-negative"):
                cls(obs, calc, w0)
        with pytest.raises(SSException, match="non-negative"):
            BMECustom(np.array([23.0, 60.0, 11.5]), calc, initial_weights=w0)
        with pytest.raises(SSException, match="all be zero"):
            BME(obs, calc, np.zeros(300))
        with pytest.raises(SSException, match="finite"):
            iCOPER(obs, calc, np.full(300, np.nan))

    def test_two_observable_offset_fit_warns(self):
        from soursop.ssbme import iBME
        from soursop.sscoper import iCOPER

        obs, calc = self._data()
        with pytest.warns(UserWarning, match="exactly determined"):
            iBME(obs[:2], calc[:, :2]).fit(theta=1.0, verbose=False)
        with pytest.warns(UserWarning, match="exactly determined"):
            iCOPER(obs[:2], calc[:, :2]).fit(chi2_limit=1.0, verbose=False)

    def test_icoper_reports_zero_prior_frames(self):
        from soursop.sscoper import iCOPER

        obs, calc = self._data()
        w0 = np.ones(300)
        w0[:10] = 0.0
        res = iCOPER(obs, calc, w0).fit(chi2_limit=1.0, verbose=False)
        assert res.metadata["n_zero_prior_frames"] == 10
        assert np.all(res.weights[:10] == 0.0)

    def test_theta_scan_ignores_failed_fits(self):
        from soursop.ssbme import theta_scan

        obs, calc = self._data()
        # a single optimizer iteration cannot converge: every fit fails
        with pytest.raises(SSException, match="every fit failed"):
            theta_scan(obs, calc, n_points=4, fit_kwargs={"max_iterations": 1})
        scan = theta_scan(obs, calc, n_points=4)
        assert all(r.success for r in scan.results)

    def test_chi2_limit_scan_prefers_feasible_limits(self):
        from soursop.sscoper import chi2_limit_scan

        obs, calc = self._data()
        # a target far outside the sampled range: the small limits cannot be met
        obs[0] = ExperimentalObservable(50.0, 1.0)
        with pytest.warns(UserWarning, match="infeasible"):
            scan = chi2_limit_scan(
                obs, calc, chi2_limits=np.array([0.25, 1.0, 1e3, 1e4])
            )
        assert scan.feasible_mask.tolist() == [False, False, True, True]
        assert scan.feasible_mask[scan.optimal_idx]
        # nothing feasible at all: still returns, with a warning
        with pytest.warns(UserWarning, match="no scanned"):
            scan = chi2_limit_scan(obs, calc, chi2_limits=np.array([0.25, 1.0]))
        assert not scan.feasible_mask.any()


def test_sampling_quality_rejects_mismatched_reference_chain(tmp_path):
    from soursop.sssampling import SamplingQuality
    from soursop.sstools import find_trajectory_files

    wt, wp = find_trajectory_files(
        os.path.join(test_data_dir, "sampling_quality/WT"), 3
    )
    ntl9_xtc = os.path.join(test_data_dir, "ntl9_AA.xtc")
    ntl9_pdb = os.path.join(test_data_dir, "ntl9_AA.pdb")
    with pytest.raises(SSException, match="same chain"):
        SamplingQuality(
            wt[:1],
            [ntl9_xtc],
            top_file=wp[0],
            ref_top=ntl9_pdb,
            truncate=True,
            force_sequential=True,
            verbose=False,
        )


# ---------------------------------------------------------------------------
# 16. contact map: plain distances, never minimum-image
# ---------------------------------------------------------------------------
def test_contact_map_ignores_the_box():
    import mdtraj as md
    from soursop.sstrajectory import SSTrajectory

    # gromacs1chain spans more than half of its 7.4 nm box, so minimum-image
    # distances differ from the plain ones for far-apart pairs
    T = SSTrajectory(
        os.path.join(test_data_dir, "gromacs1chain/traj.xtc"),
        os.path.join(test_data_dir, "gromacs1chain/top.pdb"),
    )
    P = T.proteinTrajectoryList[0]
    first, last = P.resid_with_CA[0], P.resid_with_CA[-1]
    cmap, _ = P.get_contact_map(distance_thresh=45.0, mode="closest-heavy")
    d_plain = P.get_inter_residue_atomic_distance(first, last, mode="closest-heavy")
    d_image = (
        10
        * md.compute_contacts(
            P.traj, [[first, last]], scheme="closest-heavy", periodic=True
        )[0][:, 0]
    )
    assert cmap[0, -1] == pytest.approx(np.mean(d_plain < 45.0))
    assert cmap[0, -1] != pytest.approx(np.mean(d_image < 45.0))
    # and the whole map is what a box-free copy of the trajectory gives
    nobox = P.traj[:]
    nobox.unitcell_vectors = None
    from soursop.ssprotein import SSProtein

    cmap_nobox, _ = SSProtein(nobox).get_contact_map(
        distance_thresh=45.0, mode="closest-heavy"
    )
    np.testing.assert_allclose(cmap, cmap_nobox)
