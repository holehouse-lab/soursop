# Import package, test suite, and other packages as needed
import numpy as np
import pytest

import soursop.ssnmr as nmr
from soursop.ssexceptions import SSException


#
# All test code written by Alex Keeley
#


def test_simple_test():
    output = nmr.compute_random_coil_chemical_shifts(
        "ACDE",
        temperature=15,
        pH=6.5,
        use_ggxgg=True,
        use_perdeuteration=False,
        asFloat=True,
    )
    control = [
        {
            "Res": "A",
            "Index": 0,
            "CA": 52.695,
            "CB": 18.967,
            "CO": 177.864,
            "N": 126.057,
            "HN": 8.519,
            "HA": 4.357,
        },
        {
            "Res": "C",
            "Index": 1,
            "CA": 58.584,
            "CB": 29.752,
            "CO": 174.436,
            "N": 118.754,
            "HN": 8.466,
            "HA": 4.500,
        },
        {
            "Res": "D",
            "Index": 2,
            "CA": 54.701,
            "CB": 40.940,
            "CO": 176.465,
            "N": 122.592,
            "HN": 8.476,
            "HA": 4.609,
        },
        {
            "Res": "E",
            "Index": 3,
            "CA": 56.941,
            "CB": 30.167,
            "CO": 176.755,
            "N": 121.337,
            "HN": 8.426,
            "HA": 4.263,
        },
    ]

    for i in range(4):
        amino_acid = output[i]
        assert amino_acid == control[i]


def test_simple_types():
    """Tests the ability of the function randCoilChemShifts to accurately predict random coil chemical shifts for all 20 amino acid types."""

    output = nmr.compute_random_coil_chemical_shifts(
        "ACDEFGHIKLMNPQRSTVWY",
        temperature=25,
        pH=6.5,
        use_ggxgg=True,
        use_perdeuteration=False,
        asFloat=True,
    )
    control = [
        {
            "Res": "A",
            "Index": 0,
            "CA": 52.673,
            "CB": 19.015,
            "CO": 177.793,
            "N": 126.005,
            "HN": 8.429,
            "HA": 4.364,
        },
        {
            "Res": "C",
            "Index": 1,
            "CA": 58.575,
            "CB": 29.764,
            "CO": 174.410,
            "N": 118.672,
            "HN": 8.396,
            "HA": 4.500,
        },
        {
            "Res": "D",
            "Index": 2,
            "CA": 54.896,
            "CB": 40.919,
            "CO": 176.273,
            "N": 122.765,
            "HN": 8.403,
            "HA": 4.540,
        },
        {
            "Res": "E",
            "Index": 3,
            "CA": 56.873,
            "CB": 30.051,
            "CO": 176.341,
            "N": 120.831,
            "HN": 8.335,
            "HA": 4.225,
        },
        {
            "Res": "F",
            "Index": 4,
            "CA": 57.900,
            "CB": 39.333,
            "CO": 176.236,
            "N": 119.924,
            "HN": 8.186,
            "HA": 4.637,
        },
        {
            "Res": "G",
            "Index": 5,
            "CA": 45.364,
            "CB": "**.***",
            "CO": 173.150,
            "N": 112.203,
            "HN": 8.145,
            "HA": 3.909,
        },
        {
            "Res": "H",
            "Index": 6,
            "CA": 55.775,
            "CB": 29.978,
            "CO": 174.752,
            "N": 118.609,
            "HN": 8.137,
            "HA": 4.631,
        },
        {
            "Res": "I",
            "Index": 7,
            "CA": 61.170,
            "CB": 38.698,
            "CO": 175.946,
            "N": 122.883,
            "HN": 8.143,
            "HA": 4.135,
        },
        {
            "Res": "K",
            "Index": 8,
            "CA": 56.199,
            "CB": 32.866,
            "CO": 176.184,
            "N": 126.438,
            "HN": 8.474,
            "HA": 4.309,
        },
        {
            "Res": "L",
            "Index": 9,
            "CA": 55.223,
            "CB": 42.321,
            "CO": 177.352,
            "N": 124.735,
            "HN": 8.114,
            "HA": 4.339,
        },
        {
            "Res": "M",
            "Index": 10,
            "CA": 55.397,
            "CB": 32.878,
            "CO": 175.653,
            "N": 121.485,
            "HN": 8.400,
            "HA": 4.502,
        },
        {
            "Res": "N",
            "Index": 11,
            "CA": 51.159,
            "CB": 38.132,
            "CO": 173.163,
            "N": 121.161,
            "HN": 8.416,
            "HA": 4.957,
        },
        {
            "Res": "P",
            "Index": 12,
            "CA": 63.418,
            "CB": 32.272,
            "CO": 176.884,
            "N": "***.***",
            "HN": "*.***",
            "HA": 4.423,
        },
        {
            "Res": "Q",
            "Index": 13,
            "CA": 55.989,
            "CB": 29.081,
            "CO": 176.161,
            "N": 120.570,
            "HN": 8.503,
            "HA": 4.305,
        },
        {
            "Res": "R",
            "Index": 14,
            "CA": 56.205,
            "CB": 30.960,
            "CO": 176.358,
            "N": 123.141,
            "HN": 8.392,
            "HA": 4.407,
        },
        {
            "Res": "S",
            "Index": 15,
            "CA": 58.326,
            "CB": 63.811,
            "CO": 174.642,
            "N": 117.860,
            "HN": 8.452,
            "HA": 4.549,
        },
        {
            "Res": "T",
            "Index": 16,
            "CA": 62.004,
            "CB": 69.725,
            "CO": 174.464,
            "N": 116.881,
            "HN": 8.182,
            "HA": 4.350,
        },
        {
            "Res": "V",
            "Index": 17,
            "CA": 62.918,
            "CB": 32.654,
            "CO": 175.762,
            "N": 122.829,
            "HN": 8.103,
            "HA": 3.984,
        },
        {
            "Res": "W",
            "Index": 18,
            "CA": 57.244,
            "CB": 29.237,
            "CO": 175.644,
            "N": 124.340,
            "HN": 8.098,
            "HA": 4.638,
        },
        {
            "Res": "Y",
            "Index": 19,
            "CA": 57.632,
            "CB": 38.892,
            "CO": 174.936,
            "N": 122.232,
            "HN": 7.732,
            "HA": 4.438,
        },
    ]

    for i in range(20):
        amino_acid = output[i]
        assert amino_acid == control[i]


def test_simple_deut():
    """Tests the ability of the function randCoilChemShifts to accurately predict random coil chemical shifts for perdeuterated proteins."""
    output = nmr.compute_random_coil_chemical_shifts(
        "ACDE",
        temperature=25,
        pH=6.5,
        use_ggxgg=True,
        use_perdeuteration=True,
        asFloat=True,
    )
    control = [
        {
            "Res": "A",
            "Index": 0,
            "CA": 51.993,
            "CB": 18.015,
            "CO": 177.793,
            "N": 126.005,
            "HN": 8.429,
            "HA": "*.***",
        },
        {
            "Res": "C",
            "Index": 1,
            "CA": 58.025,
            "CB": 29.054,
            "CO": 174.410,
            "N": 118.672,
            "HN": 8.396,
            "HA": "*.***",
        },
        {
            "Res": "D",
            "Index": 2,
            "CA": 54.178,
            "CB": 40.295,
            "CO": 176.417,
            "N": 122.553,
            "HN": 8.414,
            "HA": "*.***",
        },
        {
            "Res": "E",
            "Index": 3,
            "CA": 56.260,
            "CB": 29.243,
            "CO": 176.706,
            "N": 121.300,
            "HN": 8.362,
            "HA": "*.***",
        },
    ]

    for i in range(4):
        amino_acid = output[i]
        assert amino_acid == control[i]


def test_simple_undeut():
    """Tests the ability of the function randCoilChemShifts to accurately predict random coil chemical shifts for unperdeuterated proteins"""
    output = nmr.compute_random_coil_chemical_shifts(
        "ACDE",
        temperature=25,
        pH=6.5,
        use_ggxgg=True,
        use_perdeuteration=False,
        asFloat=True,
    )

    control = [
        {
            "Res": "A",
            "Index": 0,
            "CA": 52.673,
            "CB": 19.015,
            "CO": 177.793,
            "N": 126.005,
            "HN": 8.429,
            "HA": 4.364,
        },
        {
            "Res": "C",
            "Index": 1,
            "CA": 58.575,
            "CB": 29.764,
            "CO": 174.410,
            "N": 118.672,
            "HN": 8.396,
            "HA": 4.500,
        },
        {
            "Res": "D",
            "Index": 2,
            "CA": 54.728,
            "CB": 41.005,
            "CO": 176.417,
            "N": 122.553,
            "HN": 8.414,
            "HA": 4.608,
        },
        {
            "Res": "E",
            "Index": 3,
            "CA": 56.950,
            "CB": 30.213,
            "CO": 176.706,
            "N": 121.300,
            "HN": 8.362,
            "HA": 4.266,
        },
    ]

    for i in range(4):
        amino_acid = output[i]
        assert amino_acid == control[i]


def test_simple_phos():
    """Tests the ability of the function randCoilChemShifts to accurately predict random coil chemical shifts for phosphorylated amino acids."""

    output = nmr.compute_random_coil_chemical_shifts(
        "(sep)(tpo)(ptr)",
        temperature=25,
        pH=6.5,
        use_ggxgg=True,
        use_perdeuteration=False,
        asFloat=True,
    )

    control = [
        {
            "Res": "SEP",
            "Index": 0,
            "CA": 58.483,
            "CB": 65.876,
            "CO": 174.762,
            "N": 119.143,
            "HN": 9.005,
            "HA": 4.440,
        },
        {
            "Res": "TPO",
            "Index": 1,
            "CA": 62.641,
            "CB": 71.984,
            "CO": 174.201,
            "N": 117.185,
            "HN": 8.913,
            "HA": 4.261,
        },
        {
            "Res": "PTR",
            "Index": 2,
            "CA": 57.802,
            "CB": 38.635,
            "CO": 175.690,
            "N": 122.211,
            "HN": 8.087,
            "HA": 4.585,
        },
    ]

    for i in range(3):
        amino_acid = output[i]
        assert amino_acid == control[i]


def test_simple_ph():
    """Tests the ability of the function randCoilChemShifts to accurately predict random coil chemical shifts in different pH conditions."""

    output = nmr.compute_random_coil_chemical_shifts(
        "ACDE",
        temperature=25,
        pH=4,
        use_ggxgg=True,
        use_perdeuteration=False,
        asFloat=True,
    )

    control = [
        {
            "Res": "A",
            "Index": 0,
            "CA": 52.673,
            "CB": 19.015,
            "CO": 177.793,
            "N": 126.005,
            "HN": 8.429,
            "HA": 4.364,
        },
        {
            "Res": "C",
            "Index": 1,
            "CA": 58.575,
            "CB": 29.764,
            "CO": 174.410,
            "N": 118.672,
            "HN": 8.396,
            "HA": 4.500,
        },
        {
            "Res": "D",
            "Index": 2,
            "CA": 53.868,
            "CB": 39.252,
            "CO": 175.702,
            "N": 121.615,
            "HN": 8.476,
            "HA": 4.675,
        },
        {
            "Res": "E",
            "Index": 3,
            "CA": 56.281,
            "CB": 29.080,
            "CO": 176.324,
            "N": 120.696,
            "HN": 8.281,
            "HA": 4.340,
        },
    ]

    for i in range(4):
        amino_acid = output[i]
        assert amino_acid == control[i]


# ----------
def test_simple_temp():
    """Tests the ability of the function randCoilChemShifts to accurately predict random coil chemical shifts in different temperature conditions."""

    output = nmr.compute_random_coil_chemical_shifts(
        "ACDE",
        temperature=15,
        pH=6.5,
        use_ggxgg=True,
        use_perdeuteration=False,
        asFloat=True,
    )

    control = [
        {
            "Res": "A",
            "Index": 0,
            "CA": 52.695,
            "CB": 18.967,
            "CO": 177.864,
            "N": 126.057,
            "HN": 8.519,
            "HA": 4.357,
        },
        {
            "Res": "C",
            "Index": 1,
            "CA": 58.584,
            "CB": 29.752,
            "CO": 174.436,
            "N": 118.754,
            "HN": 8.466,
            "HA": 4.500,
        },
        {
            "Res": "D",
            "Index": 2,
            "CA": 54.701,
            "CB": 40.940,
            "CO": 176.465,
            "N": 122.592,
            "HN": 8.476,
            "HA": 4.609,
        },
        {
            "Res": "E",
            "Index": 3,
            "CA": 56.941,
            "CB": 30.167,
            "CO": 176.755,
            "N": 121.337,
            "HN": 8.426,
            "HA": 4.263,
        },
    ]

    for i in range(4):
        amino_acid = output[i]
        assert amino_acid == control[i]


# ----------
def test_unit_round3():
    """Tests the ability of the round3 function to accurately and consistently round input numbers to exactly 3 decimal places."""
    raw = [
        110.7000,
        103.9111,
        112.7382,
        114.8413,
        112.5684,
        107.9585,
        105.9196,
        108.2177,
        112.5628,
        112.0159,
    ]

    control = [
        110.700,
        103.911,
        112.738,
        114.841,
        112.568,
        107.959,
        105.920,
        108.218,
        112.563,
        112.016,
    ]

    for i in range(10):
        assert nmr.__round3(raw[i], asFloat=True) == control[i]


# ----------
def test_unit_setseq():
    """Tests the ability of the set_sequence function to accurately and consistently translate a sequence of amino acid abbreviations into a corresponding sequence of representative numbers."""
    key_aa1 = [
        0,
        -1,
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        -1,
        8,
        9,
        10,
        11,
        -1,
        12,
        13,
        14,
        15,
        16,
        -1,
        17,
        18,
        -1,
        19,
        -1,
    ]
    key_aa3 = {
        "ALA": 0,
        "CYS": 1,
        "ASP": 2,
        "GLU": 3,
        "PHE": 4,
        "GLY": 5,
        "HIS": 6,
        "ILE": 7,
        "LYS": 8,
        "LEU": 9,
        "MET": 10,
        "ASN": 11,
        "PRO": 12,
        "GLN": 13,
        "ARG": 14,
        "SER": 15,
        "THR": 16,
        "VAL": 17,
        "TRP": 18,
        "TYR": 19,
        "PSER": 20,
        "SEP": 20,
        "PS": 20,
        "PTHR": 21,
        "TPO": 21,
        "PT": 21,
        "PTYR": 22,
        "PTR": 22,
        "PY": 22,
    }

    raw = [
        "KEQEERL",
        "EGHVTPC",
        "GQGDHQD",
        "MHDTMII",
        "FTWFSGH",
        "KLD(pT)GQK",
        "L(pY)LQATYG",
        "QYQES(pS)RP",
        "WWGN(pS)TGI",
        "KHRDLAL(pY)",
    ]

    control = [
        [23, 23, 8, 3, 13, 3, 3, 14, 9, 23, 23],
        [23, 23, 3, 5, 6, 17, 16, 12, 1, 23, 23],
        [23, 23, 5, 13, 5, 2, 6, 13, 2, 23, 23],
        [23, 23, 10, 6, 2, 16, 10, 7, 7, 23, 23],
        [23, 23, 4, 16, 18, 4, 15, 5, 6, 23, 23],
        [23, 23, 8, 9, 2, 21, 5, 13, 8, 23, 23],
        [23, 23, 9, 22, 9, 13, 0, 16, 19, 5, 23, 23],
        [23, 23, 13, 19, 13, 3, 15, 20, 14, 12, 23, 23],
        [23, 23, 18, 18, 5, 11, 20, 16, 5, 7, 23, 23],
        [23, 23, 8, 6, 14, 2, 9, 0, 9, 22, 23, 23],
    ]

    for i in range(10):
        seq = nmr.__set_sequence(raw[i], key_aa1, key_aa3)
        assert seq[0] == control[i]


# --------------------------------------------------------------------------
# Scalar J-couplings (Karplus relation)
# --------------------------------------------------------------------------
class TestKarplusEvaluator:
    """The generic ``karplus`` Karplus-relation evaluator."""

    def test_phi60_bax2007_closed_form(self):
        """At phi=60 with Bax2007 the cosine argument is exactly 0, so
        J = A + B + C = 8.40 - 1.36 + 0.33 = 7.37 Hz exactly."""
        coef = nmr.KARPLUS_HN_HA_COEFFICIENTS["Bax2007"]
        assert nmr.karplus(60.0, **coef) == pytest.approx(7.37, abs=1e-9)

    def test_array_broadcast(self):
        """Array input preserves shape and gives the right closed-form values."""
        coef = nmr.KARPLUS_HN_HA_COEFFICIENTS["Bax2007"]
        phi = np.array([-60.0, 0.0, 60.0, 120.0, 180.0])
        # cos(phi - 60 deg) at these points is {-0.5, 0.5, 1, 0.5, -0.5}
        # so J = 8.4*c^2 - 1.36*c + 0.33
        expected = np.array(
            [
                8.4 * 0.25 - 1.36 * (-0.5) + 0.33,  # -60 -> 3.11
                8.4 * 0.25 - 1.36 * (0.5) + 0.33,  #   0 -> 1.75
                8.4 * 1.00 - 1.36 * (1.0) + 0.33,  #  60 -> 7.37
                8.4 * 0.25 - 1.36 * (0.5) + 0.33,  # 120 -> 1.75
                8.4 * 0.25 - 1.36 * (-0.5) + 0.33,  # 180 -> 3.11
            ]
        )
        out = nmr.karplus(phi, **coef)
        assert out.shape == phi.shape
        assert np.allclose(out, expected)

    def test_default_phi0_is_zero(self):
        """With phi0=0, A=1, B=0, C=0, karplus reduces to cos^2."""
        v = nmr.karplus(np.array([0.0, 90.0, 180.0]), A=1.0, B=0.0, C=0.0)
        assert np.allclose(v, np.array([1.0, 0.0, 1.0]))


class TestKarplusTables:
    """Coefficient and uncertainty tables ported from biceps."""

    def test_all_six_models_present(self):
        expected = {"Ruterjans1999", "Bax2007", "Bax1997", "Habeck", "Vuister", "Pardi"}
        assert set(nmr.KARPLUS_HN_HA_COEFFICIENTS.keys()) == expected
        assert set(nmr.KARPLUS_HN_HA_UNCERTAINTIES.keys()) == expected

    def test_default_model_is_in_table(self):
        assert nmr.KARPLUS_HN_HA_DEFAULT_MODEL in nmr.KARPLUS_HN_HA_COEFFICIENTS

    def test_each_entry_has_full_coefficient_set(self):
        for name, coef in nmr.KARPLUS_HN_HA_COEFFICIENTS.items():
            assert set(coef.keys()) == {"phi0", "A", "B", "C"}, name
            for k, v in coef.items():
                assert isinstance(v, float), f"{name}.{k} is not float"


class TestComputeJ3HNHA:
    """``compute_J3_HN_HA`` operating on a real SSProtein fixture."""

    def test_shape_and_range(self, GS6_CP):
        atoms, J = nmr.compute_J3_HN_HA(GS6_CP)
        n_frames = GS6_CP.n_frames
        n_phi = len(atoms)
        # GS6 is a short protein, so n_phi = n_residues - 1 (no phi at
        # the N-terminus). Confirm the shape lines up and the values
        # are finite + in the plausible HN-Ha range.
        assert J.shape == (n_frames, n_phi)
        assert np.all(np.isfinite(J))
        # Physically: 3J(HN, Ha) lies in roughly 0-12 Hz for protein
        # backbone phi; allow a small numerical margin either side.
        assert J.min() > -0.5
        assert J.max() < 13.0

    def test_unknown_model_raises(self, GS6_CP):
        with pytest.raises(SSException):
            nmr.compute_J3_HN_HA(GS6_CP, model="not-a-real-model")

    def test_uniform_weights_match_unweighted_mean(self, GS6_CP):
        """With uniform weights the result must equal the unweighted mean."""
        _, J = nmr.compute_J3_HN_HA(GS6_CP)
        n = GS6_CP.n_frames
        uniform = np.full(n, 1.0 / n)
        _, J_mean = nmr.compute_J3_HN_HA(GS6_CP, weights=uniform)
        assert J_mean.shape == (J.shape[1],)
        assert np.allclose(J_mean, J.mean(axis=0))

    def test_one_hot_weights_pick_a_single_frame(self, GS6_CP):
        _, J = nmr.compute_J3_HN_HA(GS6_CP)
        n = GS6_CP.n_frames
        w = np.zeros(n)
        w[3] = 1.0
        _, J_frame = nmr.compute_J3_HN_HA(GS6_CP, weights=w)
        assert np.allclose(J_frame, J[3])

    def test_invalid_weights_raise(self, GS6_CP):
        n = GS6_CP.n_frames
        bad = np.full(n, 2.0 / n)  # doesn't sum to 1
        with pytest.raises(SSException):
            nmr.compute_J3_HN_HA(GS6_CP, weights=bad)

    def test_stride_matches_subsample(self, GS6_CP):
        _, J = nmr.compute_J3_HN_HA(GS6_CP)
        _, J2 = nmr.compute_J3_HN_HA(GS6_CP, stride=2)
        # The strided result is the every-other-frame view of the full one.
        assert J2.shape[1] == J.shape[1]
        assert np.allclose(J2, J[::2])

    def test_return_uncertainty(self, GS6_CP):
        atoms, J, sigma = nmr.compute_J3_HN_HA(
            GS6_CP, model="Bax2007", return_uncertainty=True
        )
        assert sigma == pytest.approx(nmr.KARPLUS_HN_HA_UNCERTAINTIES["Bax2007"])
        assert J.shape == (GS6_CP.n_frames, len(atoms))

    def test_matches_mdtraj_reference(self, GS6_CP):
        """Cross-check against mdtraj's reference implementation."""
        try:
            from mdtraj.nmr import compute_J3_HN_HA as md_compute
        except ImportError:
            pytest.skip("mdtraj.nmr.compute_J3_HN_HA not available")
        # Underlying mdtraj trajectory lives on the SSProtein wrapper.
        traj = GS6_CP.traj
        _, J_md = md_compute(traj, model="Bax2007")
        _, J_ss = nmr.compute_J3_HN_HA(GS6_CP, model="Bax2007")
        # mdtraj returns (n_frames, n_phi) just like us.
        assert J_ss.shape == J_md.shape
        assert np.allclose(J_ss, J_md, atol=1e-6)


# --------------------------------------------------------------------------
# NOE distances
# --------------------------------------------------------------------------
class TestNOEEnsembleAverage:
    """The pure-numpy ``noe_ensemble_average`` collapse rule."""

    def test_uniform_closed_form(self):
        """Hand-computed NOE r^-6 average for a 3x2 distance array."""
        d = np.array([[2.0, 3.0], [3.0, 4.0], [4.0, 5.0]])
        out = nmr.noe_ensemble_average(d, power=6)
        # pair 0: (2^-6 + 3^-6 + 4^-6)/3, then ^(-1/6)
        expected_0 = ((2.0**-6 + 3.0**-6 + 4.0**-6) / 3.0) ** (-1.0 / 6.0)
        expected_1 = ((3.0**-6 + 4.0**-6 + 5.0**-6) / 3.0) ** (-1.0 / 6.0)
        assert out[0] == pytest.approx(expected_0)
        assert out[1] == pytest.approx(expected_1)

    def test_power_3_supported(self):
        d = np.array([[2.0, 3.0], [3.0, 4.0]])
        out6 = nmr.noe_ensemble_average(d, power=6)
        out3 = nmr.noe_ensemble_average(d, power=3)
        # power=3 collapse is the same form with a different exponent
        expected_3_0 = ((2.0**-3 + 3.0**-3) / 2.0) ** (-1.0 / 3.0)
        assert out3[0] == pytest.approx(expected_3_0)
        # The two power values should disagree for non-degenerate input
        assert not np.isclose(out6[0], out3[0])

    def test_one_hot_weights(self):
        d = np.array([[2.0, 3.0], [3.0, 4.0], [4.0, 5.0]])
        w = np.array([0.0, 1.0, 0.0])
        out = nmr.noe_ensemble_average(d, power=6, weights=w)
        # one-hot weights pick frame 1: distances [3, 4]
        assert np.allclose(out, np.array([3.0, 4.0]))

    def test_uniform_weights_match_unweighted(self):
        d = np.array([[2.0, 3.0], [3.0, 4.0], [4.0, 5.0]])
        n = d.shape[0]
        out_w = nmr.noe_ensemble_average(d, power=6, weights=np.full(n, 1.0 / n))
        out = nmr.noe_ensemble_average(d, power=6)
        assert np.allclose(out_w, out)

    def test_non_positive_distance_raises(self):
        d = np.array([[2.0, 0.0], [3.0, 4.0]])  # zero distance -> r^-p undefined
        with pytest.raises(SSException):
            nmr.noe_ensemble_average(d, power=6)

    def test_invalid_weights_raise(self):
        d = np.array([[2.0, 3.0], [3.0, 4.0]])
        bad = np.array([0.7, 0.7])  # doesn't sum to 1
        with pytest.raises(SSException):
            nmr.noe_ensemble_average(d, weights=bad)


class TestComputeNOEDistances:
    """``compute_NOE_distances`` operating on a real SSProtein fixture."""

    def test_shape_and_units(self, GS6_CP):
        # Two arbitrary atom-pair indices into the topology.
        pairs = np.array([[0, 10], [0, 20], [5, 15]])
        d = nmr.compute_NOE_distances(GS6_CP, pairs)
        assert d.shape == (GS6_CP.n_frames, 3)
        assert np.all(np.isfinite(d))
        # Distances should be in Angstroms - none should be negative,
        # nor absurdly large for a hexapeptide.
        assert d.min() > 0.0
        assert d.max() < 100.0

    def test_matches_mdtraj_directly(self, GS6_CP):
        import mdtraj as md
        pairs = np.array([[0, 10], [0, 20], [5, 15]])
        d_ss = nmr.compute_NOE_distances(GS6_CP, pairs)
        d_md = md.compute_distances(GS6_CP.traj, pairs) * 10.0  # nm -> A
        assert np.allclose(d_ss, d_md)

    def test_stride(self, GS6_CP):
        pairs = np.array([[0, 10]])
        d = nmr.compute_NOE_distances(GS6_CP, pairs)
        d2 = nmr.compute_NOE_distances(GS6_CP, pairs, stride=2)
        assert np.allclose(d2, d[::2])

    def test_bad_pairs_shape_raises(self, GS6_CP):
        with pytest.raises(SSException):
            nmr.compute_NOE_distances(GS6_CP, np.array([0, 10]))           # 1D
        with pytest.raises(SSException):
            nmr.compute_NOE_distances(GS6_CP, np.array([[0, 10, 20]]))     # n_cols != 2
