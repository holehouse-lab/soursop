# Import package, test suite, and other packages as needed
import numpy as np
import soursop
import pytest
import sys
import soursop.ssnmr as nmr


#
# All test code written by Alex Keeley 
#

def test_simple_test():
    output = nmr.compute_random_coil_chemical_shifts("ACDE", temperature=15, pH=6.5, use_ggxgg=True, use_perdeuteration=False, asFloat=True)
    control = [
        {"Res": "A", "Index": 0, "CA": 52.695, "CB": 18.967, "CO": 177.864, "N": 126.057, "HN": 8.519, "HA": 4.357},
        {"Res": "C", "Index": 1, "CA": 58.584, "CB": 29.752, "CO": 174.436, "N": 118.754, "HN": 8.466, "HA": 4.500},
        {"Res": "D", "Index": 2, "CA": 54.701, "CB": 40.940, "CO": 176.465, "N": 122.592, "HN": 8.476, "HA": 4.609},
        {"Res": "E", "Index": 3, "CA": 56.941, "CB": 30.167, "CO": 176.755, "N": 121.337, "HN": 8.426, "HA": 4.263}]

    for i in range(4):
        amino_acid = output[i]
        assert amino_acid == control[i]


def test_simple_types():
    """Tests the ability of the function randCoilChemShifts to accurately predict random coil chemical shifts for all 20 amino acid types.
    
    """

    output = nmr.compute_random_coil_chemical_shifts("ACDEFGHIKLMNPQRSTVWY", temperature=25, pH=6.5, use_ggxgg=True, use_perdeuteration=False, asFloat=True)
    control = [
        {"Res": "A", "Index": 0, "CA": 52.673, "CB": 19.015, "CO": 177.793, "N": 126.005, "HN": 8.429, "HA": 4.364},
        {"Res": "C", "Index": 1, "CA": 58.575, "CB": 29.764, "CO": 174.410, "N": 118.672, "HN": 8.396, "HA": 4.500},
        {"Res": "D", "Index": 2, "CA": 54.896, "CB": 40.919, "CO": 176.273, "N": 122.765, "HN": 8.403, "HA": 4.540},
        {"Res": "E", "Index": 3, "CA": 56.873, "CB": 30.051, "CO": 176.341, "N": 120.831, "HN": 8.335, "HA": 4.225},
        {"Res": "F", "Index": 4, "CA": 57.900, "CB": 39.333, "CO": 176.236, "N": 119.924, "HN": 8.186, "HA": 4.637},
        {"Res": "G", "Index": 5, "CA": 45.364, "CB": "**.***", "CO": 173.150, "N": 112.203, "HN": 8.145, "HA": 3.909},
        {"Res": "H", "Index": 6, "CA": 55.775, "CB": 29.978, "CO": 174.752, "N": 118.609, "HN": 8.137, "HA": 4.631},
        {"Res": "I", "Index": 7, "CA": 61.170, "CB": 38.698, "CO": 175.946, "N": 122.883, "HN": 8.143, "HA": 4.135},
        {"Res": "K", "Index": 8, "CA": 56.199, "CB": 32.866, "CO": 176.184, "N": 126.438, "HN": 8.474, "HA": 4.309},
        {"Res": "L", "Index": 9, "CA": 55.223, "CB": 42.321, "CO": 177.352, "N": 124.735, "HN": 8.114, "HA": 4.339},
        {"Res": "M", "Index": 10, "CA": 55.397, "CB": 32.878, "CO": 175.653, "N": 121.485, "HN": 8.400, "HA": 4.502},
        {"Res": "N", "Index": 11, "CA": 51.159, "CB": 38.132, "CO": 173.163, "N": 121.161, "HN": 8.416, "HA": 4.957},
        {"Res": "P", "Index": 12, "CA": 63.418, "CB": 32.272, "CO": 176.884, "N": "***.***", "HN": "*.***",
         "HA": 4.423},
        {"Res": "Q", "Index": 13, "CA": 55.989, "CB": 29.081, "CO": 176.161, "N": 120.570, "HN": 8.503, "HA": 4.305},
        {"Res": "R", "Index": 14, "CA": 56.205, "CB": 30.960, "CO": 176.358, "N": 123.141, "HN": 8.392, "HA": 4.407},
        {"Res": "S", "Index": 15, "CA": 58.326, "CB": 63.811, "CO": 174.642, "N": 117.860, "HN": 8.452, "HA": 4.549},
        {"Res": "T", "Index": 16, "CA": 62.004, "CB": 69.725, "CO": 174.464, "N": 116.881, "HN": 8.182, "HA": 4.350},
        {"Res": "V", "Index": 17, "CA": 62.918, "CB": 32.654, "CO": 175.762, "N": 122.829, "HN": 8.103, "HA": 3.984},
        {"Res": "W", "Index": 18, "CA": 57.244, "CB": 29.237, "CO": 175.644, "N": 124.340, "HN": 8.098, "HA": 4.638},
        {"Res": "Y", "Index": 19, "CA": 57.632, "CB": 38.892, "CO": 174.936, "N": 122.232, "HN": 7.732, "HA": 4.438}]

    for i in range(20):
        amino_acid = output[i]
        assert amino_acid == control[i]


def test_simple_deut():
    """Tests the ability of the function randCoilChemShifts to accurately predict random coil chemical shifts for perdeuterated proteins.
    
    """
    output = nmr.compute_random_coil_chemical_shifts("ACDE", temperature=25, pH=6.5, use_ggxgg=True, use_perdeuteration=True, asFloat=True)
    control = [
        {"Res": "A", "Index": 0, "CA": 51.993, "CB": 18.015, "CO": 177.793, "N": 126.005, "HN": 8.429, "HA": "*.***"},
        {"Res": "C", "Index": 1, "CA": 58.025, "CB": 29.054, "CO": 174.410, "N": 118.672, "HN": 8.396, "HA": "*.***"},
        {"Res": "D", "Index": 2, "CA": 54.178, "CB": 40.295, "CO": 176.417, "N": 122.553, "HN": 8.414, "HA": "*.***"},
        {"Res": "E", "Index": 3, "CA": 56.260, "CB": 29.243, "CO": 176.706, "N": 121.300, "HN": 8.362, "HA": "*.***"}]

    for i in range(4):
        amino_acid = output[i]
        assert amino_acid == control[i]


def test_simple_undeut():
    """Tests the ability of the function randCoilChemShifts to accurately predict random coil chemical shifts for unperdeuterated proteins

    """
    output = nmr.compute_random_coil_chemical_shifts("ACDE", temperature=25, pH=6.5, use_ggxgg=True, use_perdeuteration=False, asFloat=True)

    control = [
        {"Res": "A", "Index": 0, "CA": 52.673, "CB": 19.015, "CO": 177.793, "N": 126.005, "HN": 8.429, "HA": 4.364},
        {"Res": "C", "Index": 1, "CA": 58.575, "CB": 29.764, "CO": 174.410, "N": 118.672, "HN": 8.396, "HA": 4.500},
        {"Res": "D", "Index": 2, "CA": 54.728, "CB": 41.005, "CO": 176.417, "N": 122.553, "HN": 8.414, "HA": 4.608},
        {"Res": "E", "Index": 3, "CA": 56.950, "CB": 30.213, "CO": 176.706, "N": 121.300, "HN": 8.362, "HA": 4.266}]

    for i in range(4):
        amino_acid = output[i]
        assert amino_acid == control[i]


def test_simple_phos():
    """Tests the ability of the function randCoilChemShifts to accurately predict random coil chemical shifts for phosphorylated amino acids.

    """

    output = nmr.compute_random_coil_chemical_shifts("(sep)(tpo)(ptr)", temperature=25, pH=6.5, use_ggxgg=True, use_perdeuteration=False, asFloat=True)

    control = [
        {"Res": "SEP", "Index": 0, "CA": 58.483, "CB": 65.876, "CO": 174.762, "N": 119.143, "HN": 9.005, "HA": 4.440},
        {"Res": "TPO", "Index": 1, "CA": 62.641, "CB": 71.984, "CO": 174.201, "N": 117.185, "HN": 8.913, "HA": 4.261},
        {"Res": "PTR", "Index": 2, "CA": 57.802, "CB": 38.635, "CO": 175.690, "N": 122.211, "HN": 8.087, "HA": 4.585}]

    for i in range(3):
        amino_acid = output[i]
        assert amino_acid == control[i]


def test_simple_ph():
    """Tests the ability of the function randCoilChemShifts to accurately predict random coil chemical shifts in different pH conditions.

    """

    output = nmr.compute_random_coil_chemical_shifts("ACDE", temperature=25, pH=4, use_ggxgg=True, use_perdeuteration=False, asFloat=True)

    control = [
        {"Res": "A", "Index": 0, "CA": 52.673, "CB": 19.015, "CO": 177.793, "N": 126.005, "HN": 8.429, "HA": 4.364},
        {"Res": "C", "Index": 1, "CA": 58.575, "CB": 29.764, "CO": 174.410, "N": 118.672, "HN": 8.396, "HA": 4.500},
        {"Res": "D", "Index": 2, "CA": 53.868, "CB": 39.252, "CO": 175.702, "N": 121.615, "HN": 8.476, "HA": 4.675},
        {"Res": "E", "Index": 3, "CA": 56.281, "CB": 29.080, "CO": 176.324, "N": 120.696, "HN": 8.281, "HA": 4.340}]

    for i in range(4):
        amino_acid = output[i]
        assert amino_acid == control[i]


# ----------
def test_simple_temp():
    """Tests the ability of the function randCoilChemShifts to accurately predict random coil chemical shifts in different temperature conditions.

    """

    output = nmr.compute_random_coil_chemical_shifts("ACDE", temperature=15, pH=6.5, use_ggxgg=True, use_perdeuteration=False, asFloat=True)

    control = [
        {"Res": "A", "Index": 0, "CA": 52.695, "CB": 18.967, "CO": 177.864, "N": 126.057, "HN": 8.519, "HA": 4.357},
        {"Res": "C", "Index": 1, "CA": 58.584, "CB": 29.752, "CO": 174.436, "N": 118.754, "HN": 8.466, "HA": 4.500},
        {"Res": "D", "Index": 2, "CA": 54.701, "CB": 40.940, "CO": 176.465, "N": 122.592, "HN": 8.476, "HA": 4.609},
        {"Res": "E", "Index": 3, "CA": 56.941, "CB": 30.167, "CO": 176.755, "N": 121.337, "HN": 8.426, "HA": 4.263}]

    for i in range(4):
        amino_acid = output[i]
        assert amino_acid == control[i]


# ----------
def test_unit_round3():
    """Tests the ability of the round3 function to accurately and consistently round input numbers to exactly 3 decimal places.

    """
    raw = [110.7000, 103.9111, 112.7382, 114.8413, 112.5684, 107.9585, 105.9196, 108.2177, 112.5628, 112.0159]

    control = [110.700, 103.911, 112.738, 114.841, 112.568, 107.959, 105.920, 108.218, 112.563, 112.016]

    for i in range(10):
        assert nmr.__round3(raw[i], asFloat=True) == control[i]


# ----------
def test_unit_setseq():
    """Tests the ability of the set_sequence function to accurately and consistently translate a sequence of amino acid abbreviations into a corresponding sequence of representative numbers.

    """
    key_aa1 = [0, -1, 1, 2, 3, 4, 5, 6, 7, -1, 8, 9, 10, 11, -1, 12, 13, 14, 15, 16, -1, 17, 18, -1, 19, -1]
    key_aa3 = {"ALA": 0, "CYS": 1, "ASP": 2, "GLU": 3, "PHE": 4, "GLY": 5, "HIS": 6, "ILE": 7, "LYS": 8, "LEU": 9,
               "MET": 10, "ASN": 11, "PRO": 12, "GLN": 13, "ARG": 14, "SER": 15, "THR": 16, "VAL": 17, "TRP": 18,
               "TYR": 19,
               "PSER": 20, "SEP": 20, "PS": 20, "PTHR": 21, "TPO": 21, "PT": 21, "PTYR": 22, "PTR": 22, "PY": 22}

    raw = ["KEQEERL", "EGHVTPC", "GQGDHQD", "MHDTMII", "FTWFSGH", "KLD(pT)GQK", "L(pY)LQATYG", "QYQES(pS)RP",
           "WWGN(pS)TGI", "KHRDLAL(pY)"]

    control = [[23, 23, 8, 3, 13, 3, 3, 14, 9, 23, 23],
               [23, 23, 3, 5, 6, 17, 16, 12, 1, 23, 23],
               [23, 23, 5, 13, 5, 2, 6, 13, 2, 23, 23],
               [23, 23, 10, 6, 2, 16, 10, 7, 7, 23, 23],
               [23, 23, 4, 16, 18, 4, 15, 5, 6, 23, 23],
               [23, 23, 8, 9, 2, 21, 5, 13, 8, 23, 23],
               [23, 23, 9, 22, 9, 13, 0, 16, 19, 5, 23, 23],
               [23, 23, 13, 19, 13, 3, 15, 20, 14, 12, 23, 23],
               [23, 23, 18, 18, 5, 11, 20, 16, 5, 7, 23, 23],
               [23, 23, 8, 6, 14, 2, 9, 0, 9, 22, 23, 23]]

    for i in range(10):
        seq = nmr.__set_sequence(raw[i], key_aa1, key_aa3)
        assert seq[0] == control[i]
