##     _____  ____  _    _ _____   _____  ____  _____  
##   / ____|/ __ \| |  | |  __ \ / ____|/ __ \|  __ \ 
##  | (___ | |  | | |  | | |__) | (___ | |  | | |__) |
##   \___ \| |  | | |  | |  _  / \___ \| |  | |  ___/ 
##   ____) | |__| | |__| | | \ \ ____) | |__| | |     
##  |_____/ \____/ \____/|_|  \_\_____/ \____/|_|     

## Alex Holehouse (Pappu Lab and Holehouse Lab) and Jared Lalmansing (Pappu lab)
## ssnmr was largely written by Alex Keeley
## Simulation analysis package
## Copyright 2014 - 2022
##


"""

**Author(s):** Alex Keeley (with Alex Holehouse)

"""
from .ssexceptions import SSException
import re

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------



# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
def compute_random_coil_chemical_shifts(protein_sequence, temperature=25, pH=7.4, use_ggxgg=True, use_perdeuteration=False, asFloat=True):
    """

    Function that predicts the random coil chemical shifts for user-provided amino 
    acid sequence corrected for user-provided conditions. Specifically, chemical 
    shift and general sequence correction factors are from [1], temperature corrections     
    and glycine corrections are from [2] and the underlying methods associated with 
    correction-factor calculations are in [3]. The correction factors for pertdeuteration 
    are from [4].
    
    Input sequence can be a standard one-letter sequence, but phosphoresidues can also 
    be included (see examples). Code is based on JavaScript written by Alex Maltsev 
    at the NIH and can be accessed here [5] 
    

    Parameters
    ----------

    sequence : str
        The sequence of representative abbreviations for the sequence of amino acids

    temperature : float or int    
        Experiment temperature of the sample of amino acids for use in corrected chemical 
        shift calculations (note units are in degrees celcius). Default value is 25 and 
        should be between 0 and 100.

    pH : float or int
        pH of the sample of amino acids for use in corrected chemical shift calculations. 
        Default value us 7.4 and should be between 0 and 14.

    use_ggxgg : bool 
        Whether to use GGXGG-based neighbor correction for glycines. Default is True.

    use_perdeuteration : bool (default = False)
        Whether perdeuterated correction factors should be used. Note this cannot 
        work with phosphoresidues.
          

    asFloat : bool (default = True)
          Whether to populate output dictionaries with float or string variables containing 
          chemical shift numbers. False (strings) by default.
          

    Returns
    -------
    output : list of dict
           List containing a dictionary for each amino acid in the provided sequence detailing 
           abbreviation and chemical shifts for the main six different atoms. 

    Examples
    --------

    
    References
    ----------

    [1] Kjaergaard, M. and Poulsen, F.M. (2011) Sequence correction of random coil chemical shifts: 
    correlation between neighbor correction factors and changes in the Ramachandran distribution 
    J. Biomol. NMR 50(2):157-165        
    
    [2] Kjaergaard, M., Brander, S. and Poulsen, F.M. (2011) Random coil chemical shifts for 
    intrinsically disordered proteins: Effects of temperature and pH J. Biomol. NMR 49(2):139-49.
        
    [3] Schwarzinger, S., Kroon, G.J., Foss, T.R., Chung. J., Wright, P.E., Dyson, H.J. (2001) 
    Sequence-dependent correction of random coil NMR chemical shifts. JACS 123(13):2970-8.        

    [4] Cavanagh, J., Fairbrother, W.J., Palmer, A.G., Rance, M. and Skelton, N.J. (2007) 
    Protein NMR Spectroscopy - Principles and practice. 2nd edition. Academic Press

    [5] https://www1.bio.ku.dk/english/research/bms/research/sbinlab/randomchemicalshifts/
    
    

    """
    # sanity check temperature
    if temperature > 100 or temperature < 0:
        raise SSException('Temperature provided (%i) was non-physiological. Remember temperature should be in *celcius*.' %(temperature))

    # pH sanity check
    if pH < 0 or pH > 14:
        raise SSException('pH provided (%i) was non-physiological. Remember pH should be in between 0 and 14.' %(pH))



    # SETUP
    # The array 'key' is used to translate amino acid letter code into
    # numerical index. Value is -1 when there is no such amino acid letter
    key_aa1 = [0, -1, 1, 2, 3, 4, 5, 6, 7, -1, 8, 9, 10, 11, -1, 12, 13, 14, 15, 16, -1, 17, 18, -1, 19, -1]
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
        "PY": 22}

    # The array 'sequence' keeps the set of amino acid indices for the given
    # protein. The array size is not fixed and will accommodate a given
    # sequence.
    sequence = []

    # The XX_av arrays contain uncorrected random coil values for atom type XX at 5C and pH 6.5
    ca_av = [
        52.747,
        58.639,
        54.586,
        56.805,
        57.804,
        45.317,
        55.933,
        61.279,
        56.542,
        55.359,
        55.565,
        53.422,
        63.236,
        56.041,
        56.287,
        58.604,
        62.151,
        62.574,
        57.464,
        57.968,
        58.235,
        62.985,
        57.977]
    cb_av = [
        19.048,
        29.693,
        40.942,
        30.181,
        39.482,
        0,
        29.974,
        38.640,
        33.024,
        42.262,
        32.677,
        38.643,
        32.155,
        29.322,
        30.786,
        63.706,
        69.877,
        32.784,
        29.262,
        38.654,
        65.725,
        72.449,
        38.698]
    co_av = [
        177.967,
        174.736,
        176.527,
        176.707,
        175.721,
        174.363,
        175.044,
        176.451,
        176.738,
        177.580,
        176.423,
        175.436,
        177.023,
        176.118,
        176.442,
        174.798,
        174.647,
        176.330,
        176.276,
        175.773,
        174.769,
        174.566,
        175.740]
    n_av = [125.922, 121.101, 121.778, 122.87, 121.691, 110.613, 120.581, 123.509, 123.571, 124.07, 122.428,
            119.966, 0, 122.304, 123.182, 117.629, 116.375, 122.879, 122.099, 121.787, 118.870, 119.095, 121.935]
    hn_av = [8.575, 8.627, 8.572, 8.692, 8.438, 8.662, 8.653, 8.441, 8.592, 8.485, 8.63,
             8.672, 0, 8.67, 8.61, 8.585, 8.413, 8.44, 8.267, 8.396, 9.134, 9.091, 8.336]
    ha_av = [4.293, 4.498, 4.588, 4.266, 4.644, 3.977, 4.669, 4.134, 4.287, 4.338, 4.478,
             4.696, 4.435, 4.321, 4.319, 4.437, 4.343, 4.076, 4.672, 4.588, 4.435, 4.253, 4.596]

    # The XX_t arrays contain temperature corrections for atom type XX
    ca_t = [-2.2, -0.9, 2.78, 0.94, -4.74, 3.28, 7.76, -1.98, -0.76, 1.73, 4.09,
            2.78, 1.12, 2.26, -1.37, -1.7, -0.03, -2.79, -2.69, -4.99, 1.38, -5.27, -3.68]
    cb_t = [4.73, 1.27, 6.53, 4.6, 2.42, 0, 15.54, 4.6, 2.41, 4.92, 9.37,
            5.08, -0.19, 3.62, 3.54, 4.4, 2.15, 2.49, 3.07, 2.92, 4.82, 3.66, 0.01]
    co_t = [-7.09, -2.55, -4.81, -4.9, -6.9, -3.21, -8.3, -8.73, -7.12, -8.18, -8.17, -
            6.11, -4.01, -5.7, -6.9, -4.67, -5.24, -8.09, -7.88, -7.73, -1.71, -3.00, -3.00]
    n_t = [-5.25, -8.2, -3.91, -3.7, -11.15, -6.15, 3.3, -12.73, -7.6, -2.85, -6.2, -
           3.25, 0, -6.45, -5.3, -3.8, -6.7, -14.16, -10.1, -12, -22.69, -31.26, -25.64]
    hn_t = [-8.95, -7, -6.2, -6.46, -7.5, -9.1, -8.3, -7.78, -7.5, -7.45, -7.05, -
            6.95, 0, -7.2, -7.05, -7.6, -7.25, -7.64, -7.8, -7.73, -6.47, -8.34, -8.26]
    ha_t = [0.69, 0, -0.06, 0.31, 0.4, -0.02, -0.93, 0.37, 0.38, 0.05, -0.48, -
            2.92, -0.02, 0.26, 0.4, 0.05, -0.04, 0.47, 0.38, 0.53, 0.48, 1.36, 0.42]

    # Neighbor correction factors for Ca. Notice extra zero for 'no Residue' (last element).
    ca_a = [-0.007, -0.043, -0.019, -0.025, 0.168, -0.059, -0.036, -0.076, -0.003, 0.084, 0.004,
            0.01, -0.206, 0, 0, 0.003, -0.007, -0.095, -0.003, 0.14, -0.304, -0.303, 0.264, 0]
    ca_b = [-0.156, -0.011, 0.055, -0.013, -0.018, 0.132, -0.041, -0.27, -0.148, -0.205, -0.096,
            0.047, -2.249, 0, -0.128, -0.024, -0.087, -0.242, 0.24, -0.004, 0.180, -0.041, 0.040, 0]
    ca_c = [-0.076, 0.095, 0.086, -0.036, -0.301, 0.007, -0.066, -0.217, -0.082, -0.139, -0.039,
            0.157, -0.072, 0, -0.062, 0.025, -0.053, -0.212, -0.226, -0.344, -0.160, -0.100, -0.403, 0]
    ca_d = [0.007, 0.044, 0.131, 0.007, 0.105, 0.061, 0.09, -0.003, 0.048, -0.031, 0.003,
            0.1, -0.024, 0, 0.074, 0.073, 0.05, -0.01, 0.112, 0.123, -0.009, 0.058, 0.131, 0]

    # Neighbor correction factors for Cb. Notice extra zero for 'no Residue' (last element).
    cb_a = [0.017, 0.07, -0.011, -0.002, -0.086, -0.011, -0.012, 0.097, 0.033, 0.058, 0.036, -
            0.006, 0.149, 0, 0.11, -0.243, 0.039, 0.066, -0.159, -0.102, 0.013, 0.085, -0.154, 0]
    cb_b = [-0.049, -0.117, -0.03, -0.096, -0.151, -0.128, -0.08, -0.087, -0.077, -0.231, -0.153, -
            0.138, -0.701, 0, -0.098, -0.043, -0.044, -0.039, -0.343, -0.147, 0.076, 0.210, -0.511, 0]
    cb_c = [0.078, -0.025, -0.104, 0.054, 0.113, 0.031, 0.052, 0, 0.021, -0.08, -0.031, -
            0.069, 0.102, 0, -0.005, -0.028, 0.259, 0.034, 0.102, 0.14, -0.175, 0.002, 0.081, 0]
    cb_d = [0.058, 0.049, -0.111, -0.004, -0.022, -0.067, -0.011, 0.099, 0.083, 0.119, 0.08, -
            0.074, 0.107, 0, 0.031, 0.006, 0.027, 0.078, -0.145, -0.055, -0.065, -0.031, -0.016, 0]

    # Neighbor correction factors for Co. Notice extra zero for 'no Residue' (last element).
    co_a = [-0.043, -0.046, 0.043, -0.086, -0.144, -0.008, -0.07, -0.249, -0.055, -0.163, -0.138,
            0.021, -0.228, 0, -0.066, 0.033, -0.025, -0.21, 0.05, -0.135, -0.136, -0.367, -0.031, 0]
    co_b = [-0.181, -0.075, -0.251, 0.065, -0.357, 0.546, -0.195, -0.113, -0.042, -0.083, 0.02, -
            0.23, -2.01, 0, -0.047, 0.148, 0.193, -0.06, -0.183, -0.338, 0.245, 0.074, -0.251, 0]
    co_c = [0.062, -0.082, 0.143, 0, -0.546, 0.212, -0.194, -0.201, -0.068, -0.109, -0.123,
            0.022, 0.096, 0, -0.046, -0.017, -0.112, -0.146, -0.665, -0.569, 0.002, -0.129, -0.605, 0]
    co_d = [0.005, -0.044, 0.177, 0.068, -0.007, 0.069, 0.01, -0.037, -0.04, -0.018, -0.015,
            0.075, -0.069, 0, -0.051, 0.024, 0.01, -0.017, -0.002, 0.012, 0.123, 0.096, 0.052, 0]

    # Neighbor correction factors for N. Notice extra zero for 'no Residue' (last element).
    n_a = [0.003, -0.05, 0.041, 0.018, 0.212, -0.024, -0.066, -0.069, -0.036, 0.039, 0.06, -
           0.035, -0.091, 0, -0.042, 0.009, -0.02, -0.079, -0.073, 0.203, -0.018, -0.053, 0.247, 0]
    n_b = [0.141, 0.147, -0.033, -0.035, -0.445, -0.043, -0.316, 0.249, 0.068, -0.158, -0.044, -
           0.256, 1.433, 0, 0.17, 0.167, 0.229, 0.318, -0.231, -0.401, 0.721, 0.485, -0.532, 0]
    n_c = [-2.25, 1.102, -1.407, -0.432, 0.187, -2.086, -0.033, 3.166, 0.152, -0.642, -0.053, -
           1.25, -1.027, 0, 0.157, 0.212, 0.881, 2.83, 0.051, 0.395, -0.758, 1.619, 0.232, 0]
    n_d = [-0.209, -0.072, -1.003, -0.139, -0.089, -0.445, -0.049, 0.649, 0.17, -0.12, 0.015, -
           0.757, -0.082, 0, 0.183, -0.62, 0.014, 0.634, -0.764, -0.127, -0.821, -0.193, 0.081, 0]

    # Neighbor correction factors for Hn. Notice extra zero for 'no Residue' (last element).
    hn_a = [-0.002, -0.007, 0.001, 0, -0.011, 0.007, -0.009, -0.031, -0.005, -0.004, -0.009,
            0.004, -0.024, 0, -0.004, 0.006, -0.001, -0.029, -0.044, -0.013, -0.010, -0.019, 0.043, 0]
    hn_b = [-0.036, 0.032, 0.01, -0.012, -0.034, 0.032, -0.005, -0.027, -0.024, -0.004, -0.156,
            0.031, -0.007, 0, -0.013, -0.046, 0.031, -0.028, -0.023, -0.02, -0.013, -0.044, 0.109, 0]
    hn_c = [-0.101, 0.052, -0.148, -0.021, -0.265, -0.206, -0.115, 0.028, -0.078, -0.093, -0.041,
            0.002, 0.079, 0, 0.017, -0.023, -0.016, 0.042, -0.53, -0.305, -0.123, 0.091, -0.344, 0]
    hn_d = [-0.072, -0.05, -0.104, -0.034, -0.121, 0.001, 0.017, 0.008, -0.003, -0.069, -0.021, -
            0.095, -0.03, 0, 0.009, -0.132, -0.035, 0.021, -0.358, -0.151, -0.173, -0.056, -0.137, 0]

    # Neighbor correction factors for Ha. Notice extra zero for 'no Residue' (last element).
    ha_a = [-0.002, 0.006, 0.015, 0.01, -0.068, 0.014, -0.025, -0.007, -0.012, -0.019, -0.015,
            0.001, 0.014, 0, -0.013, 0.006, 0.006, -0.003, -0.075, -0.063, 0.030, 0.027, -0.080, 0]
    ha_b = [-0.002, 0.042, 0.005, 0.009, -0.055, 0.032, -0.029, 0.035, 0.002, 0.016, 0.008,
            0.016, 0.305, 0, 0.008, 0.069, 0.106, 0.047, -0.067, -0.05, 0.033, 0.073, -0.059, 0]
    ha_c = [-0.013, 0.022, -0.007, 0, -0.027, 0.014, 0.006, 0.025, 0.002, 0.012, 0.013,
            0.003, -0.027, 0, 0.008, 0.043, 0.036, 0.02, -0.149, -0.037, 0.033, 0.008, -0.055, 0]
    ha_d = [-0.002, 0.006, -0.022, -0.005, -0.057, 0.005, -0.012, -0.011, -0.008, 0.001, -0.002, -
            0.008, 0.005, 0, -0.007, -0.007, -0.012, -0.012, -0.166, -0.063, -0.028, -0.039, -0.069, 0]

    # Neighbor correction factors for Ca for Gly.
    gly_ca_a = [-0.011, -0.001, -0.085, -0.001, -0.048, 0, -0.05, -0.054, 0.016, -0.057,
                0.019, -0.05, -0.221, 0, 0, -0.006, -0.003, -0.042, -0.082, -0.034, 0, 0, 0, 0]
    gly_ca_b = [-0.149, -0.046, -0.045, -0.14, -0.244, 0, -0.09, -0.19, -0.01, -0.188, -
                0.03, -0.015, -0.8, -0.002, -0.07, -0.051, -0.041, -0.174, -0.193, -0.245, 0, 0, 0, 0]
    gly_ca_c = [0.076, 0.117, 0.329, 0.141, 0.072, 0, 0.02, 0.029, -0.061, 0.052,
                0.114, 0.238, 0.033, 0.028, -0.01, 0.129, 0.134, 0.024, 0.11, 0.068, 0, 0, 0, 0]
    gly_ca_d = [0.017, 0.007, 0.066, 0.053, 0.006, 0, 0.01, 0.024, 0.013, -0.053, 0.011,
                0.018, 0.048, 0.017, 0.02, -0.003, 0.001, 0.016, -0.08, -0.004, 0, 0, 0, 0]

    # Neighbor correction factors for Co for Gly.
    gly_co_a = [-0.113, -0.076, -0.112, -0.113, -0.297, 0, -0.137, -0.226, -0.075, -0.153, -
                0.092, -0.085, -0.558, -0.039, -0.038, -0.073, -0.084, -0.232, -0.296, -0.303, 0, 0, 0, 0]
    gly_co_b = [-0.754, -0.46, -0.753, -0.44, -0.832, 0, -0.706, -0.558, -0.494, -0.491, -
                0.424, -0.714, -2.799, -0.477, -0.444, -0.4, -0.203, -0.526, -0.565, -0.863, 0, 0, 0, 0]
    gly_co_c = [-0.042, -0.24, 0.165, -0.024, -0.23, 0, -0.174, -0.152, -0.187, -0.104, -
                0.165, -0.079, -0.058, -0.148, -0.186, -0.13, -0.127, -0.154, -0.31, -0.23, 0, 0, 0, 0]
    gly_co_d = [-0.018, -0.056, 0.049, 0.013, -0.097, 0, -0.039, -0.029, -0.034, -0.006, -
                0.024, -0.025, -0.024, -0.027, -0.021, -0.054, -0.047, -0.024, -0.185, -0.156, 0, 0, 0, 0]

    # Neighbor correction factors for N for Gly.
    gly_n_a = [-0.044, -0.081, -0.104, 0.094, -0.142, 0, -0.087, -0.176, -0.056, -0.094, -
               0.084, -0.104, -0.115, -0.074, -0.039, -0.045, -0.055, -0.15, -0.234, -0.163, 0, 0, 0, 0]
    gly_n_b = [-0.03, -0.078, 0.045, -0.023, -0.415, 0, -0.427, -0.111, -0.143, -0.145, -
               0.161, -0.149, -0.134, -0.097, -0.117, 0.032, 0, -0.054, -0.317, -0.381, 0, 0, 0, 0]
    gly_n_c = [-0.535, 2.6, 0.764, 1.235, 2.41, 0, 0.913, 4.157, 1.248, 0.809, 1.288,
               0.692, 0.756, 1.304, 1.296, 2.114, 2.444, 3.754, 2.703, 2.709, 0, 0, 0, 0]
    gly_n_d = [-0.158, -0.003, -0.027, -0.094, -0.408, 0, 0.729, -0.03, -0.061, -0.143, -
               0.03, -0.16, -0.157, -0.056, -0.031, -0.154, -0.024, -0.057, -0.821, -0.495, 0, 0, 0, 0]

    # Set pSER, pTHR, and pTYR Gly corrections to their SER, THR and TYR equivalents
    gly_ca_a[key_aa3["PSER"]] = gly_ca_a[key_aa3["SER"]]
    gly_ca_a[key_aa3["PTHR"]] = gly_ca_a[key_aa3["THR"]]
    gly_ca_a[key_aa3["PTYR"]] = gly_ca_a[key_aa3["TYR"]]
    gly_ca_b[key_aa3["PSER"]] = gly_ca_b[key_aa3["SER"]]
    gly_ca_b[key_aa3["PTHR"]] = gly_ca_b[key_aa3["THR"]]
    gly_ca_b[key_aa3["PTYR"]] = gly_ca_b[key_aa3["TYR"]]
    gly_ca_c[key_aa3["PSER"]] = gly_ca_a[key_aa3["SER"]]
    gly_ca_a[key_aa3["PTHR"]] = gly_ca_a[key_aa3["THR"]]
    gly_ca_a[key_aa3["PTYR"]] = gly_ca_a[key_aa3["TYR"]]
    gly_ca_d[key_aa3["PSER"]] = gly_ca_d[key_aa3["SER"]]
    gly_ca_d[key_aa3["PTHR"]] = gly_ca_d[key_aa3["THR"]]
    gly_ca_d[key_aa3["PTYR"]] = gly_ca_d[key_aa3["TYR"]]

    gly_co_a[key_aa3["PSER"]] = gly_co_a[key_aa3["SER"]]
    gly_co_a[key_aa3["PTHR"]] = gly_co_a[key_aa3["THR"]]
    gly_co_a[key_aa3["PTYR"]] = gly_co_a[key_aa3["TYR"]]
    gly_co_b[key_aa3["PSER"]] = gly_co_b[key_aa3["SER"]]
    gly_co_b[key_aa3["PTHR"]] = gly_co_b[key_aa3["THR"]]
    gly_co_b[key_aa3["PTYR"]] = gly_co_b[key_aa3["TYR"]]
    gly_co_c[key_aa3["PSER"]] = gly_co_a[key_aa3["SER"]]
    gly_co_a[key_aa3["PTHR"]] = gly_co_a[key_aa3["THR"]]
    gly_co_a[key_aa3["PTYR"]] = gly_co_a[key_aa3["TYR"]]
    gly_co_d[key_aa3["PSER"]] = gly_co_d[key_aa3["SER"]]
    gly_co_d[key_aa3["PTHR"]] = gly_co_d[key_aa3["THR"]]
    gly_co_d[key_aa3["PTYR"]] = gly_co_d[key_aa3["TYR"]]

    gly_n_a[key_aa3["PSER"]] = gly_n_a[key_aa3["SER"]]
    gly_n_a[key_aa3["PTHR"]] = gly_n_a[key_aa3["THR"]]
    gly_n_a[key_aa3["PTYR"]] = gly_n_a[key_aa3["TYR"]]
    gly_n_b[key_aa3["PSER"]] = gly_n_b[key_aa3["SER"]]
    gly_n_b[key_aa3["PTHR"]] = gly_n_b[key_aa3["THR"]]
    gly_n_b[key_aa3["PTYR"]] = gly_n_b[key_aa3["TYR"]]
    gly_n_c[key_aa3["PSER"]] = gly_n_a[key_aa3["SER"]]
    gly_n_a[key_aa3["PTHR"]] = gly_n_a[key_aa3["THR"]]
    gly_n_a[key_aa3["PTYR"]] = gly_n_a[key_aa3["TYR"]]
    gly_n_d[key_aa3["PSER"]] = gly_n_d[key_aa3["SER"]]
    gly_n_d[key_aa3["PTHR"]] = gly_n_d[key_aa3["THR"]]
    gly_n_d[key_aa3["PTYR"]] = gly_n_d[key_aa3["TYR"]]

    # Arrays for calculation of pH corrected shifts
    asp_ph_0 = [[53.05, 37.81, 175.25, 120.10, 8.68, 4.70],
                [54.59, 40.95, 176.53, 121.78, 8.57, 4.58]]

    glu_ph_0 = [[55.90, 28.65, 176.19, 122.04, 8.58, 4.36],
                [56.81, 30.19, 176.71, 122.86, 8.69, 4.26]]

    his_ph_0 = [[55.30, 28.89, 174.37, 120.00, 8.78, 4.69],
                [56.69, 31.27, 175.85, 121.46, 8.44, 4.63]]

    sep_ph_0 = [[57.613, 66.527, 174.062, 116.940, 8.805, 4.525],
                [58.433, 65.463, 174.977, 119.488, 9.242, 4.409]]

    tpo_ph_0 = [[61.827, 73.986, 174.063, 115.936, 8.615, 4.407],
                [63.524, 71.718, 174.793, 121.101, 9.397, 4.167]]

    ptr_ph_0 = [[57.905, 38.760, 175.605, 121.797, 8.403, 4.604],
                [58.002, 38.685, 175.788, 121.953, 8.320, 4.595]]

    # Arrays for keeping the chemical shifts calculated for the entered pH value
    asp_ph_corr = [0, 0, 0, 0, 0, 0]
    glu_ph_corr = [0, 0, 0, 0, 0, 0]
    his_ph_corr = [0, 0, 0, 0, 0, 0]
    sep_ph_corr = [0, 0, 0, 0, 0, 0]
    tpo_ph_corr = [0, 0, 0, 0, 0, 0]
    ptr_ph_corr = [0, 0, 0, 0, 0, 0]

    # Arrays for CS corrections for deuterated proteins
    ca_deut = [-0.68, -0.55, -0.55, -0.69, -0.55, -0.39, -0.55, -0.77, -0.69, -
               0.62, -0.69, -0.55, -0.69, -0.69, -0.69, -0.55, -0.55, -0.84, -0.55, -0.55]
    cb_deut = [-1.00, -0.71, -0.71, -0.97, -0.71, 0.00, -0.71, -1.28, -1.11, -
               1.26, -0.97, -0.71, -1.11, -0.97, -1.11, -0.71, -0.71, -1.20, -0.71, -0.71]

    # RUN
    cur = 0
    delta_T = temperature - 5  # difference between the given temperature and 5 degrees C
    ca_pred = 0
    cb_pred = 0
    co_pred = 0
    n_pred = 0
    hn_pred = 0
    ha_pred = 0

    output = []

    # Calculate deprotonated fractions of Asp, Glu and His at the given pH
    asp_deprot_frac = 7.78 * (10**-5) / (7.78 * (10**-5) + (10**(-pH)))
    glu_deprot_frac = 3.43 * (10**-5) / (3.43 * (10**-5) + (10**(-pH)))
    his_deprot_frac = 1.67 * (10**-7) / (1.67 * (10**-7) + (10**(-pH)))

    sep_deprot_frac = 9.76 * (10**-7) / (9.76 * (10**-7) + (10**(-pH)))
    tpo_deprot_frac = 5.00 * (10**-7) / (5.00 * (10**-7) + (10**(-pH)))
    ptr_deprot_frac = 1.47 * (10**-6) / (1.47 * (10**-6) + (10**(-pH)))

    # Calculate pH corrected chemical shifts
    for i in range(6):
        asp_ph_corr[i] = asp_deprot_frac * asp_ph_0[1][i] + (1 - asp_deprot_frac) * asp_ph_0[0][i]
        glu_ph_corr[i] = glu_deprot_frac * glu_ph_0[1][i] + (1 - glu_deprot_frac) * glu_ph_0[0][i]
        his_ph_corr[i] = his_deprot_frac * his_ph_0[1][i] + (1 - his_deprot_frac) * his_ph_0[0][i]
        sep_ph_corr[i] = sep_deprot_frac * sep_ph_0[1][i] + (1 - sep_deprot_frac) * sep_ph_0[0][i]
        tpo_ph_corr[i] = tpo_deprot_frac * tpo_ph_0[1][i] + (1 - tpo_deprot_frac) * tpo_ph_0[0][i]
        ptr_ph_corr[i] = ptr_deprot_frac * ptr_ph_0[1][i] + (1 - ptr_deprot_frac) * ptr_ph_0[0][i]

    print(protein_sequence)
    sequences = __set_sequence(protein_sequence, key_aa1, key_aa3)
    sequence = sequences[0]
    aminos = sequences[1]

    for j in range(len(sequence) - 4):
        output.append({"Res":   aminos[j], "Index": j})

    # deuterated parameters not available for phosphorylated Residues
    if ((22 in sequence) or (25 in sequence) or (28 in sequence)) and use_perdeuteration:
        raise SSException('Phosphorylated amino acids not supported in deuterated proteins')

    while cur < (len(sequence) - 4):
        if sequence[cur + 2] == 2:  # Aspartate
            ca_pred = asp_ph_corr[0]
            cb_pred = asp_ph_corr[1]
            co_pred = asp_ph_corr[2]
            n_pred = asp_ph_corr[3]
            hn_pred = asp_ph_corr[4]
            ha_pred = asp_ph_corr[5]

        elif sequence[cur + 2] == 3:  # Glutamate
            ca_pred = glu_ph_corr[0]
            cb_pred = glu_ph_corr[1]
            co_pred = glu_ph_corr[2]
            n_pred = glu_ph_corr[3]
            hn_pred = glu_ph_corr[4]
            ha_pred = glu_ph_corr[5]

        elif sequence[cur + 2] == 6:  # Histidine
            ca_pred = his_ph_corr[0]
            cb_pred = his_ph_corr[1]
            co_pred = his_ph_corr[2]
            n_pred = his_ph_corr[3]
            hn_pred = his_ph_corr[4]
            ha_pred = his_ph_corr[5]

        elif sequence[cur + 2] == key_aa3["SEP"]:  # phospho-SER
            ca_pred = sep_ph_corr[0]
            cb_pred = sep_ph_corr[1]
            co_pred = sep_ph_corr[2]
            n_pred = sep_ph_corr[3]
            hn_pred = sep_ph_corr[4]
            ha_pred = sep_ph_corr[5]

        elif sequence[cur + 2] == key_aa3["TPO"]:  # phospho-THR
            ca_pred = tpo_ph_corr[0]
            cb_pred = tpo_ph_corr[1]
            co_pred = tpo_ph_corr[2]
            n_pred = tpo_ph_corr[3]
            hn_pred = tpo_ph_corr[4]
            ha_pred = tpo_ph_corr[5]

        elif sequence[cur + 2] == key_aa3["PTR"]:  # phospho-TYR
            ca_pred = ptr_ph_corr[0]
            cb_pred = ptr_ph_corr[1]
            co_pred = ptr_ph_corr[2]
            n_pred = ptr_ph_corr[3]
            hn_pred = ptr_ph_corr[4]
            ha_pred = ptr_ph_corr[5]

        else:  # any other amino acid
            ca_pred = ca_av[sequence[cur + 2]]
            cb_pred = cb_av[sequence[cur + 2]]
            co_pred = co_av[sequence[cur + 2]]
            n_pred = n_av[sequence[cur + 2]]
            hn_pred = hn_av[sequence[cur + 2]]
            ha_pred = ha_av[sequence[cur + 2]]

        # Apply the neighbor and temperature corrections
        if sequence[cur + 2] == 5 and use_ggxgg:  # special case of glycine
            ca_pred += gly_ca_a[sequence[cur + 4]] + gly_ca_b[sequence[cur + 3]] + \
                gly_ca_c[sequence[cur + 1]] + gly_ca_d[sequence[cur]] + (delta_T * ca_t[sequence[cur + 2]] / 1000)
            co_pred += gly_co_a[sequence[cur + 4]] + gly_co_b[sequence[cur + 3]] + \
                gly_co_c[sequence[cur + 1]] + gly_co_d[sequence[cur]] + (delta_T * co_t[sequence[cur + 2]] / 1000)
            n_pred += gly_n_a[sequence[cur + 4]] + gly_n_b[sequence[cur + 3]] + \
                gly_n_c[sequence[cur + 1]] + gly_n_d[sequence[cur]] + (delta_T * n_t[sequence[cur + 2]] / 1000)

        else:  # all other Residues
            ca_pred += ca_a[sequence[cur + 4]] + ca_b[sequence[cur + 3]] + ca_c[sequence[cur + 1]] + \
                ca_d[sequence[cur]] + (delta_T * ca_t[sequence[cur + 2]] / 1000)
            co_pred += co_a[sequence[cur + 4]] + co_b[sequence[cur + 3]] + co_c[sequence[cur + 1]] + \
                co_d[sequence[cur]] + (delta_T * co_t[sequence[cur + 2]] / 1000)
            n_pred += n_a[sequence[cur + 4]] + n_b[sequence[cur + 3]] + n_c[sequence[cur + 1]] + \
                n_d[sequence[cur]] + (delta_T * n_t[sequence[cur + 2]] / 1000)

        cb_pred += cb_a[sequence[cur + 4]] + cb_b[sequence[cur + 3]] + cb_c[sequence[cur + 1]] + \
            cb_d[sequence[cur]] + (delta_T * cb_t[sequence[cur + 2]] / 1000)
        hn_pred += hn_a[sequence[cur + 4]] + hn_b[sequence[cur + 3]] + hn_c[sequence[cur + 1]] + \
            hn_d[sequence[cur]] + (delta_T * hn_t[sequence[cur + 2]] / 1000)
        ha_pred += ha_a[sequence[cur + 4]] + ha_b[sequence[cur + 3]] + ha_c[sequence[cur + 1]] + \
            ha_d[sequence[cur]] + (delta_T * ha_t[sequence[cur + 2]] / 1000)

        if use_perdeuteration:
            ca_pred += ca_deut[sequence[cur + 2]]
            cb_pred += cb_deut[sequence[cur + 2]]

        # write to output
        if sequence[cur + 2] == 5:  # special output for gly
            output[cur].update({"CA": __round3(ca_pred, asFloat)})
            output[cur].update({"CB": "**.***"})
            output[cur].update({"CO": __round3(co_pred, asFloat)})
        else:
            output[cur].update({"CA": __round3(ca_pred, asFloat)})
            output[cur].update({"CB": __round3(cb_pred, asFloat)})
            output[cur].update({"CO": __round3(co_pred, asFloat)})

        if sequence[cur + 2] == 12:  # special output for pro
            output[cur].update({"N": "***.***"})
            output[cur].update({"HN": "*.***"})
        else:
            output[cur].update({"N": __round3(n_pred, asFloat)})
            output[cur].update({"HN": __round3(hn_pred, asFloat)})

        if use_perdeuteration:
            output[cur].update({"HA": "*.***"})
        else:
            output[cur].update({"HA": __round3(ha_pred, asFloat)})

        cur += 1

    return output


def __set_sequence(sequence, key1, key3):
    """
    Translates the amino acid input string into a list of integers 
    representing the same set of amino acids.

    Parameters
    ----------
    sequence : str
           The alphabetical list of amino acid abbreviations input by the user.

    key1 : list
         The list of numeric keys used to translate single-letter amino acid 
         abbreviations into representative numbers

    key3 : dictionary
         The dictionary of alphabetical keys used to translate 2/3-letter amino 
         acid abbreviations into representative numbers

    Returns
    -------
    tuple
        Returns a tuple of length two, element one is a list of numeric 
        representatives of input amino acids and element two is a list of
        the abbreviations for those amino acids.

    """
    key_aa1 = key1
    key_aa3 = key3

    i = 0
    code = 0
    inp = sequence

    sequence = []
    aminos = []

    sequence.append(23)
    sequence.append(23)

    # Strip white space at beginning and end
    inp = inp.strip()

    regex = re.findall(r"\(([^)]+)\)|(.)", inp)
    for i in range(len(regex)):
        set = regex[i]
        if set[0] == '':
            aa1 = set[1]
            aminos.append(aa1)
            code = ord(aa1[0]) - 65
            if (code < 0 or code > 35) or (key_aa1[code] == -1):
                continue
            sequence.append(key_aa1[code])
        else:
            aa3 = set[0].upper()
            aminos.append(aa3)
            if aa3 in key_aa3:
                sequence.append(key_aa3[aa3])
            else:
                continue

    sequence.append(23)
    sequence.append(23)

    return (sequence, aminos)


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------
def __round3(num, asFloat=False):
    """

    Performs statistics-consistent rounding of float variables to precisely 
    three decimal places which can be returned as either a string or a float.

    Parameters
    ----------
    num : float
        Float number to be rounded

    asFloat : bool (default=False)
        Determines whether to return rounded number as string or float. False (returns string type) by default

    Returns
    -------
    strng : str
          Rounded number as string

    """
    strng = "" + str(round(num * 1000 + 10**(-len(str(num * 1000)) - 1)) / 1000)
    strng2 = "" + str(round(num + 10**(-len(str(num)) - 1)))
    delta = len(strng) - len(strng2)

    if delta == 0:
        strng += ".000"
    if delta == 2:
        strng += "00"
    if delta == 3:
        strng += "0"

    if asFloat:
        return float(strng)
    else:
        return strng


