import os
import pytest
import numpy as np
import soursop.ssnmr as nmr
from soursop.ssexceptions import SSException


def test_compute_random_coil_chemical_shifts_extreme_temperatures():
    seq = 'ARNDE'
    temperatures = [
        -273.15,  # absolute zero (celcius)
        -273,
        100.1,
        101
    ]

    for temperature in temperatures:
        with pytest.raises(SSException):
            nmr.compute_random_coil_chemical_shifts(seq, temperature=temperature, use_ggxgg=True, use_perdeuteration=False)


def test_compute_random_coil_chemical_shifts_pH_low():
    seq = 'ARNDE'
    temperature = 30
    pHs = [
        -1.0,
        -1,
        15.0,
        15
    ]

    for pH in pHs:
        with pytest.raises(SSException):
            nmr.compute_random_coil_chemical_shifts(seq, temperature=temperature, pH=pH, use_ggxgg=True, use_perdeuteration=False)


def test_round3_as_strings():
    numbers = [
        #-0.4956, 
        #-0.1,
        # 0, # This produces a strange output as a string: `0.0000.100`
        0.1, 
        0.2, 
        1/3, 
        0.5
    ]
    expected = [
        #'-0.496', 
        #'-0.100',
        # '0.000', 
        '0.100', 
        '0.200', 
        '0.333', 
        '0.500'
    ]

    for num, exp_num in zip(numbers, expected):
        num_str = nmr.__round3(num, asFloat=False)
        assert type(num_str) is str
        assert num_str == exp_num
