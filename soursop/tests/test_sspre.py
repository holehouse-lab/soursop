# Import package, test suite, and other packages as needed
import os

import numpy as np
import soursop
import pytest
from contextlib import contextmanager
from soursop.sspre import SSPRE


test_data_dir = soursop.get_data('test_data')

def test_pre_profiles_ctl9(CTL9_CP):

    # check we can correctly compute PRE profiles

    # PRE settings for Stenzowski et al 2020 (note in original calulations
    # we used a tau_c of 5 ns (ultimately used 4 ns for the paper but the saved
    # profiles here used 5)
    PRE = SSPRE(CTL9_CP, 5, 12, 14, 850000000)


    # for each label position
    for P in [61, 74, 96, 109, 119, 149]:
        corrected_val = P-58

        # calculate PRE profile
        vals = PRE.generate_PRE_profile(corrected_val)

        # get previously-calculated profile
        file_path = os.path.join(test_data_dir, f'pre_data/PRE_{P}_unweighted.csv')        
        pre_test_data = np.loadtxt(file_path)

        # ensure the largest differences is less than 0.02. There is some numerical inprecision
        # and some small changes in how the calculations are done between 2019 and the current
        # version of SOURSOP but these do not materially influence the profiles
        assert max(np.abs(pre_test_data - vals[0])) < 0.02
    
