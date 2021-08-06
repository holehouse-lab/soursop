import soursop
from soursop import sstrajectory
import pytest
import sys

test_data_dir = soursop.get_data('test_data')

CAPNAMES = ['ACE_NH2', 'ACE_NME', 'ACE_UCAP', 'FOR_NH2',  'FOR_NME',  'FOR_UCAP', 'UCAP_NH2', 'UCAP_NME',  'UCAP_UCAP', 'ACE_NME_start_at_5', 'ACE_NME_multichain']

#CAPNAMES = ['ACE_NME_multichain']

CCAP = ['ACE_NH2', 'ACE_NME', 'FOR_NH2',  'FOR_NME', 'UCAP_NH2', 'UCAP_NME', 'ACE_NME_start_at_5', 'ACE_NME_multichain']
NCAP = ['ACE_NH2', 'ACE_NME', 'ACE_UCAP', 'FOR_NH2',  'FOR_NME',  'FOR_UCAP', 'ACE_NME_start_at_5', 'ACE_NME_multichain']




def test_cap_DSSP():
    """
    Tests DSSP (and hence R1/R2 selection) works
    """
            
    for CN in CAPNAMES:
        CT = sstrajectory.SSTrajectory("%s/cap_tests/%s.pdb"%(test_data_dir, CN), "%s/cap_tests/%s.pdb"%(test_data_dir, CN))
        for CP in CT.proteinTrajectoryList:  
            print("")
            print(CP.get_amino_acid_sequence())
            assert len(CP.get_secondary_structure_DSSP()[0]) == 6
            
                

def test_cap_BBSEG():
    """
    Tests DSSP (and hence R1/R2 selection) works
    """
            
    for CN in CAPNAMES:
        CT = sstrajectory.SSTrajectory("%s/cap_tests/%s.pdb"%(test_data_dir, CN), "%s/cap_tests/%s.pdb"%(test_data_dir, CN))
        for CP in CT.proteinTrajectoryList:
            assert len(CP.get_secondary_structure_DSSP()[0]) == 6
                


def test_cap():
    """
    Test caps assignment works
    """
            
    for CN in CAPNAMES:
        CT = sstrajectory.SSTrajectory("%s/cap_tests/%s.pdb"%(test_data_dir, CN), "%s/cap_tests/%s.pdb"%(test_data_dir, CN))
        for CP in CT.proteinTrajectoryList:
            if CP.ccap:
                assert CN in CCAP
            if CP.ncap:
                assert CN in NCAP

