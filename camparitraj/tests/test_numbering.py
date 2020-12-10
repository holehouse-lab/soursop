import camparitraj
from camparitraj import cttrajectory
import pytest
import sys


test_data_dir = camparitraj.get_data('test_data')

#NAMES = ['ACE_NH2', 'ACE_NME', 'ACE_UCAP', 'FOR_NH2',  'FOR_NME',  'FOR_UCAP', 'UCAP_NH2', 'UCAP_NME',  'UCAP_UCAP', 'ACE_NME_start_at_5', 'ACE_NME_multichain']

#CAPNAMES = ['ACE_NME_multichain']

#CCAP = ['ACE_NH2', 'ACE_NME', 'FOR_NH2',  'FOR_NME', 'UCAP_NH2', 'UCAP_NME', 'ACE_NME_start_at_5', 'ACE_NME_multichain']
#NCAP = ['ACE_NH2', 'ACE_NME', 'ACE_UCAP', 'FOR_NH2',  'FOR_NME',  'FOR_UCAP', 'ACE_NME_start_at_5', 'ACE_NME_multichain']

NAMES = ['ACE_NME_multichain', 'ACE_NME_start_at_5']


def test_cap_DSSP():
    """
    Tests DSSP (and hence R1/R2 selection) works
    """
            
    for CN in NAMES:
        CT = cttrajectory.CTTrajectory("%s/cap_tests/%s.pdb"%(test_data_dir, CN), "%s/cap_tests/%s.pdb"%(test_data_dir, CN))
        for CP in CT.proteinTrajectoryList:  

            res_string = ''
            for i in CP.topology.residues:
                res_string = res_string + " %s " %(str(i.index))

            #atom_string = ''
            #for i in CP.topology.atoms:
            #    atom_string = atom_string + " " + str(i)

            print(CP)
            print(res_string)
            print("Residue 1: %s" %(CP.topology.select('residue 1')))
            print("Resid 1: %s" %(CP.topology.select('resid 1')))
            
            # print(CP.atom_offset)
            # print(atom_string)
            
            print('')

    return CP
                        

                

 



CP = test_cap_DSSP()
