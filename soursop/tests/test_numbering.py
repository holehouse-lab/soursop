import os
import sys
import pytest
import soursop
from soursop import sstrajectory


#NAMES = ['ACE_NH2', 'ACE_NME', 'ACE_UCAP', 'FOR_NH2',  'FOR_NME',  'FOR_UCAP', 'UCAP_NH2', 'UCAP_NME',  'UCAP_UCAP', 'ACE_NME_start_at_5', 'ACE_NME_multichain']

#CAPNAMES = ['ACE_NME_multichain']

#CCAP = ['ACE_NH2', 'ACE_NME', 'FOR_NH2',  'FOR_NME', 'UCAP_NH2', 'UCAP_NME', 'ACE_NME_start_at_5', 'ACE_NME_multichain']
#NCAP = ['ACE_NH2', 'ACE_NME', 'ACE_UCAP', 'FOR_NH2',  'FOR_NME',  'FOR_UCAP', 'ACE_NME_start_at_5', 'ACE_NME_multichain']




def test_cap_DSSP():
    """
    Tests DSSP (and hence R1/R2 selection) works
    """

    test_data_dir = soursop.get_data('test_data')
    pdb_files = ['ACE_NME_multichain', 'ACE_NME_start_at_5']
    for cap_name in pdb_files:
        cap_path = os.path.join(test_data_dir, 'cap_tests', '{}.pdb'.format(cap_name))
        cap_trajectory = sstrajectory.SSTrajectory(cap_path, cap_path)
        # print(cap_trajectory)
        # for cap_protein in cap_trajectory.proteinTrajectoryList:
        #     res_string = ''
        #     for i in cap_protein.topology.residues:
        #         #res_string = res_string + " %s " %(str(i.index))
        #         print(i)

        #     #atom_string = ''
        #     #for i in CP.topology.atoms:
        #     #    atom_string = atom_string + " " + str(i)

        #     print(cap_protein)
        #     print(res_string)
        #     print("Residue 1: %s" %(cap_protein.topology.select('residue 1')))
        #     print("Resid 1: %s" %(cap_protein.topology.select('resid 1')))
            
        #     # print(CP.atom_offset)
        #     # print(atom_string)
            

