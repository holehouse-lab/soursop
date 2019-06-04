"""
Unit and regression test for the camparitraj package.
"""

# Import package, test suite, and other packages as needed
import numpy as np
import camparitraj
import pytest
import sys
from camparitraj import cttrajectory
"""
NTL9_FILES=['ntl9.pdb','ntl9.xtc']
GS6_FILES=['gs6.pdb','gs6.xtc']

test_data_dir = camparitraj.get_data('test_data')

NTL9_CO = cttrajectory.CTTrajectory("%s/%s"%(test_data_dir, NTL9_FILES[1]), "%s/%s"%(test_data_dir, NTL9_FILES[0])) 
NTL9_CP = NTL9_CO.proteinTrajectoryList[0]    

GS6_CO = cttrajectory.CTTrajectory("%s/%s"%(test_data_dir, GS6_FILES[1]), "%s/%s"%(test_data_dir, GS6_FILES[0])) 
GS6_CP = GS6_CO.proteinTrajectoryList[0]    
"""
#@pytest.fixture()
#def ntl9_setup(request):
    

"""
def test_code_coverage():

    NTL9_CP.print_residues()
    a = NTL9_CP.get_residue_index_list()
    a = NTL9_CP.get_numberOfResidues()
    a = NTL9_CP.get_numberOfFrames()
    a = NTL9_CP.get_aminoAcidSequence()
    a = NTL9_CP.get_CAindex(10)
    a = NTL9_CP.get_CAindex(10, correctOffset=False)
    a = NTL9_CP.get_multiple_CAindex()
    a = NTL9_CP.get_multiple_CAindex([3,4,5])
    a = NTL9_CP.get_multiple_CAindex([3,4,5], correctOffset=False)

    a = NTL9_CP.calculateAllCAdistances(10)
    a = NTL9_CP.calculateAllCAdistances(10, correctOffset=False)
    a = NTL9_CP.calculateAllCAdistances(10,stride=2)
    a = NTL9_CP.calculateAllCAdistances(10,stride=2)
    a = NTL9_CP.calculateAllCAdistances(10,stride=2, mode='COM')
    a = NTL9_CP.calculateAllCAdistances(10,stride=2, mode='COM', onlyCterminalResidues=False)

    a = NTL9_CP.get_distanceMap()
    a = NTL9_CP.get_distanceMap(verbose=True)
    a = NTL9_CP.get_distanceMap(verbose=True,mode='COM')
    a = NTL9_CP.get_distanceMap(verbose=True,mode='COM',RMS=True)

    a = NTL9_CP.get_polymer_scaled_distance_map()
    a = NTL9_CP.get_polymer_scaled_distance_map(nu=0.54,A0=6)
    a = NTL9_CP.get_polymer_scaled_distance_map(nu=0.54,A0=6,min_separation=20)
    a = NTL9_CP.get_polymer_scaled_distance_map(nu=0.54,A0=6,min_separation=20, mode='signed-fractional-change')
    a = NTL9_CP.get_polymer_scaled_distance_map(nu=0.54,A0=6,min_separation=20, mode='signed-absolute-change')
    a = NTL9_CP.get_polymer_scaled_distance_map(nu=0.54,A0=6,min_separation=20, mode='scaled')

    
    a = NTL9_CP.get_local_heterogeneity(stride=1)
    a = NTL9_CP.get_local_heterogeneity(fragmentSize=5, stride=1)
    a = NTL9_CP.get_local_heterogeneity(fragmentSize=20, stride=2)
    a = NTL9_CP.get_local_heterogeneity(fragmentSize=20, stride=1)
    a = NTL9_CP.get_DVector(stride=1)

    a = NTL9_CP.get_RMSD(1)
    a = NTL9_CP.get_RMSD(1,frame2=3)
    a = NTL9_CP.get_RMSD(1,frame2=3, stride=2)
    a = NTL9_CP.get_RMSD(1,frame2=3, stride=2, region=[10,20])
    a = NTL9_CP.get_RMSD(1,frame2=3, stride=2, backbone=False)
    a = NTL9_CP.get_RMSD(1,frame2=3, stride=2, correctOffset=False)

    a = NTL9_CP.get_Q()
    a = NTL9_CP.get_Q(stride=2)
    a = NTL9_CP.get_Q(stride=2, protein_average=False)
    a = NTL9_CP.get_Q(stride=2, protein_average=False, region=[10,20])
    a = NTL9_CP.get_Q(stride=2, protein_average=False, region=[10,20], correctOffset=False)

    a = NTL9_CP.get_contact_map()
    a = NTL9_CP.get_contact_map(distance_thresh=2)
    a = NTL9_CP.get_contact_map(distance_thresh=2, stride=2)

    a = NTL9_CP.get_clusters(stride=1)
    a = NTL9_CP.get_clusters(stride=1, n_clusters=3)
    a = NTL9_CP.get_clusters(stride=1, n_clusters=3, backbone=False)
    a = NTL9_CP.get_clusters(stride=1, n_clusters=3, backbone=False, correctOffset=False)

    a = NTL9_CP.get_interResidueCOMDistance(10,30)
    a = NTL9_CP.get_interResidueCOMDistance(10,30, stride=2)
    a = NTL9_CP.get_interResidueCOMDistance(10,30, stride=2, correctOffset=False)

    a = NTL9_CP.get_interResidueCOMVector(10,30)
    a = NTL9_CP.get_interResidueCOMVector(10,30, stride=2)
    a = NTL9_CP.get_interResidueCOMVector(10,30, stride=2, correctOffset=False)

    a = NTL9_CP.get_interResidue_atomic_distance(2,10)
    a = NTL9_CP.get_interResidue_atomic_distance(2,10, A1='CB')
"""

def test_get_radius_of_gyration(GS6_CP):    

    assert len(GS6_CP.get_radius_of_gyration()) == 5
    assert abs(GS6_CP.get_radius_of_gyration()[0] - 5.728453763896514) < 0.0001
    assert abs(GS6_CP.get_radius_of_gyration(R1=1,R2=3)[0] - 2.8815625777553278) < 0.0001

def test_get_t(GS6_CP):    

    assert len(GS6_CP.get_t()) == 5
    assert abs(GS6_CP.get_t()[0] - 0.30286034750848245) < 0.0001
    assert abs(GS6_CP.get_t(R1=1,R2=3)[0] - 0.0766270900235029) < 0.0001

def test_get_internal_scaling(GS6_CP):    

    assert len(GS6_CP.get_internal_scaling()[0]) == 6
    assert abs(np.mean(GS6_CP.get_internal_scaling()[1][1]) - 3.784701967415692) < 0.0001
    assert abs(np.mean(GS6_CP.get_internal_scaling()[1][2]) - 6.524909331438218) < 0.0001
    assert abs(np.mean(GS6_CP.get_internal_scaling(R1=1,R2=4)[1][1]) - 3.7353248) < 0.001


def test_get_internal_scaling_RMS(GS6_CP):    

    assert len(GS6_CP.get_internal_scaling_RMS()[0]) == 6
    assert abs(np.mean(GS6_CP.get_internal_scaling_RMS()[1][1]) - 3.793829294427209) < 0.0001
    assert abs(np.mean(GS6_CP.get_internal_scaling_RMS()[1][2]) - 6.578195001335268) < 0.0001
    assert abs(np.mean(GS6_CP.get_internal_scaling_RMS(R1=1,R2=4)[1][1]) - 3.7455646317574534) < 0.001


def test_get_internal_scaling_RMS(GS6_CP):    

    assert len(GS6_CP.get_internal_scaling_RMS()[0]) == 6
    assert abs(np.mean(GS6_CP.get_internal_scaling_RMS()[1][1]) - 3.793829294427209) < 0.0001
    assert abs(np.mean(GS6_CP.get_internal_scaling_RMS()[1][2]) - 6.578195001335268) < 0.0001
    assert abs(np.mean(GS6_CP.get_internal_scaling_RMS(R1=1,R2=4)[1][1]) - 3.7455646317574534) < 0.001


def test_camparitraj_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "camparitraj" in sys.modules


def test_DSSP(NTL9_CP):
    SS = NTL9_CP.get_secondary_structure_DSSP()
    print(SS[0])
    
    
