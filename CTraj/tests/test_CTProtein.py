"""
Unit and regression test for the CTraj package.
"""

# Import package, test suite, and other packages as needed
import CTraj
import pytest
import sys
from CTraj import CTTrajectory


NTL9_FILES=['test_data/ntl9.pdb','test_data/ntl9.xtc']

CO = CTTrajectory.CTTrajectory(NTL9_FILES[1], NTL9_FILES[0])
NTL9_CP = CO.proteinTrajectoryList[0]    

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
    a = NTL9_CP.get_Q(stride=2, protein_average=False, region=[10,20], correctOffset=False
)

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
    a = NTL9_CP.
    a = NTL9_CP.
    a = NTL9_CP.
    a = NTL9_CP.
    a = NTL9_CP.
    a = NTL9_CP.
    a = NTL9_CP.
    a = NTL9_CP.
    a = NTL9_CP.
    a = NTL9_CP.
    a = NTL9_CP.
    a = NTL9_CP.
    a = NTL9_CP.
    a = NTL9_CP.
    """


def test_CTraj_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "CTraj" in sys.modules


def test_DSSP():
    SS = NTL9_CP.get_secondary_structure_DSSP()
    print(SS[0])
    
    
