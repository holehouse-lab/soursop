
import camparitraj
from camparitraj import cttrajectory
import pytest
import sys

GS6_FILES=['gs6.pdb','gs6.xtc']
NTL9_FILES=['ntl9.pdb','ntl9.xtc']
GROMACS_2_CHAINS=['gromacs2chains/top.pdb','gromacs2chains/traj.xtc']


test_data_dir = camparitraj.get_data('test_data')

@pytest.fixture(scope='session', autouse=True)
def GS6_CO(request):    
    GS6_CO = cttrajectory.CTTrajectory("%s/%s"%(test_data_dir, GS6_FILES[1]), "%s/%s"%(test_data_dir, GS6_FILES[0])) 
    return GS6_CO


@pytest.fixture(scope='session', autouse=True)
def GS6_CP(request):
    GS6_CO = cttrajectory.CTTrajectory("%s/%s"%(test_data_dir, GS6_FILES[1]), "%s/%s"%(test_data_dir, GS6_FILES[0])) 
    GS6_CP = GS6_CO.proteinTrajectoryList[0]    
    return GS6_CP


@pytest.fixture(scope='session', autouse=True)
def NTL9_CO(request):    
    NTL9_CO = cttrajectory.CTTrajectory("%s/%s"%(test_data_dir, NTL9_FILES[1]), "%s/%s"%(test_data_dir, NTL9_FILES[0])) 
    return NTL9_CO

@pytest.fixture(scope='session', autouse=True)
def NTL9_CP(request):    
    NTL9_CO = cttrajectory.CTTrajectory("%s/%s"%(test_data_dir, NTL9_FILES[1]), "%s/%s"%(test_data_dir, NTL9_FILES[0])) 
    NTL9_CP = NTL9_CO.proteinTrajectoryList[0]    
    return NTL9_CP

@pytest.fixture(scope='session', autouse=True)
def GMX_2CHAINS(request):    
    GMX_2CHAINS = cttrajectory.CTTrajectory("%s/%s"%(test_data_dir, GROMACS_2_CHAINS[1]), "%s/%s"%(test_data_dir, GROMACS_2_CHAINS[0])) 
    return GMX_2CHAINS
