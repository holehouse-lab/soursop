import os
import sys
import pytest
import camparitraj
from camparitraj import cttrajectory


GS6_FILES=['gs6.pdb','gs6.xtc']
NTL9_FILES=['ntl9.pdb','ntl9.xtc']
GROMACS_2_CHAINS=['gromacs2chains/top.pdb','gromacs2chains/traj.xtc']


test_data_dir = camparitraj.get_data('test_data')


@pytest.fixture(scope='session', autouse=True)
def GS6_CO(request):
    topology_path = os.path.join(test_data_dir, GS6_FILES[0])
    trajectory_path = os.path.join(test_data_dir, GS6_FILES[1])
    GS6_CO = cttrajectory.CTTrajectory(trajectory_path, topology_path)
    return GS6_CO


@pytest.fixture(scope='session', autouse=True)
def GS6_CP(request):
    topology_path = os.path.join(test_data_dir, GS6_FILES[0])
    trajectory_path = os.path.join(test_data_dir, GS6_FILES[1])
    GS6_CO = cttrajectory.CTTrajectory(trajectory_path, topology_path)
    GS6_CP = GS6_CO.proteinTrajectoryList[0]
    return GS6_CP


@pytest.fixture(scope='session', autouse=True)
def NTL9_CO(request):
    topology_path = os.path.join(test_data_dir, NTL9_FILES[0])
    trajectory_path = os.path.join(test_data_dir, NTL9_FILES[1])
    NTL9_CO = cttrajectory.CTTrajectory(trajectory_path, topology_path)
    return NTL9_CO


@pytest.fixture(scope='session', autouse=True)
def NTL9_CP(request):
    topology_path = os.path.join(test_data_dir, NTL9_FILES[0])
    trajectory_path = os.path.join(test_data_dir, NTL9_FILES[1])
    NTL9_CO = cttrajectory.CTTrajectory(trajectory_path, topology_path)
    NTL9_CP = NTL9_CO.proteinTrajectoryList[0]
    return NTL9_CP

@pytest.fixture(scope='session', autouse=True)
def GMX_2CHAINS(request):    
    GMX_2CHAINS = cttrajectory.CTTrajectory("%s/%s"%(test_data_dir, GROMACS_2_CHAINS[1]), "%s/%s"%(test_data_dir, GROMACS_2_CHAINS[0])) 
    return GMX_2CHAINS
