import os
import sys
import numpy as np
import pytest
import soursop
from soursop import sstrajectory
from soursop import ssprotein


GS6_FILES=['gs6.pdb','gs6.xtc']
CTL9_FILES=['ctl9.pdb','ctl9.xtc']
NTL9_FILES=['ntl9.pdb','ntl9.xtc']
GROMACS_2_CHAINS=['gromacs2chains/top.pdb','gromacs2chains/traj.xtc']


test_data_dir = soursop.get_data('test_data')


@pytest.fixture(scope='session', autouse=True)
def CTL9_CP(request):
    topology_path = os.path.join(test_data_dir, CTL9_FILES[0])
    trajectory_path = os.path.join(test_data_dir, CTL9_FILES[1])
    
    CTL9_CP = sstrajectory.SSTrajectory(trajectory_path, topology_path).proteinTrajectoryList[0]
    
    return CTL9_CP


@pytest.fixture(scope='session', autouse=True)
def GS6_CO(request):
    topology_path = os.path.join(test_data_dir, GS6_FILES[0])
    trajectory_path = os.path.join(test_data_dir, GS6_FILES[1])
    GS6_CO = sstrajectory.SSTrajectory(trajectory_path, topology_path)
    return GS6_CO


@pytest.fixture(scope='session', autouse=True)
def GS6_CP(request):
    topology_path = os.path.join(test_data_dir, GS6_FILES[0])
    trajectory_path = os.path.join(test_data_dir, GS6_FILES[1])
    GS6_CO = sstrajectory.SSTrajectory(trajectory_path, topology_path)
    GS6_CP = GS6_CO.proteinTrajectoryList[0]
    return GS6_CP


@pytest.fixture(scope='session', autouse=True)
def NTL9_CO(request):
    topology_path = os.path.join(test_data_dir, NTL9_FILES[0])
    trajectory_path = os.path.join(test_data_dir, NTL9_FILES[1])
    NTL9_CO = sstrajectory.SSTrajectory(trajectory_path, topology_path)
    return NTL9_CO


@pytest.fixture(scope='session', autouse=True)
def NTL9_CP(request):
    topology_path = os.path.join(test_data_dir, NTL9_FILES[0])
    trajectory_path = os.path.join(test_data_dir, NTL9_FILES[1])
    NTL9_CO = sstrajectory.SSTrajectory(trajectory_path, topology_path)
    NTL9_CP = NTL9_CO.proteinTrajectoryList[0]
    return NTL9_CP

@pytest.fixture(scope='session', autouse=True)
def GMX_2CHAINS(request):
    GMX_2CHAINS = sstrajectory.SSTrajectory("%s/%s"%(test_data_dir, GROMACS_2_CHAINS[1]), "%s/%s"%(test_data_dir, GROMACS_2_CHAINS[0]))
    return GMX_2CHAINS


# This is implemented for use in unittests for `ssanalyzer`.
# Adapted from:
# https://stackoverflow.com/questions/33508060/create-and-import-helper-functions-in-tests-without-creating-packages-in-test-di
class ProteinHelper:
    @staticmethod
    def determine_filenames(prefix, names, extension):
        expected_filenames = list()
        for name in names:
            name_parts = [prefix]
            if len(name) > 0:
                name_parts.append(name)
            savename = '%s.%s' % ('_'.join(name_parts), extension)
            expected_filenames.append(savename)
        return expected_filenames

    @staticmethod
    def lengthen_protein_trajectory(protein_traj, num_copies):
        trajectories = [protein_traj.traj for i in range(num_copies)]
        traj = protein_traj.traj.join(trajectories)
        protein = ssprotein.SSProtein(traj)
        return protein

    @staticmethod
    def validate_exported_csv_data(savename):
        assert os.path.getsize(savename) > 0  # implicit check for file existence
        with open(savename) as f:
            data = np.loadtxt(f, delimiter=',')  # numpy is operating on a buffer b/c of the temporaryfile.
            d = np.frombuffer(data).reshape(data.shape)  # convert the buffer back to an array, and compare.
            assert np.allclose(data, d)


@pytest.fixture(scope='session', autouse=True)
def cta_protein_helper():
    return ProteinHelper
