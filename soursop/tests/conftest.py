import os
import sys
import numpy as np
import pytest
import soursop
from soursop import sstrajectory
from soursop import ssprotein


GS6_FILES = ["gs6_AA.pdb", "gs6_AA.xtc"]
CTL9_FILES = ["ctl9_AA.pdb", "ctl9_AA.xtc"]
NTL9_FILES = ["ntl9_AA.pdb", "ntl9_AA.xtc"]
SIGA_CG_FILES = ["sigA_CG.pdb", "sigA_CG.xtc"]
SYNTH_1_CG_FILES = ["synth_1_CG.pdb", "synth_1_CG.xtc"]
ALL_RESIDUES_FILES = ["all_residues_AA.pdb", "all_residues_AA.xtc"]
GROMACS_1_CHAIN = ["gromacs1chain/top.pdb", "gromacs1chain/traj.xtc"]
GROMACS_2_CHAINS = ["gromacs2chains/top.pdb", "gromacs2chains/traj.xtc"]

# SWAN 2-bead (CA/CB) coarse-grained trajectories
SWAN_HELIX_FILES = ["swan_trajectories/helix.pdb", "swan_trajectories/helix.xtc"]
SWAN_ASH1_FILES = ["swan_trajectories/Ash1.pdb", "swan_trajectories/Ash1.xtc"]
SWAN_ASYN_FILES = ["swan_trajectories/alpha_syn.pdb", "swan_trajectories/alpha_syn.xtc"]


test_data_dir = soursop.get_data("test_data")


def _load_traj(files):
    return sstrajectory.SSTrajectory(
        os.path.join(test_data_dir, files[1]),
        os.path.join(test_data_dir, files[0]),
    )


@pytest.fixture(scope="session", autouse=True)
def CTL9_CO(request):
    return _load_traj(CTL9_FILES)


@pytest.fixture(scope="session", autouse=True)
def CTL9_CP(CTL9_CO):
    return CTL9_CO.proteinTrajectoryList[0]


@pytest.fixture(scope="session", autouse=True)
def GS6_CO(request):
    topology_path = os.path.join(test_data_dir, GS6_FILES[0])
    trajectory_path = os.path.join(test_data_dir, GS6_FILES[1])
    GS6_CO = sstrajectory.SSTrajectory(trajectory_path, topology_path)
    return GS6_CO


@pytest.fixture(scope="session", autouse=True)
def GS6_CP(request):
    topology_path = os.path.join(test_data_dir, GS6_FILES[0])
    trajectory_path = os.path.join(test_data_dir, GS6_FILES[1])
    GS6_CO = sstrajectory.SSTrajectory(trajectory_path, topology_path)
    GS6_CP = GS6_CO.proteinTrajectoryList[0]
    return GS6_CP


@pytest.fixture(scope="session", autouse=True)
def NTL9_CO(request):
    topology_path = os.path.join(test_data_dir, NTL9_FILES[0])
    trajectory_path = os.path.join(test_data_dir, NTL9_FILES[1])
    NTL9_CO = sstrajectory.SSTrajectory(trajectory_path, topology_path)
    return NTL9_CO


@pytest.fixture(scope="session", autouse=True)
def NTL9_CP(request):
    topology_path = os.path.join(test_data_dir, NTL9_FILES[0])
    trajectory_path = os.path.join(test_data_dir, NTL9_FILES[1])
    NTL9_CO = sstrajectory.SSTrajectory(trajectory_path, topology_path)
    NTL9_CP = NTL9_CO.proteinTrajectoryList[0]
    return NTL9_CP


@pytest.fixture(scope="session", autouse=True)
def GMX_2CHAINS(request):
    return _load_traj(GROMACS_2_CHAINS)


@pytest.fixture(scope="session", autouse=True)
def GMX_1CHAIN_CO(request):
    return _load_traj(GROMACS_1_CHAIN)


@pytest.fixture(scope="session", autouse=True)
def SIGA_CG_CO(request):
    return _load_traj(SIGA_CG_FILES)


@pytest.fixture(scope="session", autouse=True)
def SYNTH_1_CG_CO(request):
    return _load_traj(SYNTH_1_CG_FILES)


@pytest.fixture(scope="session", autouse=True)
def ALL_RESIDUES_CO(request):
    return _load_traj(ALL_RESIDUES_FILES)


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
            savename = "%s.%s" % ("_".join(name_parts), extension)
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
            data = np.loadtxt(
                f, delimiter=","
            )  # numpy is operating on a buffer b/c of the temporaryfile.
            d = np.frombuffer(data).reshape(
                data.shape
            )  # convert the buffer back to an array, and compare.
            assert np.allclose(data, d)


@pytest.fixture(scope="session", autouse=True)
def cta_protein_helper():
    return ProteinHelper


# ----------------------------------------------------------------------------
# SWAN 2-bead (CA/CB) coarse-grained fixtures
# ----------------------------------------------------------------------------
@pytest.fixture(scope="session")
def SWAN_HELIX_CO(request):
    return _load_traj(SWAN_HELIX_FILES)


@pytest.fixture(scope="session")
def SWAN_HELIX_CP(SWAN_HELIX_CO):
    return SWAN_HELIX_CO.proteinTrajectoryList[0]


@pytest.fixture(scope="session")
def SWAN_ASH1_CO(request):
    return _load_traj(SWAN_ASH1_FILES)


@pytest.fixture(scope="session")
def SWAN_ASH1_CP(SWAN_ASH1_CO):
    return SWAN_ASH1_CO.proteinTrajectoryList[0]


@pytest.fixture(scope="session")
def SWAN_ASYN_CO(request):
    return _load_traj(SWAN_ASYN_FILES)


@pytest.fixture(scope="session")
def SWAN_ASYN_CP(SWAN_ASYN_CO):
    return SWAN_ASYN_CO.proteinTrajectoryList[0]
