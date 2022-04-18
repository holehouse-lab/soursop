"""
Unit and regression test for the sstrajectory module.
"""
# Import package, test suite, and other packages as needed
import soursop
import hashlib
from soursop import sstrajectory
from soursop.ssexceptions import SSException
from pathlib import Path
import pytest
import sys
import random
import itertools
import platform
import tempfile
import numpy as np
import os
from soursop.configs import TMP_DIR


from . import conftest
#import conftest

# -------------------------------------------------------------------------------------------------
# SSTrajectory `__init__` tests.
#
# Also implicitly tests:
#       `__readTrajectory`
#       `__get_proteins`
#       `__get_proteins_by_residue`

def test_trajectory_repr_string(GS6_CO):
    repr_string = (hex(id(GS6_CO)), GS6_CO.n_proteins, GS6_CO.n_frames)
    return repr_string == repr(GS6_CO)


def test_trajectory_len(GS6_CO):
    assert len(GS6_CO) == len(GS6_CO.traj)
    assert len(GS6_CO) == GS6_CO.n_frames


def test_read_in_no_trajectory_and_no_topology():
    with pytest.raises(SSException):
        sstrajectory.SSTrajectory()


def test_read_in_no_trajectory_and_no_topology_but_use_custom_trajectory(GS6_CO):
    trajectory = GS6_CO.traj
    sstrajectory.SSTrajectory(TRJ=trajectory)


def test_read_in_no_trajectory():
    pdb_filename = os.path.join(soursop.get_data('test_data'), 'gs6.pdb')
    with pytest.raises(SSException):
        sstrajectory.SSTrajectory(pdb_filename=pdb_filename)


def test_read_in_trajectory_pdblead(GS6_CO):
    pdb_filename = os.path.join(soursop.get_data('test_data'), 'gs6.pdb')
    traj_filename = os.path.join(soursop.get_data('test_data'), 'gs6.xtc')
    traj = sstrajectory.SSTrajectory(pdb_filename=pdb_filename, trajectory_filename=traj_filename, pdblead=True)

    # Since we're using the PDB as an initial frame, the number of frames should have increased by 1.
    assert len(traj) == len(GS6_CO) + 1


def test_trajectory_initialization_debug():
    pdb_filename = os.path.join(soursop.get_data('test_data'), 'gs6.pdb')
    traj_filename = os.path.join(soursop.get_data('test_data'), 'gs6.xtc')
    sstrajectory.SSTrajectory(pdb_filename=pdb_filename, trajectory_filename=traj_filename, debug=True)


def test_read_in_no_topology_xtc():
    traj_filename = os.path.join(soursop.get_data('test_data'), 'gs6.xtc')
    with pytest.raises(SSException):
        sstrajectory.SSTrajectory(trajectory_filename=traj_filename)


def test_read_in_no_topology_dcd():
    """The `SSTrajectory.__init__` references that `.dcd` files are supported too."""
    traj_filename = os.path.join(soursop.get_data('test_data'), 'gs6.dcd')
    with pytest.raises(SSException):
        sstrajectory.SSTrajectory(trajectory_filename=traj_filename)


def test_read_in_compare_trajectories(GS6_CO):
    pdb_filename = os.path.join(soursop.get_data('test_data'), 'gs6.pdb')
    traj_filename = os.path.join(soursop.get_data('test_data'), 'gs6.xtc')
    trajectory = sstrajectory.SSTrajectory(trajectory_filename=traj_filename, pdb_filename=pdb_filename)

    assert np.allclose(trajectory.traj.xyz, GS6_CO.traj.xyz)


def test_read_in_protein_grouping_simple():
    protein_groups = [[0]]
    pdb_filename = os.path.join(soursop.get_data('test_data'), 'gs6.pdb')
    traj_filename = os.path.join(soursop.get_data('test_data'), 'gs6.xtc')
    trajectory = sstrajectory.SSTrajectory(trajectory_filename=traj_filename,
                                           pdb_filename=pdb_filename,
                                           protein_grouping=protein_groups)

    # verify that we have loaded the number of proteins expected
    assert trajectory.n_proteins == len(protein_groups)


def test_read_in_protein_grouping_multiple():
    protein_groups = [[0, 1, 2], [3, 4, 5], [6, 7, 8]]
    pdb_filename = os.path.join(soursop.get_data('test_data'), 'ntl9.pdb')
    traj_filename = os.path.join(soursop.get_data('test_data'), 'ntl9.xtc')  # ntl9 has 56 residues
    trajectory = sstrajectory.SSTrajectory(trajectory_filename=traj_filename,
                                           pdb_filename=pdb_filename,
                                           protein_grouping=protein_groups)

    # verify that we have loaded the number of proteins expected
    assert trajectory.n_proteins == len(protein_groups)


def test_read_in_protein_grouping_invalid_residues():
    protein_groups = [[1000, 1001, 1002], [1003, 1004, 1005], [1006, 1007, 1008]]
    pdb_filename = os.path.join(soursop.get_data('test_data'), 'ntl9.pdb')
    traj_filename = os.path.join(soursop.get_data('test_data'), 'ntl9.xtc')  # ntl9 has 56 residues
    trajectory = sstrajectory.SSTrajectory(trajectory_filename=traj_filename,
                                           pdb_filename=pdb_filename,
                                           protein_grouping=protein_groups)

    # Since this is a failing but non-disruptive test (i.e. no Exceptions), we
    # check the number of proteins which should be 0.
    assert trajectory.n_proteins == 0


def test_read_in_protein_grouping_multiple_mixed_order():
    # TODO: This test passes - it shouldn't. Look into this further.
    protein_groups = [[1, 2, 0], [10, 3, 7], [11, 1, 3]]
    pdb_filename = os.path.join(soursop.get_data('test_data'), 'ntl9.pdb')
    traj_filename = os.path.join(soursop.get_data('test_data'), 'ntl9.xtc')  # ntl9 has 56 residues
    trajectory = sstrajectory.SSTrajectory(trajectory_filename=traj_filename,
                                           pdb_filename=pdb_filename,
                                           protein_grouping=protein_groups)

    # verify that we have loaded the number of proteins expected
    assert trajectory.n_proteins == len(protein_groups)


# -------------------------------------------------------------------------------------------------


def test_get_intra_chain_distance_map_zero(GS6_CO):
    # TODO: This returns an upper triangle copy of the original distance map. Investigate further.
    distance_map_zero, stddev_map_zero = GS6_CO.get_interchain_distance_map(0, 0)
    distance_map, stddev_map = GS6_CO.proteinTrajectoryList[0].get_distance_map()

    assert np.allclose(np.triu(distance_map_zero), distance_map)
    assert np.allclose(np.triu(stddev_map_zero), stddev_map, atol=1e-6)


def test_get_intra_chain_distance_map_zero_using_center_of_mass(GS6_CO):
    # TODO: This returns an upper triangle copy of the original distance map. Investigate further.
    distance_map_zero, stddev_map_zero = GS6_CO.get_interchain_distance_map(0, 0, mode='COM')
    distance_map, stddev_map = GS6_CO.proteinTrajectoryList[0].get_distance_map(mode='COM')

    assert np.allclose(np.triu(distance_map_zero), distance_map)
    assert np.allclose(np.triu(stddev_map_zero), stddev_map, atol=1e-6)


def test_get_intra_chain_distance_map_protein_groups():
    protein_groups = [[0, 1, 2, 3, 4], [5, 6, 7, 8, 9], [10, 11, 12, 13, 14]]
    pdb_filename = os.path.join(soursop.get_data('test_data'), 'ntl9.pdb')
    traj_filename = os.path.join(soursop.get_data('test_data'), 'ntl9.xtc')  # ntl9 has 56 residues
    trajectory = sstrajectory.SSTrajectory(trajectory_filename=traj_filename,
                                           pdb_filename=pdb_filename,
                                           protein_grouping=protein_groups)

    for protein_group in itertools.combinations(range(len(protein_groups)), r=2):
        protein_a, protein_b = protein_group
        distance_map, stddev_map = trajectory.get_interchain_distance_map(protein_a, protein_b)
        assert distance_map.shape == stddev_map.shape
        assert np.count_nonzero(distance_map) == len(distance_map.flatten())  # since the residue indices are unique


def test_get_intra_chain_distance_map_protein_groups_with_residue_indices():
    # Note that the residues for the resID1 and resID2 indices are with respect to the residues of the protein chain
    # and not the full protein itself.
    protein_groups_residues = [[0, 1, 2, 3, 4, 5, 6, 7, 8],
                               [10, 11, 12, 13, 14, 15, 16, 17, 18],
                               [20, 21, 22, 23, 24, 25, 26, 27, 28],
                               [30, 31, 32, 33, 34, 35, 36, 37, 38],
                               [40, 41, 42, 43, 44, 45, 46, 47, 48]]
    pdb_filename = os.path.join(soursop.get_data('test_data'), 'ntl9.pdb')
    traj_filename = os.path.join(soursop.get_data('test_data'), 'ntl9.xtc')  # ntl9 has 56 residues
    trajectory = sstrajectory.SSTrajectory(trajectory_filename=traj_filename,
                                           pdb_filename=pdb_filename,
                                           protein_grouping=protein_groups_residues)

    # select pairwise protein ids (no repetitions).
    for protein_group in itertools.combinations(range(len(protein_groups_residues)), r=2):
        protein_a, protein_b = protein_group

        # now, let's randomly choose a set of residues to select for `protein_a` and `protein_b`
        num_residues_a = random.randint(1, len(protein_groups_residues[protein_a]))
        num_residues_b = random.randint(1, len(protein_groups_residues[protein_b]))

        protein_a_residues = protein_groups_residues[protein_a]
        a_residues = list(sorted(random.sample(protein_a_residues, num_residues_a)))
        a_residues = [(r - protein_a_residues[0]) for r in a_residues]  # offset by the starting residue number

        protein_b_residues = protein_groups_residues[protein_b]
        b_residues = list(sorted(random.sample(protein_b_residues, num_residues_b)))
        b_residues = [(r - protein_b_residues[0]) for r in b_residues]  # offset by the starting residue number

        distance_map, stddev_map = trajectory.get_interchain_distance_map(protein_a, protein_b)
        assert distance_map.shape == stddev_map.shape


def test_export_intra_chain_distance_map_protein_groups_with_residue_indices():
    # Note that the residues for the resID1 and resID2 indices are with respect to the residues of the protein chain
    # and not the full protein itself.
    protein_groups_residues = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                               [10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
                               [20, 21, 22, 23, 24, 25, 26, 27, 28, 29],
                               [30, 31, 32, 33, 34, 35, 36, 37, 38, 39],
                               [40, 41, 42, 43, 44, 45, 46, 47, 48, 49]]
    pdb_filename = os.path.join(soursop.get_data('test_data'), 'ntl9.pdb')
    traj_filename = os.path.join(soursop.get_data('test_data'), 'ntl9.xtc')  # ntl9 has 56 residues
    trajectory = sstrajectory.SSTrajectory(trajectory_filename=traj_filename,
                                           pdb_filename=pdb_filename,
                                           protein_grouping=protein_groups_residues)

    # select pairwise protein ids (no repetitions).
    for protein_index, protein_group in enumerate(itertools.combinations(range(len(protein_groups_residues)), r=2)):
        protein_a, protein_b = protein_group

        # now, let's randomly choose a set of residues to select for `protein_a` and `protein_b`
        num_residues_a = random.randint(1, len(protein_groups_residues[protein_a]))
        num_residues_b = random.randint(1, len(protein_groups_residues[protein_b]))

        protein_a_residues = protein_groups_residues[protein_a]
        a_residues = list(sorted(random.sample(protein_a_residues, num_residues_a)))
        a_residues = [(r - protein_a_residues[0]) for r in a_residues]  # offset by the starting residue number

        protein_b_residues = protein_groups_residues[protein_b]
        b_residues = list(sorted(random.sample(protein_b_residues, num_residues_b)))
        b_residues = [(r - protein_b_residues[0]) for r in b_residues]  # offset by the starting residue number

        distance_map, stddev_map = trajectory.get_interchain_distance_map(protein_a, protein_b)
        assert distance_map.shape == stddev_map.shape

        temp_save_filepath = os.path.join(TMP_DIR, 'gs6_map_{index}_{protein}'.format(index=protein_index,
                                                                                      protein=protein_index))
        temp_filename = '{filename}.csv'.format(filename=temp_save_filepath)


def test_intrachain_inter_residue_atomic_distance():

    # A custom version of NTL9 is needed, since we want to have multiple protein chains
    protein_groups_residues = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                               [10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
                               [20, 21, 22, 23, 24, 25, 26, 27, 28, 29],
                               [30, 31, 32, 33, 34, 35, 36, 37, 38, 39],
                               [40, 41, 42, 43, 44, 45, 46, 47, 48, 49]]
    pdb_filename = os.path.join(soursop.get_data('test_data'), 'ntl9.pdb')
    traj_filename = os.path.join(soursop.get_data('test_data'), 'ntl9.xtc')  # ntl9 has 56 residues
    trajectory = sstrajectory.SSTrajectory(trajectory_filename=traj_filename,
                                           pdb_filename=pdb_filename,
                                           protein_grouping=protein_groups_residues)

    for protein_group in itertools.combinations(range(len(protein_groups_residues)), r=2):
        protein_a, protein_b = protein_group
        a_residues = protein_groups_residues[protein_a]
        b_residues = protein_groups_residues[protein_b]

        residue_a = random.randint(0, len(a_residues) - 1)
        residue_b = random.randint(0, len(b_residues) - 1)

        distances = trajectory.get_interchain_distance(protein_a, protein_b, residue_a, residue_b)
        assert len(distances) == len(a_residues)


def test_intrachain_inter_residue_atomic_distance_by_name():
    """
    # TODO: Revise.
    # A custom version of NTL9 is needed, since we want to have multiple protein chains
    protein_groups_residues = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                               [10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
                               [20, 21, 22, 23, 24, 25, 26, 27, 28, 29],
                               [30, 31, 32, 33, 34, 35, 36, 37, 38, 39],
                               [40, 41, 42, 43, 44, 45, 46, 47, 48, 49]]
    pdb_filename = os.path.join(soursop.get_data('test_data'), 'ntl9.pdb')
    traj_filename = os.path.join(soursop.get_data('test_data'), 'ntl9.xtc')  # ntl9 has 56 residues
    trajectory = sstrajectory.SSTrajectory(trajectory_filename=traj_filename,
                                           pdb_filename=pdb_filename,
                                           protein_grouping=protein_groups_residues)

    atom_names = dict()
    for residue in trajectory.traj.topology.residues:
        atom_names[residue.index] = [atom for atom in residue.atoms]#[atom.name for atom in residue.atoms]

    for protein_group in itertools.combinations(range(len(protein_groups_residues)), r=2):
        protein_a, protein_b = protein_group
        a_residues = protein_groups_residues[protein_a]
        b_residues = protein_groups_residues[protein_b]

        residue_a = random.randint(0, len(a_residues) - 1)
        residue_a_atoms = atom_names[protein_groups_residues[protein_a][residue_a]]
        residue_a_atom = random.choice(residue_a_atoms)

        residue_b = random.randint(0, len(b_residues) - 1)
        residue_b_atoms = atom_names[protein_groups_residues[protein_b][residue_b]]
        residue_b_atom = random.choice(residue_b_atoms)

        distances = trajectory.get_intrachain_interResidue_atomic_distance(protein_a, protein_b,
                                                                           residue_a, residue_b,
                                                                           A1=residue_a_atom, A2=residue_a_atom)
        assert len(distances) == len(a_residues)
    """
    pass


def test_get_interchain_distance():
    protein_groups = [[0, 1, 2, 3, 4], [5, 6, 7, 8, 9], [10, 11, 12, 13, 14]]
    pdb_filename = os.path.join(soursop.get_data('test_data'), 'ntl9.pdb')
    traj_filename = os.path.join(soursop.get_data('test_data'), 'ntl9.xtc')  # ntl9 has 56 residues
    trajectory = sstrajectory.SSTrajectory(trajectory_filename=traj_filename,
                                           pdb_filename=pdb_filename,
                                           protein_grouping=protein_groups)

    for protein_group in itertools.combinations(range(len(protein_groups)), r=2):
        protein_a, protein_b = protein_group

        # There appears to be an issue where R1 > R2 - the calculation crashes for sidechain-heavy.
        R1 = 0 #random.choice(range(len(protein_groups[protein_a])))
        R2 = 1 #random.choice(range(len(protein_groups[protein_b])))
        A1 = 'CA'
        A2 = 'CA'
        modes = 'atom,ca,closest,closest-heavy,sidechain,sidechain-heavy'.split(',')
        for mode in modes:
            distances = trajectory.get_interchain_distance(protein_a, protein_b, R1, R2, A1, A2, mode)
            assert len(distances) > 0


def test_get_interchain_distance_invalid_mode():
    protein_groups = [[0, 1, 2, 3, 4], [5, 6, 7, 8, 9], [10, 11, 12, 13, 14]]
    pdb_filename = os.path.join(soursop.get_data('test_data'), 'ntl9.pdb')
    traj_filename = os.path.join(soursop.get_data('test_data'), 'ntl9.xtc')  # ntl9 has 56 residues
    trajectory = sstrajectory.SSTrajectory(trajectory_filename=traj_filename,
                                           pdb_filename=pdb_filename,
                                           protein_grouping=protein_groups)

    mode = 'unknown'
    for protein_group in itertools.combinations(range(len(protein_groups)), r=2):
        protein_a, protein_b = protein_group

        # There appears to be an issue where R1 > R2 - the calculation crashes for sidechain-heavy.
        R1 = 0 #random.choice(range(len(protein_groups[protein_a])))
        R2 = 1 #random.choice(range(len(protein_groups[protein_b])))
        A1 = 'CA'
        A2 = 'CA'
        
        with pytest.raises(SSException):
            trajectory.get_interchain_distance(protein_a, protein_b, R1, R2, A1, A2, mode)


def test_get_interchain_distance_failing_atom():
    protein_groups = [[0, 1, 2, 3, 4], [5, 6, 7, 8, 9], [10, 11, 12, 13, 14]]
    pdb_filename = os.path.join(soursop.get_data('test_data'), 'ntl9.pdb')
    traj_filename = os.path.join(soursop.get_data('test_data'), 'ntl9.xtc')  # ntl9 has 56 residues
    trajectory = sstrajectory.SSTrajectory(trajectory_filename=traj_filename,
                                           pdb_filename=pdb_filename,
                                           protein_grouping=protein_groups)

    mode = 'atom'
    atoms1 = ['X', 'CA']
    atoms2 = atoms1[::-1]
    for protein_group in itertools.combinations(range(len(protein_groups)), r=2):
        protein_a, protein_b = protein_group

        # There appears to be an issue where R1 > R2 - the calculation crashes for sidechain-heavy.
        # See `test_get_interchain_distance`
        R1 = 0 #random.choice(range(len(protein_groups[protein_a])))
        R2 = 1 #random.choice(range(len(protein_groups[protein_b])))

        for A1, A2 in zip(atoms1, atoms2):
            with pytest.raises(SSException):
                trajectory.get_interchain_distance(protein_a, protein_b, R1, R2, A1, A2, mode)


def test_get_interchain_distance_failing_residues():
    protein_groups = [[0, 1, 2, 3, 4], [5, 6, 7, 8, 9], [10, 11, 12, 13, 14]]
    pdb_filename = os.path.join(soursop.get_data('test_data'), 'ntl9.pdb')
    traj_filename = os.path.join(soursop.get_data('test_data'), 'ntl9.xtc')  # ntl9 has 56 residues
    trajectory = sstrajectory.SSTrajectory(trajectory_filename=traj_filename,
                                           pdb_filename=pdb_filename,
                                           protein_grouping=protein_groups)

    residues1 = [0, 100]
    residues2 = residues1[::-1]
    A1 = 'CA'
    A2 = 'CA'
    modes = 'atom,ca,closest,closest-heavy,sidechain,sidechain-heavy'.split(',')
    
    for protein_group in itertools.combinations(range(len(protein_groups)), r=2):
        protein_a, protein_b = protein_group
        for mode in modes:
            for R1, R2 in zip(residues1, residues2):
                with pytest.raises(SSException):
                    trajectory.get_interchain_distance(protein_a, protein_b, R1, R2, A1, A2, mode)


def test_get_interchain_distance_invalid_protein_ids():
    protein_groups = [[0, 1, 2, 3, 4], [5, 6, 7, 8, 9], [10, 11, 12, 13, 14]]
    pdb_filename = os.path.join(soursop.get_data('test_data'), 'ntl9.pdb')
    traj_filename = os.path.join(soursop.get_data('test_data'), 'ntl9.xtc')  # ntl9 has 56 residues
    trajectory = sstrajectory.SSTrajectory(trajectory_filename=traj_filename,
                                           pdb_filename=pdb_filename,
                                           protein_grouping=protein_groups)

    R1 = 0
    R2 = 1
    A1 = 'CA'
    A2 = 'CA'
    modes = 'atom,ca,closest,closest-heavy,sidechain,sidechain-heavy'.split(',')
    
    protein_a = 100
    protein_b = 101
    for mode in modes:
        with pytest.raises(SSException):
            trajectory.get_interchain_distance(protein_a, protein_b, R1, R2, A1, A2, mode)


# -------------------------------------------------------------------------------------------------


def test_get_overall_hydrodynamic_radius(GS6_CO, NTL9_CO):
    gs6_radius = GS6_CO.get_overall_hydrodynamic_radius()
    ntl9_radius = NTL9_CO.get_overall_hydrodynamic_radius()
    assert len(gs6_radius) == GS6_CO.n_frames
    assert len(ntl9_radius) == NTL9_CO.n_frames


def test_get_overall_radius_of_gyration(GS6_CO, NTL9_CO):
    gs6_gyration = GS6_CO.get_overall_radius_of_gyration()
    ntl9_gyration = NTL9_CO.get_overall_radius_of_gyration()
    assert len(gs6_gyration) == GS6_CO.n_frames
    assert len(ntl9_gyration) == NTL9_CO.n_frames


def test_get_overall_asphericity(GS6_CO, NTL9_CO):
    gs6_asphericity = GS6_CO.get_overall_asphericity()
    ntl9_asphericity = NTL9_CO.get_overall_asphericity()
    assert len(gs6_asphericity) == GS6_CO.n_frames
    assert len(ntl9_asphericity) == NTL9_CO.n_frames


# -------------------------------------------------------------------------------------------------

def test_read_in_trajectory(GS6_CO):        

    # ensure 
    assert len(GS6_CO.traj) == 5
    assert len(GS6_CO.proteinTrajectoryList) == 1

"""
def test_protein_identification():
    pdb_filename = os.path.join(soursop.get_data('test_data'), 'gs6_invalid_r1.pdb')
    traj_filename = os.path.join(soursop.get_data('test_data'), 'gs6.xtc')
    trajectory = sstrajectory.SSTrajectory(trajectory_filename=traj_filename, pdb_filename=pdb_filename, debug=True)
    #print(trajectory)





"""


def test_get_overal_rh(GS6_CO):
    assert np.allclose(GS6_CO.get_overall_hydrodynamic_radius(), GS6_CO.proteinTrajectoryList[0].get_hydrodynamic_radius())

def test_get_overal_asphericity(GS6_CO):
    assert np.allclose(GS6_CO.get_overall_asphericity(), GS6_CO.proteinTrajectoryList[0].get_asphericity())

def test_get_overal_rg(GS6_CO):
    assert np.allclose(GS6_CO.get_overall_radius_of_gyration(), GS6_CO.proteinTrajectoryList[0].get_radius_of_gyration())

def test_get_overal_rg(GMX_2CHAINS):
    # compares rg of 2 chains vs. same values calculated by VMD

    # rg calculated in vmd
    vmd_rg = np.array([24.497053146362305, 24.938467025756836, 25.010461807250977, 27.36289405822754, 26.659685134887695, 24.56108283996582, 24.705211639404297, 27.144346237182617, 25.043672561645508, 27.67377471923828, 48.036869049072266, 23.930212020874023, 24.449384689331055, 24.491701126098633, 25.248193740844727, 27.898420333862305, 47.11823272705078, 47.213069915771484, 49.21323776245117, 27.429519653320312])

    assert np.allclose(GMX_2CHAINS.get_overall_radius_of_gyration(), vmd_rg)


def test_get_overal_asphericity(GMX_2CHAINS):
    # compares rg of 2 chains vs. same values calculated by VMD

    asph = np.array([0.36274249, 0.4079101, 0.5837932, 0.3413606, 0.30022731, 0.25667767, 0.24720234, 0.32076792, 0.26974903, 0.21184996, 0.62855842, 0.23021391, 0.16662495, 0.45061943, 0.54406795, 0.4482196, 0.53423135, 0.73352191, 0.75562579, 0.32473148])

    assert np.allclose(GMX_2CHAINS.get_overall_asphericity(), asph)

def test_get_overal_rh(GMX_2CHAINS):
    # compares rg of 2 chains vs. same values calculated by VMD

    asph = np.array([26.99809304, 27.25590015, 27.29754026, 28.59875282, 28.22148641, 27.03576021, 27.12021072, 28.48253016, 27.31671616, 28.76251356, 36.57647347, 26.66061797, 26.96999394, 26.99493729, 27.43425109, 28.87969746, 36.31506845, 36.3423472, 36.90213999, 28.63399989])

    assert np.allclose(GMX_2CHAINS.get_overall_hydrodynamic_radius(), asph)
