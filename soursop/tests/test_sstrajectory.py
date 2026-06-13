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
# import conftest

# -------------------------------------------------------------------------------------------------
# SSTrajectory `__init__` tests.
#
# Also implicitly tests:
#       `__readTrajectory`
#       `__get_proteins`
#       `__get_proteins_by_residue`


def test_trajectory_repr_string(GS6_CO):
    repr_string = "SSTrajectory (%s): %i proteins and %i frames" % (
        hex(id(GS6_CO)),
        GS6_CO.n_proteins,
        GS6_CO.n_frames,
    )
    assert repr_string == repr(GS6_CO)


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
    pdb_filename = os.path.join(soursop.get_data("test_data"), "gs6_AA.pdb")
    with pytest.raises(SSException):
        sstrajectory.SSTrajectory(pdb_filename=pdb_filename)


def test_read_in_trajectory_pdblead(GS6_CO):
    pdb_filename = os.path.join(soursop.get_data("test_data"), "gs6_AA.pdb")
    traj_filename = os.path.join(soursop.get_data("test_data"), "gs6_AA.xtc")
    traj = sstrajectory.SSTrajectory(
        pdb_filename=pdb_filename, trajectory_filename=traj_filename, pdblead=True
    )

    # Since we're using the PDB as an initial frame, the number of frames should have increased by 1.
    assert len(traj) == len(GS6_CO) + 1


def test_trajectory_initialization_debug():
    pdb_filename = os.path.join(soursop.get_data("test_data"), "gs6_AA.pdb")
    traj_filename = os.path.join(soursop.get_data("test_data"), "gs6_AA.xtc")
    sstrajectory.SSTrajectory(
        pdb_filename=pdb_filename, trajectory_filename=traj_filename, debug=True
    )


def test_read_in_no_topology_xtc():
    traj_filename = os.path.join(soursop.get_data("test_data"), "gs6_AA.xtc")
    with pytest.raises(SSException):
        sstrajectory.SSTrajectory(trajectory_filename=traj_filename)


def test_read_in_no_topology_dcd():
    """The `SSTrajectory.__init__` references that `.dcd` files are supported too."""
    traj_filename = os.path.join(soursop.get_data("test_data"), "gs6_AA.dcd")
    with pytest.raises(SSException):
        sstrajectory.SSTrajectory(trajectory_filename=traj_filename)


def test_read_in_compare_trajectories(GS6_CO):
    pdb_filename = os.path.join(soursop.get_data("test_data"), "gs6_AA.pdb")
    traj_filename = os.path.join(soursop.get_data("test_data"), "gs6_AA.xtc")
    trajectory = sstrajectory.SSTrajectory(
        trajectory_filename=traj_filename, pdb_filename=pdb_filename
    )

    assert np.allclose(trajectory.traj.xyz, GS6_CO.traj.xyz)


def test_read_in_protein_grouping_simple():
    protein_groups = [[0]]
    pdb_filename = os.path.join(soursop.get_data("test_data"), "gs6_AA.pdb")
    traj_filename = os.path.join(soursop.get_data("test_data"), "gs6_AA.xtc")
    trajectory = sstrajectory.SSTrajectory(
        trajectory_filename=traj_filename,
        pdb_filename=pdb_filename,
        protein_grouping=protein_groups,
    )

    # verify that we have loaded the number of proteins expected
    assert trajectory.n_proteins == len(protein_groups)


def test_read_in_protein_grouping_multiple():
    protein_groups = [[0, 1, 2], [3, 4, 5], [6, 7, 8]]
    pdb_filename = os.path.join(soursop.get_data("test_data"), "ntl9_AA.pdb")
    traj_filename = os.path.join(
        soursop.get_data("test_data"), "ntl9_AA.xtc"
    )  # ntl9 has 56 residues
    trajectory = sstrajectory.SSTrajectory(
        trajectory_filename=traj_filename,
        pdb_filename=pdb_filename,
        protein_grouping=protein_groups,
    )

    # verify that we have loaded the number of proteins expected
    assert trajectory.n_proteins == len(protein_groups)


def test_read_in_protein_grouping_invalid_residues():
    protein_groups = [[1000, 1001, 1002], [1003, 1004, 1005], [1006, 1007, 1008]]
    pdb_filename = os.path.join(soursop.get_data("test_data"), "ntl9_AA.pdb")
    traj_filename = os.path.join(
        soursop.get_data("test_data"), "ntl9_AA.xtc"
    )  # ntl9 has 56 residues
    trajectory = sstrajectory.SSTrajectory(
        trajectory_filename=traj_filename,
        pdb_filename=pdb_filename,
        protein_grouping=protein_groups,
    )

    # Since this is a failing but non-disruptive test (i.e. no Exceptions), we
    # check the number of proteins which should be 0.
    assert trajectory.n_proteins == 0


def test_read_in_protein_grouping_multiple_mixed_order():
    # TODO: This test passes - it shouldn't. Look into this further.
    protein_groups = [[1, 2, 0], [10, 3, 7], [11, 1, 3]]
    pdb_filename = os.path.join(soursop.get_data("test_data"), "ntl9_AA.pdb")
    traj_filename = os.path.join(
        soursop.get_data("test_data"), "ntl9_AA.xtc"
    )  # ntl9 has 56 residues
    trajectory = sstrajectory.SSTrajectory(
        trajectory_filename=traj_filename,
        pdb_filename=pdb_filename,
        protein_grouping=protein_groups,
    )

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
    distance_map_zero, stddev_map_zero = GS6_CO.get_interchain_distance_map(
        0, 0, mode="COM"
    )
    distance_map, stddev_map = GS6_CO.proteinTrajectoryList[0].get_distance_map(
        mode="COM"
    )

    assert np.allclose(np.triu(distance_map_zero), distance_map)
    assert np.allclose(np.triu(stddev_map_zero), stddev_map, atol=1e-6)


def test_get_intra_chain_distance_map_protein_groups():
    protein_groups = [[0, 1, 2, 3, 4], [5, 6, 7, 8, 9], [10, 11, 12, 13, 14]]
    pdb_filename = os.path.join(soursop.get_data("test_data"), "ntl9_AA.pdb")
    traj_filename = os.path.join(
        soursop.get_data("test_data"), "ntl9_AA.xtc"
    )  # ntl9 has 56 residues
    trajectory = sstrajectory.SSTrajectory(
        trajectory_filename=traj_filename,
        pdb_filename=pdb_filename,
        protein_grouping=protein_groups,
    )

    for protein_group in itertools.combinations(range(len(protein_groups)), r=2):
        protein_a, protein_b = protein_group
        distance_map, stddev_map = trajectory.get_interchain_distance_map(
            protein_a, protein_b
        )
        assert distance_map.shape == stddev_map.shape
        assert np.count_nonzero(distance_map) == len(
            distance_map.flatten()
        )  # since the residue indices are unique


def test_get_intra_chain_distance_map_protein_groups_with_residue_indices():
    # Note that the residues for the resID1 and resID2 indices are with respect to the residues of the protein chain
    # and not the full protein itself.
    protein_groups_residues = [
        [0, 1, 2, 3, 4, 5, 6, 7, 8],
        [10, 11, 12, 13, 14, 15, 16, 17, 18],
        [20, 21, 22, 23, 24, 25, 26, 27, 28],
        [30, 31, 32, 33, 34, 35, 36, 37, 38],
        [40, 41, 42, 43, 44, 45, 46, 47, 48],
    ]
    pdb_filename = os.path.join(soursop.get_data("test_data"), "ntl9_AA.pdb")
    traj_filename = os.path.join(
        soursop.get_data("test_data"), "ntl9_AA.xtc"
    )  # ntl9 has 56 residues
    trajectory = sstrajectory.SSTrajectory(
        trajectory_filename=traj_filename,
        pdb_filename=pdb_filename,
        protein_grouping=protein_groups_residues,
    )

    # select pairwise protein ids (no repetitions).
    for protein_group in itertools.combinations(
        range(len(protein_groups_residues)), r=2
    ):
        protein_a, protein_b = protein_group

        # now, let's randomly choose a set of residues to select for `protein_a` and `protein_b`
        num_residues_a = random.randint(1, len(protein_groups_residues[protein_a]))
        num_residues_b = random.randint(1, len(protein_groups_residues[protein_b]))

        protein_a_residues = protein_groups_residues[protein_a]
        a_residues = list(sorted(random.sample(protein_a_residues, num_residues_a)))
        a_residues = [
            (r - protein_a_residues[0]) for r in a_residues
        ]  # offset by the starting residue number

        protein_b_residues = protein_groups_residues[protein_b]
        b_residues = list(sorted(random.sample(protein_b_residues, num_residues_b)))
        b_residues = [
            (r - protein_b_residues[0]) for r in b_residues
        ]  # offset by the starting residue number

        distance_map, stddev_map = trajectory.get_interchain_distance_map(
            protein_a, protein_b
        )
        assert distance_map.shape == stddev_map.shape


def test_export_intra_chain_distance_map_protein_groups_with_residue_indices():
    # Note that the residues for the resID1 and resID2 indices are with respect to the residues of the protein chain
    # and not the full protein itself.
    protein_groups_residues = [
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
        [10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
        [20, 21, 22, 23, 24, 25, 26, 27, 28, 29],
        [30, 31, 32, 33, 34, 35, 36, 37, 38, 39],
        [40, 41, 42, 43, 44, 45, 46, 47, 48, 49],
    ]
    pdb_filename = os.path.join(soursop.get_data("test_data"), "ntl9_AA.pdb")
    traj_filename = os.path.join(
        soursop.get_data("test_data"), "ntl9_AA.xtc"
    )  # ntl9 has 56 residues
    trajectory = sstrajectory.SSTrajectory(
        trajectory_filename=traj_filename,
        pdb_filename=pdb_filename,
        protein_grouping=protein_groups_residues,
    )

    # select pairwise protein ids (no repetitions).
    for protein_index, protein_group in enumerate(
        itertools.combinations(range(len(protein_groups_residues)), r=2)
    ):
        protein_a, protein_b = protein_group

        # now, let's randomly choose a set of residues to select for `protein_a` and `protein_b`
        num_residues_a = random.randint(1, len(protein_groups_residues[protein_a]))
        num_residues_b = random.randint(1, len(protein_groups_residues[protein_b]))

        protein_a_residues = protein_groups_residues[protein_a]
        a_residues = list(sorted(random.sample(protein_a_residues, num_residues_a)))
        a_residues = [
            (r - protein_a_residues[0]) for r in a_residues
        ]  # offset by the starting residue number

        protein_b_residues = protein_groups_residues[protein_b]
        b_residues = list(sorted(random.sample(protein_b_residues, num_residues_b)))
        b_residues = [
            (r - protein_b_residues[0]) for r in b_residues
        ]  # offset by the starting residue number

        distance_map, stddev_map = trajectory.get_interchain_distance_map(
            protein_a, protein_b
        )
        assert distance_map.shape == stddev_map.shape

        temp_save_filepath = os.path.join(
            TMP_DIR,
            "gs6_map_{index}_{protein}".format(
                index=protein_index, protein=protein_index
            ),
        )
        temp_filename = "{filename}.csv".format(filename=temp_save_filepath)


def test_intrachain_inter_residue_atomic_distance():

    # A custom version of NTL9 is needed, since we want to have multiple protein chains
    protein_groups_residues = [
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
        [10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
        [20, 21, 22, 23, 24, 25, 26, 27, 28, 29],
        [30, 31, 32, 33, 34, 35, 36, 37, 38, 39],
        [40, 41, 42, 43, 44, 45, 46, 47, 48, 49],
    ]
    pdb_filename = os.path.join(soursop.get_data("test_data"), "ntl9_AA.pdb")
    traj_filename = os.path.join(
        soursop.get_data("test_data"), "ntl9_AA.xtc"
    )  # ntl9 has 56 residues
    trajectory = sstrajectory.SSTrajectory(
        trajectory_filename=traj_filename,
        pdb_filename=pdb_filename,
        protein_grouping=protein_groups_residues,
    )

    for protein_group in itertools.combinations(
        range(len(protein_groups_residues)), r=2
    ):
        protein_a, protein_b = protein_group
        a_residues = protein_groups_residues[protein_a]
        b_residues = protein_groups_residues[protein_b]

        residue_a = random.randint(0, len(a_residues) - 1)
        residue_b = random.randint(0, len(b_residues) - 1)

        distances = trajectory.get_interchain_distance(
            protein_a, protein_b, residue_a, residue_b
        )
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
    pdb_filename = os.path.join(soursop.get_data('test_data'), 'ntl9_AA.pdb')
    traj_filename = os.path.join(soursop.get_data('test_data'), 'ntl9_AA.xtc')  # ntl9 has 56 residues
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
    pdb_filename = os.path.join(soursop.get_data("test_data"), "ntl9_AA.pdb")
    traj_filename = os.path.join(
        soursop.get_data("test_data"), "ntl9_AA.xtc"
    )  # ntl9 has 56 residues
    trajectory = sstrajectory.SSTrajectory(
        trajectory_filename=traj_filename,
        pdb_filename=pdb_filename,
        protein_grouping=protein_groups,
    )

    for protein_group in itertools.combinations(range(len(protein_groups)), r=2):
        protein_a, protein_b = protein_group

        # There appears to be an issue where R1 > R2 - the calculation crashes for sidechain-heavy.
        R1 = 0  # random.choice(range(len(protein_groups[protein_a])))
        R2 = 1  # random.choice(range(len(protein_groups[protein_b])))
        A1 = "CA"
        A2 = "CA"
        modes = "atom,ca,closest,closest-heavy,sidechain,sidechain-heavy".split(",")
        for mode in modes:
            distances = trajectory.get_interchain_distance(
                protein_a, protein_b, R1, R2, A1, A2, mode
            )
            assert len(distances) > 0


def test_get_interchain_distance_invalid_mode():
    protein_groups = [[0, 1, 2, 3, 4], [5, 6, 7, 8, 9], [10, 11, 12, 13, 14]]
    pdb_filename = os.path.join(soursop.get_data("test_data"), "ntl9_AA.pdb")
    traj_filename = os.path.join(
        soursop.get_data("test_data"), "ntl9_AA.xtc"
    )  # ntl9 has 56 residues
    trajectory = sstrajectory.SSTrajectory(
        trajectory_filename=traj_filename,
        pdb_filename=pdb_filename,
        protein_grouping=protein_groups,
    )

    mode = "unknown"
    for protein_group in itertools.combinations(range(len(protein_groups)), r=2):
        protein_a, protein_b = protein_group

        # There appears to be an issue where R1 > R2 - the calculation crashes for sidechain-heavy.
        R1 = 0  # random.choice(range(len(protein_groups[protein_a])))
        R2 = 1  # random.choice(range(len(protein_groups[protein_b])))
        A1 = "CA"
        A2 = "CA"

        with pytest.raises(SSException):
            trajectory.get_interchain_distance(
                protein_a, protein_b, R1, R2, A1, A2, mode
            )


def test_get_interchain_distance_failing_atom():
    protein_groups = [[0, 1, 2, 3, 4], [5, 6, 7, 8, 9], [10, 11, 12, 13, 14]]
    pdb_filename = os.path.join(soursop.get_data("test_data"), "ntl9_AA.pdb")
    traj_filename = os.path.join(
        soursop.get_data("test_data"), "ntl9_AA.xtc"
    )  # ntl9 has 56 residues
    trajectory = sstrajectory.SSTrajectory(
        trajectory_filename=traj_filename,
        pdb_filename=pdb_filename,
        protein_grouping=protein_groups,
    )

    mode = "atom"
    atoms1 = ["X", "CA"]
    atoms2 = atoms1[::-1]
    for protein_group in itertools.combinations(range(len(protein_groups)), r=2):
        protein_a, protein_b = protein_group

        # There appears to be an issue where R1 > R2 - the calculation crashes for sidechain-heavy.
        # See `test_get_interchain_distance`
        R1 = 0  # random.choice(range(len(protein_groups[protein_a])))
        R2 = 1  # random.choice(range(len(protein_groups[protein_b])))

        for A1, A2 in zip(atoms1, atoms2):
            with pytest.raises(SSException):
                trajectory.get_interchain_distance(
                    protein_a, protein_b, R1, R2, A1, A2, mode
                )


def test_get_interchain_distance_failing_residues():
    protein_groups = [[0, 1, 2, 3, 4], [5, 6, 7, 8, 9], [10, 11, 12, 13, 14]]
    pdb_filename = os.path.join(soursop.get_data("test_data"), "ntl9_AA.pdb")
    traj_filename = os.path.join(
        soursop.get_data("test_data"), "ntl9_AA.xtc"
    )  # ntl9 has 56 residues
    trajectory = sstrajectory.SSTrajectory(
        trajectory_filename=traj_filename,
        pdb_filename=pdb_filename,
        protein_grouping=protein_groups,
    )

    residues1 = [0, 100]
    residues2 = residues1[::-1]
    A1 = "CA"
    A2 = "CA"
    modes = "atom,ca,closest,closest-heavy,sidechain,sidechain-heavy".split(",")

    for protein_group in itertools.combinations(range(len(protein_groups)), r=2):
        protein_a, protein_b = protein_group
        for mode in modes:
            for R1, R2 in zip(residues1, residues2):
                with pytest.raises(SSException):
                    trajectory.get_interchain_distance(
                        protein_a, protein_b, R1, R2, A1, A2, mode
                    )


def test_get_interchain_distance_invalid_protein_ids():
    protein_groups = [[0, 1, 2, 3, 4], [5, 6, 7, 8, 9], [10, 11, 12, 13, 14]]
    pdb_filename = os.path.join(soursop.get_data("test_data"), "ntl9_AA.pdb")
    traj_filename = os.path.join(
        soursop.get_data("test_data"), "ntl9_AA.xtc"
    )  # ntl9 has 56 residues
    trajectory = sstrajectory.SSTrajectory(
        trajectory_filename=traj_filename,
        pdb_filename=pdb_filename,
        protein_grouping=protein_groups,
    )

    R1 = 0
    R2 = 1
    A1 = "CA"
    A2 = "CA"
    modes = "atom,ca,closest,closest-heavy,sidechain,sidechain-heavy".split(",")

    protein_a = 100
    protein_b = 101
    for mode in modes:
        with pytest.raises(SSException):
            trajectory.get_interchain_distance(
                protein_a, protein_b, R1, R2, A1, A2, mode
            )


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
    traj_filename = os.path.join(soursop.get_data('test_data'), 'gs6_AA.xtc')
    trajectory = sstrajectory.SSTrajectory(trajectory_filename=traj_filename, pdb_filename=pdb_filename, debug=True)
    #print(trajectory)





"""


def test_get_overal_rh(GS6_CO):
    assert np.allclose(
        GS6_CO.get_overall_hydrodynamic_radius(),
        GS6_CO.proteinTrajectoryList[0].get_hydrodynamic_radius(),
    )


def test_get_overal_asphericity(GS6_CO):
    assert np.allclose(
        GS6_CO.get_overall_asphericity(),
        GS6_CO.proteinTrajectoryList[0].get_asphericity(),
    )


def test_get_overal_rg(GS6_CO):
    assert np.allclose(
        GS6_CO.get_overall_radius_of_gyration(),
        GS6_CO.proteinTrajectoryList[0].get_radius_of_gyration(),
    )


def test_get_overal_rg(GMX_2CHAINS):
    # compares rg of 2 chains vs. same values calculated by VMD

    # rg calculated in vmd
    vmd_rg = np.array(
        [
            24.497053146362305,
            24.938467025756836,
            25.010461807250977,
            27.36289405822754,
            26.659685134887695,
            24.56108283996582,
            24.705211639404297,
            27.144346237182617,
            25.043672561645508,
            27.67377471923828,
            48.036869049072266,
            23.930212020874023,
            24.449384689331055,
            24.491701126098633,
            25.248193740844727,
            27.898420333862305,
            47.11823272705078,
            47.213069915771484,
            49.21323776245117,
            27.429519653320312,
        ]
    )

    assert np.allclose(GMX_2CHAINS.get_overall_radius_of_gyration(), vmd_rg)


def test_get_overal_asphericity(GMX_2CHAINS):
    # compares rg of 2 chains vs. same values calculated by VMD

    asph = np.array(
        [
            0.36274249,
            0.4079101,
            0.5837932,
            0.3413606,
            0.30022731,
            0.25667767,
            0.24720234,
            0.32076792,
            0.26974903,
            0.21184996,
            0.62855842,
            0.23021391,
            0.16662495,
            0.45061943,
            0.54406795,
            0.4482196,
            0.53423135,
            0.73352191,
            0.75562579,
            0.32473148,
        ]
    )

    assert np.allclose(GMX_2CHAINS.get_overall_asphericity(), asph)


def test_get_overal_rh(GMX_2CHAINS):
    # compares rg of 2 chains vs. same values calculated by VMD

    asph = np.array(
        [
            26.99809304,
            27.25590015,
            27.29754026,
            28.59875282,
            28.22148641,
            27.03576021,
            27.12021072,
            28.48253016,
            27.31671616,
            28.76251356,
            36.57647347,
            26.66061797,
            26.96999394,
            26.99493729,
            27.43425109,
            28.87969746,
            36.31506845,
            36.3423472,
            36.90213999,
            28.63399989,
        ]
    )

    assert np.allclose(GMX_2CHAINS.get_overall_hydrodynamic_radius(), asph)


# -------------------------------------------------------------------------------------------------
# Parameterized smoke tests across every available trajectory fixture.
#
# These exercise the SSTrajectory loader and properties without asserting any
# specific numerical values — they ensure each file format / resolution loads
# cleanly and the basic properties are internally consistent.
# -------------------------------------------------------------------------------------------------

# Every loaded trajectory fixture in conftest.py. Used for property smoke
# tests that should hold regardless of resolution or chain count.
ALL_TRAJ_FIXTURES = [
    "GS6_CO",
    "CTL9_CO",
    "NTL9_CO",
    "SIGA_CG_CO",
    "SYNTH_1_CG_CO",
    "ALL_RESIDUES_CO",
    "GMX_1CHAIN_CO",
    "GMX_2CHAINS",
]

# Single-chain trajectories only — for tests that compare the SSTrajectory
# "overall" aggregator against the underlying SSProtein's own value.
SINGLE_CHAIN_TRAJ_FIXTURES = [
    "GS6_CO",
    "CTL9_CO",
    "NTL9_CO",
    "SIGA_CG_CO",
    "SYNTH_1_CG_CO",
    "ALL_RESIDUES_CO",
    "GMX_1CHAIN_CO",
]

# Modes accepted by get_interchain_distance / get_interchain_contact_map.
INTERCHAIN_DISTANCE_MODES = [
    "atom",
    "ca",
    "closest",
    "closest-heavy",
    "sidechain",
    "sidechain-heavy",
]


@pytest.mark.parametrize("traj_fixture", ALL_TRAJ_FIXTURES)
def test_n_frames_matches_mdtraj(traj_fixture, request):
    traj = request.getfixturevalue(traj_fixture)
    assert traj.n_frames == len(traj.traj)
    assert traj.n_frames > 0


@pytest.mark.parametrize("traj_fixture", ALL_TRAJ_FIXTURES)
def test_len_matches_n_frames(traj_fixture, request):
    traj = request.getfixturevalue(traj_fixture)
    assert len(traj) == traj.n_frames


@pytest.mark.parametrize("traj_fixture", ALL_TRAJ_FIXTURES)
def test_length_property_returns_tuple(traj_fixture, request):
    traj = request.getfixturevalue(traj_fixture)
    n_proteins, n_frames = traj.length
    assert n_proteins == traj.n_proteins
    assert n_frames == traj.n_frames


@pytest.mark.parametrize("traj_fixture", ALL_TRAJ_FIXTURES)
def test_n_proteins_matches_protein_list(traj_fixture, request):
    traj = request.getfixturevalue(traj_fixture)
    assert traj.n_proteins == len(traj.proteinTrajectoryList)
    assert traj.n_proteins >= 1


@pytest.mark.parametrize("traj_fixture", ALL_TRAJ_FIXTURES)
def test_repr_format(traj_fixture, request):
    traj = request.getfixturevalue(traj_fixture)
    expected = "SSTrajectory (%s): %i proteins and %i frames" % (
        hex(id(traj)),
        traj.n_proteins,
        traj.n_frames,
    )
    assert repr(traj) == expected


@pytest.mark.parametrize("traj_fixture", ALL_TRAJ_FIXTURES)
def test_unitcell_or_typeerror(traj_fixture, request):
    """unitcell returns an (3,)-shape Å array when a periodic box exists, or
    raises TypeError when the trajectory has no box (some CG runs)."""
    traj = request.getfixturevalue(traj_fixture)
    try:
        cell = traj.unitcell
    except TypeError:
        # No periodic box recorded — acceptable for CG trajectories
        return
    cell = np.asarray(cell)
    assert cell.shape == (3,)
    assert np.all(np.isfinite(cell))
    assert np.all(cell > 0)


# -------------------------------------------------------------------------------------------------
# Overall-observable consistency: for single-chain trajectories, the
# SSTrajectory.get_overall_X() must equal the underlying SSProtein's
# corresponding per-chain method (since there is only one chain to
# "aggregate" over).
# -------------------------------------------------------------------------------------------------


@pytest.mark.parametrize("traj_fixture", SINGLE_CHAIN_TRAJ_FIXTURES)
def test_overall_radius_of_gyration_matches_protein(traj_fixture, request):
    traj = request.getfixturevalue(traj_fixture)
    overall = traj.get_overall_radius_of_gyration()
    per_chain = traj.proteinTrajectoryList[0].get_radius_of_gyration()
    assert overall.shape == per_chain.shape == (traj.n_frames,)
    assert np.allclose(overall, per_chain)


@pytest.mark.parametrize("traj_fixture", SINGLE_CHAIN_TRAJ_FIXTURES)
def test_overall_hydrodynamic_radius_matches_protein(traj_fixture, request):
    traj = request.getfixturevalue(traj_fixture)
    overall = traj.get_overall_hydrodynamic_radius()
    per_chain = traj.proteinTrajectoryList[0].get_hydrodynamic_radius()
    assert overall.shape == per_chain.shape == (traj.n_frames,)
    assert np.allclose(overall, per_chain)


@pytest.mark.parametrize("traj_fixture", SINGLE_CHAIN_TRAJ_FIXTURES)
def test_overall_asphericity_matches_protein(traj_fixture, request):
    traj = request.getfixturevalue(traj_fixture)
    overall = traj.get_overall_asphericity()
    per_chain = traj.proteinTrajectoryList[0].get_asphericity(verbose=False)
    assert overall.shape == per_chain.shape == (traj.n_frames,)
    assert np.allclose(overall, per_chain)


# -------------------------------------------------------------------------------------------------
# get_interchain_contact_map coverage — completely untested before.
# -------------------------------------------------------------------------------------------------


# Subsample frames so contact-map tests stay fast even on long trajectories
# (CTL9 has 1000 frames; full evaluation would do 94 * 94 * 1000 distance
# computations per test). Used by every contact-map test in this file.
def _cmap_stride(traj):
    """Pick a stride so at most ~100 frames are evaluated per residue pair."""
    return max(1, traj.n_frames // 100)


@pytest.mark.parametrize("traj_fixture", ["GS6_CO", "NTL9_CO", "CTL9_CO"])
def test_interchain_contact_map_self_shape(traj_fixture, request):
    """Contact map of a protein with itself must be (n_res, n_res) and finite.

    GS6 and CTL9 have ACE/NME caps; the fixed implementation skips cap
    residues silently (leaving 0 in those rows/cols) so the full shape is
    preserved.
    """
    traj = request.getfixturevalue(traj_fixture)
    n_res = traj.proteinTrajectoryList[0].n_residues
    cmap = traj.get_interchain_contact_map(0, 0, stride=_cmap_stride(traj))
    assert cmap.shape == (n_res, n_res)
    assert np.all(np.isfinite(cmap))
    # Each entry is a fraction of frames in contact — must be in [0, 1].
    assert np.all(cmap >= 0.0) and np.all(cmap <= 1.0)


@pytest.mark.parametrize("traj_fixture", ["GS6_CO", "CTL9_CO"])
def test_interchain_contact_map_self_cap_rows_zero(traj_fixture, request):
    """For trajectories with ACE/NME caps, the cap rows and columns of the
    contact map must be exactly zero (the fixed implementation can't compute
    a CA-based contact for a residue with no CA)."""
    traj = request.getfixturevalue(traj_fixture)
    protein = traj.proteinTrajectoryList[0]
    cmap = traj.get_interchain_contact_map(0, 0, stride=_cmap_stride(traj))

    cap_resids = set(range(protein.n_residues)) - set(protein.resid_with_CA)
    assert len(cap_resids) > 0, "fixture is supposed to have caps"
    for cap in cap_resids:
        assert np.all(cmap[cap, :] == 0.0), f"cap row {cap} should be zero"
        assert np.all(cmap[:, cap] == 0.0), f"cap col {cap} should be zero"


def test_interchain_contact_map_gmx_2chains_shape(GMX_2CHAINS):
    """Contact map between the two real chains has shape (n_res_0, n_res_1)."""
    n_res_0 = GMX_2CHAINS.proteinTrajectoryList[0].n_residues
    n_res_1 = GMX_2CHAINS.proteinTrajectoryList[1].n_residues
    cmap = GMX_2CHAINS.get_interchain_contact_map(
        0, 1, stride=_cmap_stride(GMX_2CHAINS)
    )
    assert cmap.shape == (n_res_0, n_res_1)
    assert np.all(np.isfinite(cmap))
    assert np.all(cmap >= 0.0) and np.all(cmap <= 1.0)


@pytest.mark.parametrize("mode", INTERCHAIN_DISTANCE_MODES)
def test_interchain_contact_map_gmx_2chains_all_modes(GMX_2CHAINS, mode):
    """Every documented mode must produce a finite (n_res_0, n_res_1) contact map.

    sidechain / sidechain-heavy used to fail on glycine-containing chains
    because mdtraj's sidechain selector returns an empty set for glycine;
    get_interchain_contact_map now catches that per-residue failure and
    leaves the affected entries at 0.
    """
    n_res_0 = GMX_2CHAINS.proteinTrajectoryList[0].n_residues
    n_res_1 = GMX_2CHAINS.proteinTrajectoryList[1].n_residues
    cmap = GMX_2CHAINS.get_interchain_contact_map(
        0,
        1,
        mode=mode,
        stride=_cmap_stride(GMX_2CHAINS),
    )
    assert cmap.shape == (n_res_0, n_res_1)
    assert np.all(np.isfinite(cmap))
    assert np.all(cmap >= 0.0) and np.all(cmap <= 1.0)


def test_interchain_contact_map_sidechain_zeros_glycine_rows(GMX_2CHAINS):
    """In sidechain mode, rows/cols for glycine residues must be exactly 0
    (glycine has no sidechain atom, so the contact is undefined)."""
    p1 = GMX_2CHAINS.proteinTrajectoryList[0]
    p2 = GMX_2CHAINS.proteinTrajectoryList[1]
    p1_gly = [
        i
        for i, name in enumerate(p1.get_amino_acid_sequence(numbered=False))
        if name.split("-")[0] == "GLY"
    ]
    p2_gly = [
        i
        for i, name in enumerate(p2.get_amino_acid_sequence(numbered=False))
        if name.split("-")[0] == "GLY"
    ]
    assert len(p1_gly) > 0 or len(p2_gly) > 0, "fixture is supposed to contain glycine"

    cmap = GMX_2CHAINS.get_interchain_contact_map(
        0,
        1,
        mode="sidechain",
        stride=_cmap_stride(GMX_2CHAINS),
    )
    for r in p1_gly:
        assert np.all(cmap[r, :] == 0.0), f"glycine row {r} (chain 0) should be zero"
    for r in p2_gly:
        assert np.all(cmap[:, r] == 0.0), f"glycine col {r} (chain 1) should be zero"


def test_interchain_contact_map_stride_default_equals_explicit_1(GMX_2CHAINS):
    """stride=1 must match the default (no stride) — backwards-compatible."""
    cmap_default = GMX_2CHAINS.get_interchain_contact_map(0, 1, mode="ca")
    cmap_stride1 = GMX_2CHAINS.get_interchain_contact_map(0, 1, mode="ca", stride=1)
    assert np.array_equal(cmap_default, cmap_stride1)


def test_interchain_contact_map_stride_subsamples(GMX_2CHAINS):
    """Larger strides produce the same matrix shape and bounded values."""
    cmap = GMX_2CHAINS.get_interchain_contact_map(
        0,
        1,
        mode="ca",
        stride=GMX_2CHAINS.n_frames,
    )
    n_res_0 = GMX_2CHAINS.proteinTrajectoryList[0].n_residues
    n_res_1 = GMX_2CHAINS.proteinTrajectoryList[1].n_residues
    assert cmap.shape == (n_res_0, n_res_1)
    assert np.all(np.isfinite(cmap))
    # With stride==n_frames only frame 0 is used, so every entry is 0 or 1.
    assert set(np.unique(cmap.astype(int))).issubset({0, 1})


@pytest.mark.parametrize("bad_stride", [0, -1, 1.5, "two"])
def test_interchain_contact_map_invalid_stride(GMX_2CHAINS, bad_stride):
    """stride must be a positive integer; other values raise SSException."""
    with pytest.raises(SSException):
        GMX_2CHAINS.get_interchain_contact_map(0, 1, stride=bad_stride)


def test_interchain_contact_map_invalid_mode(GMX_2CHAINS):
    with pytest.raises(SSException):
        GMX_2CHAINS.get_interchain_contact_map(0, 1, mode="not-a-real-mode")


def test_interchain_contact_map_invalid_protein_ids(GMX_2CHAINS):
    # Out-of-range proteinID1 or proteinID2 must raise SSException up front
    # (the implementation now validates IDs before any dereference).
    with pytest.raises(SSException):
        GMX_2CHAINS.get_interchain_contact_map(99, 1)
    with pytest.raises(SSException):
        GMX_2CHAINS.get_interchain_contact_map(0, 99)
    with pytest.raises(SSException):
        GMX_2CHAINS.get_interchain_contact_map(-1, 0)


def test_interchain_contact_map_invalid_atom_name(GMX_2CHAINS):
    with pytest.raises(SSException):
        GMX_2CHAINS.get_interchain_contact_map(0, 1, mode="atom", A1="XXX", A2="CA")


# -------------------------------------------------------------------------------------------------
# Inter-chain methods on GMX_2CHAINS (real two-chain trajectory).
# -------------------------------------------------------------------------------------------------


@pytest.mark.parametrize("mode", ["CA", "COM"])
def test_gmx_2chains_interchain_distance_map_shape(GMX_2CHAINS, mode):
    n_res_0 = GMX_2CHAINS.proteinTrajectoryList[0].n_residues
    n_res_1 = GMX_2CHAINS.proteinTrajectoryList[1].n_residues
    dmap, sdmap = GMX_2CHAINS.get_interchain_distance_map(0, 1, mode=mode)
    assert dmap.shape == sdmap.shape == (n_res_0, n_res_1)
    assert np.all(np.isfinite(dmap))
    assert np.all(np.isfinite(sdmap))
    # All distances must be positive; std must be non-negative.
    assert np.all(dmap > 0.0)
    assert np.all(sdmap >= 0.0)


def test_gmx_2chains_interchain_distance_map_invalid_mode(GMX_2CHAINS):
    with pytest.raises(SSException):
        GMX_2CHAINS.get_interchain_distance_map(0, 1, mode="not-a-real-mode")


@pytest.mark.parametrize("mode", INTERCHAIN_DISTANCE_MODES)
def test_gmx_2chains_interchain_distance_all_modes(GMX_2CHAINS, mode):
    """get_interchain_distance returns one positive finite value per frame."""
    distances = GMX_2CHAINS.get_interchain_distance(
        0, 1, R1=0, R2=0, A1="CA", A2="CA", mode=mode
    )
    assert len(distances) == GMX_2CHAINS.n_frames
    distances = np.asarray(distances)
    assert np.all(np.isfinite(distances))
    assert np.all(distances > 0.0)


# -------------------------------------------------------------------------------------------------
# Error-path coverage: parameterized error checks for get_interchain_distance.
# (The existing pre-parameterized tests cover get_interchain_distance in
# aggregate; these split per-mode for granular reporting.)
# -------------------------------------------------------------------------------------------------


@pytest.mark.parametrize("mode", INTERCHAIN_DISTANCE_MODES)
def test_interchain_distance_invalid_protein_id_per_mode(GMX_2CHAINS, mode):
    with pytest.raises(SSException):
        GMX_2CHAINS.get_interchain_distance(
            99, 1, R1=0, R2=0, A1="CA", A2="CA", mode=mode
        )


@pytest.mark.parametrize("mode", INTERCHAIN_DISTANCE_MODES)
def test_interchain_distance_out_of_bounds_residue_per_mode(GMX_2CHAINS, mode):
    with pytest.raises(SSException):
        GMX_2CHAINS.get_interchain_distance(
            0, 1, R1=99999, R2=0, A1="CA", A2="CA", mode=mode
        )
