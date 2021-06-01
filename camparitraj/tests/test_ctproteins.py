"""
Unit and regression test for the camparitraj package.
"""

# Import package, test suite, and other packages as needed
import numpy as np
import camparitraj
import pytest
import sys
import random
import itertools
from camparitraj import cttrajectory, ctprotein
from camparitraj.ctexceptions import CTException
from camparitraj.configs import DEBUGGING


def test_contacts_special(NTL9_CP):

    with pytest.raises(CTException) as error:
        a = NTL9_CP.get_contact_map(distance_thresh=2, stride=1, mode='sidechain-heavy')

def test_code_coverage(NTL9_CP):

    NTL9_CP.print_residues()

    # property
    a = NTL9_CP.residue_index_list
    a = NTL9_CP.n_residues
    a = NTL9_CP.n_frames



    a = NTL9_CP.get_amino_acid_sequence()
    a = NTL9_CP.get_CA_index(10)
    a = NTL9_CP.get_multiple_CA_index()
    a = NTL9_CP.get_multiple_CA_index([3,4,5])


    a = NTL9_CP.calculate_all_CA_distances(10)
    a = NTL9_CP.calculate_all_CA_distances(10, stride=2)
    a = NTL9_CP.calculate_all_CA_distances(10, stride=2)
    a = NTL9_CP.calculate_all_CA_distances(10, stride=2, mode='COM')
    a = NTL9_CP.calculate_all_CA_distances(10, stride=2, mode='COM', only_C_terminal_residues=False)

    a = NTL9_CP.get_distance_map()
    a = NTL9_CP.get_distance_map(verbose=True)
    a = NTL9_CP.get_distance_map(verbose=True, mode='COM')
    a = NTL9_CP.get_distance_map(verbose=True, mode='COM', RMS=True)

    a = NTL9_CP.get_polymer_scaled_distance_map()
    a = NTL9_CP.get_polymer_scaled_distance_map(nu=0.54,A0=6)
    a = NTL9_CP.get_polymer_scaled_distance_map(nu=0.54,A0=6,min_separation=20)
    a = NTL9_CP.get_polymer_scaled_distance_map(nu=0.54,A0=6,min_separation=20, mode='signed-fractional-change')
    a = NTL9_CP.get_polymer_scaled_distance_map(nu=0.54,A0=6,min_separation=20, mode='signed-absolute-change')
    a = NTL9_CP.get_polymer_scaled_distance_map(nu=0.54,A0=6,min_separation=20, mode='scaled')

    a = NTL9_CP.get_local_heterogeneity(stride=1)
    a = NTL9_CP.get_local_heterogeneity(fragment_size=5, stride=1)
    a = NTL9_CP.get_local_heterogeneity(fragment_size=20, stride=2)
    a = NTL9_CP.get_local_heterogeneity(fragment_size=20, stride=1)
    a = NTL9_CP.get_D_vector(stride=1)

    a = NTL9_CP.get_RMSD(1)
    a = NTL9_CP.get_RMSD(1,frame2=3)
    a = NTL9_CP.get_RMSD(1,frame2=3, stride=2)
    a = NTL9_CP.get_RMSD(1,frame2=3, stride=2, region=[10,20])
    a = NTL9_CP.get_RMSD(1,frame2=3, stride=2, backbone=False)


    a = NTL9_CP.get_Q()
    a = NTL9_CP.get_Q(stride=2)
    a = NTL9_CP.get_Q(stride=2, protein_average=False)
    a = NTL9_CP.get_Q(stride=2, protein_average=False, region=[10,20])


    a = NTL9_CP.get_contact_map()
    a = NTL9_CP.get_contact_map(distance_thresh=2)
    a = NTL9_CP.get_contact_map(distance_thresh=2, stride=2)
    a = NTL9_CP.get_contact_map(distance_thresh=2, stride=1, mode='ca')
    a = NTL9_CP.get_contact_map(distance_thresh=2, stride=1, mode='closest-heavy')
    a = NTL9_CP.get_contact_map(distance_thresh=2, stride=1, mode='closest')
    a = NTL9_CP.get_contact_map(distance_thresh=2, stride=1, mode='sidechain')

    # This test fails. Will revisit.
    # a = NTL9_CP.get_contact_map(distance_thresh=2, stride=1, mode='sidechain-heavy')

    a = NTL9_CP.get_clusters(stride=1)
    a = NTL9_CP.get_clusters(stride=1, n_clusters=3)
    a = NTL9_CP.get_clusters(stride=1, n_clusters=3, backbone=False)


    a = NTL9_CP.get_inter_residue_COM_distance(10,30)
    a = NTL9_CP.get_inter_residue_COM_distance(10,30, stride=2)


    a = NTL9_CP.get_inter_residue_COM_vector(10,30)
    a = NTL9_CP.get_inter_residue_COM_vector(10,30, stride=2)


    a = NTL9_CP.get_inter_residue_atomic_distance(2,10)
    a = NTL9_CP.get_inter_residue_atomic_distance(2,10, A1='CB')
    a = NTL9_CP.get_angle_decay()


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


def test_BBSEG(NTL9_CP):
    SS = NTL9_CP.get_secondary_structure_BBSEG()
    print(SS[0])



def test_phi_psi_bbseg(NTL9_CP):
    """
    Test that ensures our representation of the BBSEG2 matrix
    reproduces the version from CAMPARI

    """

    # read in BBSEG file
    bbseg_file = camparitraj.get_data('bbseg2.dat')
    with open(bbseg_file,'r') as fh:
        content = fh.readlines()

    idx=0
    for psi in range(179,-180, -10):

        local = ''
        for phi in range(-179,180, 10):
            local = local + '%i ' %(NTL9_CP._CTProtein__phi_psi_bbseg([phi],[psi])[0])

        assert(content[idx].strip() == local.strip())
        idx=idx+1



def test_get_distance_map(GS6_CO):
    for protein_index in range(len(GS6_CO.proteinTrajectoryList)):
        distance_map, stddev_map = GS6_CO.proteinTrajectoryList[0].get_distance_map()

        # verify that the shapes match
        assert distance_map.shape == stddev_map.shape

        # verify that we obtain an upper triangular matrix
        assert np.allclose(distance_map, np.triu(distance_map)) is True



def test_get_hydrodynamic_radius(GS6_CO):

    CP = GS6_CO.proteinTrajectoryList[0]
    rh = CP.get_hydrodynamic_radius()
    assert (11.840569781179006 - rh[0]) < 0.001


# ====
def test_init_from_trajectory(GS6_CO, NTL9_CO, GS6_CP, NTL9_CP):
    trajs = [GS6_CO, NTL9_CO]
    proteins = [GS6_CP, NTL9_CP]
    for ct_traj, ct_protein in zip(trajs, proteins):
        protein_from_md = ctprotein.CTProtein(ct_traj)
        assert protein_from_md.n_frames > 0
        assert protein_from_md.n_residues > 0
        assert protein_from_md.length() == ct_protein.length()

        protein_from_ct = ctprotein.CTProtein(ct_traj.traj)
        assert protein_from_ct.n_frames > 0
        assert protein_from_ct.n_residues > 0
        assert protein_from_ct.length() == ct_protein.length()


def test_init_from_trajectory_invalid():
    args = [1, 1.0, 'test', None]
    for arg in args:
        with pytest.raises(RuntimeError):
            protein = ctprotein.CTProtein(arg)


def test_init_from_trajectory_debug(GS6_CO, NTL9_CO, GS6_CP, NTL9_CP):
    trajs = [GS6_CO, NTL9_CO]
    proteins = [GS6_CP, NTL9_CP]
    for ct_traj, ct_protein in zip(trajs, proteins):
        protein_from_md = ctprotein.CTProtein(ct_traj.traj, debug=True)
        assert protein_from_md.n_frames > 0
        assert protein_from_md.n_residues > 0
        assert protein_from_md.length() == ct_protein.length()

        protein_from_ct = ctprotein.CTProtein(ct_traj, debug=True)
        assert protein_from_ct.n_frames > 0
        assert protein_from_ct.n_residues > 0
        assert protein_from_ct.length() == ct_protein.length()


def test_properties(GS6_CO, NTL9_CO, GS6_CP, NTL9_CP):
    properties = 'resid_with_CA,ncap,ccap,n_frames,n_residues,residue_index_list'.split(',')
    trajs = [GS6_CO, NTL9_CO]
    proteins = [GS6_CP, NTL9_CP]
    for ct_traj, ct_protein in zip(trajs, proteins):
        protein_from_ct = ctprotein.CTProtein(ct_traj)

        for prop in properties:
            assert getattr(protein_from_ct, prop) == getattr(ct_protein, prop)


def test_repr(GS6_CP, NTL9_CP):
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        repr_string = repr(protein)
        rebuilt_repr = "CTProtein (%s): %i res and %i frames" % (hex(id(protein)), protein.n_residues, protein.n_frames)
        assert rebuilt_repr == repr_string


def test_len(GS6_CP, NTL9_CP):
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        length = len(protein)
        attributes = protein.length()
        assert length == protein.n_frames
        assert attributes == (protein.n_residues, protein.n_frames)


def test_check_weights_invalid_weights_type(GS6_CP, NTL9_CP):
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        weights = 'abcdefgh'
        with pytest.raises(ValueError):
            protein._CTProtein__check_weights(weights=weights)


def test_check_weights_valid_uniform_weights_length(GS6_CP, NTL9_CP):
    default_weight = 1.0
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        weights = [default_weight for frame in range(protein.n_frames)]
        normalized_weights = [w/protein.n_frames for w in weights]
        protein._CTProtein__check_weights(weights=normalized_weights)


def test_check_weights_valid_uniform_weights_length_low_etol(GS6_CP, NTL9_CP):
    default_weight = 1.0
    tolerance = np.finfo(np.longdouble).precision
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        weights = [default_weight for frame in range(protein.n_frames)]
        normalized_weights = [w/protein.n_frames for w in weights]
        protein._CTProtein__check_weights(weights=normalized_weights, etol=tolerance)


def test_check_weights_valid_uniform_weights_length_large_etol(GS6_CP, NTL9_CP):
    default_weight = 1.0
    tolerance = 0.1
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        weights = [default_weight for frame in range(protein.n_frames)]
        normalized_weights = [w/protein.n_frames for w in weights]
        out = protein._CTProtein__check_weights(weights=normalized_weights, etol=tolerance)
        assert np.allclose(out, normalized_weights)


def test_check_weights_valid_nonuniform_weights_length_low_etol(GS6_CP, NTL9_CP):
    tolerance = np.finfo(np.longdouble).precision
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        weights = [w for w in range(1, protein.n_frames + 1)]
        normalized_weights = [w/sum(weights) for w in weights]
        protein._CTProtein__check_weights(weights=normalized_weights, etol=tolerance)


def test_check_weights_invalid_uniform_weights_length(GS6_CP, NTL9_CP):
    default_weight = 1.0
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        weights = [default_weight for frame in range(protein.n_frames - 1)]
        normalized_weights = [w/protein.n_frames for w in weights]
        with pytest.raises(CTException):
            protein._CTProtein__check_weights(weights=normalized_weights)


def test_check_weights_valid_uniform_weights_nondefault_stride(GS6_CP, NTL9_CP):
    stride = 2  # default stride is None
    default_weight = 1.0
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        weights = [default_weight for frame in range(protein.n_frames)]
        normalized_weights = [w/protein.n_frames for w in weights]
        protein._CTProtein__check_weights(weights=normalized_weights, stride=stride)


def test_check_weights_short_uniform_weights_nondefault_stride(GS6_CP, NTL9_CP):
    stride = 2  # default stride is None
    default_weight = 1.0
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        weights = [default_weight for frame in range(0, protein.n_frames, stride)]
        normalized_weights = [w/protein.n_frames for w in weights]
        with pytest.raises(CTException):
            protein._CTProtein__check_weights(weights=normalized_weights, stride=stride)


# == CTProtein.__get_first_and_last
# Should there be a test for R1 == R2?
def test_get_first_and_last_flip_residue_indices(GS6_CP, NTL9_CP):
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        residue_indices = list(range(protein.n_residues))
        mid_point = (protein.n_residues//2) - 1
        lower = residue_indices[:mid_point]
        upper = residue_indices[mid_point+1:]

        # Shuffle to ensure that the indices chosen, whilst random
        # are always where R1 > R2.
        random.shuffle(lower)
        random.shuffle(upper)
        r1 = random.choice(upper)
        r2 = random.choice(lower)
        protein._CTProtein__get_first_and_last(r1, r2)


# == CTProtein.__check_stride
def test_check_invalid_strides(GS6_CP, NTL9_CP):
    strides = [-1, 0, 100]  # none of the proteins have more than 100 frames.
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        for stride in strides:
            with pytest.raises(CTException):
                protein._CTProtein__check_stride(stride)


# == CTProtein.__check_single_residue
def test_check_invalid_single_residue(GS6_CP, NTL9_CP):
    invalid_residues = [-1, 100]  # none of the proteins have more than 100 residues
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        for invalid_residue_index in invalid_residues:
            with pytest.raises(CTException):
                protein._CTProtein__check_single_residue(invalid_residue_index)


# == CTProtein.__check_contains_CA
'''
def test_check_contains_CA_raises_exception(GS6_CP, NTL9_CP):
    # Requires modifying the protein topology to insert one or more fake residues where CA is missing
    # or replaced by something else.
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        residues_no_CA = protein.topology.select('resid %i to %i and atom not type "CA"' % (0, protein.n_residues - 1))
        print(residues_no_CA)
        for residue in residues_no_CA:
            with pytest.raises(CTException):
                protein._CTProtein__check_contains_CA(residue)
'''

def test_check_contains_CA_successful(GS6_CP, NTL9_CP):
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        residues_CA = protein.resid_with_CA
        for residue in residues_CA:
            return_value = protein._CTProtein__check_contains_CA(residue)
            assert return_value == None


# == CTProtein.__get_selection_atoms
def test_get_selection_atoms_region_of_size_2(GS6_CP, NTL9_CP):
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        region = (0, protein.n_residues - 1)
        choices = itertools.product([True, False], repeat=2)
        for backbone, heavy in choices:
            protein._CTProtein__get_selection_atoms(region, backbone, heavy)


def test_get_selection_atoms_invalid_region_of_size_3(GS6_CP, NTL9_CP):
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        region = (0, protein.n_residues - 1, 0)  # the last index doesn't matter (will raise exception regardless)
        choices = itertools.product([True, False], repeat=2)
        for backbone, heavy in choices:
            with pytest.raises(CTException):
                protein._CTProtein__get_selection_atoms(region, backbone, heavy)


# == CTProtein.get_amino_acid_sequence
def test_get_amino_acid_sequence(GS6_CP, NTL9_CP):
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        seq1 = protein.get_amino_acid_sequence(oneletter=True, numbered=False)

        # should there be a check for caps?
        assert len(seq1) == protein.n_residues

        seq2 = protein.get_amino_acid_sequence(oneletter=False, numbered=True)
        assert len(seq2) == protein.n_residues

        start_residue_number    = int(seq2[0].split('-')[-1])
        end_residue_number      = int(seq2[-1].split('-')[-1])

        assert start_residue_number == 1
        assert end_residue_number == len(seq2)
        assert (end_residue_number - start_residue_number) == len(seq2) - 1


def test_get_multiple_CA_index_existing(GS6_CP, NTL9_CP):
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        for residue_id in protein.resid_with_CA:
            atom_index = protein.get_multiple_CA_index(residue_id)

            # Validate the results and ensure they are ints
            assert len(atom_index) == 1
            assert residue_id in protein.resid_with_CA
            assert residue_id in protein._CTProtein__resid_with_CA
            assert type(atom_index[0]) in [np.int16, np.int32, np.int64]
            assert atom_index[0] > 0


def test_get_multiple_CA_index_invalid_residue_number(GS6_CP, NTL9_CP):
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        max_residue = protein.n_residues
        for residue_index in range(max_residue + 1, max_residue + protein.n_residues):
            with pytest.raises(CTException):
                protein.get_multiple_CA_index(residue_index)


def test_get_multiple_CA_index_invalid_residue_number_list(GS6_CP, NTL9_CP):
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        max_residue = protein.n_residues
        residue_list = list(range(max_residue + 1, max_residue + protein.n_residues))
        atoms_with_CA = protein.get_multiple_CA_index(resID_list=residue_list)

        assert type(atoms_with_CA) == list
        assert len(atoms_with_CA) == 0


def test_calculate_all_CA_distances_invalid_residue_number(GS6_CP, NTL9_CP):
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        max_residue = protein.n_residues
        for residue_index in range(max_residue + 1, max_residue + protein.n_residues):
            index = protein.calculate_all_CA_distances(residue_index)

            # Invalid indices return -1
            assert index == -1
