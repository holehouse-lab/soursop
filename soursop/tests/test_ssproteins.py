"""
Unit and regression test for the soursop package.
"""

# Import package, test suite, and other packages as needed
import numpy as np
import soursop
import pytest
import sys
import random
import itertools
from copy import deepcopy
from soursop import sstrajectory, ssprotein
from soursop.ssexceptions import SSException
from soursop.configs import DEBUGGING
from contextlib import contextmanager


def test_contacts_special(NTL9_CP):

    with pytest.raises(SSException) as error:
        a = NTL9_CP.get_contact_map(distance_thresh=2, stride=1, mode='sidechain-invalid')

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
    a = NTL9_CP.get_inter_residue_COM_distance(10,30)


    a = NTL9_CP.get_inter_residue_COM_vector(10,30)
    a = NTL9_CP.get_inter_residue_COM_vector(10,30)


    a = NTL9_CP.get_inter_residue_atomic_distance(2,10)
    a = NTL9_CP.get_inter_residue_atomic_distance(2,10, A1='CB')
    a = NTL9_CP.get_angle_decay()


def test_get_radius_of_gyration(GS6_CP):

    assert len(GS6_CP.get_radius_of_gyration()) == 5
    assert abs(GS6_CP.get_radius_of_gyration()[0] - 5.728453763896514) < 0.0001
    assert abs(GS6_CP.get_radius_of_gyration(R1=1,R2=3)[0] - 2.8815625777553278) < 0.0001
    assert GS6_CP.get_radius_of_gyration(R1=0,R2=7)[0]  == GS6_CP.get_radius_of_gyration()[0]


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


def test_soursop_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "soursop" in sys.modules


def test_DSSP(NTL9_CP):
    SS = NTL9_CP.get_secondary_structure_DSSP()
    # TO DO - complete


def test_BBSEG(NTL9_CP):
    SS = NTL9_CP.get_secondary_structure_BBSEG()
    # TO DO - complete




def test_phi_psi_bbseg(NTL9_CP):
    """
    Test that ensures our representation of the BBSEG2 matrix
    reproduces the version from CAMPARI

    """

    # read in BBSEG file
    bbseg_file = soursop.get_data('bbseg2.dat')
    with open(bbseg_file,'r') as fh:
        content = fh.readlines()

    idx=0
    for psi in range(179,-180, -10):

        local = ''
        for phi in range(-179,180, 10):
            local = local + '%i ' %(NTL9_CP._SSProtein__phi_psi_bbseg([phi],[psi])[0])

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

    # check nygaard mode implicit
    rh = CP.get_hydrodynamic_radius()
    assert abs(11.840569781179006 - rh[0]) < 0.001
    assert abs(np.mean(rh) - 11.3777) < 0.001

    # check nygaard mode explicit
    rh = CP.get_hydrodynamic_radius(mode='nygaard')
    assert abs(11.840569781179006 - rh[0]) < 0.001

    # check kirkwood-riseman mode explicit (implicit CA)
    rh = CP.get_hydrodynamic_radius(mode='kr')
    assert abs((5.9223 - rh[0])) < 0.001
    assert abs(np.mean(rh) - 5.884) < 0.001

    # check kirkwood-riseman mode explicit (explicit CA)
    rh = CP.get_hydrodynamic_radius(mode='kr', distance_mode='CA')
    assert abs((5.9223 - rh[0])) < 0.001
    assert abs(np.mean(rh) - 5.884) < 0.001

    # check kirkwood-riseman mode explicit with COM mode
    rh = CP.get_hydrodynamic_radius(mode='kr', distance_mode='COM')
    assert abs((5.8634 - rh[0])) < 0.001
    assert abs(np.mean(rh) - 5.8536) < 0.001

    # check it correctly raises an SSException if an invalid mode is passed
    with pytest.raises(SSException):
        rh = CP.get_hydrodynamic_radius(mode='nygaard_bad')
    with pytest.raises(SSException):
        rh = CP.get_hydrodynamic_radius(mode='nygaard_bad', distance_mode='CA_bad')

def test_get_molecular_volume(NTL9_CO):

    CP = NTL9_CO.proteinTrajectoryList[0]
    mol_vol = CP.get_molecular_volume()
    assert len(mol_vol) == 10
    assert mol_vol[0] - 20839.78797511 < 0.0000001
    assert np.mean(mol_vol) - 24064.633461433674 < 0.0000001




# ====
def test_init_from_trajectory(GS6_CO, NTL9_CO, GS6_CP, NTL9_CP):
    trajs = [GS6_CO, NTL9_CO]
    proteins = [GS6_CP, NTL9_CP]
    for ct_traj, ct_protein in zip(trajs, proteins):
        protein_from_md = ssprotein.SSProtein(ct_traj)
        assert protein_from_md.n_frames > 0
        assert protein_from_md.n_residues > 0
        assert protein_from_md.length() == ct_protein.length()

        protein_from_ct = ssprotein.SSProtein(ct_traj.traj)
        assert protein_from_ct.n_frames > 0
        assert protein_from_ct.n_residues > 0
        assert protein_from_ct.length() == ct_protein.length()


def test_init_from_trajectory_invalid():
    args = [1, 1.0, 'test', None]
    for arg in args:
        with pytest.raises(RuntimeError):
            protein = ssprotein.SSProtein(arg)


def test_init_from_trajectory_debug(GS6_CO, NTL9_CO, GS6_CP, NTL9_CP):
    trajs = [GS6_CO, NTL9_CO]
    proteins = [GS6_CP, NTL9_CP]
    for ct_traj, ct_protein in zip(trajs, proteins):
        protein_from_md = ssprotein.SSProtein(ct_traj.traj, debug=True)
        assert protein_from_md.n_frames > 0
        assert protein_from_md.n_residues > 0
        assert protein_from_md.length() == ct_protein.length()

        protein_from_ct = ssprotein.SSProtein(ct_traj, debug=True)
        assert protein_from_ct.n_frames > 0
        assert protein_from_ct.n_residues > 0
        assert protein_from_ct.length() == ct_protein.length()


def test_properties(GS6_CO, NTL9_CO, GS6_CP, NTL9_CP):
    properties = 'resid_with_CA,ncap,ccap,n_frames,n_residues,residue_index_list'.split(',')
    trajs = [GS6_CO, NTL9_CO]
    proteins = [GS6_CP, NTL9_CP]
    for ct_traj, ct_protein in zip(trajs, proteins):
        protein_from_ct = ssprotein.SSProtein(ct_traj)

        for prop in properties:
            assert getattr(protein_from_ct, prop) == getattr(ct_protein, prop)


def test_repr(GS6_CP, NTL9_CP):
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        repr_string = repr(protein)
        rebuilt_repr = "SSProtein (%s): %i res and %i frames" % (hex(id(protein)), protein.n_residues, protein.n_frames)
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
            protein._SSProtein__check_weights(weights=weights)


def test_check_weights_valid_uniform_weights_length(GS6_CP, NTL9_CP):
    default_weight = 1.0
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        weights = [default_weight for frame in range(protein.n_frames)]
        normalized_weights = [w/protein.n_frames for w in weights]
        protein._SSProtein__check_weights(weights=normalized_weights)


def test_check_weights_valid_uniform_weights_length_low_etol(GS6_CP, NTL9_CP):
    default_weight = 1.0
    tolerance = np.finfo(np.longdouble).precision
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        weights = [default_weight for frame in range(protein.n_frames)]
        normalized_weights = [w/protein.n_frames for w in weights]
        protein._SSProtein__check_weights(weights=normalized_weights, etol=tolerance)


def test_check_weights_valid_uniform_weights_length_large_etol(GS6_CP, NTL9_CP):
    default_weight = 1.0
    tolerance = 0.1
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        weights = [default_weight for frame in range(protein.n_frames)]
        normalized_weights = [w/protein.n_frames for w in weights]
        out = protein._SSProtein__check_weights(weights=normalized_weights, etol=tolerance)
        assert np.allclose(out, normalized_weights)


def test_check_weights_valid_nonuniform_weights_length_low_etol(GS6_CP, NTL9_CP):
    tolerance = np.finfo(np.longdouble).precision
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        weights = [w for w in range(1, protein.n_frames + 1)]
        normalized_weights = [w/sum(weights) for w in weights]
        protein._SSProtein__check_weights(weights=normalized_weights, etol=tolerance)


def test_check_weights_invalid_uniform_weights_length(GS6_CP, NTL9_CP):
    default_weight = 1.0
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        weights = [default_weight for frame in range(protein.n_frames - 1)]
        normalized_weights = [w/protein.n_frames for w in weights]
        with pytest.raises(SSException):
            protein._SSProtein__check_weights(weights=normalized_weights)


def test_check_weights_valid_uniform_weights_nondefault_stride(GS6_CP, NTL9_CP):
    stride = 2  # default stride is None
    default_weight = 1.0
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        weights = [default_weight for frame in range(protein.n_frames)]
        normalized_weights = [w/protein.n_frames for w in weights]
        protein._SSProtein__check_weights(weights=normalized_weights, stride=stride)


def test_check_weights_short_uniform_weights_nondefault_stride(GS6_CP, NTL9_CP):
    stride = 2  # default stride is None
    default_weight = 1.0
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        weights = [default_weight for frame in range(0, protein.n_frames, stride)]
        normalized_weights = [w/protein.n_frames for w in weights]
        with pytest.raises(SSException):
            protein._SSProtein__check_weights(weights=normalized_weights, stride=stride)


# == SSProtein.__get_first_and_last
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
        protein._SSProtein__get_first_and_last(r1, r2)


# == SSProtein.__check_stride
def test_check_invalid_strides(GS6_CP, NTL9_CP):
    strides = [-1, 0, 100]  # none of the proteins have more than 100 frames.
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        for stride in strides:
            with pytest.raises(SSException):
                protein._SSProtein__check_stride(stride)


# == SSProtein.__check_single_residue
def test_check_invalid_single_residue(GS6_CP, NTL9_CP):
    invalid_residues = [-1, 100]  # none of the proteins have more than 100 residues
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        for invalid_residue_index in invalid_residues:
            with pytest.raises(SSException):
                protein._SSProtein__check_single_residue(invalid_residue_index)


# == SSProtein.__check_contains_CA
'''
def test_check_contains_CA_raises_exception(GS6_CP, NTL9_CP):
    # Requires modifying the protein topology to insert one or more fake residues where CA is missing
    # or replaced by something else.
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        residues_no_CA = protein.topology.select('resid %i to %i and atom not type "CA"' % (0, protein.n_residues - 1))
        print(residues_no_CA)
        for residue in residues_no_CA:
            with pytest.raises(SSException):
                protein._SSProtein__check_contains_CA(residue)
'''

def test_check_contains_CA_successful(GS6_CP, NTL9_CP):
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        residues_CA = protein.resid_with_CA
        for residue in residues_CA:
            return_value = protein._SSProtein__check_contains_CA(residue)
            assert return_value == None


# == SSProtein.__get_selection_atoms
def test_get_selection_atoms_region_of_size_2(GS6_CP, NTL9_CP):
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        region = (0, protein.n_residues - 1)
        choices = itertools.product([True, False], repeat=2)
        for backbone, heavy in choices:
            protein._SSProtein__get_selection_atoms(region, backbone, heavy)


def test_get_selection_atoms_invalid_region_of_size_3(GS6_CP, NTL9_CP):
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        region = (0, protein.n_residues - 1, 0)  # the last index doesn't matter (will raise exception regardless)
        choices = itertools.product([True, False], repeat=2)
        for backbone, heavy in choices:
            with pytest.raises(SSException):
                protein._SSProtein__get_selection_atoms(region, backbone, heavy)


# == SSProtein.get_amino_acid_sequence
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
            assert residue_id in protein._SSProtein__resid_with_CA
            assert type(atom_index[0]) in [np.int16, np.int32, np.int64]
            assert atom_index[0] > 0


def test_get_multiple_CA_index_invalid_residue_number(GS6_CP, NTL9_CP):
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        max_residue = protein.n_residues
        for residue_index in range(max_residue + 1, max_residue + protein.n_residues):
            with pytest.raises(SSException):
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


# == SSProtein._SSProtein__residue_atom_index
def test_residue_atom_index_resid_atom_name_None(GS6_CO, NTL9_CO):
    trajs = [GS6_CO, NTL9_CO]
    for traj in trajs:
        # instantiate a new protein object since the previous reference is modified
        # elsewhere - hence our residue count will be off.
        protein = ssprotein.SSProtein(traj)
        unavailable_residue_index = protein.n_residues + 1
        protein.get_residue_atom_indices(unavailable_residue_index, atom_name=None)
        assert len(protein._SSProtein__residue_atom_table) == protein.n_residues + 1


def test_residue_atom_index_new_resid_atom_name_CA(GS6_CO, NTL9_CO):
    trajs = [GS6_CO, NTL9_CO]
    for traj in trajs:
        # instantiate a new protein object since the previous reference is modified
        # elsewhere - hence our residue count will be off.
        protein = ssprotein.SSProtein(traj)
        unavailable_residue_index = protein.n_residues + 1
        protein.get_residue_atom_indices(unavailable_residue_index, atom_name='CA')
        assert len(protein._SSProtein__residue_atom_table) == protein.n_residues + 1



def test_get_distance_map_weights(GS6_CP, NTL9_CP):
    default_weight = 1.0
    modes = 'CA,COM'.split(',')
    rms_options = [False, True]
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        weights = [default_weight for frame in range(protein.n_frames)]
        normalized_weights = [w/protein.n_frames for w in weights]
        for mode in modes:
            for rms_option in rms_options:
                results = protein.get_distance_map(mode=mode, RMS=rms_option, weights=normalized_weights)
                distance_map, std_dev = results

                protein_residues = len(protein.resid_with_CA)
                expected_shape = (protein_residues, protein_residues)

                assert distance_map.shape == expected_shape
                assert std_dev.shape == expected_shape

                if rms_option is True:
                    assert np.count_nonzero(np.tril(std_dev, -1)) == 0


def test_get_local_collapse_invalid_bins(GS6_CP, NTL9_CP):
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        # Wrong number of bins
        with pytest.raises(SSException):
            protein.get_local_collapse(bins=[1])

        # Wrong bins type
        with pytest.raises(SSException):
            protein.get_local_collapse(bins='a,b,c,d'.split(','))

        # Uneven bins spacing
        bins = [1,2,4,6]
        with pytest.raises(SSException):
            protein.get_local_collapse(bins=bins)


def test_get_local_collapse_invalid_window_size(GS6_CP, NTL9_CP):
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        number_of_residues = protein.n_residues + 1
        with pytest.raises(SSException):
            protein.get_local_collapse(window_size=number_of_residues)


def test_get_local_collapse_successful_run(NTL9_CP):
    # get_local_collapse(self, window_size=10, bins=None, verbose=True)
    window_size = 10
    meanData, stdData, histo, bins = NTL9_CP.get_local_collapse()

    expected_length = NTL9_CP.n_residues - (window_size - 1)
    assert len(meanData) == expected_length
    assert len(stdData) == expected_length
    assert len(histo) == expected_length
    assert type(bins) != None


# SSProtein.get_angle_decay
def test_get_angle_decay_return_all_pairs(GS6_CP, NTL9_CP):
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        result = protein.get_angle_decay(return_all_pairs=True)
        assert len(result) == 2

        return_matrix, all_pairs = result

        # check this has been calculated correctly  -we expect all_pairs
        # to be a dictionary that reports on every unique inter-residue pair
        assert len(all_pairs) == np.sum(list(range(len(return_matrix), 0, -1)))



def test_get_angle_decay_no_return_all_pairs(GS6_CP, NTL9_CP):
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        return_matrix = protein.get_angle_decay(return_all_pairs=False)


        num_caps = 0
        if protein.ncap:
            num_caps += 1
        if protein.ccap:
            num_caps += 1

        assert len(return_matrix) == protein.n_residues - num_caps

def test_get_angle_decay_consistent_value(GS6_CP, NTL9_CP):
    """
    Test that ensures the average returned in the return matrix
    matches the average that can be manually calculated from the
    individual pairs

    """
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:

        # do this because we only caculate vector
        # between res with CA, the indices here are position within the n-to-c
        # vector that always starts at 1
        # separation not index position, so we are always separation of 1-x
        min_res = 1
        max_res = protein.n_residues - (protein.n_residues - len(protein.resid_with_CA))

        (return_matrix, pair_dict) = protein.get_angle_decay(return_all_pairs=True)

        for window in range(1,8):

            if window + min_res >= max_res+1:
                continue
            
            all_pairs = []
            for i in range(min_res, (max_res+1)-window):
                j = i+window
                n = f"{i}-{j}"
                all_pairs.append(pair_dict[n])
            
            assert (np.mean(all_pairs) - return_matrix[window][1]) == 0
        

# SSProtein.get_contact_map
def test_get_contact_map_weights(GS6_CP, NTL9_CP):
    default_weight = 1.0

    # It appears that sidechain-heavy raises an error for GLY has returned in mdtraj 1.9.5+.
    # A dynamic check is made for a failing test when that mode is encountered.
    modes = 'ca,closest,closest-heavy,sidechain,sidechain-heavy'.split(',')
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        weights = [default_weight for frame in range(protein.n_frames)]
        normalized_weights = [w/protein.n_frames for w in weights]
        for mode in modes:
            if mode == 'sidechain-heavy':
                with pytest.raises(SSException):
                    results = protein.get_contact_map(mode=mode, weights=normalized_weights)  # the other options have been tested
                continue
            results = protein.get_contact_map(mode=mode, weights=normalized_weights)  # the other options have been tested
            contact_map, contact_map_order = results

            protein_residues = len(protein.resid_with_CA)
            expected_shape = (protein_residues, protein_residues)

            assert contact_map.shape == expected_shape
            assert len(contact_map_order.shape) == 1
            assert contact_map_order.shape[0] == protein_residues


def test_get_contact_map_weights_invalid_mode(GS6_CP, NTL9_CP):
    # Checking that
    #

    default_weight = 1.0
    proteins = [GS6_CP, NTL9_CP]
    # this should fail with an SSException

    mode = 'invalid-mode'
    for protein in proteins:
        weights = [default_weight for frame in range(protein.n_frames)]
        normalized_weights = [w/protein.n_frames for w in weights]


        with pytest.raises(SSException):
            protein.get_contact_map(mode=mode, weights=normalized_weights)  # the other options have been tested




# SSProtein.get_sidechain_alignment_angle
def test_get_sidechain_alignment_angle_unsupported_angle(GS6_CP, NTL9_CP):
    proteins = [GS6_CP, NTL9_CP]
    R1 = 0
    R2 = 0
    for protein in proteins:
        with pytest.raises(SSException):
            protein.get_sidechain_alignment_angle(R1, R2, sidechain_atom_1='invalid_atom')

        with pytest.raises(SSException):
            protein.get_sidechain_alignment_angle(R1, R2, sidechain_atom_2='invalid_atom')


# From: https://stackoverflow.com/a/42327075
@contextmanager
def not_raises(exception):
    try:
        yield
    except exception:
        raise pytest.fail("DID RAISE {0}".format(exception))


def test_get_sidechain_alignment_angle_successful_random(GS6_CP, NTL9_CP):
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        residue_indices = list(range(protein.n_residues))
        mid_point = (protein.n_residues//2) - 1
        lower = residue_indices[1:mid_point]
        upper = residue_indices[mid_point+1:-1]

        # Shuffle to ensure that the indices chosen, whilst random
        # are always where R1 > R2.
        random.shuffle(lower)
        random.shuffle(upper)
        r1 = random.choice(upper)
        r2 = random.choice(lower)

        # ensure that our residues do not reference GLY
        # otherwise this test will fail. Since we only
        # wish to test successful alignments, we skip over
        # these residues.
        df, indices = protein.topology.to_dataframe()
        error_residue_indices = list()
        error_residues = 'ACE,GLY,NME'.split(',')

        # convert from ids to indices (subtract 1)
        for residue_name in error_residues:
            error_residue_indices += list(sorted(set([i-1 for i in df[df['resName'] == residue_name]['resSeq']])))

        reset_loop = False
        for resid in [r1, r2]:
            if resid in error_residue_indices:
                reset_loop = True
                break

        if reset_loop:
            continue

        with not_raises(SSException):
            alignment = protein.get_sidechain_alignment_angle(r1, r2)


def test_get_sidechain_alignment_angle_failing_random(GS6_CP, NTL9_CP):
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        # Identify the indices of the residues that will fail.
        df, indices = protein.topology.to_dataframe()
        error_residue_indices = list()
        error_residues = 'ACE,GLY,NME'.split(',')

        # convert from ids to indices (subtract 1)
        for residue_name in error_residues:
            error_residue_indices += list(sorted(set([i-1 for i in df[df['resName'] == residue_name]['resSeq']])))

        r1 = random.choice(error_residue_indices)
        r2 = random.choice(error_residue_indices)

        with pytest.raises(SSException):
            alignment = protein.get_sidechain_alignment_angle(r1, r2)


def test_get_sidechain_alignment_angle_invalid_r1(GS6_CP, NTL9_CP):
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        max_residue = protein.n_residues
        r1 = max_residue + protein.n_residues
        r2 = 0

        with pytest.raises(SSException):
            protein.get_sidechain_alignment_angle(r1, r2)


def test_get_sidechain_alignment_angle_invalid_r2(GS6_CP, NTL9_CP):
    proteins = [GS6_CP, NTL9_CP]
    for protein in proteins:
        max_residue = protein.n_residues
        r1 = 0
        r2 = max_residue + protein.n_residues

        with pytest.raises(SSException):
            protein.get_sidechain_alignment_angle(r1, r2)


#get_local_to_global_correlation(self, mode='COM', n_cycles=100, max_num_pairs=10, stride=20, weights=False, verbose=True)
def test_get_local_to_global_correlation_unsupported_modes(GS6_CP, NTL9_CP):
    proteins = [GS6_CP, NTL9_CP]
    invalid_modes = 'unknown,unsupported'.split(',')
    for protein in proteins:
        for mode in invalid_modes:
            with pytest.raises(SSException):
                protein.get_local_to_global_correlation(mode=mode)


"""
# Requires additional work. Exception raised:
# "raise SSException('Something when wrong when comparing stride-derived Rg and internal distances, this is a bug in the code...)')"
def test_get_local_to_global_correlation_defaults(GS6_CP, NTL9_CP, cta_protein_helper):
    num_copies = 5
    proteins = [GS6_CP, NTL9_CP]
    modes = 'CA,COM'.split(',')
    for protein_obj in proteins:
        protein = cta_protein_helper.lengthen_protein_trajectory(protein_obj, num_copies)
        for mode in modes:
            results = protein.get_local_to_global_correlation(mode=mode)

            assert len(results) == 4
            for result in results:
                assert len(result) > 0
"""


def test_get_local_to_global_correlation_valid_weights(GS6_CP, NTL9_CP, cta_protein_helper):
    num_copies = 5
    stride = 1
    default_weight = 1.0
    proteins = [GS6_CP, NTL9_CP]
    modes = 'CA,COM'.split(',')
    for protein_obj in proteins:
        protein = cta_protein_helper.lengthen_protein_trajectory(protein_obj, num_copies)
        weights = [default_weight for frame in range(protein.n_frames)]
        normalized_weights = [w/protein.n_frames for w in weights]
        for mode in modes:
            results = protein.get_local_to_global_correlation(mode=mode, stride=stride, weights=normalized_weights)

            assert len(results) == 4
            for result in results:
                assert len(result) > 0


def test_get_regional_SASA(GS6_CP, NTL9_CP, cta_protein_helper):
    num_copies = 5
    proteins = [GS6_CP, NTL9_CP]
    for protein_obj in proteins:
        protein = cta_protein_helper.lengthen_protein_trajectory(protein_obj, num_copies)

        residue_indices = list(range(protein.n_residues))
        mid_point = (protein.n_residues//2) - 1
        lower = residue_indices[1:mid_point]
        upper = residue_indices[mid_point+1:-1]

        # Shuffle to ensure that the indices chosen, whilst random
        # are always where R1 > R2.
        random.shuffle(lower)
        random.shuffle(upper)
        r1 = random.choice(upper)
        r2 = random.choice(lower)

        rsasa = protein.get_regional_SASA(r1, r2)
        assert rsasa != None


def test_get_all_SASA(GS6_CP, NTL9_CP):

    proteins = [GS6_CP, NTL9_CP]

    # check default
    assert np.isclose(np.min(GS6_CP.get_all_SASA(stride=1)), 56.75124)
    assert np.isclose(np.max(GS6_CP.get_all_SASA(stride=1)), 144.38452)
    assert np.isclose(np.mean(GS6_CP.get_all_SASA(stride=1)), 108.992676)
    
    # check residue mode
    assert np.isclose(np.min(GS6_CP.get_all_SASA(stride=1, mode='residue')), 56.75124)
    assert np.isclose(np.max(GS6_CP.get_all_SASA(stride=1, mode='residue')), 144.38452)
    assert np.isclose(np.mean(GS6_CP.get_all_SASA(stride=1, mode='residue')), 108.992676)

    # check variable stride works
    assert np.isclose(np.mean(GS6_CP.get_all_SASA(stride=3, mode='residue')), 108.51016)

    # check shapes are ok
    assert GS6_CP.get_all_SASA(stride=1, mode='residue').shape == (5,8)
    assert GS6_CP.get_all_SASA(stride=1, mode='atom').shape == (5,66)

    # check atom values work
    assert np.min(GS6_CP.get_all_SASA(stride=1, mode='atom')) == 0.0
    assert np.isclose(np.max(GS6_CP.get_all_SASA(stride=1, mode='atom')), 38.101517) 
    assert np.isclose(np.mean(GS6_CP.get_all_SASA(stride=1, mode='atom')), 13.211235)

    # check backbone values work
    assert np.isclose(np.sum(GS6_CP.get_all_SASA(stride=1, mode='backbone')), 1894.9926)
    assert np.isclose(np.sum(GS6_CP.get_all_SASA(stride=1, mode='sidechain')), 1309.465)
    assert np.isclose(np.min(GS6_CP.get_all_SASA(stride=1, mode='backbone')), 27.768456)
    assert np.isclose(np.min(GS6_CP.get_all_SASA(stride=1, mode='sidechain')), 0.0)

    # check dimesions
    GS6_CP.get_all_SASA(stride=1, mode='backbone').shape == (5,6)
    GS6_CP.get_all_SASA(stride=1, mode='sidechain').shape == (5,6)
    GS6_CP.get_all_SASA(stride=1, mode='backbone').shape == (5,6)
    GS6_CP.get_all_SASA(stride=1, mode='sidechain').shape == (5,6)

    # check that 'all' works
    assert len(GS6_CP.get_all_SASA(stride=1, mode='all')) == 3

    # assert all works as expected
    assert np.sum(GS6_CP.get_all_SASA(stride=1, mode='all')[0] == GS6_CP.get_all_SASA(stride=1)) == 40
    assert np.sum(GS6_CP.get_all_SASA(stride=1, mode='all')[1] == GS6_CP.get_all_SASA(stride=1, mode='sidechain')) == 30
    assert np.sum(GS6_CP.get_all_SASA(stride=1, mode='all')[2] == GS6_CP.get_all_SASA(stride=1, mode='backbone')) == 30
        

def test_get_site_accessibility_resid(GS6_CP, NTL9_CP, cta_protein_helper):
    num_copies = 5
    proteins = [GS6_CP, NTL9_CP]
    for protein_obj in proteins:
        protein = cta_protein_helper.lengthen_protein_trajectory(protein_obj, num_copies)

        residue_indices = list(range(protein.n_residues))
        protein.get_site_accessibility(residue_indices, mode='resid')


def test_get_site_accessibility_residue_names(GS6_CP, NTL9_CP, cta_protein_helper):
    num_copies = 5
    proteins = [GS6_CP, NTL9_CP]
    residue_names = ['GLY']
    for protein_obj in proteins:
        protein = cta_protein_helper.lengthen_protein_trajectory(protein_obj, num_copies)

        protein.get_site_accessibility(residue_names, mode='residue_type')


def test_get_angles_1(GS6_CP, NTL9_CP):

    assert len(GS6_CP.get_angles('phi')[0]) == 6
    assert GS6_CP.get_angles('phi')[0][0][0] == 'ACE1-C'
    assert GS6_CP.get_angles('phi')[0][5][3] == 'SER7-C'
    assert np.isclose(GS6_CP.get_angles('phi')[1][5][3], -147.67358)
    assert np.isclose(GS6_CP.get_angles('phi')[1][5][0], -126.707436)
    assert np.isclose(GS6_CP.get_angles('phi')[1][0][0], -157.64288)
    assert np.isclose(GS6_CP.get_angles('psi')[1][0][0] , -41.66139)
    assert np.isclose(GS6_CP.get_angles('psi')[1][2][0] , 150.39041)
    assert np.isclose(GS6_CP.get_angles('psi')[1][3][3], -38.81799)
    assert np.isclose(GS6_CP.get_angles('psi')[1][4][3] , -89.94062)

    assert GS6_CP.get_angles('chi1')[0][0][3] == 'SER3-OG'
    assert GS6_CP.get_angles('chi1')[0][0][2] == 'SER3-CB'
    assert GS6_CP.get_angles('chi1')[0][0][1] == 'SER3-CA'
    assert GS6_CP.get_angles('chi1')[0][0][0] == 'SER3-N'

    assert len(NTL9_CP.get_angles('phi')[0]) == 55
    assert len(NTL9_CP.get_angles('psi')[0]) == 55
    assert len(NTL9_CP.get_angles('omega')[0]) == 55
    assert len(NTL9_CP.get_angles('chi1')[0]) == 44
    assert len(NTL9_CP.get_angles('chi2')[0]) == 40
    assert len(NTL9_CP.get_angles('chi3')[0]) == 21
    assert len(NTL9_CP.get_angles('chi4')[0]) == 12
    assert len(NTL9_CP.get_angles('chi5')[0]) == 1



def test_get_inter_residue_COM_vector(GS6_CP, NTL9_CP):
    """
    Test that the inter-residue COM vector matches the distances
    calculated using inter_residue_COM_distance.

    """


    com_distances = GS6_CP.get_inter_residue_COM_distance(0,5)
    com_distances_manually = np.sqrt(np.sum(np.square(GS6_CP.get_inter_residue_COM_vector(0,5)),1))
    assert np.sum(com_distances - com_distances_manually) < 0.001


    com_distances = NTL9_CP.get_inter_residue_COM_distance(0,25)
    com_distances_manually = np.sqrt(np.sum(np.square(NTL9_CP.get_inter_residue_COM_vector(0, 25)),1))
    assert np.sum(com_distances - com_distances_manually) < 0.001



def test_get_center_of_mass(GS6_CP, NTL9_CP):


    molecular_com = GS6_CP.get_center_of_mass()
    assert(len(molecular_com)) == GS6_CP.n_frames

    # compute COM of residues 0 and 1
    R0_com = GS6_CP.get_center_of_mass(0,0)
    R5_com = GS6_CP.get_center_of_mass(5,5)

    com_distances_from_center_of_mass = np.sqrt(np.sum(np.square(R0_com - R5_com),1))
    com_distances_manually = np.sqrt(np.sum(np.square(GS6_CP.get_inter_residue_COM_vector(0,5)),1))
    com_distances = GS6_CP.get_inter_residue_COM_distance(0,5)

    assert np.sum(com_distances_manually - com_distances_from_center_of_mass) == 0
    assert np.sum(com_distances - com_distances_from_center_of_mass) == 0
