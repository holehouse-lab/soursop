import os
import shutil
import numpy as np
import mdtraj as md
import pytest
import tempfile
import camparitraj.ctanalyzer.analyzer_analysis as cta_aa
from camparitraj.ctprotein import CTProtein


def test_run_CG(GS6_CP, NTL9_CP, cta_output_files_helper):
    ext = 'csv'
    prefix = 'RG'
    names = ',mean,std'.split(',')  # Note the lead `,`; this will be empty when split.
    expected_filenames = cta_output_files_helper.determine_filenames(prefix, names, ext)

    for protein in [GS6_CP, NTL9_CP]:
        # Create a temporary folder (which will autoclean) to save the analysis.
        with tempfile.TemporaryDirectory() as outdir:
            cta_aa.run_RG(protein, outdir)

            # Verify that the files exist, and are non-zero.
            for filename in expected_filenames:
                savename = os.path.join(outdir, filename)

                assert os.path.getsize(savename) > 0  # implicit check for file existence
                with open(savename) as f:
                    data = np.loadtxt(f)  # numpy is operating on a buffer b/c of the temporaryfile.
                    d = np.frombuffer(data)  # convert the buffer back to an array, and compare.
                    assert np.allclose(data, d)


def test_run_t_inst(GS6_CP, NTL9_CP):
    prefix = 't'
    names = 'inst,mean,std'.split(',')  # Note the lead `,`; this will be empty when split.
    expected_filenames = ['%s.csv' % '_'.join([prefix, name]) for name in names]

    for protein in [GS6_CP, NTL9_CP]:
        # Create a temporary folder (which will autoclean) to save the analysis.
        with tempfile.TemporaryDirectory() as outdir:
            cta_aa.run_t_inst(protein, outdir)

            # Verify that the files exist, and are non-zero.
            for filename in expected_filenames:
                savename = os.path.join(outdir, filename)

                assert os.path.getsize(savename) > 0  # implicit check for file existence
                with open(savename) as f:
                    data = np.loadtxt(f)  # numpy is operating on a buffer b/c of the temporaryfile.
                    d = np.frombuffer(data)  # convert the buffer back to an array, and compare.
                    assert np.allclose(data, d)


def test_run_RH(GS6_CP, NTL9_CP, cta_output_files_helper):
    ext = 'csv'
    prefix = 'RH'
    names = ',mean,std'.split(',')  # Note the lead `,`; this will be empty when split.
    expected_filenames = cta_output_files_helper.determine_filenames(prefix, names, ext)

    for protein in [GS6_CP, NTL9_CP]:
        # Create a temporary folder (which will autoclean) to save the analysis.
        with tempfile.TemporaryDirectory() as outdir:
            cta_aa.run_RH(protein, outdir)

            # Verify that the files exist, and are non-zero.
            for filename in expected_filenames:
                savename = os.path.join(outdir, filename)

                assert os.path.getsize(savename) > 0  # implicit check for file existence
                with open(savename) as f:
                    data = np.loadtxt(f)  # numpy is operating on a buffer b/c of the temporaryfile.
                    d = np.frombuffer(data)  # convert the buffer back to an array, and compare.
                    assert np.allclose(data, d)


def test_run_end_to_end(GS6_CP, NTL9_CP, cta_output_files_helper):
    ext = 'csv'
    prefix = 'end_to_end'
    names = ',mean,std'.split(',')  # Note the lead `,`; this will be empty when split.
    expected_filenames = cta_output_files_helper.determine_filenames(prefix, names, ext)

    for protein in [GS6_CP, NTL9_CP]:
        # Create a temporary folder (which will autoclean) to save the analysis.
        with tempfile.TemporaryDirectory() as outdir:
            cta_aa.run_end_to_end(protein, outdir)

            # Verify that the files exist, and are non-zero.
            for filename in expected_filenames:
                savename = os.path.join(outdir, filename)

                assert os.path.getsize(savename) > 0  # implicit check for file existence
                with open(savename) as f:
                    data = np.loadtxt(f)  # numpy is operating on a buffer b/c of the temporaryfile.
                    d = np.frombuffer(data)  # convert the buffer back to an array, and compare.
                    assert np.allclose(data, d)


def test_run_asphericity(GS6_CP, NTL9_CP, cta_output_files_helper):
    ext = 'csv'
    prefix = 'ASPH'
    names = ',mean,std'.split(',')  # Note the lead `,`; this will be empty when split.
    expected_filenames = cta_output_files_helper.determine_filenames(prefix, names, ext)

    for protein in [GS6_CP, NTL9_CP]:
        # Create a temporary folder (which will autoclean) to save the analysis.
        with tempfile.TemporaryDirectory() as outdir:
            cta_aa.run_asphericity(protein, outdir)

            # Verify that the files exist, and are non-zero.
            for filename in expected_filenames:
                savename = os.path.join(outdir, filename)

                assert os.path.getsize(savename) > 0  # implicit check for file existence
                with open(savename) as f:
                    data = np.loadtxt(f)  # numpy is operating on a buffer b/c of the temporaryfile.
                    d = np.frombuffer(data)  # convert the buffer back to an array, and compare.
                    assert np.allclose(data, d)


def test_run_distanceMap(GS6_CP, NTL9_CP, cta_output_files_helper):
    ext = 'csv'
    prefix = 'distance_map'
    names = ',std'.split(',')  # Note the lead `,`; this will be empty when split.
    expected_filenames = cta_output_files_helper.determine_filenames(prefix, names, ext)

    for protein in [GS6_CP, NTL9_CP]:
        # Create a temporary folder (which will autoclean) to save the analysis.
        with tempfile.TemporaryDirectory() as outdir:
            cta_aa.run_distanceMap(protein, outdir)

            # Verify that the files exist, and are non-zero.
            for filename in expected_filenames:
                savename = os.path.join(outdir, filename)

                assert os.path.getsize(savename) > 0  # implicit check for file existence
                with open(savename) as f:
                    data = np.loadtxt(f, delimiter=',')  # numpy is operating on a buffer b/c of the temporaryfile.
                    d = np.frombuffer(data).reshape(data.shape)  # convert the buffer back to an array, and compare.
                    assert np.allclose(data, d)


def test_run_polymer_scaling_map(NTL9_CP, cta_output_files_helper):
    ext = 'csv'
    prefix = 'polymer_deviation_map'
    names = 'absolute,fractional,params'.split(',')  # Note the lead `,`; this will be empty when split.
    expected_filenames = cta_output_files_helper.determine_filenames(prefix, names, ext)

    # Create a temporary folder (which will autoclean) to save the analysis.
    with tempfile.TemporaryDirectory() as outdir:
        cta_aa.run_polymer_scaling_map(NTL9_CP, outdir)

        # Verify that the files exist, and are non-zero.
        for filename in expected_filenames:
            savename = os.path.join(outdir, filename)

            assert os.path.getsize(savename) > 0  # implicit check for file existence
            with open(savename) as f:
                data = np.loadtxt(f, delimiter=',')  # numpy is operating on a buffer b/c of the temporaryfile.
                d = np.frombuffer(data).reshape(data.shape)  # convert the buffer back to an array, and compare.
                assert np.allclose(data, d)


def test_run_internal_scaling(GS6_CP, NTL9_CP):
    expected_filename = 'INTSCAL.csv'
    for protein in [GS6_CP, NTL9_CP]:
        with tempfile.TemporaryDirectory() as outdir:
            cta_aa.run_internal_scaling(protein, outdir)
            savename = os.path.join(outdir, expected_filename)

            assert os.path.getsize(savename) > 0  # implicit check for file existence
            with open(savename) as f:
                data = np.loadtxt(f, delimiter=',')  # numpy is operating on a buffer b/c of the temporaryfile.
                d = np.frombuffer(data).reshape(data.shape)  # convert the buffer back to an array, and compare.
                assert np.allclose(data, d)


def test_run_RMS_internal_scaling(GS6_CP, NTL9_CP):
    expected_filename = 'RMS_INTSCAL.csv'
    for protein in [GS6_CP, NTL9_CP]:
        with tempfile.TemporaryDirectory() as outdir:
            cta_aa.run_RMS_internal_scaling(protein, outdir)
            savename = os.path.join(outdir, expected_filename)

            assert os.path.getsize(savename) > 0  # implicit check for file existence
            with open(savename) as f:
                data = np.loadtxt(f, delimiter=',')  # numpy is operating on a buffer b/c of the temporaryfile.
                d = np.frombuffer(data).reshape(data.shape)  # convert the buffer back to an array, and compare.
                assert np.allclose(data, d)


def test_run_fractal_deviation(GS6_CP, NTL9_CP):
    expected_filename = 'fractional_deviation.csv'
    for protein in [GS6_CP, NTL9_CP]:
        # The ideal stride is ~ 20; 10 appears to be a reasonable low limit test.
        for stride in range(10, protein.n_frames):
            with tempfile.TemporaryDirectory() as outdir:
                cta_aa.run_fractal_deviation(protein, outdir, stride)
                savename = os.path.join(outdir, expected_filename)

                assert os.path.getsize(savename) > 0  # implicit check for file existence
                with open(savename) as f:
                    data = np.loadtxt(f, delimiter=',')  # numpy is operating on a buffer b/c of the temporaryfile.
                    d = np.frombuffer(data).reshape(data.shape)  # convert the buffer back to an array, and compare.
                    assert np.allclose(data, d)


def test_run_Q_analysis(NTL9_CP):
    # This only appears to work on NTL9; GS6 has problems as it creates an empty file.
    # Moreover, NTL9
    expected_filename = 'Q_res_by_res.txt'
    with tempfile.TemporaryDirectory() as outdir:
        cta_aa.run_Q_analysis(NTL9_CP, outdir)
        savename = os.path.join(outdir, expected_filename)

        assert os.path.getsize(savename) > 0  # implicit check for file existence
        with open(savename) as f:
            data = np.loadtxt(f, delimiter=',')  # numpy is operating on a buffer b/c of the temporaryfile.
            d = np.frombuffer(data).reshape(data.shape)  # convert the buffer back to an array, and compare.
            assert np.allclose(data, d)


def test_run_rg_re_correlation(GS6_CP, NTL9_CP):
    expected_filename = 'rg_re_corr.csv'
    for protein in [GS6_CP, NTL9_CP]:
        with tempfile.TemporaryDirectory() as outdir:
            cta_aa.run_rg_re_correlation(protein, outdir)
            savename = os.path.join(outdir, expected_filename)

            assert os.path.getsize(savename) > 0  # implicit check for file existence
            with open(savename) as f:
                data = np.loadtxt(f, delimiter=',')  # numpy is operating on a buffer b/c of the temporaryfile.
                d = np.frombuffer(data).reshape(data.shape)  # convert the buffer back to an array, and compare.
                assert np.allclose(data, d)


def test_run_rij_analysis(GS6_CP, NTL9_CP, cta_output_files_helper):
    ext = 'csv'
    for protein in [GS6_CP, NTL9_CP]:
        for start_residue_index in range(protein.n_residues - 1):
            end_residue_index = protein.n_residues - 1
            prefix = 'r_%d_%d' % (start_residue_index, end_residue_index)
            names = ',mean,std'.split(',')  # Note the lead `,`; this will be empty when split.
            expected_filenames = cta_output_files_helper.determine_filenames(prefix, names, ext)

            with tempfile.TemporaryDirectory() as outdir:
                cta_aa.run_rij_analysis(protein, outdir, start_residue_index, end_residue_index)
                for expected_filename in expected_filenames:
                    savename = os.path.join(outdir, expected_filename)

                    assert os.path.getsize(savename) > 0  # implicit check for file existence
                    with open(savename) as f:
                        data = np.loadtxt(f, delimiter=',')  # numpy is operating on a buffer b/c of the temporaryfile.
                        d = np.frombuffer(data).reshape(data.shape)  # convert the buffer back to an array, and compare.
                        assert np.allclose(data, d)


def test_run_scaling_exponent_power_default_end_effect(NTL9_CP, cta_output_files_helper):
    # GS6_CP is too short for this analysis
    ext = 'csv'
    prefix = 'scaling_exp'
    end_effect = 5  # 5 is the default
    names = 'analysis_power,idx_used_power,fit_power'.split(',')
    expected_filenames = cta_output_files_helper.determine_filenames(prefix, names, ext)
    with tempfile.TemporaryDirectory() as outdir:
        cta_aa.run_scaling_exponent_power(NTL9_CP, outdir, end_effect)
        for expected_filename in expected_filenames:
            savename = os.path.join(outdir, expected_filename)

            assert os.path.getsize(savename) > 0  # implicit check for file existence
            with open(savename) as f:
                data = np.loadtxt(f, delimiter=',')  # numpy is operating on a buffer b/c of the temporaryfile.
                d = np.frombuffer(data).reshape(data.shape)  # convert the buffer back to an array, and compare.
                assert np.allclose(data, d)


def test_run_scaling_exponent_power_nondefault_end_effect(NTL9_CP, cta_output_files_helper):
    # GS6_CP is too short for this analysis
    ext = 'csv'
    prefix = 'scaling_exp'
    end_effect = 4  # 5 is the default
    names = [n % end_effect for n in 'analysis_ee%i_power,idx_used_ee%i_power,fit_ee%i_power'.split(',')]
    expected_filenames = cta_output_files_helper.determine_filenames(prefix, names, ext)
    with tempfile.TemporaryDirectory() as outdir:
        cta_aa.run_scaling_exponent_power(NTL9_CP, outdir, end_effect)
        for expected_filename in expected_filenames:
            savename = os.path.join(outdir, expected_filename)

            assert os.path.getsize(savename) > 0  # implicit check for file existence
            with open(savename) as f:
                data = np.loadtxt(f, delimiter=',')  # numpy is operating on a buffer b/c of the temporaryfile.
                d = np.frombuffer(data).reshape(data.shape)  # convert the buffer back to an array, and compare.
                assert np.allclose(data, d)

def test_run_scaling_exponent_power_CA_default_end_effect(NTL9_CP, cta_output_files_helper):
    # GS6_CP is too short for this analysis
    ext = 'csv'
    prefix = 'scaling_exp'
    end_effect = 5  # 5 is the default
    names = 'analysis_power_CA,idx_used_power_CA,fit_power_CA'.split(',')
    expected_filenames = cta_output_files_helper.determine_filenames(prefix, names, ext)
    with tempfile.TemporaryDirectory() as outdir:
        cta_aa.run_scaling_exponent_power_CA(NTL9_CP, outdir, end_effect)
        for expected_filename in expected_filenames:
            savename = os.path.join(outdir, expected_filename)

            assert os.path.getsize(savename) > 0  # implicit check for file existence
            with open(savename) as f:
                data = np.loadtxt(f, delimiter=',')  # numpy is operating on a buffer b/c of the temporaryfile.
                d = np.frombuffer(data).reshape(data.shape)  # convert the buffer back to an array, and compare.
                assert np.allclose(data, d)


def test_run_scaling_exponent_power_CA_nondefault_end_effect(NTL9_CP, cta_output_files_helper):
    # GS6_CP is too short for this analysis
    ext = 'csv'
    prefix = 'scaling_exp'
    end_effect = 4  # 5 is the default
    names = [n % end_effect for n in 'analysis_ee%i_power_CA,idx_used_ee%i_power_CA,fit_ee%i_power_CA'.split(',')]
    expected_filenames = cta_output_files_helper.determine_filenames(prefix, names, ext)
    with tempfile.TemporaryDirectory() as outdir:
        cta_aa.run_scaling_exponent_power_CA(NTL9_CP, outdir, end_effect)
        for expected_filename in expected_filenames:
            savename = os.path.join(outdir, expected_filename)

            assert os.path.getsize(savename) > 0  # implicit check for file existence
            with open(savename) as f:
                data = np.loadtxt(f, delimiter=',')  # numpy is operating on a buffer b/c of the temporaryfile.
                d = np.frombuffer(data).reshape(data.shape)  # convert the buffer back to an array, and compare.
                assert np.allclose(data, d)


def test_run_motif_RG(GS6_CP, NTL9_CP, cta_output_files_helper):
    ext = 'csv'
    for protein in [GS6_CP, NTL9_CP]:
        for start_residue_index in range(protein.n_residues - 1):
            end_residue_index = protein.n_residues - 1
            prefix = 'motif_%d_%d_RG' % (start_residue_index, end_residue_index)
            names = ',mean,std'.split(',')  # Note the lead `,`; this will be empty when split.
            expected_filenames = cta_output_files_helper.determine_filenames(prefix, names, ext)

            with tempfile.TemporaryDirectory() as outdir:
                cta_aa.run_motif_RG(protein, outdir, start_residue_index, end_residue_index)
                for expected_filename in expected_filenames:
                    savename = os.path.join(outdir, expected_filename)

                    assert os.path.getsize(savename) > 0  # implicit check for file existence
                    with open(savename) as f:
                        data = np.loadtxt(f, delimiter=',')  # numpy is operating on a buffer b/c of the temporaryfile.
                        d = np.frombuffer(data).reshape(data.shape)  # convert the buffer back to an array, and compare.
                        assert np.allclose(data, d)


def test_run_DSSP_analysis(GS6_CP, NTL9_CP, cta_output_files_helper):
    ext = 'csv'
    prefix = 'DSSP'
    names = 'H,E,C'.split(',')
    expected_filenames = cta_output_files_helper.determine_filenames(prefix, names, ext)
    for protein in [GS6_CP, NTL9_CP]:
        with tempfile.TemporaryDirectory() as outdir:
            cta_aa.run_DSSP_analysis(protein, outdir)
            for expected_filename in expected_filenames:
                savename = os.path.join(outdir, expected_filename)

                assert os.path.getsize(savename) > 0  # implicit check for file existence
                with open(savename) as f:
                    data = np.loadtxt(f, delimiter=',')  # numpy is operating on a buffer b/c of the temporaryfile.
                    d = np.frombuffer(data).reshape(data.shape)  # convert the buffer back to an array, and compare.
                    assert np.allclose(data, d)


def test_run_old_DSSP_analysis(GS6_CP, NTL9_CP, cta_output_files_helper):
    ext = 'csv'
    prefix = 'DSSP'
    names = 'H,E,C'.split(',')
    expected_filenames = cta_output_files_helper.determine_filenames(prefix, names, ext)
    for protein in [GS6_CP, NTL9_CP]:
        with tempfile.TemporaryDirectory() as outdir:
            cta_aa.run_DSSP_analysis_OLD(protein, outdir)
            for expected_filename in expected_filenames:
                savename = os.path.join(outdir, expected_filename)
                
                assert os.path.getsize(savename) > 0  # implicit check for file existence
                with open(savename) as f:
                    data = np.loadtxt(f, delimiter=',')  # numpy is operating on a buffer b/c of the temporaryfile.
                    d = np.frombuffer(data).reshape(data.shape)  # convert the buffer back to an array, and compare.
                    assert np.allclose(data, d)
                    

def test_run_BBSEG_analysis(GS6_CP, NTL9_CP):
    expected_filenames = ['BBSEG_%d.csv' % num for num in range(9)]
    for protein in [GS6_CP, NTL9_CP]:
        with tempfile.TemporaryDirectory() as outdir:
            cta_aa.run_BBSEG_analysis(protein, outdir)
            for expected_filename in expected_filenames:
                savename = os.path.join(outdir, expected_filename)

                assert os.path.getsize(savename) > 0  # implicit check for file existence
                with open(savename) as f:
                    data = np.loadtxt(f, delimiter=',')  # numpy is operating on a buffer b/c of the temporaryfile.
                    d = np.frombuffer(data).reshape(data.shape)  # convert the buffer back to an array, and compare.
                    assert np.allclose(data, d)


# Note: this requires a protein with at least 20 frames.
def test_run_linear_heterogeneity(GS6_CP, NTL9_CP, cta_output_files_helper):
    ext = 'csv'
    prefix = 'linear_heterogeneity'
    names = 'mean,std'.split(',')
    expected_filenames = cta_output_files_helper.determine_filenames(prefix, names, ext)

    trajectories = [NTL9_CP.traj, NTL9_CP.traj, NTL9_CP.traj, NTL9_CP.traj, NTL9_CP.traj]
    traj = NTL9_CP.traj.join(trajectories)
    protein = CTProtein(traj)

    with tempfile.TemporaryDirectory() as outdir:
        cta_aa.run_linear_heterogeneity(protein, outdir)
        for expected_filename in expected_filenames:
            savename = os.path.join(outdir, expected_filename)

            assert os.path.getsize(savename) > 0  # implicit check for file existence
            with open(savename) as f:
                data = np.loadtxt(f, delimiter=',')  # numpy is operating on a buffer b/c of the temporaryfile.
                d = np.frombuffer(data).reshape(data.shape)  # convert the buffer back to an array, and compare.
                assert np.allclose(data, d)


def test_run_heterogeneity_analysis(GS6_CP, NTL9_CP, cta_output_files_helper):
    initial_stride_val = 20
    ext = 'csv'
    prefix = 'D'
    names = 'vector,mean,std'.split(',')
    expected_filenames = cta_output_files_helper.determine_filenames(prefix, names, ext)

    for protein_traj in [GS6_CP, NTL9_CP]:
        trajectories = [protein_traj.traj, protein_traj.traj, protein_traj.traj, protein_traj.traj, protein_traj.traj]
        traj = protein_traj.traj.join(trajectories)
        protein = CTProtein(traj)

        with tempfile.TemporaryDirectory() as outdir:
            for strideval in range(initial_stride_val, protein.n_frames - 1):
                cta_aa.run_heterogeneity_analysis(protein, strideval, outdir)
                for expected_filename in expected_filenames:
                    savename = os.path.join(outdir, expected_filename)

                    assert os.path.getsize(savename) > 0  # implicit check for file existence
                    with open(savename) as f:
                        data = np.loadtxt(f, delimiter=',')  # numpy is operating on a buffer b/c of the temporaryfile.
                        d = np.frombuffer(data).reshape(data.shape)  # convert the buffer back to an array, and compare.
                        assert np.allclose(data, d)


'''
def test_run_cluster_analysis(GS6_CP, NTL9_CP):
    initial_stride_val = 20
    expected_filenames = 'cluster_size.csv,cluster_centroid_traj.xtc,cluster_centroid_traj.pdb'.split(',')
    for protein_traj in [NTL9_CP]:
        trajectories = [protein_traj.traj, protein_traj.traj, protein_traj.traj, protein_traj.traj, protein_traj.traj]
        traj = protein_traj.traj.join(trajectories)
        protein = CTProtein(traj)

        with tempfile.TemporaryDirectory() as outdir:
            for strideval in range(initial_stride_val, protein.n_frames - 1):
                cta_aa.run_cluster_analysis(protein, strideval, outdir)
                for expected_filename in expected_filenames:
                    savename = os.path.join(outdir, expected_filename)

                    assert os.path.getsize(savename) > 0  # implicit check for file existence
                    with open(savename) as f:
                        data = np.loadtxt(f, delimiter=',')  # numpy is operating on a buffer b/c of the temporaryfile.
                        d = np.frombuffer(data).reshape(data.shape)  # convert the buffer back to an array, and compare.
                        assert np.allclose(data, d)
'''


def test_run_dihedral_extraction(GS6_CP, NTL9_CP):
    expected_filenames = 'PHI_matrix.csv,PSI_matrix.csv,OMEGA_matrix.csv'.split(',')
    for protein in [GS6_CP, NTL9_CP]:
        with tempfile.TemporaryDirectory() as outdir:
            cta_aa.run_dihedral_extraction(protein, outdir)
            for expected_filename in expected_filenames:
                savename = os.path.join(outdir, expected_filename)

                assert os.path.getsize(savename) > 0  # implicit check for file existence
                with open(savename) as f:
                    data = np.loadtxt(f, delimiter=',')  # numpy is operating on a buffer b/c of the temporaryfile.
                    d = np.frombuffer(data).reshape(data.shape)  # convert the buffer back to an array, and compare.
                    assert np.allclose(data, d)


def test_run_angle_mutual_information(GS6_CP, NTL9_CP):
    base_filename = '%s_mutual_information.csv'
    angles = 'chi1,phi,psi,omega'.split(',')
    
    for protein in [GS6_CP, NTL9_CP]:
        with tempfile.TemporaryDirectory() as outdir:
            for angle in angles:
                expected_filename = os.path.join(outdir, base_filename % angle)
                cta_aa.run_angle_mutual_information(protein, outdir, angle)
                savename = os.path.join(outdir, expected_filename)

                assert os.path.getsize(savename) > 0  # implicit check for file existence
                with open(savename) as f:
                    data = np.loadtxt(f, delimiter=',')  # numpy is operating on a buffer b/c of the temporaryfile.
                    d = np.frombuffer(data).reshape(data.shape)  # convert the buffer back to an array, and compare.
                    assert np.allclose(data, d)
