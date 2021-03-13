import os
import numpy as np
from camparitraj.configs import TMP_DIR
from . import conftest


def test_export_sampling_goodness_no_bins(NTL9_CO):
    base_filenames = 'mean_S_vector.csv,std_S_vector.csv'.split(',')
    for protein_index in range(len(NTL9_CO.proteinTrajectoryList)):
        temp_save_filepath = os.path.join(TMP_DIR, 'ntl9_{protein}'.format(protein=protein_index))
        mean, std_dev = NTL9_CO.export_samplingGoodness(protein_index, temp_save_filepath,
                                                        fragmentSize=10,
                                                        stride=5,
                                                        bins=None)

        # check that the files were created
        filepaths = list()
        for filename in base_filenames:
            path = os.path.join(TMP_DIR, 'ntl9_{protein}_{fname}'.format(protein=protein_index, fname=filename))
            filepaths.append(path)
            assert os.path.exists(path)

        # check that they are non-zero
        read_mean = np.loadtxt(os.path.join(TMP_DIR, filepaths[0]))
        read_stddev = np.loadtxt(os.path.join(TMP_DIR, filepaths[1]))

        assert np.allclose(mean, read_mean)
        assert np.allclose(std_dev, read_stddev)

        # clean up
        for filepath in filepaths:
            os.remove(filepath)


def test_export_sampling_goodness_with_bins(NTL9_CO):
    base_filenames = 'mean_S_vector.csv,std_S_vector.csv,distribution_S_vector.csv'.split(',')
    for protein_index in range(len(NTL9_CO.proteinTrajectoryList)):
        temp_save_filepath = os.path.join(TMP_DIR, 'ntl9_{protein}'.format(protein=protein_index))
        mean, std_dev, hist1d = NTL9_CO.export_samplingGoodness(protein_index, temp_save_filepath,
                                                                fragmentSize=10,
                                                                stride=5,
                                                                bins=np.arange(0, 1, 0.01))

        # check that the files were created
        filepaths = list()
        for filename in base_filenames:
            path = os.path.join(TMP_DIR, 'ntl9_{protein}_{fname}'.format(protein=protein_index, fname=filename))
            filepaths.append(path)
            assert os.path.exists(path)

        # check that they are non-zero
        read_mean = np.loadtxt(os.path.join(TMP_DIR, filepaths[0]), delimiter=',')
        read_stddev = np.loadtxt(os.path.join(TMP_DIR, filepaths[1]), delimiter=',')
        read_hist1d = np.loadtxt(os.path.join(TMP_DIR, filepaths[2]), delimiter=',')

        assert np.allclose(mean, read_mean)
        assert np.allclose(std_dev, read_stddev)
        assert np.allclose(hist1d, read_hist1d)

        # clean up
        for filepath in filepaths:
            os.remove(filepath)

