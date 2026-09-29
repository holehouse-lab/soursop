import numpy as np
from soursop import ssmutualinformation


def test_shan_entropy():
    assert 0.693 == np.round(ssmutualinformation.shan_entropy(np.repeat(1, 2)), 3)
    assert 1.099 == np.round(ssmutualinformation.shan_entropy(np.repeat(1, 3)), 3)
    assert 1.386 == np.round(ssmutualinformation.shan_entropy(np.repeat(1, 4)), 3)


def _two_bin_samples():
    # N(100, 0.1) split by a bin edge placed exactly at the mean, so each
    # variable occupies two bins ~50/50 and H(X) = ln 2. The edge at 100 is
    # deliberate: the previous version used np.arange(0, 300, 2) and only
    # happened to work because 100 fell on an edge.
    rng = np.random.default_rng(0)
    X = rng.normal(100, 0.1, size=10000)
    Y = rng.normal(100, 0.1, size=10000)
    bins = np.array([0.0, 100.0, 300.0])
    return X, Y, bins


def test_calc_MI():
    X, Y, bins = _two_bin_samples()
    # I(X; X) = H(X) = ln 2 for a 50/50 two-bin split
    assert 0.693 == np.round(ssmutualinformation.calc_MI(X, X, bins), 3)
    assert np.isclose(
        ssmutualinformation.calc_MI(X, X, bins),
        ssmutualinformation.shan_entropy(np.histogram(X, bins)[0]),
    )
    # independent samples share (almost) no information
    assert ssmutualinformation.calc_MI(X, Y, bins) < 1e-3


def test_calc_NMI():
    X, Y, bins = _two_bin_samples()
    assert 0.0 == np.round(ssmutualinformation.calc_MI(X, Y, bins, normalize=True), 3)
    assert 1.0 == np.round(ssmutualinformation.calc_MI(X, X, bins, normalize=True), 3)


def test_calc_MI_weights_and_range_checks():
    from soursop.ssexceptions import SSException
    import pytest

    X, Y, bins = _two_bin_samples()
    uniform = np.full(len(X), 1.0 / len(X))
    assert np.isclose(
        ssmutualinformation.calc_MI(X, Y, bins, weights=uniform),
        ssmutualinformation.calc_MI(X, Y, bins),
    )
    with pytest.raises(SSException):
        ssmutualinformation.calc_MI(X, Y[:-1], bins)
    with pytest.raises(SSException):
        ssmutualinformation.calc_MI(X, Y, np.array([0.0, 50.0]))
