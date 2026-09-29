"""
Regression tests for the smaller correctness fixes made in SOURSOP 2.0.6.

Covered here:

* ``SSTrajectory.get_interchain_distance(periodic=True)`` returned a list rather
  than an array, which made ``get_interchain_contact_map(mode='atom',
  periodic=True)`` crash.
* ``SSTrajectory(extra_valid_residue_names='XXX')`` silently split the string
  into single characters.
* ``sshdx`` treated the NME cap and the free N-terminal NH3+ as exchangeable
  backbone amides.
* ``ssnmr`` parsed a sequence with an unbalanced parenthesis letter by letter.
* ``ssnmr`` perdeuteration tables were the superseded values that the Poulsen
  web server has commented out.
* The ``calc_MI`` docstring example used bins that did not span the data.
"""

import numpy as np
import pytest

from soursop import sshdx, ssnmr, sstrajectory
from soursop.ssexceptions import SSException
from soursop.ssmutualinformation import calc_MI


# --------------------------------------------------------------------------------
# 1. periodic interchain distance / contact map
# --------------------------------------------------------------------------------
def test_interchain_distance_periodic_returns_ndarray(GMX_2CHAINS):
    d = GMX_2CHAINS.get_interchain_distance(0, 1, 3, 3, periodic=True)
    assert isinstance(d, np.ndarray)
    assert d.dtype == float
    assert d.shape == (GMX_2CHAINS.n_frames,)

    d_np = GMX_2CHAINS.get_interchain_distance(0, 1, 3, 3, periodic=False)
    assert isinstance(d_np, np.ndarray)
    assert d_np.shape == d.shape


def test_interchain_contact_map_atom_periodic(GMX_2CHAINS):
    # rows/columns are the CA-bearing residues (2.0.6)
    n1 = len(GMX_2CHAINS.proteinTrajectoryList[0].resid_with_CA)
    n2 = len(GMX_2CHAINS.proteinTrajectoryList[1].resid_with_CA)

    cmap = GMX_2CHAINS.get_interchain_contact_map(0, 1, mode="atom", periodic=True)
    assert isinstance(cmap, np.ndarray)
    assert cmap.shape == (n1, n2)
    assert np.all(cmap >= 0.0) and np.all(cmap <= 1.0)


# --------------------------------------------------------------------------------
# 2. extra_valid_residue_names
# --------------------------------------------------------------------------------
def _gmx_files():
    import os

    import soursop

    d = os.path.join(soursop.get_data("test_data"), "gromacs2chains")
    return os.path.join(d, "traj.xtc"), os.path.join(d, "top.pdb")


def test_extra_valid_residue_names_accepts_string():
    xtc, pdb = _gmx_files()
    T = sstrajectory.SSTrajectory(xtc, pdb, extra_valid_residue_names="XXX")
    assert "XXX" in T.valid_residue_names
    # the string must not have been split into characters
    for ch in "XXX":
        assert ch not in T.valid_residue_names


def test_extra_valid_residue_names_accepts_list():
    xtc, pdb = _gmx_files()
    T = sstrajectory.SSTrajectory(xtc, pdb, extra_valid_residue_names=["XXX", "YYY"])
    assert "XXX" in T.valid_residue_names
    assert "YYY" in T.valid_residue_names


@pytest.mark.parametrize("bad", [5, 3.2, ["ALA", 7]])
def test_extra_valid_residue_names_rejects_bad_input(bad):
    xtc, pdb = _gmx_files()
    with pytest.raises(SSException):
        sstrajectory.SSTrajectory(xtc, pdb, extra_valid_residue_names=bad)


# --------------------------------------------------------------------------------
# 3. sshdx: caps and the free N-terminus
# --------------------------------------------------------------------------------
def test_hdx_capped_excludes_nme(GS6_CP):
    # GS6 is ACE-GSGSGS-NME; residue 0 is ACE, residue 7 is NME
    # the pinned values were recorded with the backbone-only Wernet-Nilsson
    # H-bond definition, which is no longer the default
    res, lnP = sshdx.compute_protection_factors(
        GS6_CP,
        weights=np.ones(GS6_CP.n_frames) / GS6_CP.n_frames,
        hbond_method="wernet-nilsson",
    )
    assert list(res) == [1, 2, 3, 4, 5, 6]
    assert 7 not in res
    # residue 1 follows an ACE cap so is a genuine amide and is retained;
    # values pinned from the pre-fix code for the retained residues
    np.testing.assert_allclose(lnP, [0.07, 0.0, 0.21, 0.42, 0.28, 0.07], atol=1e-6)
    # the default (Best-Vendruscolo) definition reports the same residues
    res_default, _ = sshdx.compute_protection_factors(GS6_CP)
    assert list(res_default) == list(res)


def test_hdx_uncapped_excludes_free_n_terminus(NTL9_CP):
    res, lnP = sshdx.compute_protection_factors(
        NTL9_CP,
        weights=np.ones(NTL9_CP.n_frames) / NTL9_CP.n_frames,
        hbond_method="wernet-nilsson",
    )
    assert 0 not in res
    assert res[0] == 1
    # proline (residue 40) still excluded
    assert 40 not in res
    assert len(res) == NTL9_CP.n_residues - 2

    # values for the retained residues are unchanged (pinned pre-fix, with the
    # Wernet-Nilsson H-bond definition those values were recorded with)
    np.testing.assert_allclose(lnP[:4], [0.07, 0.805, 1.82, 1.68], atol=1e-6)
    np.testing.assert_allclose(lnP[-3:], [0.595, 0.735, 0.07], atol=1e-6)
    res_default, _ = sshdx.compute_protection_factors(NTL9_CP)
    assert list(res_default) == list(res)


def test_hdx_map_skips_cap_and_free_n_terminus(GS6_CP, NTL9_CP):
    res, _, _ = sshdx._backbone_nh_map(GS6_CP.traj.topology)
    names = [GS6_CP.traj.topology.residue(int(i)).name for i in res]
    assert "NME" not in names and "ACE" not in names

    res, _, _ = sshdx._backbone_nh_map(NTL9_CP.traj.topology)
    assert 0 not in res


# --------------------------------------------------------------------------------
# 4. ssnmr: unbalanced parentheses
# --------------------------------------------------------------------------------
@pytest.mark.parametrize("seq", ["A(SEP", "ASEP)A", "A(SEP))A", "A((SEP)A"])
def test_nmr_unbalanced_parentheses_raise(seq):
    with pytest.raises(SSException):
        ssnmr.compute_random_coil_chemical_shifts(seq)


def test_nmr_balanced_parentheses_still_work():
    out = ssnmr.compute_random_coil_chemical_shifts("A(SEP)A")
    assert [o["Res"] for o in out] == ["A", "SEP", "A"]


# --------------------------------------------------------------------------------
# 5. ssnmr: perdeuteration tables match the live Poulsen server
# --------------------------------------------------------------------------------
def test_nmr_perdeuteration_matches_server():
    out = ssnmr.compute_random_coil_chemical_shifts(
        "ACDEFGHIKLMNPQRSTVWY", temperature=5, pH=6.5, use_perdeuteration=True
    )
    ala = out[0]
    assert ala["Res"] == "A"
    # server value for Ala CA with the current (Maltsev 2012) ca_deut table
    assert ala["CA"] == pytest.approx(52.247, abs=1e-3)

    # the deuteration shift is the difference vs the protonated prediction
    prot = ssnmr.compute_random_coil_chemical_shifts(
        "ACDEFGHIKLMNPQRSTVWY", temperature=5, pH=6.5, use_perdeuteration=False
    )
    assert ala["CA"] - prot[0]["CA"] == pytest.approx(-0.47, abs=1e-3)
    assert ala["CB"] - prot[0]["CB"] == pytest.approx(-0.88, abs=1e-3)
    # cysteine (second residue) and tyrosine (last)
    assert out[1]["CA"] - prot[1]["CA"] == pytest.approx(-0.45, abs=1e-3)
    assert out[1]["CB"] - prot[1]["CB"] == pytest.approx(-0.71, abs=1e-3)
    assert out[-1]["CA"] - prot[-1]["CA"] == pytest.approx(-0.43, abs=1e-3)
    assert out[-1]["CB"] - prot[-1]["CB"] == pytest.approx(-0.86, abs=1e-3)


# --------------------------------------------------------------------------------
# 6. calc_MI docstring example runs
# --------------------------------------------------------------------------------
def test_calc_mi_docstring_example_runs():
    rng = np.random.default_rng(0)
    X = rng.uniform(-1, 1, 1000)
    Y = X + 0.05 * rng.standard_normal(1000)
    bins = np.linspace(-1.5, 1.5, 31)
    mi = calc_MI(X, Y, bins)
    assert mi > 1.0
    mi_ind = calc_MI(X, rng.uniform(-1, 1, 1000), bins)
    assert mi_ind < 0.2


# --------------------------------------------------------------------------------
# 7. interchain distance map (i, i) vs SSProtein.get_distance_map
# --------------------------------------------------------------------------------
def test_interchain_distance_map_self_matches_upper_triangle(GMX_2CHAINS):
    full, _ = GMX_2CHAINS.get_interchain_distance_map(0, 0)
    tri, _ = GMX_2CHAINS.proteinTrajectoryList[0].get_distance_map()
    assert np.allclose(full, full.T)
    iu = np.triu_indices_from(full)
    np.testing.assert_allclose(full[iu], tri[iu], atol=1e-5)
