"""
Tests for the one-bead-per-residue coarse-grained SASA path.

Coarse-grained (CG) models write every residue out as a single carbon ``CA``
bead, so mdtraj's Shrake-Rupley would give a glycine and a tryptophan bead the
same 1.7 Angstrom radius. SOURSOP instead computes CG SASA against the force
field's own bead sizes (radius = sigma/2, where sigma is the model's diagonal
pair sigma, i.e. the bead diameter).

The mechanism matters here: mdtraj's ``change_radii`` is keyed on
``atom.element.symbol`` and a name-keyed dictionary is *silently ignored*
rather than rejected. ``test_isolated_bead_areas_are_analytic`` is the test
that would catch that failure mode - it pins the computed areas against the
closed-form isolated-sphere result.
"""

import os

import mdtraj as md
import numpy as np
import pytest

import soursop
from soursop import ssdata, sstrajectory
from soursop.ssdata import CG_FORCEFIELDS, get_cg_bead_sigmas
from soursop.ssexceptions import SSException
from soursop.ssprotein import SSProtein


test_data_dir = soursop.get_data("test_data")

# the force fields whose sigma tables are bundled with SOURSOP
EXPECTED_FORCEFIELDS = (
    "cation-pi-1",
    "cation-pi-2",
    "fb-hps",
    "hps-kr",
    "hps-urry",
    "kh",
    "mpipi",
    "mpipi-gg",
    "mpipi-recharged",
)

THE_20 = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
]  # fmt: skip

RNA_BEADS = ["RPA", "RPC", "RPG", "RPU"]


# ........................................................................
# fixtures
#
# NOTE these are function-scoped and load their own trajectory rather than
# reusing the session-scoped conftest fixtures, because several tests set the
# `cg_forcefield` property and that state must not leak between tests.


@pytest.fixture
def CG_protein():
    """A fresh one-bead-per-residue coarse-grained SSProtein (sigA, 204 beads)."""
    traj = sstrajectory.SSTrajectory(
        os.path.join(test_data_dir, "sigA_CG.xtc"),
        os.path.join(test_data_dir, "sigA_CG.pdb"),
    )
    return traj.proteinTrajectoryList[0]


@pytest.fixture
def AA_protein():
    """A fresh all-atom SSProtein (GS6)."""
    traj = sstrajectory.SSTrajectory(
        os.path.join(test_data_dir, "gs6_AA.xtc"),
        os.path.join(test_data_dir, "gs6_AA.pdb"),
    )
    return traj.proteinTrajectoryList[0]


def _isolated_bead_protein(residue_names, separation_nm=10.0):
    """Build an in-memory CG chain of beads far enough apart not to occlude.

    Parameters
    ----------
    residue_names : list of str
        Three-letter residue names, one bead each.
    separation_nm : float, optional
        Spacing along x in nanometers. Default 10.0 (100 Angstroms), which is
        an order of magnitude beyond any bead radius + probe.

    Returns
    -------
    SSProtein
    """

    topology = md.Topology()
    chain = topology.add_chain()
    for idx, name in enumerate(residue_names):
        residue = topology.add_residue(name, chain, resSeq=idx)
        topology.add_atom("CA", md.core.element.carbon, residue)

    xyz = np.zeros((1, len(residue_names), 3), dtype=np.float32)
    xyz[0, :, 0] = np.arange(len(residue_names)) * separation_nm

    return SSProtein(md.Trajectory(xyz, topology))


# ........................................................................
# the bundled sigma tables (ssdata)
#


def test_forcefield_registry():
    assert set(CG_FORCEFIELDS) == set(EXPECTED_FORCEFIELDS)

    # the registry is derived from the data directory, so every slug must load
    for ff in CG_FORCEFIELDS:
        assert os.path.isfile(os.path.join(ssdata.CG_SIGMA_DIR, f"{ff}-sigma.csv"))


def test_every_forcefield_covers_the_twenty_amino_acids():
    for ff in CG_FORCEFIELDS:
        sigmas = get_cg_bead_sigmas(ff)
        missing = [r for r in THE_20 if r not in sigmas]
        assert missing == [], f"{ff} is missing {missing}"

        # sigmas are bead diameters in Angstroms - sanity-bound them
        for resname in THE_20:
            assert 3.0 < sigmas[resname] < 12.0


def test_mpipi_variants_carry_rna_beads():
    for ff in ["mpipi", "mpipi-gg"]:
        sigmas = get_cg_bead_sigmas(ff)
        assert all(r in sigmas for r in RNA_BEADS)

    # mpipi-recharged parameterises poly-U only
    assert "RPU" in get_cg_bead_sigmas("mpipi-recharged")

    # the HPS-family models are protein-only
    for ff in ["hps-urry", "kh", "fb-hps"]:
        sigmas = get_cg_bead_sigmas(ff)
        assert not any(r in sigmas for r in RNA_BEADS)


def test_buried_sigmas():
    exposed = get_cg_bead_sigmas("mpipi")
    buried = get_cg_bead_sigmas("mpipi", buried=True)

    # RNA has no buried variant, so those rows are dropped
    assert not any(r in buried for r in RNA_BEADS)
    assert all(r in buried for r in THE_20)

    # buried and exposed sigmas are close but not identical for mpipi
    assert buried != exposed
    for resname in THE_20:
        assert np.isclose(buried[resname], exposed[resname], rtol=0.05)


def test_hps_family_share_a_sigma_table():
    # every HPS-family model uses the same Kim-Hummer van der Waals diameters
    reference = get_cg_bead_sigmas("hps-urry")
    for ff in ["kh", "hps-kr", "fb-hps", "cation-pi-1", "cation-pi-2"]:
        assert get_cg_bead_sigmas(ff) == reference

    # ... and mpipi does not
    assert get_cg_bead_sigmas("mpipi") != reference


def test_unknown_forcefield_raises():
    with pytest.raises(SSException):
        get_cg_bead_sigmas("not-a-forcefield")


def test_sigma_lookup_is_memoised():
    # the loader caches by (forcefield, buried), so repeated calls return the
    # same object rather than re-reading the CSV
    assert get_cg_bead_sigmas("mpipi") is get_cg_bead_sigmas("mpipi")


# ........................................................................
# resolution detection and the force-field guard rails
#


def test_is_coarse_grained(CG_protein, AA_protein):
    assert CG_protein.is_coarse_grained is True
    assert AA_protein.is_coarse_grained is False


def test_cg_without_forcefield_raises(CG_protein):
    with pytest.raises(SSException):
        CG_protein.get_all_SASA(stride=1)

    with pytest.raises(SSException):
        CG_protein.get_regional_SASA(1, 5, stride=1)

    with pytest.raises(SSException):
        CG_protein.get_site_accessibility(["TRP"], stride=1)


def test_all_atom_with_forcefield_raises(AA_protein):
    # a bead-size table is meaningless for an all-atom chain; silently ignoring
    # the argument would be the worst outcome
    with pytest.raises(SSException):
        AA_protein.get_all_SASA(stride=1, forcefield="mpipi")

    AA_protein.cg_forcefield = "mpipi"
    with pytest.raises(SSException):
        AA_protein.get_all_SASA(stride=1)


def test_invalid_forcefield_raises(CG_protein):
    with pytest.raises(SSException):
        CG_protein.get_all_SASA(stride=1, forcefield="not-a-forcefield")

    with pytest.raises(SSException):
        CG_protein.cg_forcefield = "not-a-forcefield"


def test_decomposed_modes_raise_for_cg(CG_protein):
    # one bead is the whole residue, so there is no backbone/sidechain split
    for mode in ["sidechain", "backbone", "all"]:
        with pytest.raises(SSException):
            CG_protein.get_all_SASA(stride=1, mode=mode, forcefield="mpipi")


# ........................................................................
# the coarse-grained calculation itself
#


def test_cg_sasa_shape_and_values(CG_protein):
    sasa = CG_protein.get_all_SASA(stride=1, forcefield="mpipi")

    assert sasa.shape == (CG_protein.n_frames, CG_protein.n_residues)
    assert np.all(np.isfinite(sasa))
    assert np.all(sasa >= 0)

    # no bead in a disordered CG chain should be fully buried
    assert np.all(sasa.mean(axis=0) > 0)


def test_isolated_bead_areas_are_analytic():
    """An isolated bead's SASA must be 4*pi*(sigma/2 + probe)^2.

    This is the test that pins the whole mechanism: if the per-residue radii
    never reach mdtraj (e.g. because change_radii was keyed on atom name rather
    than element), every bead silently falls back to the carbon radius and all
    three residues below return the same wrong number.
    """

    residue_names = ["GLY", "ALA", "TRP"]
    protein = _isolated_bead_protein(residue_names)
    sigmas = get_cg_bead_sigmas("mpipi")

    probe_radius = 1.4
    sasa = protein.get_all_SASA(
        stride=1, probe_radius=probe_radius, forcefield="mpipi"
    )[0]

    for resname, observed in zip(residue_names, sasa):
        expected = 4 * np.pi * (sigmas[resname] / 2.0 + probe_radius) ** 2
        assert np.isclose(observed, expected, rtol=0.01)

    # and they are genuinely distinct from one another
    assert len(set(np.round(sasa, 3))) == len(residue_names)


def test_isolated_bead_areas_track_bead_size():
    residue_names = ["GLY", "ALA", "TRP"]
    protein = _isolated_bead_protein(residue_names)
    sigmas = get_cg_bead_sigmas("mpipi")

    sasa = protein.get_all_SASA(stride=1, forcefield="mpipi")[0]

    # sigma ordering (GLY < ALA < TRP) must be reflected in the areas
    assert sigmas["GLY"] < sigmas["ALA"] < sigmas["TRP"]
    assert sasa[0] < sasa[1] < sasa[2]


def test_carbon_default_would_have_been_wrong():
    """Confirm the CG path actually changes the answer.

    Every bead is a carbon CA, so the old behaviour gave all three residues the
    same 1.7 Angstrom radius.
    """

    protein = _isolated_bead_protein(["GLY", "ALA", "TRP"])

    carbon_area = 4 * np.pi * (1.7 + 1.4) ** 2
    sasa = protein.get_all_SASA(stride=1, forcefield="mpipi")[0]

    assert not np.any(np.isclose(sasa, carbon_area, rtol=0.01))


def test_forcefield_changes_the_answer(CG_protein):
    mpipi = CG_protein.get_all_SASA(stride=1, forcefield="mpipi")
    hps_urry = CG_protein.get_all_SASA(stride=1, forcefield="hps-urry")
    kh = CG_protein.get_all_SASA(stride=1, forcefield="kh")

    # different bead sizes -> different SASA
    assert not np.allclose(mpipi, hps_urry)

    # shared bead sizes -> identical SASA (also proves the memo key is not
    # collapsing distinct force fields onto one another by accident)
    assert np.allclose(hps_urry, kh)


def test_atom_and_residue_modes_agree_for_cg(CG_protein):
    # one bead per residue, so per-atom and per-residue SASA are the same thing
    per_residue = CG_protein.get_all_SASA(stride=1, mode="residue", forcefield="mpipi")
    per_atom = CG_protein.get_all_SASA(stride=1, mode="atom", forcefield="mpipi")

    assert np.allclose(per_residue, per_atom)


def test_memoisation_is_keyed_on_forcefield(CG_protein):
    first = CG_protein.get_all_SASA(stride=1, forcefield="mpipi")
    second = CG_protein.get_all_SASA(stride=1, forcefield="hps-urry")
    third = CG_protein.get_all_SASA(stride=1, forcefield="mpipi")

    assert not np.allclose(first, second)
    assert np.allclose(first, third)


def test_underlying_topology_is_not_mutated(CG_protein):
    CG_protein.get_all_SASA(stride=1, forcefield="mpipi")

    # the throwaway topology handed to mdtraj must not leak back to the caller
    assert {a.element.symbol for a in CG_protein.topology.atoms} == {"C"}
    assert {a.name for a in CG_protein.topology.atoms} == {"CA"}


def test_unknown_residue_type_raises():
    # 'XXX' has no bead size in any bundled table
    protein = _isolated_bead_protein(["GLY", "XXX"])

    with pytest.raises(SSException):
        protein.get_all_SASA(stride=1, forcefield="mpipi")


def test_rna_beads_are_supported():
    protein = _isolated_bead_protein(["RPA", "RPG"])
    sigmas = get_cg_bead_sigmas("mpipi")

    sasa = protein.get_all_SASA(stride=1, forcefield="mpipi")[0]
    for resname, observed in zip(["RPA", "RPG"], sasa):
        expected = 4 * np.pi * (sigmas[resname] / 2.0 + 1.4) ** 2
        assert np.isclose(observed, expected, rtol=0.01)

    # ... but not by a protein-only model
    with pytest.raises(SSException):
        protein.get_all_SASA(stride=1, forcefield="hps-urry")


# ........................................................................
# the cg_forcefield property
#


def test_cg_forcefield_property(CG_protein):
    assert CG_protein.cg_forcefield is None

    CG_protein.cg_forcefield = "mpipi-gg"
    assert CG_protein.cg_forcefield == "mpipi-gg"

    # a bare call now works and matches the explicit call
    assert np.allclose(
        CG_protein.get_all_SASA(stride=1),
        CG_protein.get_all_SASA(stride=1, forcefield="mpipi-gg"),
    )

    # ... and a per-call argument overrides it
    assert np.allclose(
        CG_protein.get_all_SASA(stride=1, forcefield="hps-urry"),
        CG_protein.get_all_SASA(stride=1, forcefield="kh"),
    )
    assert not np.allclose(
        CG_protein.get_all_SASA(stride=1, forcefield="hps-urry"),
        CG_protein.get_all_SASA(stride=1),
    )

    CG_protein.cg_forcefield = None
    assert CG_protein.cg_forcefield is None
    with pytest.raises(SSException):
        CG_protein.get_all_SASA(stride=1)


def test_cg_forcefield_survives_reset_cache(CG_protein):
    # reset_cache() discards memoised values, but the model choice is a user
    # setting rather than a cached result
    CG_protein.cg_forcefield = "mpipi"
    reference = CG_protein.get_all_SASA(stride=1)

    CG_protein.reset_cache()

    assert CG_protein.cg_forcefield == "mpipi"
    assert np.allclose(CG_protein.get_all_SASA(stride=1), reference)


# ........................................................................
# downstream SASA summaries
#


def test_regional_and_site_accessibility_on_cg(CG_protein):
    CG_protein.cg_forcefield = "mpipi"

    per_residue = CG_protein.get_all_SASA(stride=1)

    regional = CG_protein.get_regional_SASA(1, 5, stride=1)
    assert np.isclose(regional, per_residue[:, 1:5].mean(axis=0).sum(), rtol=1e-5)

    site = CG_protein.get_site_accessibility(["TRP"], stride=1)
    assert len(site) > 0
    for key, (mean_sasa, std_sasa) in site.items():
        assert key.startswith("TRP-")
        assert mean_sasa > 0
        assert std_sasa >= 0


def test_regional_sasa_forcefield_kwarg(CG_protein):
    mpipi = CG_protein.get_regional_SASA(1, 5, stride=1, forcefield="mpipi")
    hps = CG_protein.get_regional_SASA(1, 5, stride=1, forcefield="hps-urry")
    assert not np.isclose(mpipi, hps)


def test_weights_path_on_cg(CG_protein):
    n_frames = CG_protein.n_frames
    uniform = np.repeat(1.0 / n_frames, n_frames)

    per_frame = CG_protein.get_all_SASA(stride=1, forcefield="mpipi")
    weighted = CG_protein.get_all_SASA(stride=1, forcefield="mpipi", weights=uniform)

    assert weighted.shape == (CG_protein.n_residues,)
    assert np.allclose(weighted, per_frame.mean(axis=0), atol=1e-3)


# ........................................................................
# all-atom behaviour must be untouched
#


def test_all_atom_sasa_unchanged(AA_protein):
    # these are the values pinned in test_ssproteins.py::test_get_all_SASA and
    # must not move as a result of the coarse-grained path
    assert np.isclose(np.min(AA_protein.get_all_SASA(stride=1)), 56.75124)
    assert np.isclose(np.max(AA_protein.get_all_SASA(stride=1)), 144.38452)
    assert np.isclose(np.mean(AA_protein.get_all_SASA(stride=1)), 108.992676)

    assert np.isclose(
        np.sum(AA_protein.get_all_SASA(stride=1, mode="backbone")), 1894.9926
    )
    assert np.isclose(
        np.sum(AA_protein.get_all_SASA(stride=1, mode="sidechain")), 1309.465
    )
