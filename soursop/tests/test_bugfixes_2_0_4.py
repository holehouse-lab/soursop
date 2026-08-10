"""
Regression tests for the correctness fixes made in SOURSOP 2.0.4.

Currently this covers ``get_secondary_structure_BBSEG`` on uncapped chains.
mdtraj returns one phi for every residue with a preceding residue and one psi
for every residue with a following one. On a *capped* chain the ACE/NME caps
supply both, so the two arrays line up. On an *uncapped* chain they do not:
phi is undefined at the N-terminus and psi at the C-terminus, so ``phi[k]`` and
``psi[k]`` describe two different (adjacent) residues. The old implementation
paired them positionally and separately built the residue list from the full
CA-bearing range, which meant that on an uncapped chain:

  * each classified (phi, psi) pair was drawn from two adjacent residues, and
  * ``resid_list`` came back one element longer than the class arrays.

Capped chains were - and must remain - unaffected.
"""

import os

import mdtraj as md
import numpy as np
import pytest

import soursop
from soursop import sstrajectory
from soursop.ssexceptions import SSException


test_data_dir = soursop.get_data("test_data")

# (name, is_capped) for the bundled all-atom trajectories
ALL_ATOM_SYSTEMS = [
    ("ntl9_AA", False),
    ("ctl9_AA", True),
    ("gs6_AA", True),
    ("all_residues_AA", True),
]


def _load(name):
    return sstrajectory.SSTrajectory(
        os.path.join(test_data_dir, f"{name}.xtc"),
        os.path.join(test_data_dir, f"{name}.pdb"),
    ).proteinTrajectoryList[0]


@pytest.fixture(scope="module")
def NTL9():
    """Uncapped 56-residue chain - the case the fix targets."""
    return _load("ntl9_AA")


@pytest.fixture(scope="module")
def CTL9():
    """ACE/NME-capped 94-residue chain - must be unchanged by the fix."""
    return _load("ctl9_AA")


# ........................................................................
# the returned resid list must match the classification arrays
#


@pytest.mark.parametrize("name,capped", ALL_ATOM_SYSTEMS)
def test_bbseg_resid_list_matches_class_arrays(name, capped):
    protein = _load(name)
    resids, per_class = protein.get_secondary_structure_BBSEG()

    assert len(resids) > 0
    for c in range(9):
        assert len(per_class[c]) == len(resids), (
            f"{name}: class {c} array has {len(per_class[c])} entries but "
            f"resid_list has {len(resids)}"
        )


@pytest.mark.parametrize("name,capped", ALL_ATOM_SYSTEMS)
def test_bbseg_resid_list_matches_per_frame_arrays(name, capped):
    protein = _load(name)
    resids, per_class = protein.get_secondary_structure_BBSEG(return_per_frame=True)

    for c in range(9):
        assert per_class[c].shape == (protein.n_frames, len(resids))


@pytest.mark.parametrize("name,capped", ALL_ATOM_SYSTEMS)
def test_bbseg_class_fractions_partition(name, capped):
    """The nine classes are exhaustive, so the fractions sum to 1 everywhere."""
    protein = _load(name)
    resids, per_class = protein.get_secondary_structure_BBSEG()

    total = sum(per_class[c] for c in range(9))
    assert np.allclose(total, 1.0)


# ........................................................................
# each classified (phi, psi) pair must come from a single residue
#


@pytest.mark.parametrize("name,capped", ALL_ATOM_SYSTEMS)
def test_bbseg_phi_psi_pairs_are_same_residue(name, capped):
    """The core of the bug: verify phi and psi are drawn from one residue.

    Recomputes the dihedrals independently, maps each back to its central
    residue via the atom indices mdtraj returns, and checks the residues
    SOURSOP reports are exactly those with *both* defined.
    """
    protein = _load(name)
    resids, _ = protein.get_secondary_structure_BBSEG()

    # reproduce the selection the method makes internally
    R1, R2, selection = protein._SSProtein__get_first_and_last(None, None, withCA=False)
    target = protein.traj.atom_slice(protein.topology.select(selection))

    phi_atoms, _ = md.compute_phi(target)
    psi_atoms, _ = md.compute_psi(target)

    # phi = [C(i-1), N(i), CA(i), C(i)]; psi = [N(i), CA(i), C(i), N(i+1)]
    # so atom index 1 is on the central residue in both cases
    phi_resids = {target.topology.atom(a[1]).residue.index for a in phi_atoms}
    psi_resids = {target.topology.atom(a[1]).residue.index for a in psi_atoms}

    expected = sorted(r + R1 for r in (phi_resids & psi_resids))
    assert resids == expected


def test_bbseg_uncapped_drops_both_termini(NTL9):
    """An uncapped chain loses the first and last residue, not zero or one."""
    assert NTL9.ncap is False
    assert NTL9.ccap is False

    resids, per_class = NTL9.get_secondary_structure_BBSEG()

    assert len(resids) == NTL9.n_residues - 2
    assert resids[0] == 1
    assert resids[-1] == NTL9.n_residues - 2
    assert len(per_class[1]) == len(resids)


def test_bbseg_capped_covers_every_real_residue(CTL9):
    """A capped chain classifies every non-cap residue."""
    assert CTL9.ncap is True
    assert CTL9.ccap is True

    resids, per_class = CTL9.get_secondary_structure_BBSEG()

    # 94 residues, minus the two caps
    assert len(resids) == CTL9.n_residues - 2
    assert resids == list(range(1, CTL9.n_residues - 1))
    assert len(per_class[1]) == len(resids)


# ........................................................................
# capped-chain values must not have moved
#


def test_bbseg_capped_values_unchanged(CTL9):
    """Pin the capped-chain output, which the fix leaves byte-identical."""
    resids, per_class = CTL9.get_secondary_structure_BBSEG()

    assert len(resids) == 92
    assert resids[0] == 1 and resids[-1] == 92

    assert np.isclose(np.sum(per_class[1]), 17.454)
    assert np.isclose(np.sum(per_class[2]), 22.849)
    assert np.isclose(np.sum(per_class[4]), 7.406)
    assert np.isclose(np.sum(per_class[8]), 1.164)


def test_bbseg_capped_matches_positional_pairing(CTL9):
    """On a capped chain, positional pairing was already correct.

    This is the explicit statement that the fix is a no-op for capped chains:
    the old (positional) implementation is reproduced here and must agree.
    """
    R1, R2, selection = CTL9._SSProtein__get_first_and_last(None, None, withCA=False)
    selector = CTL9._SSProtein__get_first_and_last(None, None, withCA=True)

    target = CTL9.traj.atom_slice(CTL9.topology.select(selection))
    phi = np.degrees(md.compute_phi(target)[1])
    psi = np.degrees(md.compute_psi(target)[1])

    old_classes = np.array(
        [CTL9._SSProtein__phi_psi_bbseg(phi[f], psi[f]) for f in range(phi.shape[0])]
    )
    old_resids = list(range(selector[0], selector[1] + 1))

    resids, per_class = CTL9.get_secondary_structure_BBSEG()

    assert resids == old_resids
    for c in range(9):
        old = np.sum((old_classes == c) * 1, 0) / phi.shape[0]
        assert np.allclose(per_class[c], old)


# ........................................................................
# region selection and the weights path still behave
#


def test_bbseg_subregion(NTL9):
    """R1/R2 selection still trims to the residues with both dihedrals."""
    resids, per_class = NTL9.get_secondary_structure_BBSEG(R1=10, R2=20)

    assert len(resids) == len(per_class[1])
    # the sub-selection is itself uncapped, so it loses its own two termini
    assert resids[0] == 11
    assert resids[-1] == 19


def test_bbseg_weights_uniform_matches_unweighted(NTL9):
    n = NTL9.n_frames
    uniform = np.repeat(1.0 / n, n)

    resids_u, unweighted = NTL9.get_secondary_structure_BBSEG()
    resids_w, weighted = NTL9.get_secondary_structure_BBSEG(weights=uniform)

    assert resids_u == resids_w
    for c in range(9):
        assert np.allclose(unweighted[c], weighted[c])


def test_bbseg_weights_one_hot_selects_frame(NTL9):
    n = NTL9.n_frames
    resids, per_frame = NTL9.get_secondary_structure_BBSEG(return_per_frame=True)

    for k in [0, n // 2, n - 1]:
        one_hot = np.zeros(n)
        one_hot[k] = 1.0
        _, weighted = NTL9.get_secondary_structure_BBSEG(weights=one_hot)
        for c in range(9):
            assert np.allclose(weighted[c], per_frame[c][k])


def test_bbseg_weights_rejected_with_per_frame(NTL9):
    n = NTL9.n_frames
    with pytest.raises(SSException):
        NTL9.get_secondary_structure_BBSEG(
            return_per_frame=True, weights=np.repeat(1.0 / n, n)
        )


def test_bbseg_too_short_region_raises(NTL9):
    """A region with no residue carrying both dihedrals raises clearly."""
    with pytest.raises(SSException):
        NTL9.get_secondary_structure_BBSEG(R1=10, R2=11)


def test_local_heterogeneity_uses_every_exact_size_window(NTL9, monkeypatch):
    """The sliding window count is n-k+1 and each RMSD region has k residues."""
    regions = []

    def fake_rmsd(frame1, frame2=-1, region=None, backbone=True, stride=1):
        regions.append(tuple(region))
        return np.array([0.0])

    monkeypatch.setattr(NTL9, "get_RMSD", fake_rmsd)
    fragment_size = NTL9.n_residues
    mean, std, histo, _ = NTL9.get_local_heterogeneity(
        fragment_size=fragment_size,
        stride=NTL9.n_frames,
        verbose=False,
    )

    assert len(mean) == len(std) == len(histo) == 1
    assert regions == [(0, NTL9.n_residues - 1)]


# ........................................................................
# fixes from the post-2.0.4-commit review round
#


def test_explicit_residue_checking_skips_all_solvent_chains():
    """A chain with zero protein residues must be skipped, not crash.

    Under ``explicit_residue_checking=True`` an all-solvent chain used to
    contribute an empty atom list; atom_slice([]) built an empty topology and
    the first-residue lookup raised IndexError - defeating the exact use case
    (solvent-heavy topologies) the flag exists for.
    """
    gs6 = md.load(os.path.join(test_data_dir, "gs6_AA.pdb"))

    topology = md.Topology()
    chain = topology.add_chain()
    residue = topology.add_residue("HOH", chain)
    topology.add_atom("O", md.core.element.oxygen, residue)
    water = md.Trajectory(np.zeros((1, 1, 3), dtype=np.float32), topology)

    combo = gs6[0].stack(water)

    traj = sstrajectory.SSTrajectory(TRJ=combo, explicit_residue_checking=True)
    assert traj.n_proteins == 1
    assert traj.proteinTrajectoryList[0].n_residues == gs6.n_residues


def test_chemical_shift_sequence_parsing_is_case_insensitive():
    """Lowercase one-letter codes were silently dropped, corrupting the
    nearest-neighbour context of every surrounding residue."""
    from soursop.ssnmr import compute_random_coil_chemical_shifts

    upper = compute_random_coil_chemical_shifts("ASGAS")
    lower = compute_random_coil_chemical_shifts("asgas")
    mixed = compute_random_coil_chemical_shifts("AsGaS")

    assert len(upper) == len(lower) == len(mixed) == 5
    for res_up, res_low, res_mixed in zip(upper, lower, mixed):
        assert res_up == res_low == res_mixed


def test_get_rmsd_empty_region_raises_not_segfaults(NTL9):
    """A fully out-of-range region selects zero atoms; passing that to
    mdtraj's rmsd() returned garbage and then segfaulted the interpreter."""
    with pytest.raises(SSException, match="selects no atoms"):
        NTL9.get_RMSD(0, region=[NTL9.n_residues + 5, NTL9.n_residues + 10])


def test_get_rmsd_partially_out_of_range_region_still_clips(NTL9):
    # a partially out-of-range region is clipped by mdtraj and remains
    # bitwise identical to the properly bounded call - preserve that
    a = NTL9.get_RMSD(0, region=[0, NTL9.n_residues])
    b = NTL9.get_RMSD(0, region=[0, NTL9.n_residues - 1])
    assert np.array_equal(a, b)


def test_get_rmsd_non_integer_frame2_raises(NTL9):
    """frame2=3.5 used to silently fall into the compare-vs-all-frames
    branch, answering a different question than the caller asked."""
    with pytest.raises(SSException, match="frame2 must be an integer"):
        NTL9.get_RMSD(0, frame2=3.5)


def test_overlap_concentration_uses_exact_avogadro():
    """Pin c* against the closed form with the exact (2019 SI) Avogadro
    number; the truncated 6.023e23 gave a 0.014% systematic error."""
    from soursop.sspolymer import get_overlap_concentration

    rg_angstrom = 25.0
    volume_litres = (4.0 / 3.0) * np.pi * (rg_angstrom * 1e-10) ** 3 * 1000.0
    expected = 1.0 / (volume_litres * 6.02214076e23)

    assert np.isclose(get_overlap_concentration(rg_angstrom), expected, rtol=1e-12)
