##     _____  ____  _    _ _____   _____  ____  _____
##   / ____|/ __ \| |  | |  __ \ / ____|/ __ \|  __ \
##  | (___ | |  | | |  | | |__) | (___ | |  | | |__) |
##   \___ \| |  | | |  | |  _  / \___ \| |  | |  ___/
##   ____) | |__| | |__| | | \ \ ____) | |__| | |
##  |_____/ \____/ \____/|_|  \_\_____/ \____/|_|

## Alex Holehouse (Pappu Lab and Holehouse Lab) and Jared Lalmansing (Pappu lab)
## Simulation analysis package
## Copyright 2014 - 2026
##

"""
sshdx - HDX protection factors via the Best-Vendruscolo model.

Predicts per-residue ln(protection factor) (ln P) for backbone amide
hydrogen-deuterium exchange (HDX) from a structural ensemble, using the
empirical Best-Vendruscolo relation

    ln P_i  =  beta_c * N_c(i)  +  beta_h * N_h(i)  +  beta_0

where, for each residue ``i`` with a backbone amide N-H:

* ``N_c(i)`` is the number of *heavy atoms* in the protein within
  ``contact_cutoff`` (default 6.5 A) of the amide N of residue ``i``,
  ignoring residues whose sequence separation from ``i`` is at most
  ``exclude_neighbours`` (default 2);
* ``N_h(i)`` is the number of *backbone* H-bonds the amide H of residue
  ``i`` forms as donor to a backbone carbonyl O at sequence separation
  ``> exclude_neighbours`` (per-frame H-bonds are detected with
  ``mdtraj.wernet_nilsson``).

Both counts are returned per frame per residue, so the resulting
``(n_frames, n_residues)`` arrays plug directly into the SOURSOP BME /
COPER reweighters. The optional ``weights`` argument collapses the frame
axis to a single per-residue ensemble mean (per the package-wide
:func:`soursop.ssutils.validate_weights` contract).

Proline (no backbone H) and any residue lacking a recognisable backbone
amide H are dropped from the residue list; the function returns the
residue indices it covered alongside the data array.

Public entry points
-------------------
* :func:`compute_Nc` - per-residue heavy-atom contacts (per frame).
* :func:`compute_Nh` - per-residue backbone H-bonds (per frame).
* :func:`compute_protection_factors` - per-residue ln(P) (per frame, or
  optionally collapsed via ``weights``).

References
----------
* Best, R. B. & Vendruscolo, M. *Structural Interpretation of Hydrogen
  Exchange Protection Factors in Proteins.* Structure **14**, 97-106
  (2006). doi:`10.1016/j.str.2005.09.012
  <https://doi.org/10.1016/j.str.2005.09.012>`_.
* Vendruscolo, M., Paci, E., Dobson, C. M. & Karplus, M. *Rare
  Fluctuations of Native Proteins Sampled by Equilibrium Hydrogen
  Exchange.* J. Am. Chem. Soc. **125**, 15686-15687 (2003).

**Author(s):** Alex Holehouse
"""

import mdtraj as md
import numpy as np

from .ssexceptions import SSException
from .ssutils import validate_weights, weighted_mean

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Defaults (Best-Vendruscolo 2006)
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------

#: Heavy-atom contact weight (Best & Vendruscolo, Structure 14, 2006).
DEFAULT_BETA_C = 0.35
#: H-bond weight (Best & Vendruscolo, Structure 14, 2006).
DEFAULT_BETA_H = 2.0
#: Intercept; defaults to 0 in Best-Vendruscolo's published form.
DEFAULT_BETA_0 = 0.0
#: Heavy-atom contact cutoff, in **nanometres** (6.5 A).
DEFAULT_CONTACT_CUTOFF_NM = 0.65
#: Sequence-separation exclusion: residues with ``|i - j| <= this`` are
#: excluded from both N_c and N_h.
DEFAULT_EXCLUDE_NEIGHBOURS = 2

#: Atom-name fallbacks for the backbone amide hydrogen across common
#: force fields (CHARMM/AMBER use ``H`` or ``HN``; some N-terminal
#: parameterisations use ``H1``).
_BACKBONE_H_NAMES = ("H", "HN", "H1")


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
def _backbone_nh_map(topology):
    """Resolve backbone amide N+H atom indices per residue.

    Skips residues with no backbone N (rare hetero residues), no
    recognisable backbone H name (e.g. proline), so the returned arrays
    are aligned 1:1.

    Returns
    -------
    residue_indices : numpy.ndarray (n_res,)
    n_atom_indices  : numpy.ndarray (n_res,)
    h_atom_indices  : numpy.ndarray (n_res,)
    """
    res_idx, n_idx, h_idx = [], [], []
    for residue in topology.residues:
        if residue.name == "PRO":
            continue
        try:
            n_atom = next(a for a in residue.atoms if a.name == "N")
        except StopIteration:
            continue
        h_atom = None
        for cand in _BACKBONE_H_NAMES:
            for a in residue.atoms:
                if a.name == cand:
                    h_atom = a
                    break
            if h_atom is not None:
                break
        if h_atom is None:
            continue
        res_idx.append(residue.index)
        n_idx.append(n_atom.index)
        h_idx.append(h_atom.index)
    return (
        np.asarray(res_idx, dtype=int),
        np.asarray(n_idx, dtype=int),
        np.asarray(h_idx, dtype=int),
    )


def _backbone_carbonyl_o_map(topology):
    """Atom indices of backbone carbonyl O per residue (acceptor)."""
    res_idx, o_idx = [], []
    for residue in topology.residues:
        try:
            o_atom = next(a for a in residue.atoms if a.name == "O")
        except StopIteration:
            continue
        res_idx.append(residue.index)
        o_idx.append(o_atom.index)
    return np.asarray(res_idx, dtype=int), np.asarray(o_idx, dtype=int)


def _heavy_atom_residue_map(topology):
    """Atom indices and residue indices of all non-hydrogen atoms."""
    atom_idx, res_idx = [], []
    for a in topology.atoms:
        if a.element is None or a.element.symbol == "H":
            continue
        atom_idx.append(a.index)
        res_idx.append(a.residue.index)
    return np.asarray(atom_idx, dtype=int), np.asarray(res_idx, dtype=int)


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
def compute_Nc(
    protein,
    contact_cutoff=DEFAULT_CONTACT_CUTOFF_NM,
    exclude_neighbours=DEFAULT_EXCLUDE_NEIGHBOURS,
    stride=1,
):
    """Per-residue per-frame heavy-atom contact count ``N_c(i)``.

    ``N_c(i)`` is the number of protein heavy atoms within
    ``contact_cutoff`` of the backbone amide N of residue ``i``, excluding
    atoms in residues at sequence separation ``|i - j| <=
    exclude_neighbours``. The cutoff is in **nanometres** (matching
    mdtraj's convention).

    Parameters
    ----------
    protein : soursop.ssprotein.SSProtein
    contact_cutoff : float, optional
        Distance cutoff in nm. Default ``0.65`` (= 6.5 A,
        Best-Vendruscolo).
    exclude_neighbours : int, optional
        Sequence-separation exclusion radius. Default ``2``.
    stride : int, optional
        Frame subsample. Default ``1``.

    Returns
    -------
    residue_indices : numpy.ndarray, shape (n_res,)
        Zero-based residue indices for which N_c is defined (the
        backbone-amide residues; proline and any residue lacking a
        recognisable backbone H are dropped).
    Nc : numpy.ndarray, shape (n_frames, n_res), int
        Per-frame, per-residue heavy-atom contact counts.

    Raises
    ------
    SSException
        If the protein has no recognisable backbone amide N+H residues.
    """
    top = protein.traj.topology
    res_NH, n_atom_idx, _ = _backbone_nh_map(top)
    heavy_atom_idx, heavy_res_idx = _heavy_atom_residue_map(top)
    if len(res_NH) == 0:
        raise SSException("compute_Nc: no residues with a backbone amide N+H found")

    traj = protein.traj[::stride] if stride != 1 else protein.traj
    n_frames = traj.n_frames

    Nc = np.zeros((n_frames, len(res_NH)), dtype=int)
    for k, (i_res, n_idx) in enumerate(zip(res_NH, n_atom_idx)):
        mask = np.abs(heavy_res_idx - i_res) > exclude_neighbours
        eligible = heavy_atom_idx[mask]
        if len(eligible) == 0:
            continue
        pairs = np.column_stack([np.full(len(eligible), n_idx, dtype=int), eligible])
        d = md.compute_distances(traj, pairs)  # nm, shape (n_frames, n_pairs)
        Nc[:, k] = np.sum(d < contact_cutoff, axis=1)

    return res_NH, Nc


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
def compute_Nh(
    protein,
    exclude_neighbours=DEFAULT_EXCLUDE_NEIGHBOURS,
    stride=1,
):
    """Per-residue per-frame backbone H-bond count ``N_h(i)``.

    Counts backbone H-bonds in which the **amide H of residue i** is the
    donor and a backbone carbonyl O is the acceptor, at sequence
    separation ``|i - j| > exclude_neighbours``. Per-frame H-bonds are
    detected by :func:`mdtraj.wernet_nilsson` (cone criterion on the
    Donor-H...Acceptor geometry).

    Parameters
    ----------
    protein : soursop.ssprotein.SSProtein
    exclude_neighbours : int, optional
        Sequence-separation exclusion radius. Default ``2``.
    stride : int, optional
        Frame subsample. Default ``1``.

    Returns
    -------
    residue_indices : numpy.ndarray, shape (n_res,)
        Same residue list as :func:`compute_Nc`.
    Nh : numpy.ndarray, shape (n_frames, n_res), int
    """
    top = protein.traj.topology
    res_NH, _, h_atom_idx = _backbone_nh_map(top)
    res_O, o_atom_idx = _backbone_carbonyl_o_map(top)

    h_to_res = dict(zip(h_atom_idx.tolist(), res_NH.tolist()))
    o_to_res = dict(zip(o_atom_idx.tolist(), res_O.tolist()))
    res_to_k = {int(r): k for k, r in enumerate(res_NH)}

    traj = protein.traj[::stride] if stride != 1 else protein.traj
    n_frames = traj.n_frames

    Nh = np.zeros((n_frames, len(res_NH)), dtype=int)
    hbonds_per_frame = md.wernet_nilsson(traj)
    # mdtraj returns a list of length n_frames; each entry is an
    # ``(n_hbonds, 3)`` int array: (donor_heavy_idx, h_idx, acceptor_heavy_idx).
    for f, hb in enumerate(hbonds_per_frame):
        if hb.shape[0] == 0:
            continue
        for h_idx, a_idx in zip(hb[:, 1], hb[:, 2]):
            h_idx_i = int(h_idx)
            a_idx_i = int(a_idx)
            i_res = h_to_res.get(h_idx_i)
            j_res = o_to_res.get(a_idx_i)
            if i_res is None or j_res is None:
                continue
            if abs(i_res - j_res) > exclude_neighbours:
                Nh[f, res_to_k[i_res]] += 1

    return res_NH, Nh


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
def compute_protection_factors(
    protein,
    beta_c=DEFAULT_BETA_C,
    beta_h=DEFAULT_BETA_H,
    beta_0=DEFAULT_BETA_0,
    contact_cutoff=DEFAULT_CONTACT_CUTOFF_NM,
    exclude_neighbours=DEFAULT_EXCLUDE_NEIGHBOURS,
    stride=1,
    weights=False,
    etol=1e-7,
):
    """Per-residue ln(protection factor) via the Best-Vendruscolo formula.

    ``ln P_i = beta_c * N_c(i) + beta_h * N_h(i) + beta_0`` evaluated per
    frame, then optionally collapsed to a single per-residue ensemble
    mean by a per-frame weight vector. The output shape
    ``(n_frames, n_res)`` is the natural input for
    :class:`soursop.ssbme.BME` / :class:`soursop.sscoper.COPER`
    reweighting against experimental HDX protection factors.

    Parameters
    ----------
    protein : soursop.ssprotein.SSProtein
    beta_c, beta_h, beta_0 : float, optional
        Best-Vendruscolo coefficients. Defaults
        ``(0.35, 2.0, 0.0)``.
    contact_cutoff : float, optional
        Heavy-atom contact cutoff (nm). Default ``0.65``.
    exclude_neighbours : int, optional
        Sequence-separation exclusion radius. Default ``2``.
    stride : int, optional
        Frame subsample. Default ``1``.
    weights : numpy.ndarray or False, optional
        Optional per-frame weights collapsing the frame axis to a
        per-residue mean (validated by
        :func:`soursop.ssutils.validate_weights`). Default ``False``
        (per-frame array returned).
    etol : float, optional
        Tolerance on ``sum(weights) == 1``. Default ``1e-7``.

    Returns
    -------
    residue_indices : numpy.ndarray, shape (n_res,)
        Residue indices for which ln(P) is defined.
    lnP : numpy.ndarray
        ``(n_frames, n_res)`` per-frame ln(P) by default;
        ``(n_res,)`` ensemble mean when ``weights`` is supplied.

    Raises
    ------
    SSException
        If the protein has no backbone amide N+H residues, or if
        ``weights`` fails validation.
    """
    res_NH, Nc = compute_Nc(
        protein,
        contact_cutoff=contact_cutoff,
        exclude_neighbours=exclude_neighbours,
        stride=stride,
    )
    res_NH_check, Nh = compute_Nh(
        protein,
        exclude_neighbours=exclude_neighbours,
        stride=stride,
    )
    if not np.array_equal(res_NH, res_NH_check):
        raise SSException(
            "Internal inconsistency: compute_Nc and compute_Nh returned "
            "different residue lists"
        )

    lnP = beta_c * Nc.astype(np.float64) + beta_h * Nh.astype(np.float64) + beta_0

    n_frames_total = protein.traj.n_frames
    validated_weights = validate_weights(
        weights, n_frames_total, stride=stride, etol=etol
    )
    if validated_weights is not False:
        lnP = weighted_mean(lnP, validated_weights, axis=0)

    return res_NH, lnP
