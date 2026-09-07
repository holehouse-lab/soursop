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
* ``N_h(i)`` is the number of hydrogen bonds formed by the amide H of
  residue ``i``. Best & Vendruscolo (Structure 14:98) identify a hydrogen
  bond by "a cutoff of 2.4 A between the donor hydrogen and the acceptor",
  with no sequence-separation exclusion. Following the Best lab's
  reference implementation (HDXer), an H-bond is counted here for every
  protein oxygen atom - backbone or sidechain, in any residue - within
  ``hbond_cutoff`` (default 2.4 A) of the amide H. Two details are HDXer's
  conventions rather than the paper's: the acceptors are restricted to
  oxygen, and every H-bond is counted (the paper fitted ``beta_h`` with
  native H-bonds only, and reports that a non-native definition "did not
  make an appreciable difference"). The backbone-only geometric definition
  used by SOURSOP 2.0.5 (``mdtraj.wernet_nilsson`` H-bonds to a backbone
  carbonyl O at ``|i - j| > 2``) is retained as
  ``hbond_method='wernet-nilsson'``.

The ``|i - j| <= 2`` exclusion applied to ``N_c`` is likewise HDXer's
convention; the paper defines ``N_c`` simply as the heavy atoms within 6.5
A of the amide N and notes that "the choice of the threshold distance for
counting ... made little difference".

Both counts are returned per frame per residue, so the resulting
``(n_frames, n_residues)`` arrays plug directly into the SOURSOP BME /
COPER reweighters. The optional ``weights`` argument collapses the frame
axis to a single per-residue ensemble mean (per the package-wide
:func:`soursop.ssutils.validate_weights` contract).

Proline (no backbone H), capping groups (ACE, NME, NMA, NH2 - these are
not amino acids and do not carry an exchangeable backbone amide), a free
N-terminal residue (whose N carries an NH3+ group rather than an amide
N-H) and any residue lacking a recognisable backbone amide H are dropped
from the residue list; the function returns the residue indices it
covered alongside the data array. Note that an N-terminal residue
preceded by an ACE cap is a genuine amide and is retained.

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
from .ssutils import (
    validate_keyword_option,
    validate_stride,
    validate_weights,
    weighted_mean,
)

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
#: Sequence-separation exclusion for the heavy-atom contact count ``N_c``:
#: atoms in residues with ``|i - j| <= this`` are not counted. Also the
#: default exclusion for the ``'wernet-nilsson'`` H-bond method.
DEFAULT_EXCLUDE_NEIGHBOURS = 2
#: Acceptor cutoff for the Best-Vendruscolo H-bond count ``N_h``, in
#: **nanometres** (2.4 A): every protein oxygen closer than this to the
#: amide H counts as one H-bond (Best & Vendruscolo 2006; HDXer ``cut_Nh``).
DEFAULT_HBOND_CUTOFF_NM = 0.24
#: H-bond definitions accepted by :func:`compute_Nh`.
HBOND_METHODS = ("distance", "wernet-nilsson")

#: Atom-name fallbacks for the backbone amide hydrogen across common
#: force fields (CHARMM/AMBER use ``H`` or ``HN``; some N-terminal
#: parameterisations use ``H1``).
_BACKBONE_H_NAMES = ("H", "HN", "H1")

#: Residue names treated as capping groups. These carry no exchangeable
#: backbone amide and are never reported (mirrors the cap names in
#: :mod:`soursop.ssdata`).
_CAP_RESIDUE_NAMES = ("ACE", "NME", "NMA", "NH2", "FOR")


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
def _backbone_nh_map(topology):
    """Resolve backbone amide N+H atom indices per residue.

    Only residues with a genuine exchangeable backbone amide N-H are
    returned. Skipped are: proline (no amide H); capping groups
    (``ACE``, ``NME``, ``NMA``, ``NH2``, ``FOR``); a free N-terminal
    residue, whose N is an NH3+ group rather than an amide (identified
    as an N with no preceding residue that carries a backbone C - so a
    residue following an ACE cap is retained); residues with no backbone
    N (rare hetero residues); and residues with no recognisable backbone
    H name. The returned arrays are aligned 1:1.

    Parameters
    ----------
    topology : mdtraj.Topology
        Topology to scan.

    Returns
    -------
    residue_indices : numpy.ndarray (n_res,)
    n_atom_indices  : numpy.ndarray (n_res,)
    h_atom_indices  : numpy.ndarray (n_res,)
    """
    res_idx, n_idx, h_idx = [], [], []
    residues = list(topology.residues)
    for pos, residue in enumerate(residues):
        if residue.name == "PRO":
            continue
        if residue.name.upper() in _CAP_RESIDUE_NAMES:
            continue
        try:
            n_atom = next(a for a in residue.atoms if a.name == "N")
        except StopIteration:
            continue

        # a free N-terminus: the amide N must be bonded to the carbonyl C
        # of the preceding residue in the same chain, otherwise it is an
        # NH3+ group and does not exchange like a backbone amide
        if pos == 0 or residues[pos - 1].chain.index != residue.chain.index:
            continue
        if not any(a.name == "C" for a in residues[pos - 1].atoms):
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


def _oxygen_atom_map(topology):
    """Atom indices and residue indices of every oxygen atom (any residue).

    These are the H-bond acceptors of the Best-Vendruscolo ``'distance'``
    definition: backbone carbonyls, sidechain oxygens and C-terminal
    carboxylate oxygens alike.
    """
    atom_idx, res_idx = [], []
    for a in topology.atoms:
        if a.element is not None and a.element.symbol == "O":
            atom_idx.append(a.index)
            res_idx.append(a.residue.index)
    return np.asarray(atom_idx, dtype=int), np.asarray(res_idx, dtype=int)


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
    mdtraj's convention). Best & Vendruscolo define ``N_c`` as the heavy
    atoms within 6.5 A of the amide N and do not state a neighbour
    exclusion; the default ``exclude_neighbours=2`` is the convention of
    the Best lab's HDXer.

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
        backbone-amide residues; proline, capping groups, a free
        N-terminal residue and any residue lacking a recognisable
        backbone H are dropped - see :func:`_backbone_nh_map`).
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

    stride = validate_stride(stride, protein.traj.n_frames)
    traj = protein.traj[::stride] if stride != 1 else protein.traj
    n_frames = traj.n_frames

    Nc = np.zeros((n_frames, len(res_NH)), dtype=int)
    for k, (i_res, n_idx) in enumerate(zip(res_NH, n_atom_idx)):
        mask = np.abs(heavy_res_idx - i_res) > exclude_neighbours
        eligible = heavy_atom_idx[mask]
        if len(eligible) == 0:
            continue
        pairs = np.column_stack([np.full(len(eligible), n_idx, dtype=int), eligible])
        # plain distances (no minimum image), as everywhere else in soursop
        d = md.compute_distances(traj, pairs, periodic=False)  # nm, (n_frames, n_pairs)
        Nc[:, k] = np.sum(d < contact_cutoff, axis=1)

    return res_NH, Nc


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
def compute_Nh(
    protein,
    exclude_neighbours=None,
    stride=1,
    hbond_method="distance",
    hbond_cutoff=DEFAULT_HBOND_CUTOFF_NM,
):
    """Per-residue per-frame H-bond count ``N_h(i)`` of the backbone amide H.

    Two definitions are available:

    * ``'distance'`` (default) - the Best-Vendruscolo definition ("a
      cutoff of 2.4 A between the donor hydrogen and the acceptor",
      Structure 14:98) as implemented in the Best lab's HDXer: ``N_h(i)``
      is the number of protein **oxygen** atoms (backbone carbonyl,
      sidechain and carboxylate oxygens, in any residue) within
      ``hbond_cutoff`` of the amide H of residue ``i``. Best & Vendruscolo
      apply no sequence-separation exclusion, so by default none is
      applied here either; the same-residue carbonyl can therefore count
      (the C5 interaction of extended backbones). Restricting the
      acceptors to oxygen, and counting every H-bond rather than only
      native ones, follow HDXer; the paper reports that a non-native
      H-bond definition "did not make an appreciable difference".
    * ``'wernet-nilsson'`` - the definition SOURSOP used up to 2.0.5:
      backbone-to-backbone H-bonds with the amide H of residue ``i`` as
      donor and a backbone carbonyl O as acceptor, detected per frame by
      :func:`mdtraj.wernet_nilsson` (cone criterion on the
      Donor-H...Acceptor geometry), at sequence separation
      ``|i - j| > exclude_neighbours`` (default 2).

    The default ``beta_h = 2.0`` of :func:`compute_protection_factors` was
    fitted to the ``'distance'`` definition, which counts several times
    more H-bonds than the backbone-only geometric one; mixing the two
    biases ``ln P`` low.

    Parameters
    ----------
    protein : soursop.ssprotein.SSProtein
    exclude_neighbours : int or None, optional
        Acceptors in residues with ``|i - j| <= exclude_neighbours`` are
        ignored. ``None`` (default) uses the method's own convention: no
        exclusion for ``'distance'`` (Best-Vendruscolo), ``2`` for
        ``'wernet-nilsson'``.
    stride : int, optional
        Frame subsample. Default ``1``.
    hbond_method : {'distance', 'wernet-nilsson'}, optional
        H-bond definition, see above. Default ``'distance'``.
    hbond_cutoff : float, optional
        Amide-H to oxygen distance cutoff in **nanometres** for the
        ``'distance'`` method. Default ``0.24`` (= 2.4 A). Ignored by
        ``'wernet-nilsson'``.

    Returns
    -------
    residue_indices : numpy.ndarray, shape (n_res,)
        Same residue list as :func:`compute_Nc`.
    Nh : numpy.ndarray, shape (n_frames, n_res), int

    Raises
    ------
    SSException
        If the protein has no recognisable backbone amide N+H residues
        (e.g. a hydrogen-free or coarse-grained topology), if
        ``hbond_method`` is unknown, or if ``stride`` is invalid.
    """
    validate_keyword_option(hbond_method, list(HBOND_METHODS), "hbond_method")

    top = protein.traj.topology
    res_NH, _, h_atom_idx = _backbone_nh_map(top)
    if len(res_NH) == 0:
        raise SSException("compute_Nh: no residues with a backbone amide N+H found")

    stride = validate_stride(stride, protein.traj.n_frames)
    traj = protein.traj[::stride] if stride != 1 else protein.traj
    n_frames = traj.n_frames

    Nh = np.zeros((n_frames, len(res_NH)), dtype=int)

    if hbond_method == "distance":
        if hbond_cutoff <= 0:
            raise SSException(
                f"hbond_cutoff must be positive (in nm), got {hbond_cutoff}"
            )
        o_atom_idx, o_res_idx = _oxygen_atom_map(top)
        for k, (i_res, h_idx) in enumerate(zip(res_NH, h_atom_idx)):
            if exclude_neighbours is None:
                eligible = o_atom_idx
            else:
                eligible = o_atom_idx[np.abs(o_res_idx - i_res) > exclude_neighbours]
            if len(eligible) == 0:
                continue
            pairs = np.column_stack(
                [np.full(len(eligible), h_idx, dtype=int), eligible]
            )
            # plain distances (no minimum image), as everywhere else in soursop
            d = md.compute_distances(traj, pairs, periodic=False)
            Nh[:, k] = np.sum(d < hbond_cutoff, axis=1)
        return res_NH, Nh

    # 'wernet-nilsson': backbone donor H to backbone carbonyl O only
    if exclude_neighbours is None:
        exclude_neighbours = DEFAULT_EXCLUDE_NEIGHBOURS

    res_O, o_atom_idx = _backbone_carbonyl_o_map(top)
    h_to_res = dict(zip(h_atom_idx.tolist(), res_NH.tolist()))
    o_to_res = dict(zip(o_atom_idx.tolist(), res_O.tolist()))
    res_to_k = {int(r): k for k, r in enumerate(res_NH)}

    hbonds_per_frame = md.wernet_nilsson(traj, periodic=False)
    # mdtraj returns a list of length n_frames; each entry is an
    # ``(n_hbonds, 3)`` int array: (donor_heavy_idx, h_idx, acceptor_heavy_idx).
    for f, hb in enumerate(hbonds_per_frame):
        if hb.shape[0] == 0:
            continue
        for h_idx, a_idx in zip(hb[:, 1], hb[:, 2]):
            i_res = h_to_res.get(int(h_idx))
            j_res = o_to_res.get(int(a_idx))
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
    hbond_method="distance",
    hbond_cutoff=DEFAULT_HBOND_CUTOFF_NM,
    hbond_exclude_neighbours=None,
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
        Heavy-atom contact cutoff (nm) for ``N_c``. Default ``0.65``.
    exclude_neighbours : int, optional
        Sequence-separation exclusion radius for the heavy-atom contact
        count ``N_c``. Default ``2``.
    stride : int, optional
        Frame subsample. Default ``1``.
    weights : numpy.ndarray or False, optional
        Optional per-frame weights collapsing the frame axis to a
        per-residue mean (validated by
        :func:`soursop.ssutils.validate_weights`). Default ``False``
        (per-frame array returned).
    etol : float, optional
        Tolerance on ``sum(weights) == 1``. Default ``1e-7``.
    hbond_method : {'distance', 'wernet-nilsson'}, optional
        H-bond definition for ``N_h`` (see :func:`compute_Nh`). Default
        ``'distance'``, the Best-Vendruscolo definition the default
        ``beta_h`` was fitted to. Prior to 2.0.6 only the backbone-only
        ``'wernet-nilsson'`` definition existed.
    hbond_cutoff : float, optional
        Amide-H to oxygen cutoff (nm) for ``hbond_method='distance'``.
        Default ``0.24``.
    hbond_exclude_neighbours : int or None, optional
        Sequence-separation exclusion for the H-bond acceptors. ``None``
        (default) uses the method's convention: none for ``'distance'``,
        ``2`` for ``'wernet-nilsson'``.

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
        If the protein has no backbone amide N+H residues, if
        ``hbond_method`` is unknown, or if ``weights`` fails validation.
    """
    res_NH, Nc = compute_Nc(
        protein,
        contact_cutoff=contact_cutoff,
        exclude_neighbours=exclude_neighbours,
        stride=stride,
    )
    res_NH_check, Nh = compute_Nh(
        protein,
        exclude_neighbours=hbond_exclude_neighbours,
        stride=stride,
        hbond_method=hbond_method,
        hbond_cutoff=hbond_cutoff,
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
