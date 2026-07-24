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
ssdata contains all specific data used, and is disinct from configs, which are
package specific global configurations, while ssdata contains information that
is unambigious as well-defined.
"""

import csv
import os

import numpy as np
import soursop

from .ssexceptions import SSException

THREE_TO_ONE = {
    "ALA": "A",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "ASN": "N",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y",
    "ACE": "<",
    "NME": ">",
    "NAC": ">",
    "NH2": "(",
    "FOR": ")",
}

ONE_TO_THREE = {
    "A": "ALA",
    "C": "CYS",
    "D": "ASP",
    "E": "GLU",
    "F": "PHE",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "K": "LYS",
    "L": "LEU",
    "M": "MET",
    "N": "ASN",
    "P": "PRO",
    "Q": "GLN",
    "R": "ARG",
    "S": "SER",
    "T": "THR",
    "V": "VAL",
    "W": "TRP",
    "Y": "TYR",
    "<": "ACE",
    ">": "NME",
    "(": "NH2",
    ")": "FOR",
}

DEFAULT_SIDECHAIN_VECTOR_ATOMS = {
    "ALA": "CB",
    "CYS": "SG",
    "ASP": "CG",
    "ASH": "CG",
    "GLU": "CD",
    "GLH": "CD",
    "PHE": "CZ",
    "GLY": "ERROR",
    "HIS": "NE2",
    "HID": "NE2",
    "HIE": "NE2",
    "HIP": "NE2",
    "ILE": "CD1",
    "LYS": "NZ",
    "LYD": "NZ",
    "KAC": "NZ",
    "KM1": "NZ",
    "KM2": "NZ",
    "KM3": "NZ",
    "LEU": "CG",
    "MET": "CE",
    "ASN": "CG",
    "PRO": "CG",
    "GLN": "CD",
    "ARG": "CZ",
    "SER": "OG",
    "SEP": "OG",
    "THR": "CB",
    "TPO": "CB",
    "VAL": "CB",
    "TRP": "CG",
    "TYR": "CZ",
    "PTR": "CZ",
    "ACE": "ERROR",
    "NME": "ERROR",
    "FOR": "ERROR",
    "NH2": "ERROR",
}

# list of valid residue names as supported by SOURSOP. These are the resnames compared againt when assessing for if a molecule is a protein or not.
ALL_VALID_RESIDUE_NAMES = [
    "ALA",
    "CYS",
    "ASP",
    "ASH",
    "GLU",
    "GLH",
    "PHE",
    "GLY",
    "HIE",
    "HIS",
    "HID",
    "HIP",
    "ILE",
    "LEU",
    "LYS",
    "LYD",
    "MET",
    "ASN",
    "PRO",
    "GLN",
    "ARG",
    "SER",
    "THR",
    "VAL",
    "TRP",
    "TYR",
    "AIB",
    "ABA",
    "NVA",
    "NLE",
    "ORN",
    "DAB",
    "PTR",
    "TPO",
    "SEP",
    "KAC",
    "KM1",
    "KM2",
    "KM3",
    "ACE",
    "NME",
    "FOR",
    "NH2",
]


EV_RESIDUE_MAPPER = {
    # proline is special
    "PRO": "PRO",
    # approximately alanine - not really
    "GLY": "ALA",
    "ALA": "ALA",
    "CYS": "ALA",
    "ASN": "ALA",
    "GLN": "ALA",
    "SER": "ALA",
    "THR": "ALA",
    # approximately leucine - not really
    "LEU": "LEU",
    "MET": "LEU",
    "ASP": "LEU",
    "GLU": "LEU",
    "ARG": "LEU",
    "VAL": "LEU",
    "TRP": "LEU",
    "TYR": "LEU",
    "PHE": "LEU",
    "HIS": "LEU",
    "ILE": "LEU",
    "LYS": "LEU",
}

try:
    PSI_EV_ANGLES_DICT = np.load(
        soursop.get_data("psi_excluded_volume_tripeptides.pickle"), allow_pickle=True
    )

    PHI_EV_ANGLES_DICT = np.load(
        soursop.get_data("phi_excluded_volume_tripeptides.pickle"), allow_pickle=True
    )
except Exception:
    PSI_EV_ANGLES_DICT = None
    PHI_EV_ANGLES_DICT = None
    print("WARNING: No Phi/Psi data loaded for PENGUIN analysis...")


# .....................................................................................
#
# Coarse-grained (one-bead-per-residue) force-field bead sizes.
#
# Every model in this family shares the same bead typing (1-20 solvent-exposed amino
# acids, 21-40 the 'buried' copies, 41-44 RNA nucleotides in the Mpipi variants), so
# a model is fully described here by its per-residue-type sigma. Sigma is the diagonal
# (i,i) pair sigma from the model's LAMMPS parameter file - i.e. the bead *diameter*,
# the separation at which the pair potential crosses zero - in Angstroms.
#
# Each table lives in soursop/data/cg/<slug>-sigma.csv, named for the force field's
# slug, so dropping a new table into that directory registers the model automatically.

#: Directory holding the bundled per-force-field bead-size tables.
CG_SIGMA_DIR = soursop.get_data("cg")

#: Suffix every bundled bead-size table carries; the slug is the leading part.
_CG_SIGMA_SUFFIX = "-sigma.csv"


def _discover_cg_forcefields():
    """Slugs of every coarse-grained force field with a bundled sigma table.

    Derived by listing ``soursop/data/cg`` rather than hard-coded, so the set of
    supported models cannot drift away from the files actually shipped.

    Returns
    -------
    tuple of str
        Sorted force-field slugs (e.g. ``('cation-pi-1', ..., 'mpipi-gg')``).
        Empty if the data directory is missing.
    """

    try:
        names = os.listdir(CG_SIGMA_DIR)
    except OSError:
        return ()

    return tuple(
        sorted(
            n[: -len(_CG_SIGMA_SUFFIX)] for n in names if n.endswith(_CG_SIGMA_SUFFIX)
        )
    )


#: Valid ``forcefield`` slugs for the coarse-grained SASA calculation.
CG_FORCEFIELDS = _discover_cg_forcefields()

# memoisation for get_cg_bead_sigmas(), keyed by (forcefield, buried)
_CG_SIGMA_CACHE = {}


def get_cg_bead_sigmas(forcefield, buried=False):
    """Per-residue bead diameters (sigma) for a one-bead-per-residue force field.

    Reads the bundled ``soursop/data/cg/<forcefield>-sigma.csv`` table and
    returns the mapping from three-letter residue name to the model's bead
    sigma. Sigma here is the diagonal (i,i) pair sigma taken from the model's
    LAMMPS parameter file, which is the bead *diameter* - so a Shrake-Rupley
    bead radius is ``sigma/2``.

    Results are memoised, so repeated calls are free.

    Parameters
    ----------
    forcefield : str
        Force-field slug; must be one of :data:`CG_FORCEFIELDS`.
    buried : bool, optional
        If True read the ``sigma_buried_angstrom`` column - the sigmas of the
        'buried' bead types (21-40), the ~30% weaker copies these models define
        for residues sequestered inside a folded domain. Default is False (the
        solvent-exposed types, 1-20). Residues with no buried sigma defined
        (the RNA nucleotides) are omitted when this is True.

    Returns
    -------
    dict
        Three-letter residue name -> sigma in Angstroms.

    Raises
    ------
    SSException
        If ``forcefield`` is not a recognised slug, or the table cannot be read.

    Example
    -------
    >>> get_cg_bead_sigmas('mpipi')['GLY']
    4.69511
    """

    key = (forcefield, bool(buried))
    if key in _CG_SIGMA_CACHE:
        return _CG_SIGMA_CACHE[key]

    if forcefield not in CG_FORCEFIELDS:
        raise SSException(
            f"Unknown coarse-grained force field '{forcefield}'. Valid options are: "
            f"{', '.join(CG_FORCEFIELDS)}"
        )

    column = "sigma_buried_angstrom" if buried else "sigma_angstrom"
    fname = os.path.join(CG_SIGMA_DIR, f"{forcefield}{_CG_SIGMA_SUFFIX}")

    sigmas = {}
    try:
        with open(fname, newline="") as fh:
            for row in csv.DictReader(fh):
                value = (row.get(column) or "").strip()

                # rows with no value in the requested column simply do not
                # define that bead type (e.g. RNA has no 'buried' variant)
                if len(value) == 0:
                    continue

                sigmas[row["residue"].strip().upper()] = float(value)

    except OSError as e:
        raise SSException(
            f"Could not read the bead-size table for force field '{forcefield}' "
            f"[{fname}]: {e}"
        )

    if len(sigmas) == 0:
        raise SSException(
            f"The bead-size table for force field '{forcefield}' [{fname}] defined no "
            f"'{column}' values"
        )

    _CG_SIGMA_CACHE[key] = sigmas
    return sigmas
