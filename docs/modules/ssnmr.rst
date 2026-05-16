
ssnmr
=========================================================

Overview
----------------------------

``ssnmr`` provides functions for predicting nuclear magnetic resonance (NMR) observables from protein sequence data. Unlike :class:`~soursop.ssprotein.SSProtein` or :class:`~soursop.sstrajectory.SSTrajectory`, this module is **stateless** — it contains standalone functions rather than a class, and most functions take a sequence string directly rather than a trajectory object.

**Random coil chemical shifts.** The primary function, ``compute_random_coil_chemical_shifts``, predicts sequence-corrected random coil :sup:`1`\ H, :sup:`13`\ C, and :sup:`15`\ N chemical shifts for a given amino acid sequence. These are useful as a disordered-state reference baseline when interpreting experimental NMR spectra of intrinsically disordered proteins (IDPs) or unfolded proteins.

Corrections applied include:

* **Nearest-neighbour sequence effects** — shifts are adjusted for the two residues on either side of each position, using the correction factors of Kjaergaard & Poulsen (2011) and Schwarzinger et al. (2001).
* **Temperature** — linear corrections are applied relative to a 5 °C baseline.
* **pH** — charged-state populations for Asp, Glu, His, and phosphorylated residues (pSer, pThr, pTyr) are accounted for via fractional deprotonation at the given pH.
* **Perdeuteration** — optional corrections for fully deuterated protein samples.

**Supported residue types.** All 20 canonical amino acids are supported, along with three phosphorylated residues: phosphoserine (``pSer`` / ``SEP`` / ``PS``), phosphothreonine (``pThr`` / ``PTHR`` / ``PT``), and phosphotyrosine (``pTyr`` / ``PTYR`` / ``PY``). Phosphorylated residues cannot be combined with the perdeuteration corrections.

**Output format.** The function returns a list of per-residue dictionaries, one per position (excluding the two terminal padding residues), each containing keys ``Res``, ``Index``, ``CA``, ``CB``, ``CO``, ``N``, ``HN``, and ``HA``. Glycine lacks a Cβ (``CB`` is ``"**.***"``) and proline lacks a backbone amide (``N`` and ``HN`` are ``"*.***"``). Shifts are returned as floats or three-decimal-place strings depending on the ``asFloat`` flag.

**Example usage**::

    from soursop.ssnmr import compute_random_coil_chemical_shifts

    sequence = "MAEQKLISEEDL"
    shifts = compute_random_coil_chemical_shifts(sequence, temperature=25, pH=7.4)

    for residue in shifts:
        print(residue['Res'], residue['CA'], residue['N'])


.. automodule:: soursop.ssnmr
.. autofunction:: compute_random_coil_chemical_shifts
