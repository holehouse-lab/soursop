
ssnmr
=========================================================

Overview
----------------------------

``ssnmr`` covers three complementary NMR-observable predictions for IDP/IDR ensembles:

1. **Sequence-based random-coil chemical shifts** (:sup:`1`\ H, :sup:`13`\ C, :sup:`15`\ N), via :func:`~soursop.ssnmr.compute_random_coil_chemical_shifts`. Stateless — takes a sequence string, no trajectory needed.
2. **Structure-based scalar (J) couplings**: ³J(HN, Hα) computed per frame per residue from the φ dihedral via the Karplus relation, via :func:`~soursop.ssnmr.compute_J3_HN_HA` (and the generic Karplus evaluator :func:`~soursop.ssnmr.karplus` for arbitrary coefficients). Takes an :class:`~soursop.ssprotein.SSProtein`.
3. **NOE distances**: per-frame inter-atom distances via :func:`~soursop.ssnmr.compute_NOE_distances`, collapsed to an :math:`\langle r^{-p}\rangle^{-1/p}` ensemble average by :func:`~soursop.ssnmr.noe_ensemble_average`.


Random coil chemical shifts
----------------------------

The primary function, :func:`~soursop.ssnmr.compute_random_coil_chemical_shifts`, predicts sequence-corrected random coil :sup:`1`\ H, :sup:`13`\ C, and :sup:`15`\ N chemical shifts for a given amino acid sequence. These are useful as a disordered-state reference baseline when interpreting experimental NMR spectra of intrinsically disordered proteins (IDPs) or unfolded proteins.

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


Scalar (J) couplings
----------------------------

For each residue with a defined backbone φ dihedral, the three-bond ³J(HN, Hα) scalar coupling is well approximated by the Karplus relation

.. math::

   {}^{3}\!J(\mathrm{HN}, \mathrm{H}\alpha)
       = A\cos^{2}\!\big(\phi + \phi_{0}\big)
       + B\cos\!\big(\phi + \phi_{0}\big)
       + C ,

where ``A``, ``B``, ``C`` and ``φ₀`` are empirical coefficients fitted to experimental data. ``ssnmr`` ships **six** literature parameterisations of these coefficients, ported from the `biceps <https://github.com/vvoelz/biceps>`_ package (Voelz lab), which itself adapts MDTraj's ``mdtraj/nmr/scalar_couplings.py`` (Beauchamp et al.):

.. list-table::
   :header-rows: 1
   :widths: 22 12 12 12 18

   * - Model
     - A
     - B
     - C
     - σ (Hz)
   * - ``Bax2007`` (default)
     - 8.40
     - -1.36
     - 0.33
     - 0.36
   * - ``Bax1997``
     - 7.09
     - -1.42
     - 1.55
     - 0.39
   * - ``Ruterjans1999``
     - 7.90
     - -1.05
     - 0.65
     - 0.25
   * - ``Habeck``
     - 7.13
     - -1.31
     - 1.56
     - 0.34
   * - ``Vuister``
     - 6.51
     - -1.76
     - 1.60
     - 0.73
   * - ``Pardi``
     - 6.40
     - -1.40
     - 1.90
     - 0.76

All six models share ``φ₀ = −60°`` (the convention is to phase-shift the Karplus form by ``−60°`` so that ``θ = 0`` corresponds to the ideal HN–Cα–N–C′ eclipsed geometry). The per-model ``σ`` is the RMSD of the parameterisation against its training experimental dataset and is a sensible **forward-model uncertainty** to use when feeding J-couplings into the :doc:`BME <bme>` or :doc:`COPER <coper>` reweighters.

**Units.** ``ssnmr`` takes φ in **degrees** (consistent with :meth:`SSProtein.get_angles <soursop.ssprotein.SSProtein.get_angles>`) and returns J in **Hz**.

**Computing 3J(HN, Hα) from an ensemble**::

    from soursop.sstrajectory import SSTrajectory
    from soursop.ssnmr import compute_J3_HN_HA

    traj    = SSTrajectory('traj.xtc', 'start.pdb')
    protein = traj.proteinTrajectoryList[0]

    # Per-frame, per-residue J-couplings (shape: n_frames x n_phi).
    atoms, J = compute_J3_HN_HA(protein, model="Bax2007")

    # Ensemble mean + the model's forward-model uncertainty in Hz.
    atoms, J_mean, sigma = compute_J3_HN_HA(
        protein, model="Bax2007", weights=False, return_uncertainty=True)
    J_mean = J.mean(axis=0)

The ``(n_frames, n_phi)`` matrix is the natural input for the reweighters - so a typical workflow against an experimental J vector is::

    from soursop.ssbme import BME, ExperimentalObservable

    atoms, J_calc, sigma = compute_J3_HN_HA(protein, return_uncertainty=True)

    # one experimental observable per residue with a defined phi
    obs = [ExperimentalObservable(value=J_exp[k], uncertainty=sigma,
                                  name=f"3J_HN_HA_res{k}")
           for k in range(J_calc.shape[1])]

    result = BME(obs, J_calc).fit(theta=2.0, auto_theta=False)
    weights = result.weights

**Other Karplus parameterisations.** For any other Karplus-form coupling - ³J(Hα, C′), ³J(HN, Cβ), Bothner-By, Tvaroska, Aydin, ... - call :func:`~soursop.ssnmr.karplus` directly with the appropriate dihedral angle (in degrees) and your coefficients.

**Citations.** The HN-Hα coefficients are due to: Vögeli, B. *et al.* *J. Am. Chem. Soc.* **129**, 9377-9385 (2007) (Bax2007); Hu, J.-S. & Bax, A. *J. Am. Chem. Soc.* **119**, 6360-6368 (1997) (Bax1997); Schmidt, J. M. *et al.* *J. Biomol. NMR* **14**, 1-12 (1999) (Ruterjans1999); Habeck, M., Rieping, W. & Nilges, M. *J. Magn. Reson.* **177**, 160-165 (2005) (Habeck); Vuister, G. W. & Bax, A. *J. Am. Chem. Soc.* **115**, 7772-7777 (1993) (Vuister); Pardi, A., Billeter, M. & Wüthrich, K. *J. Mol. Biol.* **180**, 741-751 (1984) (Pardi).


NOE distances
----------------------------

The nuclear Overhauser effect (NOE) cross-peak between two protons depends on the inverse sixth power of the inter-proton distance averaged over the ensemble:

.. math::

   I_{ij} \propto \langle r_{ij}^{-6} \rangle .

Two helpers cover the typical workflow:

* :func:`~soursop.ssnmr.compute_NOE_distances` returns the per-frame inter-atom distance matrix in Angstroms for a list of atom pairs — the natural per-frame structural primitive.
* :func:`~soursop.ssnmr.noe_ensemble_average` collapses such a distance array via the NOE convention :math:`(\sum_i w_i\,d_i^{-p})^{-1/p}` (default :math:`p = 6`; some studies use :math:`p = 3`). It honours the package-wide ``weights`` contract: ``weights=False`` (default) gives the uniform mean, a per-frame weight vector gives a reweighted NOE distance.

**Computing ensemble-averaged NOE distances**::

    import numpy as np
    from soursop.ssnmr import compute_NOE_distances, noe_ensemble_average

    pairs = np.array([[0, 10], [0, 20], [5, 15]])   # atom indices
    d = compute_NOE_distances(protein, pairs)       # (n_frames, 3) A
    r_noe = noe_ensemble_average(d, power=6)        # (3,) A

**Feeding NOEs to BME / COPER.** The linear-additive observable for reweighting is :math:`r^{-p}`, *not* :math:`r` (because NOE intensity itself is the linear ensemble average of :math:`r^{-p}`). So the BME-ready ``calculated_values`` matrix is ``d ** -p`` and the experimental observable is ``r_exp ** -p``::

    from soursop.ssbme import BME, ExperimentalObservable

    calc = d ** -6                                  # (n_frames, n_pairs)
    obs = [ExperimentalObservable((r_exp[k]) ** -6,
                                  uncertainty=6.0 * (r_exp[k]) ** -7 * sigma_r[k],
                                  name=f"NOE_pair{k}")
           for k in range(len(r_exp))]              # uncertainty propagated from r_exp
    weights = BME(obs, calc).fit(theta=2.0, auto_theta=False).weights

(Note the experimental σ on the linear observable :math:`r^{-6}` is obtained by error-propagation from σ on :math:`r`: ``σ(r^{-6}) ≈ 6·r^{-7}·σ(r)``.)


.. automodule:: soursop.ssnmr
.. autofunction:: compute_random_coil_chemical_shifts
.. autofunction:: karplus
.. autofunction:: compute_J3_HN_HA
.. autofunction:: compute_NOE_distances
.. autofunction:: noe_ensemble_average
