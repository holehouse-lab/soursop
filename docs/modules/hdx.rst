
sshdx
=========================================================

Overview
----------------------------

``sshdx`` predicts per-residue **HDX protection factors** (PFs) for backbone amide hydrogen-deuterium exchange from a structural ensemble, using the empirical **Best-Vendruscolo** relation

.. math::

   \ln P_i \;=\; \beta_c\, N_c(i) \;+\; \beta_h\, N_h(i) \;+\; \beta_0 ,

where for each residue :math:`i` with a backbone amide N-H:

* :math:`N_c(i)` is the number of *heavy atoms* in the protein within ``contact_cutoff`` (default 6.5 Å) of the amide N of residue :math:`i`, ignoring residues at sequence separation :math:`|i - j| \le` ``exclude_neighbours`` (default 2). The paper defines :math:`N_c` without a neighbour exclusion and notes that the choice of threshold "made little difference"; the ±2 exclusion is HDXer's convention;
* :math:`N_h(i)` is the number of hydrogen bonds formed by the amide H of residue :math:`i`. Best & Vendruscolo identify a hydrogen bond by "a cutoff of 2.4 Å between the donor hydrogen and the acceptor" (*Structure* 14:98), with no sequence-separation exclusion. Following the Best lab's reference implementation, `HDXer <https://github.com/Lucy-Forrest-Lab/HDXer>`_, an H-bond is counted for every protein **oxygen** atom - backbone or sidechain, in any residue - within ``hbond_cutoff`` (default 2.4 Å) of the amide H. Restricting the acceptors to oxygen and counting every H-bond are HDXer's conventions: the paper fitted :math:`\beta_h` with native H-bonds only and reports that a non-native definition "did not make an appreciable difference". The backbone-only geometric definition SOURSOP used up to 2.0.5 (:func:`mdtraj.wernet_nilsson` H-bonds to a backbone carbonyl O at :math:`|i - j| > 2`) is still available as ``hbond_method='wernet-nilsson'``.

Proline (no backbone H) and any residue that lacks a recognisable backbone amide H are dropped from the residue list. The three exported helpers all return ``(residue_indices, array)`` so the caller knows the per-residue mapping. All distances are plain Euclidean distances (no minimum-image convention), as everywhere else in SOURSOP.

Public entry points
----------------------------

* :func:`~soursop.sshdx.compute_Nc` — per-residue per-frame heavy-atom contacts.
* :func:`~soursop.sshdx.compute_Nh` — per-residue per-frame H-bond counts of the amide H (``hbond_method='distance'`` or ``'wernet-nilsson'``).
* :func:`~soursop.sshdx.compute_protection_factors` — per-residue ln(P), optionally collapsed to an ensemble mean via the package-wide ``weights`` argument.

All three accept ``stride`` and return ``(n_frames, n_residues)`` arrays — the natural input shape for :class:`soursop.ssbme.BME` / :class:`soursop.sscoper.COPER` reweighting against experimental protection factors.

**Defaults** (Best & Vendruscolo, *Structure* **14**, 97-106, 2006):

* :data:`~soursop.sshdx.DEFAULT_BETA_C` ``= 0.35``
* :data:`~soursop.sshdx.DEFAULT_BETA_H` ``= 2.0``
* :data:`~soursop.sshdx.DEFAULT_BETA_0` ``= 0.0``
* :data:`~soursop.sshdx.DEFAULT_CONTACT_CUTOFF_NM` ``= 0.65`` (= 6.5 Å)
* :data:`~soursop.sshdx.DEFAULT_EXCLUDE_NEIGHBOURS` ``= 2`` (for :math:`N_c`)
* :data:`~soursop.sshdx.DEFAULT_HBOND_CUTOFF_NM` ``= 0.24`` (= 2.4 Å, for :math:`N_h`)


Computing ln(P) for an ensemble
---------------------------------------------------------

::

    from soursop.sstrajectory import SSTrajectory
    from soursop.sshdx import compute_protection_factors

    traj    = SSTrajectory('traj.xtc', 'start.pdb')
    protein = traj.proteinTrajectoryList[0]

    # Per-frame, per-residue ln(P)
    residues, lnP = compute_protection_factors(protein)
    # lnP.shape == (n_frames, len(residues))

    # Ensemble-averaged ln(P) per residue (uniform weights)
    import numpy as np
    w = np.full(protein.n_frames, 1.0 / protein.n_frames)
    residues, lnP_mean = compute_protection_factors(protein, weights=w)

Tune the empirical coefficients if you are using a different parameterisation::

    residues, lnP = compute_protection_factors(
        protein, beta_c=0.5, beta_h=2.5, beta_0=-0.3)


Reweighting against experimental protection factors
---------------------------------------------------------

The ``(n_frames, n_residues)`` ln(P) matrix plugs directly into the BME / COPER reweighters::

    from soursop.sshdx import compute_protection_factors
    from soursop.ssbme import BME, ExperimentalObservable

    residues, lnP_calc = compute_protection_factors(protein)
    # experimental ln(P) per residue with its associated uncertainty
    obs = [ExperimentalObservable(lnP_exp[k], uncertainty=sigma_lnP[k],
                                  name=f"PF_res{residues[k]}")
           for k in range(len(residues))]

    result = BME(obs, lnP_calc).fit(theta=2.0, auto_theta=False)
    weights = result.weights


Pitfalls and conventions
---------------------------------------------------------

* **Units.** ``contact_cutoff`` is in **nanometres** (matches mdtraj). 6.5 Å = 0.65 nm.
* **H-bond convention.** By default ``compute_Nh`` counts every protein oxygen within 2.4 Å of the amide H, with no sequence exclusion - the Best-Vendruscolo definition that :math:`\beta_h = 2.0` was fitted to (the same-residue carbonyl of an extended residue, the C5 interaction, therefore counts). This finds several times more H-bonds than the backbone-only :func:`mdtraj.wernet_nilsson` definition SOURSOP used up to 2.0.5; on the folded NTL9 fixture the ensemble-mean :math:`\ln P` moves by about 0.3 on average. Pass ``hbond_method='wernet-nilsson'`` to reproduce the old numbers, and ``hbond_exclude_neighbours`` to impose a sequence exclusion on either definition.
* **Termini and caps.** A free N-terminal residue carries an NH3+ group rather than an amide N-H, so ``sshdx`` drops it (a residue whose N has no preceding backbone C). Capping groups (``ACE``, ``NME``, ``NMA``, ``NH2``) are likewise never reported. If the first residue follows an ``ACE`` cap it is a genuine amide and is retained. Proline is always excluded.
* **Per-frame statistics.** The per-frame ln(P) array is what the reweighters consume; the function only collapses the frame axis when you supply ``weights``. Always inspect ``lnP.mean(axis=0)`` (or feed a uniform-weight vector) when comparing to experimental ln(P) directly.
* **Parameterisation.** The default :math:`\beta_c = 0.35`, :math:`\beta_h = 2.0` come from Best & Vendruscolo's original fit on globular proteins. For IDRs / partially-folded states it is reasonable to keep them as a baseline, but a forward-model uncertainty in BME / COPER (or treating :math:`\beta_c`, :math:`\beta_h` as nuisance parameters) is recommended.


Key references
----------------------------

* Best, R. B. & Vendruscolo, M. *Structural Interpretation of Hydrogen Exchange Protection Factors in Proteins.* Structure **14**, 97-106 (2006). doi:`10.1016/j.str.2005.09.012 <https://doi.org/10.1016/j.str.2005.09.012>`_.
* Vendruscolo, M., Paci, E., Dobson, C. M. & Karplus, M. *Rare Fluctuations of Native Proteins Sampled by Equilibrium Hydrogen Exchange.* J. Am. Chem. Soc. **125**, 15686-15687 (2003).


.. automodule:: soursop.sshdx
.. autofunction:: compute_Nc
.. autofunction:: compute_Nh
.. autofunction:: compute_protection_factors
