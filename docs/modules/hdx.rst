
sshdx
=========================================================

Overview
----------------------------

``sshdx`` predicts per-residue **HDX protection factors** (PFs) for backbone amide hydrogen-deuterium exchange from a structural ensemble, using the empirical **Best-Vendruscolo** relation

.. math::

   \ln P_i \;=\; \beta_c\, N_c(i) \;+\; \beta_h\, N_h(i) \;+\; \beta_0 ,

where for each residue :math:`i` with a backbone amide N-H:

* :math:`N_c(i)` is the number of *heavy atoms* in the protein within ``contact_cutoff`` (default 6.5 Å) of the amide N of residue :math:`i`, ignoring residues at sequence separation :math:`|i - j| \le` ``exclude_neighbours`` (default 2);
* :math:`N_h(i)` is the number of *backbone* H-bonds the amide H of residue :math:`i` forms (as donor) to a backbone carbonyl O at sequence separation :math:`|i - j| >` ``exclude_neighbours``. Per-frame H-bonds are detected by :func:`mdtraj.wernet_nilsson`.

Proline (no backbone H) and any residue that lacks a recognisable backbone amide H are dropped from the residue list. The three exported helpers all return ``(residue_indices, array)`` so the caller knows the per-residue mapping.

Public entry points
----------------------------

* :func:`~soursop.sshdx.compute_Nc` — per-residue per-frame heavy-atom contacts.
* :func:`~soursop.sshdx.compute_Nh` — per-residue per-frame backbone H-bond counts.
* :func:`~soursop.sshdx.compute_protection_factors` — per-residue ln(P), optionally collapsed to an ensemble mean via the package-wide ``weights`` argument.

All three accept ``stride`` and return ``(n_frames, n_residues)`` arrays — the natural input shape for :class:`soursop.ssbme.BME` / :class:`soursop.sscoper.COPER` reweighting against experimental protection factors.

**Defaults** (Best & Vendruscolo, *Structure* **14**, 97-106, 2006):

* :data:`~soursop.sshdx.DEFAULT_BETA_C` ``= 0.35``
* :data:`~soursop.sshdx.DEFAULT_BETA_H` ``= 2.0``
* :data:`~soursop.sshdx.DEFAULT_BETA_0` ``= 0.0``
* :data:`~soursop.sshdx.DEFAULT_CONTACT_CUTOFF_NM` ``= 0.65`` (= 6.5 Å)
* :data:`~soursop.sshdx.DEFAULT_EXCLUDE_NEIGHBOURS` ``= 2``


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
* **H-bond convention.** ``compute_Nh`` only counts *backbone-to-backbone* H-bonds with the amide H of residue :math:`i` as donor and a backbone carbonyl O at :math:`|i - j| >` ``exclude_neighbours`` as acceptor — sidechain contributions and short-range :math:`i \rightarrow i, i+1, i+2` interactions are dropped, in line with the Best-Vendruscolo definition.
* **N-terminus.** The N-terminal amine of most force fields is named ``H1/H2/H3``; ``sshdx`` keeps the residue if it can find any of ``H``/``HN``/``H1`` on the backbone N. If your N-terminus is parameterised differently, mask out the first residue post hoc.
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
