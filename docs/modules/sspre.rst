
sspre
=========================================================

Overview
----------------------------

``sspre`` provides tools for computing synthetic paramagnetic relaxation enhancement (PRE) profiles from MD simulation ensembles and comparing them to experimental NMR data. The module exposes a single class, ``SSPRE``, which wraps an :class:`~soursop.ssprotein.SSProtein` object together with the experimental parameters of the PRE measurement.

**What PRE profiles measure.** In a PRE experiment, a nitroxide spin label (typically a MTSL-labelled cysteine) is introduced at a specific sequence position. The unpaired electron of the spin label accelerates the relaxation of nearby backbone amide protons, reducing their observed NMR signal intensity. The key observables are:

* **Intensity ratio** :math:`I_\text{para} / I_\text{dia}` — ratio of signal intensities in the paramagnetic vs. diamagnetic sample, ranging from 0 (near the label) to 1 (far away).
* **Transverse relaxation rate** :math:`\Gamma_2` — the spin-label-induced relaxation rate enhancement (in Hz), a complementary observable.

**Physical model.** The calculation uses the Solomon–Bloembergen framework: for each residue, the :math:`r^{-6}`-weighted mean distance to the paramagnetic centre is averaged over all frames of the trajectory. This ensemble averaging correctly captures the non-linear distance dependence of PRE relaxation and is the approach validated in the references below. As of SOURSOP 2.0.2 the paramagnetic centre is, by default, a coarse-grained spin-label cloud (see below) rather than a point on the ``spin_label_atom``.

.. warning::

   **Breaking change in SOURSOP 2.0.2.** ``generate_PRE_profile`` now defaults to ``use_label=True``, so the paramagnetic centre is modelled as a DEER-PREdict-calibrated spin-label cloud rather than sitting directly on the CB atom. This yields different — and more accurate — profiles than SOURSOP ≤ 2.0.1. To reproduce results from earlier versions exactly, pass ``use_label=False``.

**Workflow.** An ``SSPRE`` object is created by supplying an ``SSProtein`` and the four experimental parameters: effective correlation time :math:`\tau_c` (ns), INEPT delay :math:`t_\text{delay}` (ms), diamagnetic transverse relaxation rate :math:`R_{2D}` (Hz), and proton Larmor frequency :math:`\omega_H` (Hz). ``generate_PRE_profile`` then computes the full-length intensity ratio and :math:`\Gamma_2` profiles for a given label position::

    from soursop.sstrajectory import SSTrajectory
    from soursop.sspre import SSPRE

    TrajOb = SSTrajectory('traj.xtc', 'start.pdb')
    protein = TrajOb.proteinTrajectoryList[0]

    # 600 MHz magnet, tau_c = 5 ns, t_delay = 16 ms, R_2D = 10 Hz
    pre = SSPRE(protein, tau_c=5, t_delay=16, R_2D=10, W_H=600000000)

    # spin label at residue 42; uses the calibrated label-cloud model by default
    intensity_ratio, gamma2 = pre.generate_PRE_profile(label_position=42)

**Coarse-grained spin-label cloud (the default since 2.0.2).** Rather than placing the paramagnetic centre directly on the ``spin_label_atom`` (``CB``), the nitroxide is modelled as a coarse-grained cloud of beads placed a fixed ``label_distance`` (default 7.0 Å) from the anchor atom — the geometry of an MTSL side chain — without requiring an explicit all-atom rotamer library. When a distinct CB is available the cloud is a cone about the CA→CB direction (so it points away from the backbone); on coarse-grained CA-only chains it falls back to a full isotropic sphere, so the same call works at both resolutions. The relaxation is averaged over the whole cloud as well as over frames, so the :math:`r^{-6}` non-linearity is preserved across both the conformer cloud and the ensemble. The cloud is generated deterministically (a Fibonacci lattice), so results are reproducible. To recover the classic point-at-CB model used by SOURSOP ≤ 2.0.1, pass ``use_label=False``::

    # explicit calibrated label-cloud model (this is now the default)
    intensity_ratio, gamma2 = pre.generate_PRE_profile(label_position=42, use_label=True)

    # classic point-at-CB model (SOURSOP <= 2.0.1 behaviour)
    intensity_ratio, gamma2 = pre.generate_PRE_profile(label_position=42, use_label=False)

Steric exclusion of the label bead against the chain is always applied under ``use_label=True``, either as a ``'hard'`` cutoff (beads within ``label_bead_radius`` of any non-label CA are dropped) or a ``'soft'`` WCA-like repulsive wall (the default). The defaults (``label_steric='soft'``, ``label_bead_radius=5.5`` Å, ``label_wall_stiffness=6.0``) were calibrated against DEER-PREdict's explicit MTSL rotamer PRE profiles across 85 disordered and ~2700 folded label sites; the soft wall is marginally more accurate than the hard cutoff and transfers better between folded and disordered chains. This stops the cloud from being placed through the protein interior, which matters most for folded/compact structures. See the ``generate_PRE_profile`` docstring below for the full list of ``label_*`` arguments.

**Key references.** The methodology is described in the supplementary information of:

* Meng, W., Lyle, N., Luan, B., Raleigh, D.P., and Pappu, R.V. (2013). *Proc. Natl. Acad. Sci. U. S. A.* 110, 2123–2128.
* Das, R.K., Huang, Y., Phillips, A.H., Kriwacki, R.W., and Pappu, R.V. (2016). *Proc. Natl. Acad. Sci. U. S. A.* 113, 5616–5621.
* Peran, I., Holehouse, A. S., Carrico, I. S., Pappu, R. V., Bilsel, O., & Raleigh, D. P. (2019). *Proc. Natl. Acad. Sci. U. S. A.* 116(25), 12301–12310.


.. autoclass:: soursop.sspre.SSPRE

        .. automethod:: __init__
        .. automethod:: generate_PRE_profile

