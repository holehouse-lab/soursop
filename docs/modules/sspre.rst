
sspre
=========================================================

Overview
----------------------------

``sspre`` provides tools for computing synthetic paramagnetic relaxation enhancement (PRE) profiles from MD simulation ensembles and comparing them to experimental NMR data. The module exposes a single class, ``SSPRE``, which wraps an :class:`~soursop.ssprotein.SSProtein` object together with the experimental parameters of the PRE measurement.

**What PRE profiles measure.** In a PRE experiment, a nitroxide spin label (typically a MTSL-labelled cysteine) is introduced at a specific sequence position. The unpaired electron of the spin label accelerates the relaxation of nearby backbone amide protons, reducing their observed NMR signal intensity. The key observables are:

* **Intensity ratio** :math:`I_\text{para} / I_\text{dia}` — ratio of signal intensities in the paramagnetic vs. diamagnetic sample, ranging from 0 (near the label) to 1 (far away).
* **Transverse relaxation rate** :math:`\Gamma_2` — the spin-label-induced relaxation rate enhancement (in Hz), a complementary observable.

**Physical model.** The calculation uses the Solomon–Bloembergen framework: for each residue, the :math:`r^{-6}`-weighted mean distance to the spin-label atom is averaged over all frames of the trajectory. This ensemble averaging correctly captures the non-linear distance dependence of PRE relaxation and is the approach validated in the references below.

**Workflow.** An ``SSPRE`` object is created by supplying an ``SSProtein`` and the four experimental parameters: effective correlation time :math:`\tau_c` (ns), INEPT delay :math:`t_\text{delay}` (ms), diamagnetic transverse relaxation rate :math:`R_{2D}` (Hz), and proton Larmor frequency :math:`\omega_H` (Hz). ``generate_PRE_profile`` then computes the full-length intensity ratio and :math:`\Gamma_2` profiles for a given label position::

    from soursop.sstrajectory import SSTrajectory
    from soursop.sspre import SSPRE

    TrajOb = SSTrajectory('traj.xtc', 'start.pdb')
    protein = TrajOb.proteinTrajectoryList[0]

    # 600 MHz magnet, tau_c = 5 ns, t_delay = 16 ms, R_2D = 10 Hz
    pre = SSPRE(protein, tau_c=5, t_delay=16, R_2D=10, W_H=600000000)

    # spin label at residue 42 (CB atom by default)
    intensity_ratio, gamma2 = pre.generate_PRE_profile(label_position=42)

**Key references.** The methodology is described in the supplementary information of:

* Meng, W., Lyle, N., Luan, B., Raleigh, D.P., and Pappu, R.V. (2013). *Proc. Natl. Acad. Sci. U. S. A.* 110, 2123–2128.
* Das, R.K., Huang, Y., Phillips, A.H., Kriwacki, R.W., and Pappu, R.V. (2016). *Proc. Natl. Acad. Sci. U. S. A.* 113, 5616–5621.
* Peran, I., Holehouse, A. S., Carrico, I. S., Pappu, R. V., Bilsel, O., & Raleigh, D. P. (2019). *Proc. Natl. Acad. Sci. U. S. A.* 116(25), 12301–12310.


.. autoclass:: soursop.sspre.SSPRE

        .. automethod:: __init__
        .. automethod:: generate_PRE_profile

