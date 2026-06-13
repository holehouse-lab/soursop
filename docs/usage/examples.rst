
Examples
=========================================================

This page walks through common IDP analysis workflows using SOURSOP, organized from basic trajectory loading through to advanced sampling-quality assessment. All examples assume you have a trajectory file (e.g. ``traj.xtc``) and a topology file (e.g. ``start.pdb``).

A larger collection of worked examples and Jupyter notebooks, including the analyses from the SOURSOP publication, can be found in the `supporting data repository <https://github.com/holehouse-lab/supportingdata/tree/master/2023/lalmansingh_2023>`_.

.. contents:: On this page
   :local:
   :depth: 1


1. Reading a trajectory
---------------------------------------------------------

All SOURSOP analyses begin by reading a trajectory with :class:`~soursop.sstrajectory.SSTrajectory`. Upon loading, SOURSOP automatically identifies every protein chain and exposes each as an :class:`~soursop.ssprotein.SSProtein` object::

    from soursop.sstrajectory import SSTrajectory
    import numpy as np

    # read in a trajectory and topology
    traj = SSTrajectory('traj.xtc', 'start.pdb')

    # retrieve the first (and often only) protein chain
    protein = traj.proteinTrajectoryList[0]

    # basic properties
    print(f"Residues : {protein.n_residues}")
    print(f"Frames   : {protein.n_frames}")
    print(f"Sequence : {protein.get_amino_acid_sequence(oneletter=True)}")

SOURSOP also accepts pre-loaded ``mdtraj`` trajectory objects via the ``TRJ`` keyword, which is useful for scripting pipelines that already have a trajectory in memory::

    import mdtraj as md
    mdtraj_traj = md.load('traj.xtc', top='start.pdb')
    traj = SSTrajectory(TRJ=mdtraj_traj)
    protein = traj.proteinTrajectoryList[0]


2. Global dimensions
---------------------------------------------------------

The most common IDP observables describe the overall size and shape of the chain. All functions return per-frame NumPy arrays, so standard statistics apply directly.

**Radius of gyration** :math:`R_g`, **hydrodynamic radius** :math:`R_h`, and **end-to-end distance** :math:`r_{ee}`::

    rg  = protein.get_radius_of_gyration()      # Angstroms, shape (n_frames,)
    rh  = protein.get_hydrodynamic_radius()     # Angstroms, shape (n_frames,)
    e2e = protein.get_end_to_end_distance()     # Angstroms, shape (n_frames,)

    print(f"Mean Rg  = {np.mean(rg):.2f} ± {np.std(rg):.2f} Å")
    print(f"Mean Rh  = {np.mean(rh):.2f} ± {np.std(rh):.2f} Å")
    print(f"Mean e2e = {np.mean(e2e):.2f} ± {np.std(e2e):.2f} Å")

**Asphericity** describes how far the chain deviates from a sphere (0 = perfectly spherical, 1 = rod-like)::

    asph = protein.get_asphericity()
    print(f"Mean asphericity = {np.mean(asph):.3f}")

**Correlation between** :math:`r_{ee}` **and** :math:`R_g` — a useful diagnostic for whether the two global size metrics are capturing consistent information::

    pearson_r, pval = protein.get_end_to_end_vs_rg_correlation()
    print(f"Pearson r(e2e, Rg) = {pearson_r:.3f}  (p = {pval:.2e})")

**Overlap concentration** :math:`c^*` estimates the concentration above which chains begin to crowd one another::

    c_star = protein.get_overlap_concentration()
    print(f"c* = {c_star:.4f} mg/mL")


3. Polymer scaling and internal structure
---------------------------------------------------------

IDR conformational behaviour is often interpreted through the lens of polymer physics. The **internal scaling profile** :math:`\langle r^2(i,j) \rangle` reports the mean-square inter-residue distance as a function of sequence separation :math:`|i - j|`.

**Internal scaling** (mean across the ensemble)::

    import matplotlib.pyplot as plt

    separation, mean_r2 = protein.get_internal_scaling(mode='CA', mean_vals=True)

    plt.loglog(separation, mean_r2)
    plt.xlabel("Sequence separation |i - j|")
    plt.ylabel(r"$\langle r^2 \rangle$ (Å²)")
    plt.title("Internal scaling profile")
    plt.show()

**Scaling exponent** :math:`\nu` — the Flory exponent extracted by fitting :math:`\sqrt{\langle r^2 \rangle} = A_0\,|i-j|^{\nu}`. ``get_scaling_exponent`` returns a 10-element tuple; the first two entries are the point estimates ``nu`` and ``A0``, entries 2–5 are the bootstrap confidence-interval bounds on each, and entries 6–7 are the reduced :math:`\chi^2` of the fit::

    result = protein.get_scaling_exponent(mode='CA')
    nu, A0 = result[0], result[1]
    nu_lo, nu_hi = result[2], result[3]          # 95% CI on nu (entries 4,5 are the A0 CI)
    redchi = result[7]                           # reduced chi^2 across all separations

    print(f"Flory exponent ν = {nu:.3f}  (95% CI {nu_lo:.3f}–{nu_hi:.3f})")
    print(f"homopolymer-fit reduced χ² = {redchi:.2f}")
    # ν ≈ 0.5  → Gaussian / theta-solvent behaviour
    # ν ≈ 0.6  → self-avoiding random coil (good solvent)
    # ν < 0.5  → compact / collapsed chain

Uncertainties on :math:`\nu` and :math:`A_0` are confidence intervals from a
frame-level bootstrap (frames resampled with replacement), so they tighten as
more — and more decorrelated — frames are supplied; tune the resampling with
``n_bootstrap`` (default 200) and ``confidence_interval`` (default 95.0). The
returned reduced :math:`\chi^2` is a genuine goodness-of-fit for the power-law
model (≈1 indicates the data are consistent with a single homopolymer scaling
law within their bootstrap errors; substantially larger values flag systematic
deviation, e.g. heteropolymeric structure).

**Polymer-scaled distance map** normalises the mean inter-residue distance matrix by the expected excluded-volume scaling, highlighting regions that are more compact or more expanded than a reference random coil. It returns the map together with the ``nu``, ``A0`` and reduced :math:`\chi^2` of the homopolymer fit it performs internally::

    dmap, nu, A0, redchi = protein.get_polymer_scaled_distance_map(mode='CA')

    plt.imshow(dmap, origin='lower', cmap='RdBu_r')
    plt.colorbar(label='Normalized distance')
    plt.xlabel('Residue index')
    plt.ylabel('Residue index')
    plt.title('Polymer-scaled distance map')
    plt.show()

**Local heterogeneity** in the scaling behaviour measures how the local Flory exponent varies along the chain, revealing compact or expanded subregions::

    local_het = protein.get_local_heterogeneity(window_size=10)


4. Secondary structure and backbone angles
---------------------------------------------------------

**DSSP secondary structure** assigns a secondary structure label to each residue in every frame::

    dssp = protein.get_secondary_structure_DSSP()
    # returns an (n_frames, n_residues) array of single-character labels

    # mean helicity per residue
    helicity = np.mean(dssp == 'H', axis=0)

    plt.bar(range(protein.n_residues), helicity)
    plt.xlabel('Residue index')
    plt.ylabel('Fractional helicity')
    plt.title('Per-residue α-helix propensity')
    plt.show()

.. note::

   For **SWAN** two-bead (``CA``/``CB``) coarse-grained trajectories,
   ``get_secondary_structure_DSSP`` automatically assigns helix/β/coil
   from the ``CA`` trace (DSSP itself needs the full backbone, which SWAN
   does not have), so the call above works unchanged. The φ/ψ-based
   functions in this section — ``get_secondary_structure_BBSEG`` and
   ``get_angles`` — are *not* available for SWAN chains and raise an
   ``SSException``.

**BBSEG backbone-torsion classification** provides an 8-state assignment based on φ/ψ backbone dihedral regions, which is particularly useful for IDRs where the DSSP labels can be sparse or noisy::

    bbseg = protein.get_secondary_structure_BBSEG()
    # returns an (n_frames, n_residues) array of integer labels (0–7)
    # 0=unassigned, 1=α-helix, 2=PPII, 3=β-strand, 4=turn, ...

**Backbone dihedral angles**::

    angles = protein.get_angles()
    # returns a dict with keys 'phi' and 'psi', each (n_frames, n_residues)

    phi = angles['phi']
    psi = angles['psi']

    # Ramachandran plot for residue 10
    residue_idx = 10
    plt.scatter(np.degrees(phi[:, residue_idx]),
                np.degrees(psi[:, residue_idx]),
                alpha=0.3, s=2)
    plt.xlabel('φ (°)')
    plt.ylabel('ψ (°)')
    plt.title(f'Ramachandran plot — residue {residue_idx}')
    plt.show()


5. Distance maps and contacts
---------------------------------------------------------

**Mean inter-residue distance map**::

    mean_dist, std_dist = protein.get_distance_map(mode='CA')

    plt.imshow(mean_dist, origin='lower', cmap='viridis_r')
    plt.colorbar(label='Mean CA–CA distance (Å)')
    plt.xlabel('Residue j')
    plt.ylabel('Residue i')
    plt.title('Mean inter-residue distance map')
    plt.show()

**Contact map** — fraction of frames in which each residue pair is within a cutoff distance (default 8 Å)::

    contact_map = protein.get_contact_map(distance_thresh=8.0, mode='CA')

    plt.imshow(contact_map, origin='lower', cmap='hot_r', vmin=0, vmax=1)
    plt.colorbar(label='Contact frequency')
    plt.xlabel('Residue j')
    plt.ylabel('Residue i')
    plt.title('Contact map (CA, 8 Å cutoff)')
    plt.show()

**RMSD** to a reference frame (here, frame 0)::

    rmsd = protein.get_RMSD(frame_index=0, region=[0, protein.n_residues])
    print(f"Mean RMSD from frame 0: {np.mean(rmsd):.2f} Å")


6. Solvent accessibility
---------------------------------------------------------

**Per-residue SASA** averaged across the trajectory::

    mean_sasa, std_sasa = protein.get_all_SASA()
    # mean_sasa: shape (n_residues,), units Å²

    plt.bar(range(protein.n_residues), mean_sasa, yerr=std_sasa, capsize=2)
    plt.xlabel('Residue index')
    plt.ylabel('SASA (Å²)')
    plt.title('Per-residue solvent accessibility')
    plt.show()

**Regional SASA** for a specific stretch of residues — useful for assessing the accessibility of a functional linear motif::

    # SASA for residues 10 to 20 (inclusive)
    mean_region, std_region = protein.get_regional_SASA(R1=10, R2=20)
    print(f"Mean SASA (residues 10–20) = {mean_region:.1f} ± {std_region:.1f} Å²")

**Site accessibility** returns the fraction of frames in which a residue's SASA exceeds a threshold, giving a per-residue accessibility score::

    accessibility = protein.get_site_accessibility()


7. Multi-chain systems
---------------------------------------------------------

For systems with more than one protein chain, system-level analyses are performed on the :class:`~soursop.sstrajectory.SSTrajectory` object rather than on individual :class:`~soursop.ssprotein.SSProtein` objects::

    traj = SSTrajectory('traj.xtc', 'start.pdb')

    # overall Rg combining all chains
    rg_total = traj.get_overall_radius_of_gyration()
    print(f"System Rg: {np.mean(rg_total):.2f} Å")

    # inter-chain distance map between chain 0 and chain 1
    mean_dist, std_dist = traj.get_interchain_distance_map(
        proteinID1=0, proteinID2=1, mode='CA'
    )

    # inter-chain contact map
    contact_map = traj.get_interchain_contact_map(
        proteinID1=0, proteinID2=1, distance_thresh=8.0
    )


8. NMR comparison
---------------------------------------------------------

**Random coil chemical shifts** (via :mod:`~soursop.ssnmr`) predict the expected NMR chemical shifts for a sequence in the disordered limit, corrected for temperature, pH, and nearest-neighbour effects. These are a useful reference baseline when interpreting experimental IDP spectra::

    from soursop.ssnmr import compute_random_coil_chemical_shifts

    sequence = protein.get_amino_acid_sequence(oneletter=True)
    shifts = compute_random_coil_chemical_shifts(
        sequence, temperature=25, pH=7.4
    )

    # extract CA shifts
    ca_shifts = [res['CA'] for res in shifts]
    print("Predicted CA chemical shifts:", ca_shifts)

**Backbone scalar (J) couplings.** ³J(HN, Hα) is computed per frame per residue from the φ dihedral via the Karplus relation, using any of the six literature parameterisations shipped in :mod:`~soursop.ssnmr` (Bax2007, Bax1997, Ruterjans1999, Habeck, Vuister, Pardi). The returned ``(n_frames, n_phi)`` matrix is the natural input for the BME / COPER reweighters::

    from soursop.ssnmr import compute_J3_HN_HA

    # per-frame, per-residue J-couplings (shape n_frames x n_phi)
    atoms, J = compute_J3_HN_HA(protein, model="Bax2007")

    # ensemble mean + the model's forward-model uncertainty in Hz
    atoms, J_mean, sigma = compute_J3_HN_HA(
        protein, weights=False, return_uncertainty=True)
    J_mean = J.mean(axis=0)

    # feeding the result to BME against an experimental J vector
    from soursop.ssbme import BME, ExperimentalObservable
    obs = [ExperimentalObservable(J_exp[k], sigma, name=f"3J_res{k}")
           for k in range(J.shape[1])]
    weights = BME(obs, J).fit(theta=2.0, auto_theta=False).weights

**NOE ⟨r⁻⁶⟩ ensemble distances.** Inter-proton distances reported by NOE cross-peaks are not linear ensemble averages but :math:`\langle r^{-6}\rangle^{-1/6}` averages. :mod:`~soursop.ssnmr` exposes the per-frame distance primitive and the NOE collapse rule separately so both BME-style reweighting and direct experimental comparison are clean::

    import numpy as np
    from soursop.ssnmr import compute_NOE_distances, noe_ensemble_average

    pairs = np.array([[0, 10], [0, 20], [5, 15]])    # atom indices
    d = compute_NOE_distances(protein, pairs)        # (n_frames, n_pairs) in Å
    r_noe = noe_ensemble_average(d, power=6)         # (n_pairs,) Å

    # For BME / COPER the linear observable is r^-p, not r:
    from soursop.ssbme import BME, ExperimentalObservable
    calc = d ** -6
    obs  = [ExperimentalObservable(r_exp[k] ** -6,
                                   uncertainty=6 * r_exp[k] ** -7 * sigma_r[k],
                                   name=f"NOE_{k}")
            for k in range(len(r_exp))]
    weights = BME(obs, calc).fit(theta=2.0, auto_theta=False).weights

**HDX protection factors (Best–Vendruscolo).** :mod:`~soursop.sshdx` predicts per-residue ln(P) from per-frame heavy-atom contacts (``N_c``) and backbone H-bond counts (``N_h``)::

    from soursop.sshdx import compute_protection_factors

    residues, lnP = compute_protection_factors(protein)
    # lnP.shape == (n_frames, len(residues))

    # ensemble-mean ln(P) per residue
    import numpy as np
    w = np.full(protein.n_frames, 1.0 / protein.n_frames)
    residues, lnP_mean = compute_protection_factors(protein, weights=w)

    # feeding to BME against experimental ln(P)
    from soursop.ssbme import BME, ExperimentalObservable
    obs = [ExperimentalObservable(lnP_exp[k], sigma_lnP[k],
                                  name=f"PF_res{residues[k]}")
           for k in range(len(residues))]
    weights = BME(obs, lnP).fit(theta=2.0, auto_theta=False).weights

**Paramagnetic relaxation enhancement (PRE)** profiles compare the ensemble to an experiment in which a nitroxide spin label at a chosen position relaxes neighbouring amide protons. The intensity ratio I_para/I_dia decays toward 0 for residues near the label and stays near 1 for distant residues::

    from soursop.sspre import SSPRE

    # 600 MHz magnet; tau_c = 5 ns, t_delay = 16 ms, R_2D = 10 Hz
    pre = SSPRE(protein, tau_c=5, t_delay=16, R_2D=10, W_H=600000000)

    # spin label at residue 20 (CB atom)
    intensity_ratio, gamma2 = pre.generate_PRE_profile(label_position=20)

    plt.plot(intensity_ratio)
    plt.axhline(0.5, linestyle='--', color='gray', label='0.5 threshold')
    plt.xlabel('Residue index')
    plt.ylabel(r'$I_\mathrm{para} / I_\mathrm{dia}$')
    plt.title('Simulated PRE profile — label at residue 20')
    plt.legend()
    plt.show()


9. Assessing sampling quality with PENGUIN
---------------------------------------------------------

For disordered proteins it is important to verify that independent replicate simulations have adequately explored conformational space and have not become trapped in local energetic minima. PENGUIN (Lotthammer & Holehouse, *J. Chem. Inf. Model.* 2025) addresses this by comparing per-residue backbone dihedral distributions to the excluded-volume (EV) polymer limit and to one another across replicates.

See the :doc:`../modules/sssampling` page for a full description of the methodology.

**Quickstart — using the precomputed EV reference**::

    from soursop.sssampling import SamplingQuality

    # list of replicate trajectory files (all share the same topology)
    replicates = ['rep0.xtc', 'rep1.xtc', 'rep2.xtc', 'rep3.xtc']

    sq = SamplingQuality(
        traj_list=replicates,
        top_file='topology.pdb',
    )

    # per-residue Hellinger distances: FH vs. precomputed EV reference
    fh_vs_ev = sq.compute_dihedral_hellingers()
    # shape: (n_replicates, n_residues)
    # near 0 → EV-like sampling; near 1 → restricted dihedrals

    # all-to-all inter-replica Hellinger distances
    fh_vs_fh = sq.get_all_to_all_2d_trj_comparison()
    # near 0 → replicates agree; large → replicates are in distinct trapped states

    # four-panel PENGUIN summary figure (saves or displays)
    sq.quality_plot()

**Interpreting the output:**

+----------------------------+--------------------------------+------------------------------------------+
| FH vs. EV distance         | FH vs. FH distance             | Interpretation                           |
+============================+================================+==========================================+
| Near 0                     | Near 0                         | Well-sampled, disordered                 |
+----------------------------+--------------------------------+------------------------------------------+
| Elevated, consistent       | Near 0                         | Locally folded / transient structure     |
+----------------------------+--------------------------------+------------------------------------------+
| Elevated, variable         | Large                          | Energetically trapped — poor sampling    |
+----------------------------+--------------------------------+------------------------------------------+

**Using your own EV trajectories** — if you have run dedicated excluded-volume simulations (e.g. with CAMPARI), supply them explicitly::

    sq = SamplingQuality(
        traj_list=['fh_rep0.xtc', 'fh_rep1.xtc'],
        reference_list=['ev_rep0.xtc', 'ev_rep1.xtc'],
        top_file='fh_topology.pdb',
        ref_top='ev_topology.pdb',
    )
    sq.quality_plot()


10. Reweighting against experiment (BME / iBME)
---------------------------------------------------------

When a simulated ensemble does not quite reproduce an experimental measurement, :mod:`~soursop.ssbme` can compute a new set of per-frame weights that reconcile the two while perturbing the prior ensemble as little as possible. The resulting weights plug straight back into any SOURSOP observable via its ``weights=`` argument. See the :doc:`../modules/bme` page for the theory and pitfalls.

**Standard BME — match an experimental** :math:`R_g` **and** :math:`r_{ee}`. We compute the per-frame observables, define the experimental targets with their uncertainties, fit, and then read back *reweighted* ensemble averages::

    import numpy as np
    from soursop.sstrajectory import SSTrajectory
    from soursop.ssbme import BME, ExperimentalObservable

    traj    = SSTrajectory('traj.xtc', 'start.pdb')
    protein = traj.proteinTrajectoryList[0]

    # per-frame calculated observables -> (n_frames, n_observables)
    rg  = protein.get_radius_of_gyration()
    e2e = protein.get_end_to_end_distance()
    calc = np.column_stack([rg, e2e])

    # experimental values ± uncertainty (same units, here Å)
    obs = [
        ExperimentalObservable(value=23.0, uncertainty=1.0, name="Rg"),
        ExperimentalObservable(value=60.0, uncertainty=2.0, name="Ree"),
    ]

    bme = BME(obs, calc)
    result = bme.fit(theta=2.0, auto_theta=False)
    result.print_diagnostics()          # chi2 before/after, phi, warnings

    w = result.weights

    print(f"Rg  : {np.mean(rg):.2f}  ->  {np.average(rg,  weights=w):.2f} Å")
    print(f"Ree : {np.mean(e2e):.2f} ->  {np.average(e2e, weights=w):.2f} Å")

    # the weights are consistent across *every* SOURSOP observable
    cmap_rew = protein.get_contact_map(weights=w)

**Letting the L-curve choose** :math:`\theta`. Rather than guessing the regularisation strength, scan it and pick the knee::

    scan = bme.scan_theta(theta_range=(0.01, 50.0), n_points=20)
    scan.print_summary()

    result = bme.fit(theta=scan.optimal_theta, auto_theta=False)
    # equivalently: result = bme.fit(auto_theta=True)

**Predicting an independent observable.** A fair quality check is to apply the fitted weights to an observable that was *not* used in the fit::

    asph = protein.get_asphericity()
    print("reweighted asphericity:", np.average(asph, weights=result.weights))

**Iterative BME — data with an unknown scale/offset (e.g. SAXS).** Here each scattering-vector point is one observable, and the calculated intensities differ from the experiment by an unknown global scale and background. ``iBME`` fits those nuisance parameters jointly with the ensemble::

    from soursop.ssbme import iBME, ExperimentalObservable

    # q, I_exp, sigma_exp : experimental SAXS curve (length n_q)
    # calc_I : calculated intensities, shape (n_frames, n_q)
    obs = [ExperimentalObservable(I_exp[k], sigma_exp[k], name=f"q{k}")
           for k in range(len(q))]

    ib = iBME(obs, calc_I)
    result = ib.fit(theta=10.0, ftol=0.01, max_ibme_iterations=50)

    print(f"fitted scale  = {result.scale:.4g}")
    print(f"fitted offset = {result.offset:.4g}")
    print(f"chi2 {result.chi_squared_initial:.2f} -> "
          f"{result.chi_squared_final:.2f}  (phi = {result.phi:.2f})")

    # per-iteration convergence log
    for it in result.ibme_iterations:
        print(it)

    saxs_weights = result.weights       # use with any SOURSOP observable


11. Reweighting with COPER (hard chi-squared constraint)
---------------------------------------------------------

:mod:`~soursop.sscoper` offers an alternative to BME: instead of a tunable penalty it *maximises the ensemble entropy subject to a hard* :math:`\chi^2 \le 1` *constraint* (Leung et al. 2016). There is no :math:`\theta`; the knob is the chi-squared limit, and the method reports whether the data can be satisfied at all. The user-facing API mirrors :mod:`~soursop.ssbme`, so the same patterns apply. See :doc:`../modules/coper` for the theory and the "COPER vs BME" comparison.

**Standard COPER — match an experimental** :math:`R_g` **and** :math:`r_{ee}`::

    import numpy as np
    from soursop.sstrajectory import SSTrajectory
    from soursop.sscoper import COPER, ExperimentalObservable

    traj    = SSTrajectory('traj.xtc', 'start.pdb')
    protein = traj.proteinTrajectoryList[0]

    rg  = protein.get_radius_of_gyration()
    e2e = protein.get_end_to_end_distance()
    calc = np.column_stack([rg, e2e])              # (n_frames, n_observables)

    obs = [
        ExperimentalObservable(value=23.0, uncertainty=1.0, name="Rg"),
        ExperimentalObservable(value=60.0, uncertainty=2.0, name="Ree"),
    ]

    coper  = COPER(obs, calc)
    result = coper.fit(chi2_limit=1.0)
    result.print_diagnostics()                     # chi2, phi, delta_S, warnings

    # ALWAYS check feasibility before using the weights
    if result.feasible:
        w = result.weights
        print(f"Rg  : {np.mean(rg):.2f}  ->  {np.average(rg,  weights=w):.2f} Å")
        print(f"Ree : {np.mean(e2e):.2f} ->  {np.average(e2e, weights=w):.2f} Å")
        cmap_rew = protein.get_contact_map(weights=w)
    else:
        print("Data infeasible: no reweighting reproduces them. "
              "Improve sampling / the force field, or loosen the limit.")

**Per-data-type chi-squared.** When fitting several kinds of data, tag each observable with a ``group`` so COPER constrains each :math:`\chi^2_\alpha` separately (here, "size" vs. "shape")::

    obs = [
        ExperimentalObservable(23.0, 1.0, name="Rg",   group="size"),
        ExperimentalObservable(60.0, 2.0, name="Ree",  group="size"),
        ExperimentalObservable(0.45, 0.05, name="Asph", group="shape"),
    ]
    calc = np.column_stack([rg, e2e, protein.get_asphericity()])
    result = COPER(obs, calc).fit(chi2_limit=1.0)

**Choosing the chi-squared limit** by scanning it (the error-scaling analogue of BME's L-curve)::

    scan = coper.scan_chi2_limit(chi2_limits=(0.25, 4.0), n_points=10)
    scan.print_summary()
    result = coper.fit(chi2_limit=scan.optimal_chi2_limit)

**Iterative COPER — data with an unknown scale/offset (e.g. SAXS)**::

    from soursop.sscoper import iCOPER, ExperimentalObservable

    obs = [ExperimentalObservable(I_exp[k], sigma_exp[k], name=f"q{k}")
           for k in range(len(q))]
    result = iCOPER(obs, calc_I).fit(chi2_limit=1.0, ftol=0.01)

    print(f"fitted scale  = {result.scale:.4g}")
    print(f"fitted offset = {result.offset:.4g}")
    if result.feasible:
        saxs_weights = result.weights
