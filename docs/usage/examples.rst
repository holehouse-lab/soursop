
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

Both ``get_hydrodynamic_radius`` and ``get_t`` accept ``R1`` / ``R2`` to restrict the calculation to a sub-region. For the hydrodynamic radius both estimators work on the CA-bearing residues of the region (caps never contribute): ``mode='nygaard'`` feeds the Nygaard equation the Rg of their CA atoms and ``N`` equal to their number, as in the paper, and ``mode='kr'`` sums the Kirkwood-Riseman ``1/r_ij`` over every pair of them. Note that Pesce et al. report an ensemble :math:`R_h` as the harmonic mean of the per-frame values, ``1 / np.mean(1 / rh)``. For ``get_t`` the chain length ``N`` is the number of residues in the region (``R2 - R1 + 1``), every residue including caps by default.

**Asphericity** describes how far the chain deviates from a sphere (0 = perfectly spherical, 1 = rod-like)::

    asph = protein.get_asphericity()
    print(f"Mean asphericity = {np.mean(asph):.3f}")

**Correlation between** :math:`r_{ee}` **and** :math:`R_g` — a useful diagnostic for whether the two global size metrics are capturing consistent information. This returns a single Pearson correlation coefficient::

    pearson_r = protein.get_end_to_end_vs_rg_correlation()
    print(f"Pearson r(e2e, Rg) = {pearson_r:.3f}")

**Overlap concentration** :math:`c^*` estimates the concentration above which chains begin to crowd one another::

    c_star = protein.get_overlap_concentration()
    print(f"c* = {c_star:.4f} M")


3. Polymer scaling and internal structure
---------------------------------------------------------

IDR conformational behaviour is often interpreted through the lens of polymer physics. The **internal scaling profile** reports how the inter-residue distance grows with sequence separation :math:`|i - j|`. ``get_internal_scaling`` gives the mean distance :math:`\langle r(i,j) \rangle` in Ångströms; ``get_internal_scaling_RMS`` gives the root-mean-square distance :math:`\sqrt{\langle r^2(i,j) \rangle}`, the order parameter the scaling-exponent fit below uses.

**Internal scaling** (mean across the ensemble)::

    import matplotlib.pyplot as plt

    separation, mean_r = protein.get_internal_scaling(mode='CA', mean_vals=True)
    separation, rms_r = protein.get_internal_scaling_RMS(mode='CA')

    plt.loglog(separation, mean_r, label=r"$\langle r \rangle$")
    plt.loglog(separation, rms_r, label=r"$\sqrt{\langle r^2 \rangle}$")
    plt.xlabel("Sequence separation |i - j|")
    plt.ylabel("Inter-residue distance (Å)")
    plt.legend()
    plt.show()

**Scaling exponent** :math:`\nu` — the Flory exponent extracted by fitting :math:`\sqrt{\langle r^2 \rangle} = A_0\,|i-j|^{\nu}`. ``get_scaling_exponent`` returns a 10-element list; the first two entries are the point estimates ``nu`` and ``A0``, entries 2–5 are the bootstrap confidence-interval bounds on each, and entries 6–7 are the reduced :math:`\chi^2` of the fit (over the log-spaced fit points and over every separation in the fit region, respectively)::

    result = protein.get_scaling_exponent(mode='CA')
    nu, A0 = result[0], result[1]
    nu_lo, nu_hi = result[2], result[3]          # 95% CI on nu (entries 4,5 are the A0 CI)
    redchi = result[7]                           # reduced chi^2 across the fit region

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

**Polymer-scaled distance map** normalises the root-mean-square inter-residue (centre-of-mass) distance matrix by the expected homopolymer scaling (so any ``nu`` / ``A0`` you pass by hand should come from an RMS fit, as returned by ``get_scaling_exponent``), highlighting regions that are more compact or more expanded than a reference random coil. It returns the map together with the ``nu``, ``A0`` and reduced :math:`\chi^2` of the homopolymer fit it performs internally. Note ``mode`` here selects how the deviation is *expressed*, not which atoms are used — the options are ``'fractional-change'`` (the default), ``'signed-fractional-change'``, ``'scaled'``, and ``'signed-absolute-change'``::

    dmap, nu, A0, redchi = protein.get_polymer_scaled_distance_map(
        mode='signed-fractional-change')

    plt.imshow(dmap, origin='lower', cmap='RdBu_r')
    plt.colorbar(label='Fractional deviation from homopolymer')
    plt.xlabel('Residue index')
    plt.ylabel('Residue index')
    plt.title(f'Polymer-scaled distance map (ν = {nu:.3f})')
    plt.show()

If you already know the scaling behaviour you want to normalise against — for example the AFRC or an independently fitted :math:`\nu` — pass it in directly rather than letting the function fit it::

    dmap, nu, A0, redchi = protein.get_polymer_scaled_distance_map(nu=0.588, A0=5.5)

**Local heterogeneity** measures how much the local conformational behaviour varies along the chain, by computing the intra-window RMSD distribution for a sliding window of ``fragment_size`` residues. It returns the per-window mean and standard deviation of that distribution, along with the underlying histograms::

    mean_rmsd, std_rmsd, histograms, bins = protein.get_local_heterogeneity(
        fragment_size=10, stride=1)

    plt.plot(mean_rmsd)
    plt.xlabel('Window start residue')
    plt.ylabel('Mean intra-window RMSD (Å)')
    plt.title('Local conformational heterogeneity')
    plt.show()


4. Secondary structure and backbone angles
---------------------------------------------------------

**DSSP secondary structure** collapses the per-frame DSSP assignment into per-residue fractional occupancies of three buckets — helix, extended (strand), and coil. It returns four arrays: the residue indices covered, and the helix / extended / coil fractions, which sum to 1 at every residue::

    resids, helix, extended, coil = protein.get_secondary_structure_DSSP()
    # each of helix/extended/coil has shape (len(resids),)

    plt.bar(resids, helix)
    plt.xlabel('Residue index')
    plt.ylabel('Fractional helicity')
    plt.title('Per-residue α-helix propensity')
    plt.show()

Pass ``return_per_frame=True`` if you need the raw per-frame assignment instead — the three arrays then become ``(n_frames, len(resids))`` binary masks::

    resids, helix, extended, coil = protein.get_secondary_structure_DSSP(
        return_per_frame=True)

    # fraction of frames helical at each residue (equivalent to the above)
    helicity = helix.mean(axis=0)

**BBSEG backbone-torsion classification** provides a 9-state assignment based on φ/ψ backbone dihedral regions, which is particularly useful for IDRs where the DSSP labels can be sparse or noisy. It returns the residue indices plus a dictionary keyed by class (0–8), where ``bbseg[k]`` is the fractional occupancy of class ``k``::

    resids, bbseg = protein.get_secondary_structure_BBSEG()

    plt.plot(resids, bbseg[4], label='alpha helix (class 4)')
    plt.plot(resids, bbseg[2], label='PPII (class 2)')
    plt.plot(resids, bbseg[1], label='beta (class 1)')
    plt.xlabel('Residue index')
    plt.ylabel('Fractional occupancy')
    plt.legend()
    plt.show()

Only residues with *both* φ and ψ defined are classified, and ``resid_list`` always matches the class arrays element-for-element. On a capped chain (ACE/NME) that is every real residue; on an uncapped chain the two terminal residues are excluded, since φ is undefined at the N-terminus and ψ at the C-terminus.

.. note::

   This changed in SOURSOP 2.0.4. Before that release, uncapped chains returned
   a ``resid_list`` one element longer than the class arrays, and the φ/ψ pair
   classified at a given index was drawn from two *adjacent* residues rather
   than one — so BBSEG assignments on uncapped chains were wrong. Capped chains
   were correct then and are unchanged now.

**Backbone dihedral angles.** ``get_angles`` takes the name of the dihedral to compute (``'phi'``, ``'psi'``, ``'omega'``, or ``'chi1'``–``'chi5'``) and returns a two-element list: the ``mdtraj.Atom`` objects defining each dihedral, and the angles themselves. Note the angle array is ``(n_dihedrals, n_frames)`` — dihedrals first — and is already **in degrees**::

    phi_atoms, phi = protein.get_angles('phi')
    psi_atoms, psi = protein.get_angles('psi')
    # phi.shape == (n_phi, n_frames)

    # Ramachandran plot for the 10th phi/psi pair. Note phi is undefined at the
    # N-terminus and psi at the C-terminus of an uncapped chain, so the two
    # arrays can be offset relative to one another - check the atom lists if
    # you need to map a dihedral back to a specific residue.
    idx = 10
    plt.scatter(phi[idx], psi[idx], alpha=0.3, s=2)
    plt.xlabel('φ (°)')
    plt.ylabel('ψ (°)')
    plt.xlim(-180, 180)
    plt.ylim(-180, 180)
    plt.title(f'Ramachandran plot — dihedral {idx}')
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

**Contact map** — fraction of frames in which each residue pair is within a cutoff distance (``distance_thresh``, default 5 Å). It returns both the map and the per-residue contact order. Valid ``mode`` values are ``'closest-heavy'`` (the default), ``'ca'``, ``'closest'``, ``'sidechain'`` and ``'sidechain-heavy'``::

    contact_map, contact_order = protein.get_contact_map(
        distance_thresh=8.0, mode='ca')

    plt.imshow(contact_map, origin='lower', cmap='hot_r', vmin=0, vmax=1)
    plt.colorbar(label='Contact frequency')
    plt.xlabel('Residue j')
    plt.ylabel('Residue i')
    plt.title('Contact map (CA, 8 Å cutoff)')
    plt.show()

    plt.plot(contact_order)
    plt.xlabel('Residue index')
    plt.ylabel('Mean fractional contacts')
    plt.show()

**RMSD** to a reference frame (here, frame 0). ``frame1`` is the reference; leaving ``frame2`` at its default compares every frame against it::

    rmsd = protein.get_RMSD(frame1=0, region=[0, protein.n_residues - 1])
    print(f"Mean RMSD from frame 0: {np.mean(rmsd):.2f} Å")

    # ... or compare two specific frames
    rmsd_0_to_5 = protein.get_RMSD(frame1=0, frame2=5)


6. Solvent accessibility
---------------------------------------------------------

**Per-residue SASA.** ``get_all_SASA`` returns the *per-frame* array, shape ``(n_frames_strided, n_residues)``, in Å². Because the calculation is expensive the default ``stride`` is 20, so reduce it if your trajectory is short::

    sasa = protein.get_all_SASA(mode='residue', stride=10)
    # sasa.shape == (n_frames_strided, n_residues)

    mean_sasa = sasa.mean(axis=0)
    std_sasa = sasa.std(axis=0)

    plt.bar(range(len(mean_sasa)), mean_sasa, yerr=std_sasa, capsize=2)
    plt.xlabel('Residue index')
    plt.ylabel('SASA (Å²)')
    plt.title('Per-residue solvent accessibility')
    plt.show()

Other granularities are available via ``mode``: ``'atom'`` for per-atom SASA, ``'sidechain'`` / ``'backbone'`` for the per-residue sum over those atoms only, and ``'all'`` for a 3-tuple of ``(residue, sidechain, backbone)`` arrays. In ``'all'`` mode the three arrays are column-aligned over the CA-bearing residues (``protein.resid_with_CA``), so on a capped chain the residue array here excludes the ACE/NME caps that ``mode='residue'`` on its own includes::

    residue_sasa, sidechain_sasa, backbone_sasa = protein.get_all_SASA(
        mode='all', stride=10)
    assert residue_sasa.shape == sidechain_sasa.shape == backbone_sasa.shape

**Regional SASA** for a specific stretch of residues — useful for assessing the accessibility of a functional linear motif. This returns a single number: the sum, over residues ``R1`` to ``R2-1``, of each residue's time-averaged SASA::

    # total SASA of residues 10-19
    region_sasa = protein.get_regional_SASA(R1=10, R2=20, stride=10)
    print(f"Mean total SASA (residues 10–19) = {region_sasa:.1f} Å²")

**Site accessibility** summarises the SASA of chosen residues, returning a dictionary keyed by ``"RESNAME-RESID"`` whose values are ``[mean_SASA, std_SASA]`` in Å². By default the input is interpreted as residue *types*, so every matching residue in the chain is reported — handy for comparing chemically identical residues at different positions::

    accessibility = protein.get_site_accessibility(['TRP', 'TYR'], stride=10)
    for site, (mean_sasa, std_sasa) in accessibility.items():
        print(f"{site}: {mean_sasa:.1f} ± {std_sasa:.1f} Å²")

    # ... or ask for specific residues by index
    accessibility = protein.get_site_accessibility(
        [10, 11, 12], mode='resid', stride=10)

**Coarse-grained ensembles.** One-bead-per-residue models represent every residue as a single ``CA`` bead, so mdtraj's atomic van der Waals radii would give a glycine and a tryptophan the same 1.7 Å radius. For these chains you must name the force field, and SOURSOP uses that model's own per-residue bead sizes (radius = :math:`\sigma/2`) instead::

    protein.is_coarse_grained        # True for a one-bead-per-residue chain

    # pass the model per call ...
    mean_sasa = protein.get_all_SASA(mode='residue', forcefield='mpipi-gg')

    # ... or set it once
    protein.cg_forcefield = 'mpipi-gg'
    mean_sasa = protein.get_all_SASA(mode='residue')

Calling SASA on a coarse-grained chain without a force field raises an ``SSException`` rather than silently returning atomic-radius numbers. See :ref:`cg-sasa` for the list of supported models and the details.


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

    # inter-chain contact map; like the distance map above it is indexed by
    # the CA-bearing residues of each chain (caps excluded), so row k is
    # residue traj.proteinTrajectoryList[0].resid_with_CA[k]
    contact_map = traj.get_interchain_contact_map(
        proteinID1=0, proteinID2=1, threshold=8.0, mode='ca'
    )

.. note::

   Two easy things to trip over here. The contact-map cutoff keyword is
   ``threshold`` (not ``distance_thresh``, which is what the single-chain
   :meth:`~soursop.ssprotein.SSProtein.get_contact_map` uses), and the two
   inter-chain functions spell their ``mode`` values differently:
   ``get_interchain_distance_map`` takes ``'CA'`` / ``'COM'``, while
   ``get_interchain_contact_map`` takes lower-case ``'atom'`` / ``'ca'`` /
   ``'closest'`` / ``'closest-heavy'`` / ``'sidechain'`` /
   ``'sidechain-heavy'``.


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

    # the same call, additionally returning the model's forward-model
    # uncertainty in Hz (a single scalar for the chosen parameterisation)
    atoms, J, sigma = compute_J3_HN_HA(protein, return_uncertainty=True)

    # ensemble mean per dihedral
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

    import numpy as np
    from soursop.sspre import SSPRE

    # 600 MHz magnet; tau_c = 5 ns, t_delay = 16 ms, R_2D = 10 Hz.
    # NOTE W_H is the ANGULAR proton Larmor frequency, omega_H = 2*pi*nu_H
    # (rad/s) - so a 600 MHz magnet is 2*pi*600e6, not 600e6.
    pre = SSPRE(protein, tau_c=5, t_delay=16, R_2D=10, W_H=2 * np.pi * 600e6)

    # spin label at residue 20; since SOURSOP 2.0.2 this uses the calibrated
    # coarse-grained spin-label cloud model by default (pass use_label=False
    # for the classic point-at-CB behaviour of SOURSOP <= 2.0.1)
    intensity_ratio, gamma2 = pre.generate_PRE_profile(label_position=20)

    plt.plot(intensity_ratio)
    plt.axhline(0.5, linestyle='--', color='gray', label='0.5 threshold')
    plt.xlabel('Residue index')
    plt.ylabel(r'$I_\mathrm{para} / I_\mathrm{dia}$')
    plt.title('Simulated PRE profile — label at residue 20')
    plt.legend()
    plt.show()

Since SOURSOP 2.0.2 the paramagnetic centre is, by default, a coarse-grained cloud of beads offset from the anchor (calibrated against DEER-PREdict), which better reproduces the geometry of an MTSL side chain and also works on coarse-grained CA-only trajectories. This is a breaking change: earlier versions placed the label directly on the CB atom, so pass ``use_label=False`` to reproduce older results exactly::

    # classic point-at-CB model (SOURSOP <= 2.0.1 behaviour)
    intensity_ratio, gamma2 = pre.generate_PRE_profile(label_position=20, use_label=False)


9. Assessing sampling quality with PENGUIN
---------------------------------------------------------

For disordered proteins it is important to verify that independent replicate simulations have adequately explored conformational space and have not become trapped in local energetic minima. PENGUIN (Lotthammer & Holehouse, *J. Chem. Inf. Model.* 2025) addresses this by comparing per-residue backbone dihedral distributions to the excluded-volume (EV) polymer limit and to one another across replicates.

See the :doc:`../modules/sssampling` page for a full description of the methodology.

**Quickstart — using the precomputed EV reference** (with several trajectories the files are loaded in parallel, so in a script put this under ``if __name__ == "__main__":``, or pass ``force_sequential=True``)::

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
