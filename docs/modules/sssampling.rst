
sssampling
=========================================================

Overview
----------------------------

``sssampling`` implements **PENGUIN** (**P**\ ipeline for **E**\ valuating co\ **N**\ formational hetero\ **G**\ eneity in **U**\ nstructured prote\ **IN**\ s), an algorithmic framework for quantifying per-residue local conformational sampling quality in molecular dynamics or Monte Carlo simulations of intrinsically disordered proteins (IDRs).

The method and its applications are described in full in:

  Lotthammer J.M. & Holehouse A.S. *"Disentangling Folding from Energetic Traps in Simulations of Disordered Proteins."* J. Chem. Inf. Model. (2025). https://doi.org/10.1021/acs.jcim.4c02005


Why sampling assessment is hard for IDRs
..........................................

Folded proteins are naturally assessed by convergence of structural metrics around a reference structure. IDRs have no such reference: they occupy a vast, nearly flat energy landscape populated by many distinct low-energy configurations. Two fundamentally different failure modes can therefore afflict IDR simulations:

* **Energetic trapping / poor sampling** — a subregion becomes stuck in one (or several independent but incompatible) local minima, yielding simulation-to-simulation variation that cannot be attributed to sequence-encoded conformational preference.
* **Local or global folding** — a subregion reproducibly collapses to the same specific conformation across all independent replicas, reflecting a genuine (transient or persistent) structural preference rather than simulation artefact.

Global metrics such as :math:`R_g` or autocorrelation convergence are insensitive to these local failures: a single poorly-sampled subregion is easily masked by a well-sampled majority. PENGUIN therefore operates at residue-level granularity and is designed to distinguish these two scenarios.


The excluded-volume reference model
.....................................

PENGUIN leverages the **excluded-volume (EV) limiting polymer model** as a theoretical upper bound on local conformational heterogeneity. In the EV limit all inter-residue attractive terms are removed from the force field; only steric repulsion, bond geometry, and short-range repulsive potentials are retained. The result is a self-avoiding random coil whose backbone dihedral distributions are maximally broad — effectively a sequence-specific ceiling on what *well-sampled* disordered behaviour looks like.

Comparing dihedral distributions from a full-Hamiltonian (FH) simulation to the EV reference yields an information-theoretic score for each residue. If the FH and EV distributions are similar, the residue is broadly sampling its available dihedral space. If they are very different, the residue has adopted a specific (or limited) set of dihedrals that deviate from random-coil statistics.

SOURSOP provides two routes to the EV reference:

1. **Run your own EV simulations** and supply the resulting trajectory files as ``reference_list`` to :class:`SamplingQuality`.
2. **Use precomputed distributions** (the default). SOURSOP ships with tabulated EV dihedral angle distributions for all 8000 possible amino-acid tripeptide contexts (all 20³ combinations of i-1, i, i+1). The :class:`PrecomputedDihedralInterface` class handles look-up and sampling from these tables, eliminating the need for bespoke EV runs.


The Hellinger distance
.......................

PENGUIN uses the **Hellinger distance** :math:`H(P, Q)` to compare dihedral probability distributions. For two discrete probability vectors :math:`P = (p_1, \ldots, p_k)` and :math:`Q = (q_1, \ldots, q_k)`:

.. math::

    H(P, Q) = \frac{1}{\sqrt{2}} \sqrt{\sum_{i=1}^{k} \left(\sqrt{p_i} - \sqrt{q_i}\right)^2}

The distance is bounded in :math:`[0, 1]`: 0 means the distributions are identical; 1 means they have completely disjoint support. The Hellinger distance is symmetric, handles empty bins naturally, and is fast to compute — important for running many replicate comparisons. A 2D variant (simultaneous :math:`\varphi`/:math:`\psi` comparison) is the default in SOURSOP and is generally preferred because it captures correlated backbone preferences that the marginal 1D distributions can miss.

Two complementary comparisons are performed for each residue:

* **FH vs. EV** — how far are the simulated dihedrals from the maximally heterogeneous reference? Values near 0 indicate EV-like sampling (a necessary and sufficient condition for good sampling); values near 1 indicate restricted dihedrals.
* **FH vs. FH (all-to-all)** — how consistent are the dihedral distributions across independent replicas? Near-zero inter-replica distances indicate reproducible behaviour (either good sampling or consistent folding). Large inter-replica scatter indicates distinct trapped states, i.e., poor sampling.

Together, these two views let you distinguish the three main scenarios:

+------------------------------+--------------------------------------+-------------------------------------------+
| FH vs. EV                    | FH vs. FH (all-to-all)               | Interpretation                            |
+==============================+======================================+===========================================+
| Near 0                       | Near 0                               | Well-sampled, disordered                  |
+------------------------------+--------------------------------------+-------------------------------------------+
| Elevated, consistent         | Near 0                               | Locally folded / transient structure      |
+------------------------------+--------------------------------------+-------------------------------------------+
| Elevated, variable           | Large                                | Energetically trapped (poor sampling)     |
+------------------------------+--------------------------------------+-------------------------------------------+


Typical workflow
.................

Most users will interact with :class:`SamplingQuality`, which accepts a list of replicate trajectory files and (optionally) a list of reference trajectories. If no reference is provided it falls back to the precomputed EV distributions via :class:`PrecomputedDihedralInterface`::

    from soursop.sssampling import SamplingQuality

    # 4 independent replicate simulations, all vs. precomputed EV reference
    sq = SamplingQuality(
        traj_list=['rep0.xtc', 'rep1.xtc', 'rep2.xtc', 'rep3.xtc'],
        top_file='topology.pdb',
    )

    # Per-residue Hellinger distances (FH vs. EV and FH vs. FH all-to-all)
    fh_vs_ev   = sq.compute_dihedral_hellingers()
    fh_vs_fh   = sq.get_all_to_all_2d_trj_comparison()

    # 4-panel PENGUIN summary figure
    sq.quality_plot()

For cases where you want to supply your own EV trajectories (e.g., from a CAMPARI run in the excluded-volume limit)::

    sq = SamplingQuality(
        traj_list=['fh_rep0.xtc', 'fh_rep1.xtc'],
        reference_list=['ev_rep0.xtc', 'ev_rep1.xtc'],
        top_file='topology.pdb',
        ref_top='ev_topology.pdb',
    )

Low-level distance functions (:func:`hellinger_distance`, :func:`compute_joint_hellinger_distance`, :func:`rel_entropy`) are also available for custom analyses.


SamplingQuality Functions
----------------------------

.. autoclass:: soursop.sssampling.SamplingQuality

        .. automethod:: compute_frac_helicity
        .. automethod:: compute_dihedral_hellingers
        .. automethod:: compute_dihedral_rel_entropy
        .. automethod:: compute_series_of_histograms_along_axis
        .. automethod:: compute_pdf
        .. automethod:: get_all_to_all_2d_trj_comparison
        .. automethod:: get_all_to_all_trj_comparisons
        .. automethod:: get_degree_bins
        .. automethod:: quality_plot
        .. automethod:: trj_pdfs
        .. automethod:: ref_pdfs
        .. automethod:: hellingers_distances
        .. automethod:: fractional_helicity
		
		
PrecomputedDihedralInterface Functions
----------------------------------------------

.. autoclass:: soursop.sssampling.PrecomputedDihedralInterface	
	
        .. automethod:: sample_angles
        .. automethod:: gather_phi_reference_dihedrals
        .. automethod:: gather_psi_reference_dihedrals


Stand-alone sssampling functions
----------------------------------------------
		
.. automodule:: soursop.sssampling
.. autofunction:: rel_entropy
.. autofunction:: hellinger_distance
.. autofunction:: compute_joint_hellinger_distance
