
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
2. **Use precomputed distributions** (the default). SOURSOP ships with tabulated EV dihedral angle distributions for each of the 20 amino acids in three coarse-grained neighbour contexts (ALA-, LEU- and PRO-like, via ``EV_RESIDUE_MAPPER``): phi is conditioned on the class of the preceding residue (i-1) and psi on the class of the following residue (i+1). The :class:`PrecomputedDihedralInterface` class handles look-up and sampling from these tables, eliminating the need for bespoke EV runs.

The precomputed route resamples each residue's tabulated distribution (4000 angles per context) by inverse-CDF sampling on the histogram grid, drawing as many synthetic frames as there are frames in your trajectories. Each draw picks a bin in proportion to its mass and then a uniform position within that bin, so the resampled histogram reproduces the tabulated one (including the bins at ±180°) up to sampling noise. Because this is a random draw, the resulting Hellinger distances vary slightly from run to run; pass ``seed=<int>`` to :class:`SamplingQuality` to make the EV reference, and everything computed from it, reproducible. The tables only cover the 20 standard amino acids, so a chain containing a non-standard residue raises an :class:`~soursop.ssexceptions.SSException` naming that residue; supply your own EV trajectories in that case.


The Hellinger distance
.......................

PENGUIN uses the **Hellinger distance** :math:`H(P, Q)` to compare dihedral probability distributions. For two discrete probability vectors :math:`P = (p_1, \ldots, p_k)` and :math:`Q = (q_1, \ldots, q_k)`:

.. math::

    H(P, Q) = \frac{1}{\sqrt{2}} \sqrt{\sum_{i=1}^{k} \left(\sqrt{p_i} - \sqrt{q_i}\right)^2}

The distance is bounded in :math:`[0, 1]`: 0 means the distributions are identical; 1 means they have completely disjoint support. The Hellinger distance is symmetric, handles empty bins naturally, and is fast to compute — important for running many replicate comparisons. A 2D variant (simultaneous :math:`\varphi`/:math:`\psi` comparison) is the default in SOURSOP and is generally preferred because it captures correlated backbone preferences that the marginal 1D distributions can miss. One caveat applies on the precomputed route: the tabulated EV distributions are :math:`\varphi` and :math:`\psi` *marginals*, resampled independently, so the synthetic reference has no :math:`\varphi`/:math:`\psi` correlation and any correlation present in the simulated ensemble adds to the 2D distance. Supply your own EV trajectories via ``reference_list`` if the joint reference distribution matters.

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

Every dihedral array is restricted to residues that carry both a phi and a psi (all non-cap residues of a capped chain, the interior residues of an uncapped one); the residue index of each column is available as ``SamplingQuality.residue_indices``.

Most users will interact with :class:`SamplingQuality`, which accepts a list of replicate trajectory files and (optionally) a list of reference trajectories. If no reference is provided it falls back to the precomputed EV distributions via :class:`PrecomputedDihedralInterface`. With more than one trajectory the files are read in parallel with ``multiprocessing``, so in a script the call must sit under an ``if __name__ == "__main__":`` guard (required on macOS and Windows, where worker processes re-import the script); alternatively pass ``force_sequential=True``::

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

Trajectories and references are paired one-to-one, so ``reference_list`` must have exactly as many entries as ``traj_list``; a mismatch raises an :class:`~soursop.ssexceptions.SSException` at construction. For multi-chain topologies, ``proteinID`` selects the chain that every downstream quantity (dihedrals, sequence, fractional helicity, :meth:`~soursop.sssampling.SamplingQuality.quality_plot`) is computed on. Note also that ``dihedral='2D'`` in ``quality_plot`` requires ``method='2D angle distributions'`` (the default), while ``dihedral='phi'`` / ``'psi'`` require ``method='1D angle distributions'``.

Trajectories are stacked into a single array, so every trajectory (and every reference) must contain the same number of frames. If they do not, construction raises an :class:`~soursop.ssexceptions.SSException`; pass ``truncate=True`` to slice everything to the shortest trajectory instead (set ``verbose=False`` to silence the accompanying message). Topologies can be given either as a single path shared by every trajectory or as a list with one entry per trajectory, for both ``top_file`` and ``ref_top``.


Bin width
..........

Every histogram in the pipeline is built on a fixed grid spanning ``[-180, 180]`` degrees. The bin width is set by ``bwidth`` (in radians; the default is ``np.deg2rad(15)``) and must divide 360 degrees into a whole number of bins, so 15°, 10°, 7.5° or 5° are all fine but 7° is rejected with an :class:`~soursop.ssexceptions.SSException` at construction. The grid is laid out with ``np.linspace`` so the outer edges land exactly on ±180° with no rounding of the width. Note that the 2D method builds an ``n_bins x n_bins`` joint histogram per residue and per trajectory, so very fine grids get expensive quickly.


Relative entropy
.................

Alongside the Hellinger distance, :meth:`~soursop.sssampling.SamplingQuality.compute_dihedral_rel_entropy` and :meth:`~soursop.sssampling.SamplingQuality.get_all_to_all_trj_comparisons` (with ``metric='relative entropy'``) report the Kullback-Leibler divergence :math:`D_{KL}(P\,\|\,Q)` between the 1D φ and ψ distributions. KL divergence is undefined (infinite) wherever the second distribution has an empty bin that the first does not, and with 15° bins on a finite ensemble that is the rule rather than the exception. Both methods therefore add a small pseudocount to every bin of both distributions and renormalise before taking the divergence (see :func:`smooth_pdf`); the ``pseudocount`` argument (default ``1e-6``) controls the amount, and ``pseudocount=0`` returns the raw, possibly infinite, value. The Hellinger distance needs no such treatment, which is one reason it is the primary metric.


The summary figure
...................

:meth:`~soursop.sssampling.SamplingQuality.quality_plot` produces the four-panel PENGUIN figure: **A** the per-residue Hellinger distance against the reference (titled for the EV limit or for a supplied reference ensemble, as appropriate) with the across-replica mean; **B** the all-to-all replica comparison; **C** the per-residue spread (max minus min) of the panel A distances across replicas; **D** the per-residue fractional helicity of each replica, with the mean reference helicity drawn as a dashed line when reference trajectories were supplied. All four panels share one residue axis, the 1-based position of each residue among the CA-bearing residues (caps are not counted). The dihedral panels only cover residues that carry both a φ and a ψ, so on an uncapped chain they run from position 2 to :math:`n-1` while the helicity panel runs from 1 to :math:`n` (prior to 2.0.6 both were drawn from position 1, which put the dihedral panels one residue to the left of the helicity panel on uncapped chains); on a capped chain every panel covers positions 1 to :math:`n`. Passing ``save_dir`` writes the figure to ``<save_dir>/<dihedral>_<figname>`` (e.g. ``2D_hellingers.pdf``).

Low-level distance functions (:func:`hellinger_distance`, :func:`compute_joint_hellinger_distance`, :func:`rel_entropy`, :func:`smooth_pdf`) are also available for custom analyses.


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
.. autofunction:: smooth_pdf
.. autofunction:: hellinger_distance
.. autofunction:: compute_joint_hellinger_distance
