Ensemble reweighting (frame weights)
=========================================================

Overview
----------------------------

By default every analysis routine in SOURSOP treats all frames of a
trajectory as contributing equally to an ensemble average. In many
situations — enhanced-sampling simulations, Markov-state-model
reweighting, maximum-entropy / experimentally-restrained reweighting,
or simple importance reweighting — each frame should instead contribute
according to a **statistical weight**.

SOURSOP functions that return an *ensemble-average value* therefore
accept an optional ``weights`` keyword: a per-frame vector that defines
how much each frame contributes to the average (the few exceptions are
listed under `What can be reweighted`_ below). The default,
``weights=False``, is an exact no-op and reproduces the original
unweighted behaviour bit-for-bit.

Such a weight vector can be supplied from any external source (an MSM,
an enhanced-sampling estimator, importance weights), or generated
*within* SOURSOP by reweighting the ensemble against experimental data.
Two complementary maximum-entropy reweighters are provided:

* :class:`~soursop.ssbme.BME` / :class:`~soursop.ssbme.iBME` — Bayesian
  Maximum Entropy, a tunable penalty balancing fit and diversity (see
  :doc:`../modules/bme`); and
* :class:`~soursop.sscoper.COPER` / :class:`~soursop.sscoper.iCOPER` —
  Convex Optimization for Ensemble Reweighting, a hard chi-squared
  constraint with no free regularisation parameter (see
  :doc:`../modules/coper`).

The weights returned by either satisfy the contract below and can be
passed straight into any reweighting-capable method.

Reweighting in SOURSOP is **deterministic**. A weighted average is the
closed-form expectation under the supplied weights — there is no
stochastic resampling, so results are exactly reproducible and do not
depend on a random seed.

The weight vector contract
----------------------------

A valid ``weights`` vector is a *probability vector over frames*:

* one entry per frame (``len(weights) == n_frames``);
* every entry in the closed interval ``[0, 1]``;
* all entries finite (no ``nan`` / ``inf``);
* the entries sum to ``1`` within a small tolerance ``etol``
  (default ``1e-7``).

These conditions are enforced centrally by
:func:`soursop.ssutils.validate_weights`. If any condition fails an
``SSException`` is raised with a descriptive message, rather than
silently producing a meaningless number.

Every reweighting-capable method also exposes an ``etol`` keyword so the
sum-to-one tolerance can be tightened or relaxed.

Interaction with ``stride``
----------------------------

When a method is called with both a frame ``stride`` and ``weights``,
the full weight vector is first validated (including the sum-to-one
check), then subsampled (``weights[::stride]``) and **renormalised** so
it still sums to ``1`` over the retained frames. A
:class:`~soursop.ssexceptions.SoursopWarning` is emitted because
per-stride reweighting is rarely what you want unless the weights were
computed for exactly those frames.

The deterministic helpers
----------------------------

All weighting is funnelled through a small set of shared helpers in
:mod:`soursop.ssutils` so that every method applies weights identically:

* :func:`~soursop.ssutils.validate_weights` — validate (and
  stride-normalise) a weight vector; the single source of truth.
* :func:`~soursop.ssutils.weighted_mean` — ``sum(w * x)``.
* :func:`~soursop.ssutils.weighted_rms` — ``sqrt(sum(w * x**2))`` (the
  polymer-physics RMS order parameter).
* :func:`~soursop.ssutils.weighted_std` — the reliability-weighted
  *population* standard deviation. Frame weights are probability weights
  with no associated sample size, so no ``ddof`` correction is
  well-defined; the population estimator is the unambiguous, reproducible
  choice.
* :func:`~soursop.ssutils.weighted_corr` — the weighted Pearson
  correlation (used by the correlation observables).

Example
----------------------------

.. code-block:: python

    import numpy as np
    from soursop.sstrajectory import SSTrajectory

    TrajOb  = SSTrajectory('traj.xtc', 'start.pdb')
    protein = TrajOb.proteinTrajectoryList[0]

    n = protein.n_frames

    # 1. default: every frame contributes equally
    rg_per_frame = protein.get_radius_of_gyration()        # array, length n

    # 2. uniform weights == the ordinary mean of the per-frame values
    w_uniform = np.full(n, 1.0 / n)
    rg_mean   = protein.get_radius_of_gyration(weights=w_uniform)
    assert np.isclose(rg_mean, rg_per_frame.mean())

    # 3. a one-hot weight isolates a single frame
    w_one      = np.zeros(n); w_one[10] = 1.0
    rg_frame10 = protein.get_radius_of_gyration(weights=w_one)
    assert np.isclose(rg_frame10, rg_per_frame[10])

    # 4. a real reweighting vector (e.g. from MSM / MaxEnt), here from BME
    from soursop.ssbme import BME, ExperimentalObservable
    obs  = [ExperimentalObservable(value=23.0, uncertainty=1.0, name="Rg")]
    calc = protein.get_radius_of_gyration().reshape(-1, 1)   # (n, 1)
    w    = BME(obs, calc).fit(theta=2.0, auto_theta=False).weights
    dmap_w, _ = protein.get_distance_map(weights=w)
    nu, A0    = protein.get_scaling_exponent(weights=w)[:2]

Behaviour at the extremes
----------------------------

* **Uniform weights** (``1/n`` everywhere) reproduce the ordinary
  (unweighted) ensemble mean to within floating-point precision.
* **One-hot weights** (a single ``1.0``) return exactly that frame's
  value.
* **Invalid weights** (wrong length, an element outside ``[0, 1]``, a
  non-finite element, or a sum that differs from ``1`` by more than
  ``etol``, checked before any stride is applied) raise an
  ``SSException``.
* **Degenerate weights** for correlations: a weighted correlation needs
  at least two frames with non-zero weight, so a one-hot weight vector
  passed to ``get_end_to_end_vs_rg_correlation`` or
  ``get_local_to_global_correlation`` raises rather than returning
  ``nan``.

What can be reweighted
----------------------------

Any method that collapses the trajectory to an ensemble value accepts
``weights``. This includes the global dimensions
(``get_radius_of_gyration``, ``get_hydrodynamic_radius``,
``get_asphericity``, ``get_end_to_end_distance``,
``get_gyration_tensor``), the per-frame getters' ensemble means, the
distance/contact maps, ``get_Q`` (with ``protein_average=False``; the
per-frame protein average is not collapsed), the polymer-scaling observables
(``get_internal_scaling`` with ``mean_vals=True``, ``get_internal_scaling_RMS``,
``get_scaling_exponent``, ``get_local_to_global_correlation``), the
SASA summaries (``get_all_SASA``, ``get_site_accessibility``,
``get_regional_SASA``), ``get_angle_decay``, ``get_local_collapse`` (whose
weighted histograms are reported in effective frame counts, so uniform
weights reproduce the unweighted counts), ``get_end_to_end_vs_rg_correlation``,
``get_secondary_structure_DSSP`` (per-residue fractions), the dihedral mutual
information, and the ``SSTrajectory`` overall observables
(``get_overall_radius_of_gyration`` etc.) and
``get_interchain_distance_map``. ``get_interchain_contact_map`` and
``get_interchain_distance`` do not take weights.

A small number of routines describe a *distribution*, a *per-frame*
quantity or a *pairwise frame-vs-frame* quantity rather than a single
ensemble average: ``get_internal_scaling(mean_vals=False)``,
``get_local_heterogeneity``, ``get_Q(protein_average=True)`` and
``get_secondary_structure_DSSP(return_per_frame=True)``. A single
per-frame probability vector is not well-defined for those, so they raise
an ``SSException`` if ``weights`` is supplied (instead of silently
returning a questionable number).

Wherever a standard deviation accompanies a weighted average
(``get_distance_map``, ``get_interchain_distance_map``,
``get_local_collapse``, ``get_site_accessibility``) it is the weighted
population standard deviation (:func:`~soursop.ssutils.weighted_std`).
Every method that takes ``weights`` also takes ``etol``.

API reference
----------------------------

.. autofunction:: soursop.ssutils.validate_weights
.. autofunction:: soursop.ssutils.weighted_mean
.. autofunction:: soursop.ssutils.weighted_rms
.. autofunction:: soursop.ssutils.weighted_std
.. autofunction:: soursop.ssutils.weighted_corr

Shared reweighting primitives
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These helpers in :mod:`soursop.ssutils` are shared by the BME
(:doc:`../modules/bme`) and COPER (:doc:`../modules/coper`) reweighters so
that both expose an identical interface.

.. autoclass:: soursop.ssutils.ExperimentalObservable

   .. automethod:: get_bounds

.. autofunction:: soursop.ssutils.relative_entropy
.. autofunction:: soursop.ssutils.weighted_linear_regression
.. autofunction:: soursop.ssutils.find_optimal_theta
.. autofunction:: soursop.ssutils.build_scan_grid
.. autofunction:: soursop.ssutils.validate_reweighting_inputs
.. autofunction:: soursop.ssutils.constraint_chi_squared
