
ssbme
=========================================================

Overview
----------------------------

``ssbme`` provides **Bayesian Maximum Entropy (BME)** and **iterative BME (iBME)** reweighting of conformational ensembles against experimental observables. Given a set of per-frame calculated quantities and the corresponding experimental measurements, these methods compute a new set of frame *weights* that brings the ensemble averages into agreement with experiment while perturbing the prior ensemble as little as possible. The resulting weight vector plugs directly into the consistent SOURSOP reweighting system — every ensemble-average observable accepts a ``weights=`` argument (see :doc:`../usage/weights`).

The module is self-contained (NumPy + SciPy only) and is adapted from the reference BME implementation of Bottaro, Bengtsen & Lindorff-Larsen with the clean array-based API used in STARLING.

It exposes three user-facing pieces:

* :class:`~soursop.ssutils.ExperimentalObservable` — a container for one experimental data point (value, uncertainty, constraint type). This class lives in :mod:`soursop.ssutils` and is shared with the COPER reweighter (:doc:`coper`); it is re-exported from ``soursop.ssbme`` so ``from soursop.ssbme import ExperimentalObservable`` continues to work.
* :class:`~soursop.ssbme.BME` — standard maximum-entropy reweighting.
* :class:`~soursop.ssbme.iBME` — iterative BME, which additionally fits an unknown global scale and offset between calculated and experimental data.

Both classes return a :class:`~soursop.ssbme.BMEResult` and share an L-curve :func:`~soursop.ssbme.theta_scan` helper for selecting the regularisation strength.

SOURSOP also provides :doc:`COPER <coper>`, an alternative maximum-entropy reweighter that imposes a *hard* :math:`\chi^2 \le 1` constraint instead of BME's tunable penalty. The two share the same ``ExperimentalObservable`` interface; see the "COPER vs BME" comparison on the :doc:`coper` page.


The reweighting problem
----------------------------

A simulation produces an ensemble of :math:`N` conformations with prior weights :math:`w_i^{0}` (uniform unless you supply enhanced-sampling weights). For each conformation we can compute :math:`M` observables (e.g. a radius of gyration, a set of NOE distances, a SAXS intensity curve). The ensemble average of observable :math:`k` is

.. math::

   \langle F_k \rangle = \sum_{i=1}^{N} w_i \, F_k(x_i).

In general these averages do not match the experimental values :math:`F_k^{\text{exp}} \pm \sigma_k`, either because the force field is imperfect or because the ensemble is undersampled. **Reweighting** asks: what is the *minimal* modification of the weights :math:`w_i^{0} \rightarrow w_i` that reconciles the ensemble with the data?

BME answers this by minimising the functional

.. math::

   \mathcal{L}(w) = \tfrac{1}{2}\chi^2(w) \;+\; \theta \, D_{\mathrm{KL}}(w \,\|\, w^{0}),

where :math:`\chi^2` measures the disagreement with experiment and :math:`D_{\mathrm{KL}}(w\|w^{0}) = \sum_i w_i \ln (w_i / w_i^{0})` is the relative entropy (information loss) between the new and prior ensembles. The hyperparameter :math:`\theta` controls the trade-off: large :math:`\theta` keeps the ensemble close to the prior (trusts the simulation), small :math:`\theta` fits the data aggressively (trusts the experiment). The maximum-entropy solution has the closed exponential form

.. math::

   w_i \;\propto\; w_i^{0}\,\exp\!\Big(\!-\!\sum_{k}\lambda_k F_k(x_i)\Big),

so in practice only the :math:`M` Lagrange multipliers :math:`\lambda_k` are optimised (via ``scipy.optimize.minimize`` with ``L-BFGS-B``), not the :math:`N` weights directly. This makes BME efficient even for very large ensembles.

A key diagnostic is the **fraction of effective frames**

.. math::

   \phi \;=\; \exp\!\big(\!-\!D_{\mathrm{KL}}(w\,\|\,w^{0})\big) \in (0, 1],

reported on every result. :math:`\phi \approx 1` means the data were already well reproduced (little reweighting); :math:`\phi \ll 1` means a small subset of frames dominates the reweighted ensemble and the result should be treated with caution.


BME
----------------------------

:class:`~soursop.ssbme.BME` performs the standard reweighting described above. Use it when the calculated observables are already on the **same scale** as the experimental values (e.g. an :math:`R_g` in Å compared to a SAXS-derived :math:`R_g` in Å, or NOE-derived distances).

You construct it with a list of observables and a ``(n_frames, n_observables)`` array of calculated values, then call :meth:`~soursop.ssbme.BME.fit`. ``fit`` either uses a manually supplied ``theta`` or, with ``auto_theta=True``, runs an L-curve scan to choose one automatically::

    from soursop.ssbme import BME, ExperimentalObservable
    import numpy as np

    # per-frame calculated values: shape (n_frames, n_observables)
    calc = np.column_stack([rg_per_frame, ree_per_frame])

    obs = [
        ExperimentalObservable(value=23.0, uncertainty=1.0, name="Rg"),
        ExperimentalObservable(value=60.0, uncertainty=2.0, name="Ree"),
    ]

    bme = BME(obs, calc)
    result = bme.fit(theta=2.0, auto_theta=False)

    print(result)                       # chi2 before/after, phi, theta
    weights = result.weights            # plug into any SOURSOP observable

Constraints. Each observable carries a ``constraint``:

* ``"equality"`` (default) — the observable should match ``value`` within ``uncertainty`` (a two-sided restraint).
* ``"upper"`` — the observable should not exceed ``value``; being below is not penalised. Useful for e.g. an experimentally determined upper bound on a distance.
* ``"lower"`` — the observable should not fall below ``value``.

Inequality constraints are enforced through bounds on the corresponding Lagrange multiplier, so an already-satisfied bound leaves the ensemble essentially untouched (:math:`\phi \approx 1`).


iBME — iterative BME
----------------------------

For several important experimental techniques the calculated and experimental quantities differ by an **unknown global scale and/or offset**. The canonical example is small-angle X-ray scattering (SAXS): the shape of the calculated intensity curve is meaningful, but its absolute scale (and a flat background offset) are nuisance parameters that must be fitted simultaneously with the ensemble.

:class:`~soursop.ssbme.iBME` (ported from the reference ``Reweight.ibme``) solves this by alternating two steps until the :math:`\chi^2` stabilises:

#. **Scale/offset fit** — perform a (by default :math:`1/\sigma^2`-weighted) linear regression of the current ensemble-averaged calculated values against the experimental values, obtaining a scale :math:`\alpha` and offset :math:`\beta`, and rescale ``calc → α·calc + β``.
#. **Maxent step** — run a standard BME fit on the rescaled data and adopt the new weights.

The scale/offset accumulates across iterations; the returned :class:`~soursop.ssbme.BMEResult` reports the **net** ``scale`` and ``offset`` relative to the original input, together with a per-iteration log (``ibme_iterations``)::

    from soursop.ssbme import iBME, ExperimentalObservable
    import numpy as np

    # e.g. SAXS: one observable per scattering-vector (q) point
    obs = [ExperimentalObservable(I_q[k], sigma_q[k], name=f"q{k}")
           for k in range(len(q))]
    calc = calculated_intensities          # shape (n_frames, len(q))

    ib = iBME(obs, calc)
    result = ib.fit(theta=10.0, ftol=0.01, max_ibme_iterations=50)

    print(result.scale, result.offset)     # fitted nuisance parameters
    print(result)                          # chi2, phi, n_iterations
    weights = result.weights

Convergence is declared when the change in :math:`\chi^2` between iterations drops below ``ftol`` (or ``max_ibme_iterations`` is reached). Set ``fit_offset=False`` to fit a scale only (no additive background), and ``lr_weights=False`` to use an unweighted regression instead of the default :math:`1/\sigma^2` weighting.

**When to use which.** Use ``BME`` when calculated and experimental observables are directly comparable. Use ``iBME`` when there is an unknown multiplicative scale and/or additive offset linking them (SAXS being the textbook case). If in doubt, ``iBME`` with ``fit_offset=True`` reduces to ``BME`` when the true scale is 1 and offset is 0, at the cost of a few extra iterations.


BMECustom — vector / matrix BME with a pluggable cost
---------------------------------------------------------

For *profile* data (a SAXS curve, a PRE intensity profile, an NMR chemical-shift vector — anything where the experiment is a length-:math:`m` vector and you can back-calculate the same observable per conformer) and for cases where you want a **non-Gaussian goodness-of-fit** (log-scale :math:`\chi^2`, a covariance :math:`\chi^2`, a Pearson-correlation cost, ...), :class:`~soursop.ssbme.BMECustom` reweights the ensemble against:

* an experimental **vector** ``experiment`` of shape ``(m,)``,
* a per-conformer **matrix** ``calculated_values`` of shape ``(n_frames, m)``,
* an optional **uncertainty** scalar or length-``m`` vector (used by the default cost), and
* an optional **cost function** ``cost(experiment, calculated_values, weights) -> float`` (lower = better fit).

It is the **penalty form** of Bayesian Ensemble Refinement: rather than BME's closed-form exponential solution (special to the Gaussian :math:`\chi^2`), it optimises the weights directly to minimise

.. math::

   \mathcal{L}(w) = \text{cost}(w) \;+\; \theta\, D_{\mathrm{KL}}(w\,\|\,w^{0}),

over the simplex (:math:`w_i \ge 0`, :math:`\sum_i w_i = 1`). Internally the simplex is enforced via a softmax reparameterisation (:math:`w = \mathrm{softmax}(z)`), which reduces the fit to an unconstrained minimisation in :math:`z \in \mathbb{R}^{n}` solved by SciPy's L-BFGS-B (and avoids boundary singularities in :math:`\log w`). With no ``cost_function`` it reproduces :class:`BME`'s reduced :math:`\chi^2` exactly; with a custom one it is the natural way to plug in any scalar fit metric.

The cost function **must take the weights** — reweighting is the act of choosing them, and the cost can only be reweighting-sensitive through the weighted prediction :math:`\langle F_k \rangle = \sum_i w_i\, \mathrm{calc}[i, k]`::

    import numpy as np
    from soursop.ssbme import BMECustom

    # SAXS-like: m points, n conformers
    bme = BMECustom(I_exp, calc_I, uncertainty=sigma_exp)
    result = bme.fit(theta=1.0)
    weights = result.weights

    # custom cost — chi-squared on a log scale
    def log_chi2(experiment, calculated_values, weights):
        avg = weights @ calculated_values
        return float(np.mean((np.log(avg) - np.log(experiment)) ** 2))

    result = BMECustom(I_exp, calc_I, cost_function=log_chi2).fit(theta=1.0)

The result is a :class:`~soursop.ssbme.BMECustomResult` carrying ``cost_initial`` / ``cost_final`` in place of :math:`\chi^2`, ``phi``, the per-frame ``reweighting_factors`` :math:`r_i = w_i / w_i^{0}`, a ``predict`` for new observables, and the usual ``diagnostics`` / ``print_diagnostics``. :meth:`BMECustom.scan_theta <soursop.ssbme.BMECustom.scan_theta>` runs the same L-curve trade-off using the cost in place of :math:`\chi^2`.

**Performance note.** The default :math:`\chi^2` path uses an analytic gradient; a custom ``cost_function`` triggers a finite-difference gradient (one cost evaluation per weight), so an expensive cost on a very large ensemble can be slow. Convexity (and hence a unique optimum) is guaranteed only when the cost is convex in ``w``; an arbitrary cost is optimised locally.


Choosing theta (the L-curve scan)
---------------------------------------------------------

:math:`\theta` is the single most important choice and should **not** be left at an arbitrary value. The principled approach is an L-curve scan: fit over a grid of :math:`\theta`, then plot :math:`\chi^2` against the relative entropy (or :math:`\phi`). The curve is L-shaped — there is a "knee" beyond which fitting the data better requires a disproportionate loss of ensemble diversity. :func:`~soursop.ssbme.theta_scan` (also available as ``BME.scan_theta`` / ``iBME.scan_theta``) automates this and selects the knee via a perpendicular-distance or Menger-curvature criterion::

    scan = bme.scan_theta(theta_range=(0.01, 50.0), n_points=20)
    scan.print_summary()
    print("recommended theta:", scan.optimal_theta)

    # refit at the recommended value
    result = bme.fit(theta=scan.optimal_theta, auto_theta=False)

Calling ``bme.fit(auto_theta=True)`` (the default when no ``theta`` is given) performs the scan and returns the knee result in one step.


Interpreting the result and common pitfalls
---------------------------------------------------------

Every fit returns a :class:`~soursop.ssbme.BMEResult`. Inspect it before trusting the weights — ``result.print_diagnostics()`` produces a formatted report, and ``result.diagnostics()`` returns the same information as a dict with explicit ``warnings``. Key things to check:

* **χ² actually decreased.** ``chi_squared_initial`` vs. ``chi_squared_final``. If it barely moved, the ensemble may already agree with the data, or the optimisation failed (check ``result.success``).
* **φ is not tiny.** A very low :math:`\phi` (e.g. :math:`< 0.1`) means a handful of frames carry essentially all the weight; the "reweighted ensemble" is then statistically a few structures, not an ensemble. This usually indicates the prior ensemble does not contain conformations consistent with the data — reweighting cannot create structures that were never sampled. Prefer increasing :math:`\theta`, loosening uncertainties, or improving sampling.
* **Two effective sample sizes.** The diagnostics report both the entropy-based :math:`N_{\text{eff}} = N\phi` (the usual BME measure) and the Rényi-2 / participation ratio :math:`1/\sum_i w_i^2`, which is more sensitive to a few dominant weights.
* **Realistic uncertainties.** :math:`\sigma_k` must include experimental error *and* forward-model error. Uncertainties that are too tight force aggressive, low-:math:`\phi` reweighting; too loose and nothing happens.
* **Don't validate on the fitting data.** :math:`\chi^2` against the observables used for reweighting will always improve. Assess quality with cross-validation (the L-curve knee) or by predicting *independent* observables with :meth:`~soursop.ssbme.BME.predict`.
* **Garbage in, garbage out.** Reweighting only redistributes weight among conformations that were already sampled. It is not a substitute for adequate sampling (see PENGUIN in :doc:`sssampling`) or a reasonable force field.


Key references
----------------------------

The implementation here is adapted from the reference BME software and protocol; please cite the primary reference if you use it:

* Bottaro, S., Bengtsen, T. & Lindorff-Larsen, K. *Integrating Molecular Simulation and Experimental Data: A Bayesian/Maximum Entropy Reweighting Approach.* Methods Mol. Biol. **2112**, 219–240 (2020). doi:`10.1007/978-1-0716-0270-6_15 <https://doi.org/10.1007/978-1-0716-0270-6_15>`_ (preprint: `bioRxiv 2018 <https://doi.org/10.1101/457952>`_). Reference code: `github.com/KULL-Centre/BME <https://github.com/KULL-Centre/BME>`_.

Foundational and closely related methodology:

* Hummer, G. & Köfinger, J. *Bayesian ensemble refinement by replica simulations and reweighting.* J. Chem. Phys. **143**, 243150 (2015) — the BioEn reweighting functional to which BME is mathematically equivalent.
* Cesari, A., Gil-Ley, A. & Bussi, G. *Combining simulations and solution experiments as a paradigm for RNA force field refinement.* J. Chem. Theory Comput. **12**, 6192–6200 (2016) — the maximum-entropy-with-error objective :math:`\Gamma(\lambda)` minimised by this module.
* Cesari, A., Reißer, S. & Bussi, G. *Using the maximum entropy principle to combine simulations and solution experiments.* Computation **6**, 15 (2018) — a review of the method and its failure modes.
* Różycki, B., Kim, Y. C. & Hummer, G. *SAXS ensemble refinement of ESCRT-III CHMP3 conformational transitions.* Structure **19**, 109–116 (2011) — EROS; origin of the global scaling parameter :math:`\theta` and the SAXS scale/offset problem that motivates iBME.
* Bottaro, S. & Lindorff-Larsen, K. *Biophysical experiments and biomolecular simulations: A perfect match?* Science **361**, 355–360 (2018) — perspective on integrating simulation and experiment.


The ``ExperimentalObservable`` container is documented under
:doc:`../usage/weights` (it is shared with COPER).

.. autoclass:: soursop.ssbme.BME

        .. automethod:: __init__
        .. automethod:: fit
        .. automethod:: scan_theta
        .. automethod:: predict

.. autoclass:: soursop.ssbme.iBME

        .. automethod:: __init__
        .. automethod:: fit
        .. automethod:: scan_theta
        .. automethod:: predict
        .. automethod:: get_ibme_weights
        .. automethod:: get_ibme_stats

.. autoclass:: soursop.ssbme.BMECustom

        .. automethod:: __init__
        .. automethod:: fit
        .. automethod:: scan_theta
        .. automethod:: predict

.. autoclass:: soursop.ssbme.BMEResult

        .. automethod:: predict
        .. automethod:: diagnostics
        .. automethod:: print_diagnostics

.. autoclass:: soursop.ssbme.BMECustomResult

        .. automethod:: predict
        .. automethod:: diagnostics
        .. automethod:: print_diagnostics

.. autoclass:: soursop.ssbme.ThetaScanResult

        .. automethod:: print_summary

.. autofunction:: soursop.ssbme.theta_scan
