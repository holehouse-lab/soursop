
sscoper
=========================================================

Overview
----------------------------

``sscoper`` provides **COPER** (Convex Optimization for Ensemble Reweighting; Leung *et al.*, 2016) and its iterative variant **iCOPER** for reweighting conformational ensembles against experimental observables. Like :doc:`bme`, it computes a new set of frame *weights* that brings the ensemble averages into agreement with experiment while perturbing the prior ensemble as little as possible, and the resulting weight vector plugs straight into the consistent SOURSOP reweighting system — every ensemble-average observable accepts a ``weights=`` argument (see :doc:`../usage/weights`).

COPER and BME answer the same question but pose it differently. BME minimises a *penalty* (``½χ² + θ·D_KL``) with a tunable regularisation :math:`\theta`; COPER instead solves a *hard-constrained* problem — maximise the ensemble entropy **subject to** :math:`\chi^2 \le 1` — with no free parameter beyond the χ² limit itself. See `COPER vs BME`_ below.

The module is self-contained (NumPy + SciPy only). Its user-facing API is deliberately analogous to :doc:`bme`, so the same code patterns work for both:

* :class:`~soursop.ssutils.ExperimentalObservable` — a container for one experimental data point (value, uncertainty, constraint type, optional data-type ``group``). This is the **same** class used by BME (shared via :mod:`soursop.ssutils`).
* :class:`~soursop.sscoper.COPER` — standard COPER reweighting.
* :class:`~soursop.sscoper.iCOPER` — iterative COPER, which additionally fits an unknown global scale and offset between calculated and experimental data.

Both classes return a :class:`~soursop.sscoper.COPERResult` and share a :func:`~soursop.sscoper.chi2_limit_scan` helper (the COPER analogue of BME's ``theta_scan``).


The reweighting problem (COPER)
---------------------------------------------------------

A simulation produces an ensemble of :math:`N` conformations with prior weights :math:`w_i^{0}` (uniform unless you supply enhanced-sampling weights). For each conformation we compute :math:`M` observables; the ensemble average of observable :math:`k` is :math:`\langle F_k \rangle = \sum_i w_i F_k(x_i)`, and the disagreement with experiment is measured by the reduced chi-squared

.. math::

   \chi^2(w) = \frac{1}{M}\sum_{k=1}^{M}
       \left(\frac{\langle F_k\rangle - F_k^{\text{exp}}}{\sigma_k}\right)^2 .

COPER seeks the weights that **maximise the entropy** (stay as close as possible to the prior) **subject to fitting the data within error**:

.. math::

   \max_{w}\; S(w) = -\sum_i w_i \ln w_i
   \qquad\text{subject to}\qquad
   \chi^2(w) \le 1,\;\; w_i \ge 0,\;\; \textstyle\sum_i w_i = 1 .

(SOURSOP generalises :math:`S` to the relative entropy :math:`-D_{\mathrm{KL}}(w\|w^0)` so a non-uniform prior is supported; for the usual uniform prior the two differ only by the constant :math:`\ln N`.) Because the entropy is concave and the χ² and simplex constraints are convex, this is a **convex optimisation problem with a unique global solution**, solved here as a *primal* problem over the :math:`N` weights with SciPy's ``trust-constr`` interior-point method, following the two-step recipe of the paper:

#. **χ² minimisation / feasibility.** Minimise :math:`\chi^2(w)` over the simplex. If the minimum satisfies :math:`\chi^2 \le 1`, that point is a feasible interior point for step 2; **if not, the problem has no solution** — no reweighting of the prior ensemble can reproduce the data, a diagnostic in its own right.
#. **Entropy maximisation.** Maximise the entropy subject to :math:`\chi^2 \le 1`. SOURSOP starts this step from the uniform prior — the unconstrained entropy maximum — and tightens toward the constraint, which is better conditioned than starting from the (low-entropy) step-1 point; the problem is convex, so the optimum is the same either way.

The χ² limit need not be exactly 1: scaling it is equivalent to scaling the experimental uncertainties (the paper varies it over ~0.25–4). Tightening it fits the data harder at the cost of ensemble diversity.

Two derived quantities are central to interpreting a COPER fit. The **entropy reduction**

.. math::

   \Delta S = S(w) - S(w^0) = -D_{\mathrm{KL}}(w\|w^0) \le 0

measures the *information content* of the experimental data relative to the prior: a large :math:`|\Delta S|` means the data force a big change. It equals (in magnitude) the mean free-energy change :math:`\langle\Delta G\rangle / kT` needed to reconcile model and experiment. The per-frame **reweighting factor** :math:`r_i = w_i / w_i^{0}` (with :math:`\Delta G_i = -kT\ln r_i`) shows which conformations are up- or down-weighted. As in BME, the fraction of effective frames :math:`\phi = \exp(-D_{\mathrm{KL}}) \in (0,1]` summarises how much diversity survived.


How COPER solves this (and why it scales)
---------------------------------------------------------

**In plain terms.** Think of the weights as a *recipe* for mixing the :math:`N`
conformers. Maximising entropy means keeping the recipe as even as possible
(closest to the prior); the χ² constraint is a *wall* you may not cross (the fit
must be at least this good). COPER simply walks "downhill in unevenness" until it
just touches that wall — and because the landscape is bowl-shaped (convex), there
is one and only one place it can come to rest.

A computer hits two practical snags doing this, and COPER removes both:

#. *The recipe must stay valid.* Every weight has to be non-negative and the
   whole set must sum to one — :math:`N+1` fiddly side-conditions that bog the
   optimiser down. COPER makes them disappear with the **softmax trick**: instead
   of tuning the weights directly, it tunes unconstrained numbers
   :math:`z_i` and *defines* :math:`w_i = e^{z_i}/\sum_j e^{z_j}`. Any choice of
   the :math:`z_i` automatically gives positive weights that sum to one, so the
   optimiser is free to roam.

#. *Taking confident steps needs the local curvature.* To jump straight to the
   bottom of a bowl (rather than inching along), the solver wants the **Hessian** —
   the table of second derivatives describing how the landscape curves. Written
   out, that table is :math:`N\times N`: for 100,000 conformers it would hold
   :math:`10^{10}` numbers (~80 GB) — hopeless to even store. COPER never builds
   it. The curvature here has a very simple structure that we know *exactly*: the
   entropy part is essentially **diagonal** (each conformer curves on its own), and
   the data part is **"low rank"** — it carries only as much complexity as there
   are data points :math:`M`, which is tiny next to :math:`N`. So instead of a
   giant table we keep a recipe for computing "curvature × direction" on the fly,
   costing about :math:`N\times M` work and only :math:`O(N)` memory.

The payoff is that COPER converges in a few dozen steps and fits a
hundred-thousand-conformer ensemble in seconds-to-tens-of-seconds, returning the
*same* maximum-entropy solution it always would.

**In detail.** With :math:`F_{ik} = F_k(x_i)` the calculated value of observable
:math:`k` for frame :math:`i`, the simplex is removed by optimising over
:math:`z\in\mathbb{R}^N` with :math:`w = \mathrm{softmax}(z)`. For any scalar
function :math:`g(w)` the chain rule through the softmax (whose Jacobian is
:math:`J = \mathrm{diag}(w) - w w^{\mathsf T}`) gives the gradient

.. math::

   (\nabla_z g)_i = w_i\big[(\nabla_w g)_i - w\cdot\nabla_w g\big],

and the Hessian-vector product used by the trust-region solver is

.. math::

   H_z\,v = u\circ(g_w - a) \;+\; w\circ(H_w u) \;-\; (b + c)\,w,

where :math:`g_w=\nabla_w g`, :math:`H_w=\nabla_w^2 g`,
:math:`u = J v = w\circ v - w\,(w\cdot v)`, and the scalars
:math:`a = w\cdot g_w`, :math:`b = u\cdot g_w`, :math:`c = w\cdot(H_w u)`
(:math:`\circ` is the elementwise product). Only the *w-space* pieces
:math:`g_w` and :math:`H_w u` are problem-specific:

* **Entropy** :math:`g = D_{\mathrm{KL}}(w\|w^0)`:
  :math:`g_w = \ln(w/w^0) + 1` and :math:`H_w = \mathrm{diag}(1/w)`. The factor of
  :math:`w` in :math:`u` cancels the :math:`1/w`, so :math:`H_w u = v - (w\cdot v)`
  — diagonal-cheap and free of any division by small weights.

* **Per-group χ²** :math:`\chi^2_\alpha` is *quadratic* in :math:`w`, so its
  w-space Hessian is the **constant, rank-≤ :math:`M_\alpha`** matrix

  .. math::

     H_w^{(\alpha)} = \frac{2}{M_\alpha}\sum_{k\in\alpha}
        \frac{f_k\, f_k^{\mathsf T}}{\sigma_k^2}
        = \frac{2}{M_\alpha}\,C_\alpha\,\mathrm{diag}(\sigma_k^{-2})\,C_\alpha^{\mathsf T},

  where :math:`f_k` is the length-:math:`N` column :math:`F_{\cdot k}` and
  :math:`C_\alpha` stacks the active columns. We never form this :math:`N\times N`
  matrix: the product :math:`H_w^{(\alpha)} u = \tfrac{2}{M_\alpha}
  C_\alpha\big(\mathrm{diag}(\sigma_k^{-2})\,(C_\alpha^{\mathsf T}u)\big)` is
  evaluated right-to-left in :math:`O(N M_\alpha)`. (For ``upper`` / ``lower``
  constraints only the currently violated observables are included, exactly as in
  the gradient.)

Both Hessians are handed to ``trust-constr`` as matrix-free
``scipy.sparse.linalg.LinearOperator``\ s, so memory stays :math:`O(N)` and the
solver uses its iterative trust-region subproblem solver — the reason COPER scales
to :math:`\sim 10^5` conformers. None of this changes the optimum: it is the same
convex problem, just posed so it can be solved cheaply.


COPER
----------------------------

:class:`~soursop.sscoper.COPER` performs the reweighting described above. You construct it with a list of observables and a ``(n_frames, n_observables)`` array of calculated values, then call :meth:`~soursop.sscoper.COPER.fit`::

    from soursop.sscoper import COPER, ExperimentalObservable
    import numpy as np

    calc = np.column_stack([rg_per_frame, ree_per_frame])   # (n_frames, 2)

    obs = [
        ExperimentalObservable(value=23.0, uncertainty=1.0, name="Rg"),
        ExperimentalObservable(value=60.0, uncertainty=2.0, name="Ree"),
    ]

    coper = COPER(obs, calc)
    result = coper.fit(chi2_limit=1.0)

    if result.feasible:
        weights = result.weights        # plug into any SOURSOP observable
    result.print_diagnostics()

Constraints (``"equality"`` / ``"upper"`` / ``"lower"``) behave exactly as in :doc:`bme`: an ``upper`` / ``lower`` bound only penalises the disallowed side, so an already-satisfied bound leaves the ensemble essentially untouched.

**Feasibility.** Unlike BME, COPER can report that the data are *infeasible*: if the smallest achievable χ² already exceeds the limit, ``result.feasible`` is ``False``, ``result.success`` is ``False``, and ``result.weights`` hold the closest the prior ensemble can come (the χ²-minimiser). This is a genuine result — it says the data cannot be reproduced by reweighting alone (a sampling or force-field problem) — so always check ``result.feasible`` before using the weights.

**Per-data-type χ².** When you fit several kinds of data at once (e.g. RDCs and J-couplings, as in the paper), a single pooled χ² can let one data type dominate. Assign each observable a ``group`` label and COPER imposes a *separate* :math:`\chi^2_\alpha \le \text{limit}` for every group::

    obs = [
        ExperimentalObservable(-5.1, 0.4, name="RDC1", group="RDC"),
        ExperimentalObservable( 2.3, 0.4, name="RDC2", group="RDC"),
        ExperimentalObservable( 6.0, 0.5, name="J1",   group="Jcoupling"),
    ]
    result = COPER(obs, calc).fit(chi2_limit=1.0)

Observables without a ``group`` are pooled into one default group.


iCOPER — iterative COPER
----------------------------

For data with an **unknown global scale and/or offset** (the canonical case being SAXS, where the calculated intensity curve has an arbitrary scale and a flat background), :class:`~soursop.sscoper.iCOPER` alternates two steps until χ² stabilises, exactly as :class:`soursop.ssbme.iBME` does for BME:

#. a (by default :math:`1/\sigma^2`-weighted) linear regression of the ensemble-averaged calculated values against the experimental values, giving a scale :math:`\alpha` and offset :math:`\beta`, after which ``calc → α·calc + β``; then
#. a standard COPER step on the rescaled data.

The returned :class:`~soursop.sscoper.COPERResult` reports the **net** ``scale`` and ``offset`` relative to the original input plus a per-iteration log (``icoper_iterations``)::

    from soursop.sscoper import iCOPER, ExperimentalObservable

    obs  = [ExperimentalObservable(I_q[k], sigma_q[k], name=f"q{k}")
            for k in range(len(q))]
    result = iCOPER(obs, calc_intensities).fit(chi2_limit=1.0, ftol=0.01)
    print(result.scale, result.offset)
    weights = result.weights

Set ``fit_offset=False`` for a scale-only fit, and ``lr_weights=False`` for an unweighted regression.


Choosing the χ² limit (the error-scaling scan)
---------------------------------------------------------

The χ² limit plays the role BME's :math:`\theta` plays, and is chosen the same way — by scanning it and looking for the knee of the trade-off between fit quality (χ²) and ensemble diversity (relative entropy / :math:`\phi`). :func:`~soursop.sscoper.chi2_limit_scan` (also ``COPER.scan_chi2_limit`` / ``iCOPER.scan_chi2_limit``) sweeps the limit over the paper's error-scaling range and selects the knee::

    scan = coper.scan_chi2_limit(chi2_limits=(0.25, 4.0), n_points=10)
    scan.print_summary()
    print("recommended limit:", scan.optimal_chi2_limit)

    result = coper.fit(chi2_limit=scan.optimal_chi2_limit)


Interpreting the result and common pitfalls
---------------------------------------------------------

Every fit returns a :class:`~soursop.sscoper.COPERResult`; ``result.print_diagnostics()`` produces a formatted report and ``result.diagnostics()`` returns the same as a dict with explicit ``warnings``. Key things to check:

* **Feasibility first.** If ``result.feasible`` is ``False`` the data cannot be fit by any reweighting at this limit. Loosen the χ² limit / uncertainties, or — more likely — improve sampling or the force field. Reweighting cannot create conformations the ensemble never sampled.
* **φ / ΔS is not extreme.** A very low :math:`\phi` (large :math:`|\Delta S|`) means a handful of frames dominate the reweighted ensemble; the result is then effectively a few structures, not an ensemble. A large :math:`\langle\Delta G\rangle/kT` (well beyond a few :math:`kT`) signals a model that disagrees strongly with the data.
* **Two effective sample sizes.** The diagnostics report both the entropy-based :math:`N_{\text{eff}} = N\phi` and the Rényi-2 / participation ratio :math:`1/\sum_i w_i^2`, which is more sensitive to a few dominant weights.
* **Information content of data types.** Because :math:`\Delta S` is well defined, comparing the entropy reduction induced by different observable groups quantifies how much each data type actually constrains the ensemble (a central use in the paper).
* **Realistic uncertainties.** :math:`\sigma_k` must include experimental *and* forward-model error; scaling the χ² limit is exactly equivalent to scaling all :math:`\sigma_k`.
* **Don't validate on the fitting data.** χ² against the fitted observables only improves. Assess quality by predicting *independent* observables with :meth:`~soursop.sscoper.COPER.predict`, or by checking the reproducibility of :math:`\Delta S` across independent subsets of the ensemble (a sampling-density test).


COPER vs BME
----------------------------

Both methods produce the least-biased (maximum-entropy) reweighting consistent with the data; they differ in formulation and in which parameter you tune:

.. list-table::
   :header-rows: 1
   :widths: 22 39 39

   * -
     - **BME** (:doc:`bme`)
     - **COPER** (this page)
   * - Formulation
     - Penalty: minimise :math:`\tfrac12\chi^2 + \theta\,D_{\mathrm{KL}}`
     - Hard constraint: maximise :math:`S` s.t. :math:`\chi^2 \le 1`
   * - Free parameter
     - :math:`\theta` (regularisation)
     - χ² limit (= error scaling)
   * - Solved over
     - :math:`M` Lagrange multipliers (the dual problem; the weights then
       follow in closed form, :math:`w_i \propto w_i^0 e^{-\sum_k \lambda_k F_k(x_i)}`)
     - :math:`N` frame weights directly (the primal problem)
   * - Infeasibility
     - not signalled
     - reported explicitly
   * - Choosing the knob
     - L-curve ``theta_scan``
     - ``chi2_limit_scan``

For most ensembles either method gives similar weights. The difference is computational: BME optimises the :math:`M`-dimensional *dual* problem (one Lagrange multiplier per experimental observable) and recovers the weights afterwards, which is cheaper when there are many frames but few observables (:math:`M \ll N`); COPER optimises the :math:`N` weights directly, but its explicit χ² constraint and feasibility test make the "can these data be fit at all?" question direct.

**Scaling note.** The faithful primal optimisation here uses ``trust-constr`` over the :math:`N` weights; it is convex and fast for ensembles up to ~10\ :sup:`3`–10\ :sup:`4` frames. The original paper used IPOPT and reached ~10\ :sup:`5`; for very large ensembles, cluster or subsample the trajectory before reweighting.


Key references
----------------------------

The implementation here is adapted from the COPER method; please cite the primary reference if you use it:

* Leung, H. T. A., Bignucolo, O., Aregger, R., Dames, S. A., Mazur, A., Bernèche, S. & Grzesiek, S. *A Rigorous and Efficient Method To Reweight Very Large Conformational Ensembles Using Average Experimental Data and To Determine Their Relative Information Content.* J. Chem. Theory Comput. **12**, 383–394 (2016). doi:`10.1021/acs.jctc.5b00759 <https://doi.org/10.1021/acs.jctc.5b00759>`_.

Foundational and closely related work:

* Kullback, S. & Leibler, R. A. *On Information and Sufficiency.* Ann. Math. Stat. **22**, 79–86 (1951) — the relative entropy COPER maximises.
* Bottaro, S., Bengtsen, T. & Lindorff-Larsen, K. *Integrating Molecular Simulation and Experimental Data: A Bayesian/Maximum Entropy Reweighting Approach.* Methods Mol. Biol. **2112**, 219–240 (2020) — the BME / iBME method in :doc:`bme`; COPER differs by using a hard χ² constraint rather than a penalty.
* Hummer, G. & Köfinger, J. *Bayesian ensemble refinement by replica simulations and reweighting.* J. Chem. Phys. **143**, 243150 (2015).
* Różycki, B., Kim, Y. C. & Hummer, G. *SAXS ensemble refinement of ESCRT-III CHMP3 conformational transitions.* Structure **19**, 109–116 (2011) — EROS; the SAXS scale/offset problem that motivates iCOPER.


The ``ExperimentalObservable`` container is documented under
:doc:`../usage/weights` (it is shared with BME).

.. autoclass:: soursop.sscoper.COPER

        .. automethod:: __init__
        .. automethod:: fit
        .. automethod:: scan_chi2_limit
        .. automethod:: predict

.. autoclass:: soursop.sscoper.iCOPER

        .. automethod:: __init__
        .. automethod:: fit
        .. automethod:: scan_chi2_limit
        .. automethod:: predict
        .. automethod:: get_icoper_weights
        .. automethod:: get_icoper_stats

.. autoclass:: soursop.sscoper.COPERResult

        .. automethod:: predict
        .. automethod:: diagnostics
        .. automethod:: print_diagnostics

.. autoclass:: soursop.sscoper.COPERScanResult

        .. automethod:: print_summary

.. autofunction:: soursop.sscoper.chi2_limit_scan
