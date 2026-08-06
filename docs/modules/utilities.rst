
Utilities and supporting modules
=========================================================

Overview
----------------------------

Alongside the main analysis modules, SOURSOP ships a handful of smaller modules that hold shared helpers, reference data, and the package's exception types. Most users will never need to import these directly — they underpin the analysis routines rather than replacing them — but they are documented here because several are genuinely useful on their own, and because knowing what lives where makes the package easier to extend.

* ``ssdata`` — reference data tables: residue-name mappings, the coarse-grained force-field bead sizes, and the precomputed excluded-volume dihedral distributions.
* ``sstools`` — small numerical and file-system helpers shared across the package.
* ``sspolymer`` — polymer-physics relations that do not need a trajectory.
* ``ssmutualinformation`` — the mutual-information machinery behind ``SSProtein.get_dihedral_mutual_information``.
* ``ssutils`` — validation and reduction helpers, including the package-wide ``weights`` validator. The reweighting-specific parts are documented on the :doc:`../usage/weights` page.
* ``ssexceptions`` — the exception and warning types SOURSOP raises.


ssdata — reference data
----------------------------

``ssdata`` holds the static lookup tables SOURSOP needs. The most useful of these for end users are the coarse-grained bead-size tables introduced for the coarse-grained SASA calculation (see :ref:`cg-sasa`).

**Coarse-grained force-field bead sizes.** SOURSOP bundles per-residue bead diameters (sigma) for nine one-bead-per-residue force fields, taken from each model's own LAMMPS parameter file. ``CG_FORCEFIELDS`` is derived by listing the bundled data directory, so it always matches what is actually shipped::

    from soursop.ssdata import CG_FORCEFIELDS, get_cg_bead_sigmas

    print(CG_FORCEFIELDS)
    # ('cation-pi-1', 'cation-pi-2', 'fb-hps', 'hps-kr', 'hps-urry',
    #  'kh', 'mpipi', 'mpipi-gg', 'mpipi-recharged')

    sigmas = get_cg_bead_sigmas('mpipi')
    print(sigmas['GLY'], sigmas['TRP'])       # 4.69511  7.066548

Sigma is the model's diagonal (i,i) pair sigma — that is, the bead **diameter** — in Angstroms, so a Shrake–Rupley bead radius is ``sigma/2``. The Mpipi variants additionally define the RNA beads ``RPA`` / ``RPC`` / ``RPG`` / ``RPU``. Passing ``buried=True`` returns the sigmas of bead types 21–40, the roughly 30% weaker copies these models define for residues sequestered inside a folded domain::

    buried = get_cg_bead_sigmas('mpipi', buried=True)

Note the six HPS-family models (``hps-urry``, ``hps-kr``, ``fb-hps``, ``kh``, ``cation-pi-1``, ``cation-pi-2``) all share the same Kim–Hummer diameters, so they return identical tables; Mpipi and its variants do not.

**Other tables.** ``ssdata`` also exposes ``THREE_TO_ONE`` / ``ONE_TO_THREE`` residue-name maps, ``ALL_VALID_RESIDUE_NAMES`` (the residue names SOURSOP recognises as protein during chain detection), ``DEFAULT_SIDECHAIN_VECTOR_ATOMS`` (the per-residue sidechain-tip atom used for sidechain-vector analyses), and the precomputed excluded-volume φ/ψ tripeptide distributions used by :doc:`sssampling`.

.. autofunction:: soursop.ssdata.get_cg_bead_sigmas


sstools — general helpers
----------------------------

Small, stateless helpers used across the package.

:func:`~soursop.sstools.find_trajectory_files` is the most immediately useful of these — it walks a directory tree and returns matched trajectory/topology pairs for a set of numbered replicate directories, which is exactly the input :class:`~soursop.sssampling.SamplingQuality` and :func:`~soursop.sstrajectory.parallel_load_trjs` expect::

    from soursop.sstools import find_trajectory_files

    trj_files, top_files = find_trajectory_files(
        root_dir='simulations/', num_replicates=5)

:func:`~soursop.sstools.powermodel` evaluates the polymer power law :math:`R_0 |i-j|^{\nu}` and is handy for plotting a fitted scaling law alongside an internal-scaling profile::

    import numpy as np
    from soursop.sstools import powermodel

    nu, A0 = protein.get_scaling_exponent(mode='CA')[:2]
    seps = np.arange(1, protein.n_residues)
    plt.loglog(seps, powermodel(seps, nu, A0), '--', label=f'fit (ν={nu:.2f})')

.. autofunction:: soursop.sstools.find_trajectory_files
.. autofunction:: soursop.sstools.powermodel
.. autofunction:: soursop.sstools.get_distance_periodic
.. autofunction:: soursop.sstools.find_nearest
.. autofunction:: soursop.sstools.fix_histadine_name
.. autofunction:: soursop.sstools.chunks


sspolymer — polymer relations
--------------------------------

Polymer-physics relations that operate on numbers rather than trajectories. :func:`~soursop.sspolymer.get_overlap_concentration` converts a radius of gyration into the overlap concentration :math:`c^*` — the concentration at which chains begin to crowd one another — and is what :meth:`SSProtein.get_overlap_concentration <soursop.ssprotein.SSProtein.get_overlap_concentration>` calls internally. Use it directly when you have an :math:`R_g` from elsewhere (an experiment, or the AFRC)::

    from soursop.sspolymer import get_overlap_concentration

    c_star = get_overlap_concentration(25.0)      # Rg in Angstroms
    print(f"c* = {c_star:.4f} mg/mL")

.. autofunction:: soursop.sspolymer.get_overlap_concentration


ssmutualinformation — mutual information
--------------------------------------------

The information-theoretic machinery behind :meth:`SSProtein.get_dihedral_mutual_information <soursop.ssprotein.SSProtein.get_dihedral_mutual_information>`. :func:`~soursop.ssmutualinformation.calc_MI` computes the mutual information :math:`I(X; Y)` between two observables by histogramming them, and honours the package-wide ``weights`` contract, so it can be applied to a reweighted ensemble::

    import numpy as np
    from soursop.ssmutualinformation import calc_MI

    _, phi = protein.get_angles('phi')
    mi = calc_MI(phi[10], phi[11], bins=np.arange(-180, 190, 10))

Pass ``normalize=True`` for the normalised mutual information (bounded in ``[0, 1]``), which is easier to compare across residue pairs.

.. autofunction:: soursop.ssmutualinformation.calc_MI
.. autofunction:: soursop.ssmutualinformation.shan_entropy


ssutils — validation and reduction helpers
----------------------------------------------

The ``weights`` validator and the deterministic weighted reducers (``validate_weights``, ``weighted_mean``, ``weighted_std``, ``weighted_rms``, ``weighted_corr``) plus the shared reweighting primitives (``ExperimentalObservable`` and friends) are documented on the :doc:`../usage/weights` page, since that is where they are explained in context.

The remaining helpers are listed here.

:func:`~soursop.ssutils.validate_keyword_option` is the validator SOURSOP uses for every string ``mode`` / ``scheme`` argument, and is worth reusing if you are writing a plugin (see :doc:`../usage/development`) so your errors look like the rest of the package::

    from soursop.ssutils import validate_keyword_option

    validate_keyword_option(mode, ['CA', 'COM'], 'mode')

:func:`~soursop.ssutils.set_numpy_threads` caps the number of threads NumPy's BLAS backend uses, which matters when running many SOURSOP analyses in parallel and BLAS would otherwise oversubscribe the machine::

    from soursop.ssutils import set_numpy_threads

    set_numpy_threads(1)      # one thread per worker process

:func:`~soursop.ssutils.kabsch_rmsd` gives the minimal RMSD between two point sets after optimal superposition, and :func:`~soursop.ssutils.ideal_helix_ca` / :func:`~soursop.ssutils.ideal_extended_ca` generate idealised CA geometries — these underpin the two-bead (CA/CB) secondary-structure assignment.

.. autofunction:: soursop.ssutils.validate_keyword_option
.. autofunction:: soursop.ssutils.set_numpy_threads
.. autofunction:: soursop.ssutils.kabsch_rmsd
.. autofunction:: soursop.ssutils.ideal_helix_ca
.. autofunction:: soursop.ssutils.ideal_extended_ca
.. autofunction:: soursop.ssutils.is_swan_topology


ssexceptions — errors and warnings
--------------------------------------

SOURSOP raises a single exception type, :class:`~soursop.ssexceptions.SSException`, for every user-facing error — invalid keyword options, malformed ``weights`` vectors, residue indices out of range, analyses that are undefined for the resolution of the chain, and so on. Catching it is therefore enough to catch anything SOURSOP itself raises::

    from soursop.ssexceptions import SSException

    try:
        sasa = protein.get_all_SASA()
    except SSException as e:
        print(f"SOURSOP could not compute this: {e}")

Note that errors originating in mdtraj, NumPy, or SciPy propagate unchanged and are *not* wrapped in an ``SSException``.

:func:`~soursop.ssexceptions.SSWarning` is a thin helper used internally to emit non-fatal advisories — for example ``sspre`` warning that a parameter looks far outside its normal range. Note it is a *function* that calls :func:`warnings.warn`, not a warning category, so it raises a plain :class:`UserWarning`. To escalate SOURSOP's advisories into errors, filter on ``UserWarning`` (optionally narrowed by module)::

    import warnings

    warnings.filterwarnings('error', category=UserWarning, module='soursop')

.. autoclass:: soursop.ssexceptions.SSException
.. autoclass:: soursop.ssexceptions.notYetImplementedException
.. autofunction:: soursop.ssexceptions.SSWarning
