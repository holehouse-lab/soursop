Overview
==========

SOURSOP is a Python-based simulation analysis package developed for the
analysis of conformational ensembles of disordered and unfolded
proteins. Its goal is to make it as easy as possible to read in an
ensemble of an intrinsically disordered region (IDR) and quickly extract
polymer-physics-aware observables. In addition to a large library of
pre-built analysis routines, SOURSOP provides easy and rapid access to
all inter-residue and inter-atomic distances, contact maps, dimensions,
secondary-structure content, and more.

SOURSOP is built on top of `MDTraj <https://mdtraj.org/>`_, which
handles trajectory I/O and the low-level atomic representation. SOURSOP
focuses on the analysis layer, with routines specifically chosen to be
useful for characterizing disordered and unfolded states through the
lens of polymer physics.


A first example
----------------------

.. code-block:: python

  # import soursop
  from soursop.sstrajectory import SSTrajectory

  # read in the simulation trajectory (trajectory file + topology file)
  TO = SSTrajectory('traj.xtc', 'start.pdb')

  # once the trajectory has been read in, individual protein chains can
  # be extracted from the proteinTrajectoryList
  protein = TO.proteinTrajectoryList[0]

  # per-residue center-of-mass distance between residues 10 and 20
  d_10_20 = protein.get_inter_residue_COM_distance(10, 20)

  # ensemble-average asphericity
  asph = protein.get_asphericity()

  # ensemble-average radius of gyration
  rg = protein.get_radius_of_gyration()

  # ensemble-average inter-residue distance map
  dm = protein.get_distance_map()


Core concepts
----------------------

Two objects underpin almost all analysis in SOURSOP:

``SSTrajectory``
~~~~~~~~~~~~~~~~~~~~~~~~~~

``SSTrajectory`` is the top-level, *system-level* object. You construct
it from a trajectory file and a topology (PDB) file. It wraps the
underlying MDTraj trajectory, identifies the distinct protein chains in
the system, and exposes system-wide and inter-chain analyses (for
example inter-chain distance and contact maps for multi-chain
simulations of, e.g., phase-separating systems or protein complexes).

``SSProtein``
~~~~~~~~~~~~~~~~~~~~~~~~~~

Each individual protein chain identified during loading is represented
as an ``SSProtein`` object, accessed through
``SSTrajectory.proteinTrajectoryList``. ``SSProtein`` is where the bulk
of single-chain analysis lives: global dimensions (Rg, Rh,
end-to-end distance, asphericity), polymer scaling, distance and contact
maps, secondary structure, solvent accessibility, dihedral angles, and
arbitrary inter-residue / inter-atomic distances. Expensive lookups are
memoised, so repeated queries on the same object are fast.

The typical pattern is therefore:

.. code-block:: python

   TO = SSTrajectory('traj.xtc', 'start.pdb')   # system-level object
   P  = TO.proteinTrajectoryList[0]             # one protein chain
   # ... run analyses on P ...

For multi-chain systems, iterate over ``proteinTrajectoryList`` (one
``SSProtein`` per chain) and use the ``SSTrajectory`` inter-chain
methods for cross-chain observables.


The SOURSOP modules
----------------------

SOURSOP is organised into a small number of focused modules:

* ``sstrajectory`` - the ``SSTrajectory`` class; trajectory loading,
  chain detection, system-level and inter-chain analysis, and helpers
  for parallel loading of many trajectories.

* ``ssprotein`` - the ``SSProtein`` class; the main single-chain
  analysis engine (dimensions, scaling, maps, secondary structure,
  SASA, angles, and inter-residue/atomic distances).

* ``ssnmr`` - NMR observables for comparison against experiment:
  sequence-based random-coil backbone chemical shifts (CA, CB, CO, N,
  HN, HA; with temperature, pH and perdeuteration corrections and
  phospho-residue support), backbone ³J(HN, Hα) scalar couplings from
  the φ dihedral via the Karplus relation, and per-frame NOE distances.

* ``sspre`` - the ``SSPRE`` class; fast calculation of synthetic
  paramagnetic relaxation enhancement (PRE) intensity ratios and gamma
  profiles for a spin label placed at an arbitrary sequence position.

* ``sssampling`` - the ``SSSampling`` class and PENGUIN support, for
  assessing the sampling quality / convergence of disordered-protein
  ensembles.

* ``ssbme`` - the ``BME``, ``iBME`` and ``BMECustom`` classes; Bayesian
  Maximum Entropy reweighting of an ensemble against experimental
  observables (``iBME`` additionally fits an unknown data scale/offset;
  ``BMECustom`` takes a profile/vector observable with an arbitrary
  user-supplied goodness-of-fit), producing per-frame ``weights`` that
  plug directly into every reweighting-capable analysis routine (see
  :doc:`weights`).

* ``sscoper`` - the ``COPER`` and ``iCOPER`` classes; Convex Optimization
  for Ensemble Reweighting, an alternative maximum-entropy reweighter that
  imposes a hard chi-squared constraint (rather than BME's tunable
  penalty) and reports whether the data can be fit at all. Shares the
  ``ExperimentalObservable`` interface with ``ssbme``.

* ``sshdx`` - per-residue HDX protection factors via the
  Best-Vendruscolo formula. The ``(n_frames, n_residues)`` output is the
  natural input for BME / COPER reweighting against experimental HDX
  protection-factor data.

* ``ssmutualinformation`` - dihedral mutual-information helpers
  (``calc_MI``) underpinning ``SSProtein.get_dihedral_mutual_information``,
  for quantifying correlated backbone/side-chain motions.

* ``ssutils`` - shared validation and reduction helpers, including the
  single ``weights`` validator/reducers used package-wide and the
  reweighting primitives (``ExperimentalObservable``, ...) shared by
  ``ssbme`` and ``sscoper``.

* ``sstools`` - miscellaneous numerical helper functions shared across
  the package (chunking, residue-name normalisation, the polymer
  power-law model, minimum-image distances, trajectory-file discovery).

In addition, SOURSOP is designed to be extended via user-contributed
plugins in ``soursop/plugins`` - see the :doc:`development` page.


Where to go next
----------------------

* :doc:`installation` - install SOURSOP with pip, uv, or conda.

* :doc:`examples` - worked, end-to-end IDP analysis examples.

* :doc:`development` - extending SOURSOP and contributing plugins.

* The per-module API references (``sstrajectory``, ``ssprotein``,
  ``ssnmr``, ``sspre``, ``sssampling``, ``ssbme``, ``sscoper``,
  ``sshdx``) for the full list of available analysis routines.

* :doc:`../modules/bme` and :doc:`../modules/coper` - reweighting an
  ensemble against experimental data (BME / iBME, or COPER / iCOPER) to
  generate frame ``weights``.

* :doc:`../modules/hdx` - HDX protection factors via the
  Best-Vendruscolo formula (per-residue, per-frame).
