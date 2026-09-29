
sstrajectory
=========================================================

Overview
----------------------------

``SSTrajectory`` is the entry point for all SOURSOP analyses. It reads a simulation trajectory from disk (or accepts a pre-loaded ``mdtraj.Trajectory`` object), identifies every protein chain present, and exposes each chain as an :class:`~soursop.ssprotein.SSProtein` object via the ``proteinTrajectoryList`` attribute.

**Relationship to SSProtein.** ``SSTrajectory`` operates at the *system* level; ``SSProtein`` operates at the *chain* level. For analyses of a single isolated chain — radius of gyration, end-to-end distance, secondary structure, internal scaling, etc. — retrieve the relevant ``SSProtein`` from ``proteinTrajectoryList`` and call methods on that object. ``SSTrajectory`` functions are reserved for analyses that span multiple chains, such as inter-chain distance maps and contact maps, or for system-wide observables that combine all chains into a single pseudo-chain (overall :math:`R_g`, :math:`R_h`, asphericity).

**Supported formats.** Any trajectory format supported by `mdtraj <https://mdtraj.org>`_ can be read in (XTC, DCD, TRR, NetCDF, etc.), paired with a compatible topology file (PDB, GRO, etc.).

**Input contract: whole molecules.** SOURSOP computes every distance, size and contact directly from the coordinates it is handed. It never applies the minimum-image convention (the only exception is the opt-in ``periodic=True`` flag on the inter-chain methods below), and it never re-images or unwraps coordinates. This is deliberate - it keeps every observable a plain function of the conformation - but it means the trajectory you load must already contain **whole molecules**. A chain that the simulation engine has wrapped back into the primary cell, so that part of it sits on the far side of the box, gives wrong radii of gyration, distance maps, contact maps and everything derived from them, with no error. Make molecules whole before analysis, for example with ``gmx trjconv -pbc mol -center`` for GROMACS output, or in Python with ``traj = traj.make_molecules_whole()`` / ``traj = traj.image_molecules()`` on an mdtraj trajectory before passing it in via ``TRJ=`` (both mdtraj methods return a new trajectory unless you pass ``inplace=True``, so calling them without assigning the result changes nothing). CAMPARI's default output is already whole.

Because a wrapped chain is easy to miss by eye, ``SSTrajectory`` checks for it on load: in a whole chain every covalent bond is 1-2 Å long, whereas a wrapped one has two bonded atoms separated by roughly a box vector. SOURSOP tests every bond within a residue or between sequence-adjacent residues of the same chain (so a split sidechain or cap is caught as well as a split backbone; topologies without bonds fall back to consecutive CA-CA distances). If any of these exceeds half the shortest box vector in any frame, a warning names the chain, the atoms and residues involved and the number of affected frames (pass ``check_whole_molecules='raise'`` to make this an ``SSException`` instead, or ``False`` to skip the check). A whole chain that is simply larger than the box is not flagged, and trajectories without a unit cell cannot be tested. The same test is available at any time as :meth:`~soursop.sstrajectory.SSTrajectory.check_molecules_whole`, which returns a per-chain report::

    from soursop.sstrajectory import SSTrajectory

    TrajOb = SSTrajectory('traj.xtc', 'start.pdb')
    for entry in TrajOb.check_molecules_whole():
        if entry['split']:
            print(entry['protein'], entry['worst_pair'], entry['n_frames_split'])

**Coarse-grained models.** In addition to all-atom trajectories, SOURSOP supports one-bead-per-residue coarse-grained trajectories (a single ``CA`` bead per residue), which are auto-detected on load. The bulk of the geometric API (dimensions, distance/contact maps, polymer scaling, etc.) works unchanged on coarse-grained ensembles, since every residue still carries a ``CA``.

**Typical workflow**::

    from soursop.sstrajectory import SSTrajectory
    import numpy as np

    # read in a trajectory and topology
    TrajOb = SSTrajectory('traj.xtc', 'start.pdb')

    # SSProtein objects for each chain
    protein = TrajOb.proteinTrajectoryList[0]

    # system-level observable (all chains combined)
    rg_overall = TrajOb.get_overall_radius_of_gyration()
    print(f"Mean overall Rg = {np.mean(rg_overall):.2f} Å")

    # chain-level observable via SSProtein
    print(f"Mean chain Rg   = {np.mean(protein.get_radius_of_gyration()):.2f} Å")

The SSProtein API is documented on the :doc:`ssprotein` page (quick-link in the left navigation bar). The remainder of this page covers trajectory initialisation, class properties, and SSTrajectory-level functions.


.. autoclass:: soursop.sstrajectory.SSTrajectory

        .. automethod:: __init__



SSTrajectory Properties
----------------------------
Properties are class variables that are dynamically calculated. They do not require parentheses at the end. For example::

  from soursop.sstrajectory import SSTrajectory

  # read in a trajectory
  TrajOb = SSTrajectory('traj.xtc', 'start.pdb')

  print(TrajOb.n_frames)

Will print the number of frames associated with this trajectory.


Properties
................

.. autoattribute:: soursop.sstrajectory.SSTrajectory.n_frames
.. autoattribute:: soursop.sstrajectory.SSTrajectory.n_proteins
.. autoattribute:: soursop.sstrajectory.SSTrajectory.length
.. autoattribute:: soursop.sstrajectory.SSTrajectory.unitcell

In addition to these functions, the Python builtin ``len()`` returns the number of frames in a trajectory.


SSTrajectory Functions
----------------------------
Functions enable operations to be performed on the entire system. Note if you wish to obtain information about a specific protein in isolation, you should use the ``SSProtein`` objects, but for information where multiple proteins are considered all operations should be performed at the level of the ``SSTrajectory`` object.

Functions
................

.. automethod:: soursop.sstrajectory.SSTrajectory.check_molecules_whole
.. automethod:: soursop.sstrajectory.SSTrajectory.get_overall_radius_of_gyration
.. automethod:: soursop.sstrajectory.SSTrajectory.get_overall_hydrodynamic_radius
.. automethod:: soursop.sstrajectory.SSTrajectory.get_overall_asphericity
.. automethod:: soursop.sstrajectory.SSTrajectory.get_interchain_distance_map
.. automethod:: soursop.sstrajectory.SSTrajectory.get_interchain_contact_map
.. automethod:: soursop.sstrajectory.SSTrajectory.get_interchain_distance

.. note::

   The ``mode`` keyword is spelled differently by the two inter-chain map
   functions: ``get_interchain_distance_map`` accepts ``'CA'`` / ``'COM'``,
   whereas ``get_interchain_contact_map`` accepts lower-case ``'atom'`` /
   ``'ca'`` / ``'closest'`` / ``'closest-heavy'`` / ``'sidechain'`` /
   ``'sidechain-heavy'``. The contact cutoff keyword is also ``threshold``
   here, while the single-chain
   :meth:`~soursop.ssprotein.SSProtein.get_contact_map` calls it
   ``distance_thresh``.

.. warning::

   **Breaking change in 2.0.6.** ``get_interchain_contact_map`` now indexes
   its rows and columns by the CA-bearing residues of each chain
   (``protein.resid_with_CA``), so caps (ACE, NME, FOR, NH2, ...) are
   excluded and the map is ``(n_CA_1, n_CA_2)``. Row ``k`` is residue
   ``protein.resid_with_CA[k]``. This is the same convention as
   ``get_interchain_distance_map`` and ``SSProtein.get_contact_map``, so the
   three maps can now be compared element by element. Earlier versions
   indexed every residue by resid, caps included (as all-zero rows and
   columns), which made the map one row/column larger at each capped end
   and offset it by one from the other maps on N-terminally capped chains.
   Code that indexes the contact map by resid on a capped chain needs
   updating; uncapped chains are unaffected.


Loading many trajectories in parallel
--------------------------------------

Reading a large set of replicate trajectories is usually I/O- and parse-bound,
so SOURSOP provides a helper that loads them across multiple processes and
returns a list of ``SSTrajectory`` objects::

    from soursop.sstrajectory import parallel_load_trjs

    trj_files = ['rep0/traj.xtc', 'rep1/traj.xtc', 'rep2/traj.xtc']
    top_files = ['rep0/start.pdb', 'rep1/start.pdb', 'rep2/start.pdb']

    trajectories = parallel_load_trjs(trj_files, top_files, n_procs=4)

    for trajectory in trajectories:
        protein = trajectory.proteinTrajectoryList[0]
        print(protein.get_radius_of_gyration().mean())

Any additional keyword arguments are passed straight through to the
``SSTrajectory`` constructor.

.. autofunction:: soursop.sstrajectory.parallel_load_trjs
			
