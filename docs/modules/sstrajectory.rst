
sstrajectory
=========================================================

Overview
----------------------------

``SSTrajectory`` is the entry point for all SOURSOP analyses. It reads a simulation trajectory from disk (or accepts a pre-loaded ``mdtraj.Trajectory`` object), identifies every protein chain present, and exposes each chain as an :class:`~soursop.ssprotein.SSProtein` object via the ``proteinTrajectoryList`` attribute.

**Relationship to SSProtein.** ``SSTrajectory`` operates at the *system* level; ``SSProtein`` operates at the *chain* level. For analyses of a single isolated chain — radius of gyration, end-to-end distance, secondary structure, internal scaling, etc. — retrieve the relevant ``SSProtein`` from ``proteinTrajectoryList`` and call methods on that object. ``SSTrajectory`` functions are reserved for analyses that span multiple chains, such as inter-chain distance maps and contact maps, or for system-wide observables that combine all chains into a single pseudo-chain (overall :math:`R_g`, :math:`R_h`, asphericity).

**Supported formats.** Any trajectory format supported by `mdtraj <https://mdtraj.org>`_ can be read in (XTC, DCD, TRR, NetCDF, etc.), paired with a compatible topology file (PDB, GRO, etc.).

**Coarse-grained and SWAN models.** SOURSOP handles all-atom trajectories, one-bead-per-residue coarse-grained trajectories (a single ``CA`` bead per residue), and **SWAN** two-bead coarse-grained trajectories. SWAN (a currently unpublished backbone-sidechain model) represents every residue with a ``CA`` backbone bead plus a ``CB`` sidechain bead, except glycine which has only a ``CA``. SWAN topologies are **auto-detected on load** and reported via the ``SSTrajectory.swan_trajectory`` attribute (and propagated to each chain's :attr:`~soursop.ssprotein.SSProtein.is_swan`). Auto-detection can be overridden with the ``swan_trajectory`` constructor keyword (``None`` = auto-detect, the default; ``True`` / ``False`` = force). For SWAN chains, sidechain-vector and secondary-structure analyses switch to ``CA``/``CB`` definitions automatically (see the :doc:`ssprotein` page); functions that fundamentally require the full backbone or sidechain heavy atoms (backbone/sidechain dihedrals and the BBSEG classification) raise an ``SSException``.

::

    TrajOb = SSTrajectory('swan_traj.xtc', 'swan_top.pdb')
    TrajOb.swan_trajectory                       # True (auto-detected)
    TrajOb.proteinTrajectoryList[0].is_swan      # True

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

In addition to these three functions, the Python builtin ``len()`` returns the number of frames in a trajectory.


SSTrajectory Functions
----------------------------
Functions enable operations to be performed on the entire system. Note if you wish to obtain information about a specific protein in isolation, you should use the ``SSProtein`` objects, but for information where multiple proteins are considered all operations should be performed at the level of the ``SSTrajectory`` object.

Functions
................

.. automethod:: soursop.sstrajectory.SSTrajectory.get_overall_radius_of_gyration
.. automethod:: soursop.sstrajectory.SSTrajectory.get_overall_hydrodynamic_radius
.. automethod:: soursop.sstrajectory.SSTrajectory.get_overall_asphericity
.. automethod:: soursop.sstrajectory.SSTrajectory.get_interchain_distance_map
.. automethod:: soursop.sstrajectory.SSTrajectory.get_interchain_contact_map
.. automethod:: soursop.sstrajectory.SSTrajectory.get_interchain_distance
			
