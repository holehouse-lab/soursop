
ssprotein
=========================================================

Overview
----------------------------

``SSProtein`` is the core analysis class in SOURSOP and the primary interface for characterizing the conformational behaviour of a single protein chain across a simulation trajectory.

``SSProtein`` objects are not created directly. Instead, they are extracted from an :class:`~soursop.sstrajectory.SSTrajectory` object after reading a trajectory from disk. Each protein chain in the system is represented by its own ``SSProtein`` instance, accessible via ``SSTrajectory.proteinTrajectoryList``.

Analyses fall into several broad categories:

* **Residue utilities** — sequence retrieval, atom/CA index lookup, residue center-of-mass positions and masses.
* **Inter-residue distances** — pairwise CA or COM distance matrices, distance maps, and polymer-scaled distance maps.
* **Global size and shape** — radius of gyration, hydrodynamic radius, end-to-end distance, asphericity, gyration tensor, and the :math:`t`-parameter.
* **Secondary structure** — per-frame DSSP assignments and BBSEG backbone-torsion-based classification.
* **Polymer scaling** — internal scaling profiles (:math:`\langle r^2 \rangle` vs sequence separation), the scaling exponent :math:`\nu`, and local heterogeneity in scaling behaviour.
* **Contact and RMSD analysis** — contact maps (with configurable threshold and mode), RMSD to a reference structure, and fraction of native contacts :math:`Q`.
* **Solvent accessibility** — per-residue and region-level SASA via ``get_all_SASA``, ``get_regional_SASA``, and ``get_site_accessibility``.
* **Local dynamics** — local collapse profiles, sidechain alignment angles, dihedral mutual information, local-to-global correlation, and angle decay.
* **Clustering and overlap** — conformational clustering via ``get_clusters``, overlap concentration :math:`c^*`.

Most functions return NumPy arrays; per-frame results have shape ``(n_frames,)`` or ``(n_frames, ...)`` so standard NumPy operations (``np.mean``, ``np.std``, etc.) apply directly.

.. note::

   **SWAN coarse-grained chains.** When a chain is a SWAN two-bead (``CA``/``CB``) model — auto-detected by :class:`~soursop.sstrajectory.SSTrajectory` and exposed here via the ``is_swan`` property — a few methods adapt automatically. ``get_sidechain_alignment_angle`` uses the ``CA``→``CB`` vector for every residue (glycine, which has no ``CB``, raises). ``get_secondary_structure_DSSP`` assigns helix/β/coil from the ``CA`` trace against idealized α-helix and extended-strand templates (DSSP itself needs the N/C/O backbone SWAN lacks), with optional ``helix_window`` / ``helix_rmsd_thresh`` / ``beta_window`` / ``beta_rmsd_thresh`` knobs. Methods that require backbone or sidechain dihedrals — ``get_angles``, ``get_secondary_structure_BBSEG`` and ``get_dihedral_mutual_information`` — raise an ``SSException`` for SWAN chains. All other (``CA``-based) analyses work unchanged.

By way of an example::

  from soursop.sstrajectory import SSTrajectory
  import numpy as np

  # read in a trajectory
  TrajOb = SSTrajectory('traj.xtc', 'start.pdb')

  # get the first protein chain (this is an SSProtein object)
  ProtObj = TrajOb.proteinTrajectoryList[0]

  # print the mean end-to-end distance
  mean_e2e = np.mean(ProtObj.get_end_to_end_distance())
  print(mean_e2e)


SSProtein Properties
----------------------------
SSProtein objects have a set of object variables associated with them.

.. autoclass:: soursop.ssprotein.SSProtein

        .. autoattribute:: resid_with_CA
        .. autoattribute:: ncap
        .. autoattribute:: ccap
        .. autoattribute:: is_swan
        .. autoattribute:: n_frames
        .. autoattribute:: n_residues
        .. autoattribute:: residue_index_list
        .. autoattribute:: length


SSProtein Functions
----------------------------

.. autoclass:: soursop.ssprotein.SSProtein
        :no-index:

        .. automethod:: reset_cache
        .. automethod:: print_residues
        .. automethod:: get_amino_acid_sequence
        .. automethod:: get_residue_atom_indices
        .. automethod:: get_CA_index
        .. automethod:: get_multiple_CA_index
        .. automethod:: get_residue_COM
        .. automethod:: get_residue_mass
        .. automethod:: calculate_all_CA_distances
        .. automethod:: get_inter_residue_COM_distance
        .. automethod:: get_inter_residue_COM_vector
        .. automethod:: get_inter_residue_atomic_distance
        .. automethod:: get_distance_map
        .. automethod:: get_polymer_scaled_distance_map
        .. automethod:: get_radius_of_gyration
        .. automethod:: get_hydrodynamic_radius
        .. automethod:: get_end_to_end_distance
        .. automethod:: get_end_to_end_vs_rg_correlation
        .. automethod:: get_asphericity
        .. automethod:: get_t
        .. automethod:: get_gyration_tensor
        .. automethod:: get_angles
        .. automethod:: get_secondary_structure_DSSP
        .. automethod:: get_secondary_structure_BBSEG
        .. automethod:: get_internal_scaling
        .. automethod:: get_internal_scaling_RMS
        .. automethod:: get_scaling_exponent
        .. automethod:: get_local_heterogeneity
        .. automethod:: get_D_vector
        .. automethod:: get_RMSD
        .. automethod:: get_Q
        .. automethod:: get_contact_map
        .. automethod:: get_clusters
        .. automethod:: get_all_SASA
        .. automethod:: get_regional_SASA
        .. automethod:: get_site_accessibility
        .. automethod:: get_local_collapse
        .. automethod:: get_sidechain_alignment_angle
        .. automethod:: get_dihedral_mutual_information
        .. automethod:: get_local_to_global_correlation
        .. automethod:: get_overlap_concentration
        .. automethod:: get_angle_decay


