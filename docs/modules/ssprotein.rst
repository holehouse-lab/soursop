
ssprotein
=========================================================

Overview
----------------------------

``SSProtein`` is the core analysis class in SOURSOP and the primary interface for characterizing the conformational behaviour of a single protein chain across a simulation trajectory.

``SSProtein`` objects are not created directly. Instead, they are extracted from an :class:`~soursop.sstrajectory.SSTrajectory` object after reading a trajectory from disk. Each protein chain in the system is represented by its own ``SSProtein`` instance, accessible via ``SSTrajectory.proteinTrajectoryList``.

Analyses fall into several broad categories:

* **Residue utilities** — sequence retrieval, atom/CA index lookup, residue center-of-mass positions and masses.
* **Inter-residue distances** — pairwise CA or COM distance matrices, distance maps, and polymer-scaled distance maps.
* **Global size and shape** — radius of gyration, hydrodynamic radius, end-to-end distance, asphericity, acylindricity, prolateness, gyration tensor, and the :math:`t`-parameter.
* **Secondary structure** — per-frame DSSP assignments and BBSEG backbone-torsion-based classification.
* **Polymer scaling** — internal scaling profiles (:math:`\langle r^2 \rangle` vs sequence separation), the scaling exponent :math:`\nu`, and local heterogeneity in scaling behaviour.
* **Contact and RMSD analysis** — contact maps (with configurable threshold and mode), RMSD to a reference structure, and fraction of native contacts :math:`Q`.
* **Solvent accessibility** — per-residue and region-level SASA via ``get_all_SASA``, ``get_regional_SASA``, and ``get_site_accessibility``. On one-bead-per-residue coarse-grained chains these use the force field's own per-residue bead sizes (see :ref:`cg-sasa` below) rather than atomic van der Waals radii.
* **Local dynamics** — local collapse profiles, sidechain alignment angles, dihedral mutual information, local-to-global correlation, and angle decay.
* **Clustering and overlap** — conformational clustering via ``get_clusters``, overlap concentration :math:`c^*`.

Most functions return NumPy arrays; per-frame results have shape ``(n_frames,)`` or ``(n_frames, ...)`` so standard NumPy operations (``np.mean``, ``np.std``, etc.) apply directly.

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


.. _cg-sasa:

SASA for coarse-grained ensembles
-------------------------------------

One-bead-per-residue coarse-grained models (Mpipi, HPS, KH and friends) write every residue out as a single ``CA`` atom, so mdtraj — which assigns van der Waals radii from the atomic *element* — would give a glycine bead and a tryptophan bead the same 1.7 Å radius. That is not a useful approximation when the actual bead diameters in these models span 4.5 to 8.5 Å.

SOURSOP therefore computes SASA on these chains using the force field's own bead sizes. The bead radius is taken as :math:`\sigma/2`, where :math:`\sigma` is the model's diagonal (i,i) pair sigma — i.e. the bead diameter. Because SOURSOP cannot infer which model produced an ensemble, the force field has to be named explicitly; a coarse-grained SASA call without one raises an ``SSException`` rather than silently returning atomic-radius numbers::

  from soursop.sstrajectory import SSTrajectory
  from soursop.ssdata import CG_FORCEFIELDS

  TrajOb = SSTrajectory('cg_traj.xtc', 'cg_start.pdb')
  ProtObj = TrajOb.proteinTrajectoryList[0]

  ProtObj.is_coarse_grained          # True for one-bead-per-residue chains

  # name the model on the call ...
  sasa = ProtObj.get_all_SASA(mode='residue', forcefield='mpipi-gg')

  # ... or set it once for this protein
  ProtObj.cg_forcefield = 'mpipi-gg'
  sasa = ProtObj.get_all_SASA(mode='residue')

  print(CG_FORCEFIELDS)
  # ('cation-pi-1', 'cation-pi-2', 'fb-hps', 'hps-kr', 'hps-urry',
  #  'kh', 'mpipi', 'mpipi-gg', 'mpipi-recharged')

The same ``forcefield`` keyword is accepted by ``get_regional_SASA`` and ``get_site_accessibility``. On a one-bead chain only ``mode='residue'`` and ``mode='atom'`` are meaningful (and are equivalent) — a single bead represents the whole residue, so there is no backbone/sidechain decomposition to make and those modes raise.

The underlying bead sizes are readable directly if you want them for something else::

  from soursop.ssdata import get_cg_bead_sigmas

  sigmas = get_cg_bead_sigmas('mpipi')     # residue name -> sigma in Angstroms
  sigmas['GLY']                            # 4.69511

Note this applies only to one-bead-per-residue models. All-atom chains — and two-bead (CA/CB) models, which carry a CB for every non-glycine residue — use mdtraj's atomic radii as before, and passing a ``forcefield`` to them raises.


SSProtein Properties
----------------------------
SSProtein objects have a set of object variables associated with them.

.. autoclass:: soursop.ssprotein.SSProtein

        .. autoattribute:: resid_with_CA
        .. autoattribute:: ncap
        .. autoattribute:: ccap
        .. autoattribute:: is_coarse_grained
        .. autoattribute:: is_swan
        .. autoattribute:: cg_forcefield
        .. autoattribute:: n_frames
        .. autoattribute:: n_residues
        .. autoattribute:: residue_index_list
        .. autoattribute:: unitcell
        .. autoattribute:: length


SSProtein Functions
----------------------------

.. autoclass:: soursop.ssprotein.SSProtein
        :no-index:

        .. automethod:: reset_cache
        .. automethod:: print_residues
        .. automethod:: get_amino_acid_sequence
        .. automethod:: get_residue_atom_indices
        .. automethod:: get_all_atomic_indices
        .. automethod:: get_CA_index
        .. automethod:: get_multiple_CA_index
        .. automethod:: get_residue_COM
        .. automethod:: get_residue_mass
        .. automethod:: get_center_of_mass
        .. automethod:: get_molecular_volume
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
        .. automethod:: get_acylindricity
        .. automethod:: get_prolateness
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


