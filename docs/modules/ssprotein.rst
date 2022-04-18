
ssprotein
=========================================================

The ssprotein module holds the SSProtein class. This is the main class used for working with simulations of individual proteins.

Once a trajectory has been read in and an SSProtein object extracted from the SSTrajectory object, a wide range of analyses are available. 

By way of an example::

  from soursop.sstrajectory import SSTrajectory
  import numpy as np

  # read in a trajectory
  TrajOb = SSTrajectory('traj.xtc', 'start.pdb')

  # get the first protein chain (this is an SSProtein object)
  ProtObj = TrajOb.proteinTrajectoryList[0]

  # print the mean radius of gyration
  print(np.mean(ProtObj.get_end_to_end_distance))


SSProtein Properties
----------------------------
SSProtein objects have a set of object variables associated with them.

.. autoclass:: soursop.ssprotein.SSProtein

        .. automethod:: resid_with_CA
        .. automethod:: ncap
        .. automethod:: ccap
        .. automethod:: n_frames
        .. automethod:: n_residues
        .. automethod:: residue_index_list
        .. automethod:: length


SSProtein Functions
----------------------------

.. autoclass:: soursop.ssprotein.SSProtein

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
        .. automethod:: get_secondary_structure_DSSP
        .. automethod:: get_secondary_structure_BBSEG
        .. automethod:: get_internal_scaling
        .. automethod:: get_internal_scaling_RMS
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


