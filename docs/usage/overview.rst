Overview
==========

SOURSOP is a Python-based simulation analysis package developed for the analysis of conformational ensembles of disordered proteins. The goal is to make it as easy as possible to quickly read in and analyze an ensemble of an IDR. In addition to the many pre-built analysis routines, SOURSOP provides easy and rapid access to all inter-residue or inter-atomic distances.

As an example:

.. code-block:: python

  # import soursop		
  from soursop.sstrajectory import SSTrajectory

  # read in the simulation trajectory 
  TO = SSTrajectory('traj.xtc', 'start.pdb')

  # once the trajectory has been read in, proteins can be extracted
  # from the proteinTrajectoryList
  
  protein = TO.proteinTrajectoryList[0]

  # calculate per-residue distance between residues 10 and 20
  d_10_20 = protein.get_inter_residue_COM_distance(10, 20)

  # calculate the ensemble asphericity
  asph = protein.get_asphericity()

  # calculate the ensemble asphericity
  rg = protein.get_radius_of_gyration()

  # calculate the ensemble distance map
  dm = protein.get_distance_map()

  




