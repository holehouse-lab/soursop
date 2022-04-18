
sstrajectory 
=========================================================

SSTrajectory Class
----------------------------
Reading in trajectories starts with the SSTrajectory class. Upon reading, SOURSOP extracts out individual protein chains which can then be subsequently analyzed. 

By way of a quickstart example::

    from soursop.sstrajectory import SSTrajectory
    import numpy as np

    # read in a trajectory
    TrajOb = SSTrajectory('traj.xtc', 'start.pdb')

    # get the first protein chain (this is an SSProtein object)
    ProtObj = TrajOb.proteinTrajectoryList[0]

    # print the mean radius of gyration
    print(np.mean(ProtObj.get_radius_of_gyration()))
  

The set of possible SSProtein functions can be found in the ssprotein documentation page (quick-link to the left in the navigation bar). The remainder of these documents walk through trajectory initialization, SSTrajectory class variables, SSTrajectory properties, and SSTrajectory functions.


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

.. autoclass:: soursop.sstrajectory.SSTrajectory

        .. automethod:: n_frames
        .. automethod:: n_proteins
        .. automethod:: length

SSTrajectory Functions
----------------------------

.. autoclass:: soursop.sstrajectory.SSTrajectory

        .. automethod:: get_overall_radius_of_gyration
        .. automethod:: get_overall_asphericity
        .. automethod:: get_overall_hydrodynamic_radius
        .. automethod:: get_interchain_distance_map
        .. automethod:: get_interchain_distance
			
