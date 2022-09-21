
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


Properties
................

.. automethod:: soursop.sstrajectory.SSTrajectory.n_frames
.. automethod:: soursop.sstrajectory.SSTrajectory.n_proteins
.. automethod:: soursop.sstrajectory.SSTrajectory.length

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
			
