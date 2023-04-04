Development
=============

SOURSOP was designed to be highly extendable. In particular, third parties (that's YOU!) are invited to provide plugins in the soursop/plugins directory. These can be simple stand-alone functions or entire classes that implement statless or statefull functionality. The bar to developing and sharing a plugin was meant to be kept as low as possible. If you have an analysis routine you think we be useful to include please consider making a pull-request, or contact Alex or Jared about integrating new code into SOURSOP.

General workflow for adding plugins
----------------------------------------

1. Fork the SOURSOP repo and add your code into the plugins directory
2. Once your code is complete, please write a set of stand-alone tests using PyTest - if this is challenging PLEASE don't hesistate to reach out to Alex and Jared about how best to do this.
3. Finally, once your code is working and tests are run, you can make a pull request to merge your fork back into the main SOURSOP branch.

Example
---------------
By way of example, we have a simple demo of a possible plugin in the `soursop/plugins
<https://github.com/holehouse-lab/soursop/tree/master/soursop/plugins/>`_ directory.

To use this code, a short example is provided below::

	from soursop.sstrajectory import SSTrajectory
	T = SSTrajectory('ntl9.xtc', 'ntl9.pdb')
	
	NTL9_CP = T.proteinTrajectoryList[0]


	# import the sparrow_plugin module
	from soursop.plugins import sparrow_plugin 
	
	print(sparrow_plugin.get_protein_net_charge(NTL9_CP))
	
	spobj = sparrow_plugin.SparrowProtein(NTL9_CP)
	print(spobj.NCPR)
	
	




