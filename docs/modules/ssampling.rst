
ssampling 
=========================================================

ssampling contains functions related to assessing sampling quality. This is where PENGUIN (Pipeline for Evaluating coNformational heteroGeneity in Unstructured proteINs) is implemented


SamplingQuality Functions
----------------------------

.. autoclass:: soursop.sssampling.SamplingQuality

        .. automethod:: compute_frac_helicity
        .. automethod:: compute_dihedral_hellingers
        .. automethod:: compute_dihedral_rel_entropy
        .. automethod:: compute_series_of_histograms_along_axis
        .. automethod:: compute_pdf
        .. automethod:: get_all_to_all_2d_trj_comparison
        .. automethod:: get_all_to_all_trj_comparisons
        .. automethod:: get_degree_bins
        .. automethod:: quality_plot
        .. automethod:: trj_pdfs
        .. automethod:: ref_pdfs
        .. automethod:: hellingers_distances
        .. automethod:: fractional_helicity
		
		
PrecomputedDihedralInterface Functions
----------------------------------------------

.. autoclass:: soursop.sssampling.PrecomputedDihedralInterface	
	
        .. automethod:: sample_angles
        .. automethod:: gather_phi_reference_dihedrals
        .. automethod:: gather_psi_reference_dihedrals


Stand-alone ssampling functions
----------------------------------------------
		
.. automodule:: soursop.sssampling
.. autofunction:: rel_entropy
.. autofunction:: hellinger_distance
.. autofunction:: compute_joint_hellinger_distance
