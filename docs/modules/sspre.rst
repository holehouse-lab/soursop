
sspre
=========================================================

sspre contains functions related to paramagnetic relaxation enhancement (PRE) analysis of IDPs.

    
SSPRE is a class that takes a SSProtein objecs and can perform PRE-related calculations.
       
This is, in many ways, a purely functional class, but given it's very specific goal it is separated out into its own class in the interest of more robust modularity. 

There are a number of internal functions, but the only public facing function is the  generate_PRE_profile function below. This gives both the intensity ratio profile and  the transverse relaxation rate (Gamma_2) profiles. For more detail on the calculation of these profiles see the help associated with generate_PRE_profile function

For more information on caculation of PREs using method see the Supplementary information in the following two papers.

Meng, W., Lyle, N., Luan, B., Raleigh, D.P., and Pappu, R.V. (2013). Experiments and simulations 
show how long-range contacts can form in expanded unfolded proteins with negligible secondary structure. 
Proc. Natl. Acad. Sci. U. S. A. 110, 2123-2128.

Das, R.K., Huang, Y., Phillips, A.H., Kriwacki, R.W., and Pappu, R.V. (2016). Cryptic sequence 
features within the disordered protein p27Kip1 regulate cell cycle signaling. 
Proc. Natl. Acad. Sci. U. S. A. 113, 5616- 5621.

Peran, I., Holehouse, A. S., Carrico, I. S., Pappu, R. V., Bilsel, O., & Raleigh, D. P. (2019). Unfolded states under 
folding conditions accommodate sequence-specific conformational preferences with random coil-like dimensions. Proceedings 
of the National Academy of Sciences of the United States of America, 116(25), 12301–12310.

.. autoclass:: soursop.sspre.SSPRE

        .. automethod:: __init__
        .. automethod:: generate_PRE_profile

