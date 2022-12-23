##     _____  ____  _    _ _____   _____  ____  _____
##   / ____|/ __ \| |  | |  __ \ / ____|/ __ \|  __ \
##  | (___ | |  | | |  | | |__) | (___ | |  | | |__) |
##   \___ \| |  | | |  | |  _  / \___ \| |  | |  ___/
##   ____) | |__| | |__| | | \ \ ____) | |__| | |
##  |_____/ \____/ \____/|_|  \_\_____/ \____/|_|

## Alex Holehouse (Pappu Lab and Holehouse Lab) and Jared Lalmansing (Pappu lab)
## Simulation analysis package
## Copyright 2014 - 2022
##

import mdtraj as md
import numpy as np
import scipy
from .ssexceptions import SSWarning, SSException
from .ssprotein import SSProtein

### SSPRE contains all the functionality associated with calculating
### PRE profiles.
###
###
###
###


original_K = 1.2300e-32       # K constant in cm6*s-2
K_IN_NM6   = original_K*1e42  # K constant in nm6 s-2
#W_H        = 267530000        # Proton Larmor frequency
#W_H_SQUARED = W_H * W_H

class SSPRE:
    """
    Class for generating synthetic paramagnetic resonance enhancement
    profiles.

    """

    # ........................................................................
    #
    def __init__(self, SSProteinObject, tau_c, t_delay, R_2D, W_H):
        """
        Initialization function for creating a SSPRE object. The resulting 
        object can then be used to calculate PRE profiles from the underlying 
        ensemble. This calculation is extremely fast.

        Parameters
        ---------------

        SSProteinObject : SSProtein-derived object
            SSProtein object extracted from a SSTrajectory objected. The 
            SSProtein object is the main object that most protein-based 
            analysis is performed over in SOURSOP.

        tau_c : float
            tau_c is the effective correlation time, measured in nanoseconds, 
            which is typically between 1 and 30.

        t_delay : float
            Total duration of the INEPT delays from the PRE experiment, as 
            measured in ms. This will depend on the pulse sequence used, 
            but is typically around 1-30 ms for HSQC.

        R_2D : float
            Is the transverse relaxation rate of the backbone amide protons in 
            the diamagnetic form of the protein, measured in Herz (i.e. 'per 
            second'). A value of around 10 might be expected.

        W_H : float
            Is the proton Larmor frequency, which is typically the "MHz" value
            associated with the magnet, given in Hz. For examle, a 600 MHz 
            magnet would use the value 600000000. Note that the proton 
            Larmor frequency at 1 Tesla = 267530000 per second per Tesla.

        """

        # NOTE: NO CHANGES SHOULD BE MADE TO THE SSPO by the object - this should be treated
        # as a read only object (no such explicit control in Python)
        self.SSPO = SSProteinObject

        if type(self.SSPO) is not SSProtein:
            raise SSException(f'SSPRE requires an SSProtein object to be passed, but instead the first argument is of type {type(self.SSPO)}') 

        # set the INEPT delay value and the backbone amide transverse relaxation rate which
        # is used explicitly later
        self.t_delay = float(t_delay)   # in ms - INDEPT delay
        self.R_2D = float(R_2D)         # in Hz - backbone amide transverse relaxation rate
        self.tau_c = tau_c              # in ns - effecive correlation time
        self.W_H   = W_H                # in Hz - Proton Larmor frequency in the magnet

        # ------------------------------------------
        # sanity checks to warn if any of the input values seem dratsically wrong. NOTE that these won't block the analysis but will
        # throw up errors
        if self.R_2D < 0.01 or self.R_2D > 100:
            SSWarning("WARNING: The value of R_2D (bacbone amide transverese relaxation rate) is far from the normal expected value of ~10 (R_2D = %4.4e) - recal this value is in units of Herz" %(self.R_2D))

        if self.t_delay < 0.01 or self.t_delay > 100:
            SSWarning("WARNING: The value of t_delay (INEPT delay) is far from the normal expected value of ~15 (t_delay = %4.4e) - recal this is in units of ms" %(self.t_delay))

        if tau_c < 0.01 or tau_c > 100:
            SSWarning("WARNING: The value of tau_c (effective correlation time) is far from the normal expected value of ~5 (t_delay = %4.4e) - recal this is in units of ns" %(self.tau_c))

        # if Larmor frequency less than 100 MhZ or above 2 GHz assume something is wrong
        if W_H < 50000000 or W_H > 2000000000:
            SSWarning(f"WARNING: The value of W_h {self.W_H} (proton Larmor frequency) is far from the normal expected value of ~600 000 000 - recal this value should be provided in Herz")

        # # convert tau_c to seconds and calculate tau_c squared
        tau_c = float(tau_c)/1000000000     # tau c in seconds
        tau_c_squared = tau_c * tau_c       #

        # compute the prefactor term which will be used when computing the PRE dependent relaxation profile
        # by the generate_PRE_profile function
        W_H_SQUARED = W_H*W_H
        PREFACTOR = (3 * tau_c)/(1 + W_H_SQUARED * tau_c_squared)
        PREFACTOR = (4*tau_c + PREFACTOR)
        self.PREFACTOR = PREFACTOR * K_IN_NM6


    # ........................................................................
    #
    def __repr__(self):
        """

        """
        return "["+hex(id(self)) + "]: SSPRE OBJ - (R_2D = %3.2f Hz, t_delay = %3.2f ms, tau_c = %3.2f ns, H1 Larmor = %3.3e Hz)" % (self.R_2D, self.t_delay, self.tau_c, self.W_H)


    # ........................................................................
    #
    def generate_PRE_profile(self, label_position, spin_label_atom='CB', target_relaxation_atom='N'):
        """
        Construct a PRE intensity profile and gamma profile based on a nitroxide
        spin label being placed at the label position position on the 
        $spin_label_atom atom. By default this is the CB and the PRE distance 
        to be used to asses relaxation comes from the CB-backbone N distances 
        (as used in work by Meng & Lyle [1], Das [2], and Peran & Holehouse [3]).

        The PRE profile describes the intensity ratio (I_paramagnetic / 
        I_diamagnetic), and typically varies between 0 and 1. When the intensity 
        ratio ~0 the spin label dominantes relaxation and suggests the residue 
        in question is near the spin label. When the intensity ratio is ~1 the 
        relaxation is primarily through non-spin label mediated mechanisms
        suggesting the label and the residue are far apart.

        The gamma profile describes the spin-label induced amide proton 
        relaxation rate in units of per-second, providing another observable 
        that can be directly compared with experiment.
        
        This function is extremely fast (sub 10 seconds on a ~6000 frame 
        ensemble).

        Parameters
        --------------

        label_position : int
            Position in the sequence at which the spin-label is located. 
            Should ideally contain a CB atom (i.e. not be be glycine), 
            else the label atom must be set to 'CA' see below.
            
        spin_label_atom : str (default = 'CB')
            Name of the atom upon which the spin label is located. Should 
            really be CB but may be changed if a residue lacks a CB atom 
            (e.g. a glycine is in the place of the Cys nitroxide spin 
            labeled residue).

        target_relaxation_atom : str (default='N')
            Name of the atom where relaxation is being performed. This 
            should be 'N' (backbone amide) as that's how this approach is 
            parameterized - highly recommended that this isn't changed. 
            If it is changed the method will look for an atom of this name 
            in every residue. Again, it is STRONGLY recommended this
            isn't changed.

        Returns
        -------
        tuple
            Returns a 2 place tuple - tuple position 0 is the PRE intensity
            profile and tuple position 1 is the PRE H1 relaxatation profile.
            

        References
        -------------
        [1] Meng, W., Lyle, N., Luan, B., Raleigh, D.P., and Pappu, R.V. 
        (2013). Experiments and simulations show how long-range
        contacts can form in expanded unfolded proteins with negligible 
        secondary structure.
        Proc. Natl. Acad. Sci. U. S. A. 110, 2123-2128.

        [2] Das, R.K., Huang, Y., Phillips, A.H., Kriwacki, R.W., and Pappu, 
        R.V. (2016). Cryptic sequence features within the disordered protein 
        p27Kip1 regulate cell cycle signaling. Proc. Natl. Acad. Sci. U. S. A. 
        113, 5616- 5621.
        
        [3] Peran, I., Holehouse, A. S., Carrico, I. S., Pappu, R. V., Bilsel, 
        O., & Raleigh, D. P. (2019). Unfolded states under folding conditions 
        accommodate sequence-specific conformational preferences with random 
        coil-like dimensions. Proceedings of the National Academy of Sciences 
        of the United States of America, 116(25), 12301–12310.

        """

        # get index value of all residues
        #residue_list = self.SSPO.get_residue_index_list()

        tmp = list(self.SSPO._SSProtein__CA_residue_atom.keys())
        residue_list = sorted(tmp)

        # first calculate mean rij distance for pair-residue distances
        gamma = []

        # finally for each residue calculate the r^6 distances associated with each frame and for EACH FRAME calculate
        # the PREFACTOR / r^6 value and then take the mean. THIS gives a different answer to if you take the mean distance
        # and calculate the PREFACTOR/<R^6> value because there is a non-linear mapping between relaxation and distance so
        # it's important the former method is used (i.e. only average at the end). This calculates the gamma coefficient for
        # each residue, which measures relaxation
        for idx in residue_list:
            r_6_nm = np.power(0.1*self.SSPO.get_inter_residue_atomic_distance(label_position, idx, A1=spin_label_atom, A2=target_relaxation_atom),6)
            gamma.append(np.mean(self.PREFACTOR/r_6_nm))

        # convert the t_delay from ms to seconds
        t_delay_in_seconds = self.t_delay/1000

        # for each gamma compute the intensity ration
        profile = []
        for g in gamma:
            profile.append((self.R_2D * np.exp(-g*t_delay_in_seconds)) / (self.R_2D + g))

        return (profile, gamma)
