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
from scipy.special import rel_entr
from typing import Type, List, Union
import pathlib

from .ssprotein import SSProtein
from .ssexceptions import SSException
from . import ssutils
from .sstrajectory import SSTrajectory, parallel_load_trjs
from glob import glob
from natsort import natsorted

def hellinger_distance(p : np.ndarray, q : np.ndarray) -> np.ndarray:
    """
    Computes the hellinger distance between a set probability distributions p and q.
    The hellinger distances is defined by:
        H(P,Q) = \frac{1}{\sqrt{2}} \times \sqrt{\sum_{i=1}^{k}(\sqrt{p_i}-\sqrt{q_i})^2}
    where k is the length of the probability vectors being compared.

    Parameters
    ----------
    p : np.ndarray
        a probability density function, or series of probabiliy density functions, to compute the hellingers distance on.
        p and q must be the same shape

    q : np.ndarray
        a probability density function, or series of probabiliy density functions, to compute the hellingers distance on
        p and q must be the same shape
    
    Returns
    -------
    np.ndarray
        The hellingers distance for each residue in the sequence
    """
    if p.ndim == 2 and q.ndim == 2:
        hellingers = np.sqrt(np.sum((np.sqrt(p) - np.sqrt(q)) ** 2, axis=1)) / np.sqrt(2)
    else:
        hellingers = np.sqrt(np.sum((np.sqrt(p) - np.sqrt(q)) ** 2)) / np.sqrt(2)
    return hellingers


def rel_entropy(p : np.ndarray, q : np.ndarray) -> np.ndarray:
    """Computes the relative entropy between two probability distributions p and q."""
    if p.ndim == 2 and q.ndim == 2:
        relative_entropy = np.sum(rel_entr(p,q), axis=1)
    else:
        relative_entropy = np.sum(rel_entr(p,q))
    return relative_entropy
    

def glob_traj_paths(filepath : Union[str, pathlib.Path],num_reps : int, traj_name="__traj.xtc",top_name="__START.pdb",mode=None):
    """
    This function just assembles the list of trajectory and topology paths for a set of simulations

    Parameters
    ----------
    filepath : Union[str, pathlib.path]
        _description_
    num_reps : int
        _description_
    mode : _type_, optional
        _description_, by default None
    traj_name : str, optional
        _description_, by default "__traj.xtc"
    top_name : str, optional
        _description_, by default "__START.pdb"
    
    """
    top_paths, traj_paths = [], []
    if mode == "mega":
        directories = glob("*start")
        for directory in directories:
            for rep in range(1,num_reps):
                traj_paths.extend(glob(f"{directory}/{rep}/{traj_name}"))
                top_paths.extend(glob(f"{directory}/{rep}/{top_name}"))
    else:
        traj_paths = glob(f"{filepath}/{traj_name}")
    
    return top_paths, traj_paths

class SamplingQuality:
    """Compare the sampling quality for a simulated trajectory relative to a polymer model limit.
    """
 
    def __init__(self, traj_list, 
                       polymer_model_traj_list, 
                       method, 
                       bwidth = np.pi / 5,
                ):

        super(SamplingQuality, self).__init__()
        self.traj_list = traj_list
        self.polymer_model_traj_list = polymer_model_traj_list
        self.method = method
        self.bwidth = bwidth

        if self.method == "rmsd" or self.method == "p_vects":
            raise NotImplementedError("This functionality has not been implemented yet")

        # ssutils.validate_keyword_option(method, ['dihedral', 'rmsd', 'p_vects'], 'method')
        if self.method == "dihedral":
            if self.bwidth > 2*np.pi or not (self.bwidth > 0):
                raise SSException(f'The bwidth parameter must be between 2*pi and greater than 0. Received {self.bwidth}')
                        
            # n_res (angle) x n_frames              
            self.psi_traj1 = self.traj1.get_angles("psi")[1]
            self.phi_traj1 = self.traj1.get_angles("phi")[1]
            self.psi_traj2 = self.traj2.get_angles("psi")[1]
            self.phi_traj2 = self.traj2.get_angles("phi")[1]

    def compute_dihedral_hellingers(self):
        """Compute the hellingers distance for both the phi and psi angles between two trajectories.

        Returns
        -------
        np.ndarray
            The hellinger distances between the probability density distributions for the phi and psi angles for a pair of trajectories.
        """
        bins = self.get_degree_bins()
        psi_traj1_pdf = self.compute_pdf(self.psi_traj1, bins=bins)
        psi_traj2_pdf = self.compute_pdf(self.psi_traj2, bins=bins)
        
        phi_traj1_pdf = self.compute_pdf(self.phi_traj1, bins=bins)
        phi_traj2_pdf = self.compute_pdf(self.phi_traj2, bins=bins)

        psi_hellingers = hellinger_distance(psi_traj1_pdf, psi_traj2_pdf)
        phi_hellingers = hellinger_distance(phi_traj1_pdf, phi_traj2_pdf)

        return np.array([phi_hellingers, psi_hellingers])

    def compute_pdf(self, arr : np.ndarray, bins : np.ndarray) -> np.ndarray:
        """
        Computes a probability density by constructing a histogram of the data array with a specified set of bins. 

        Parameters
        ----------
        arr : np.ndarray
            A vector of shape (n_res x n_frames). The
        bins : np.ndarray
            The set of bin edges that specify the range for the histogram buckets.

        Returns
        -------
        np.ndarray
            Returns a set of histograms of the probabilities densities for each residue in the amino acid sequence.
            Shape (n_res, len(bins) - 1) 
        """
        pdf = np.apply_along_axis(lambda col: np.histogram(col, bins=bins, density=True)[0], axis=1, arr=arr)
        return pdf


    def get_radian_bins(self) -> np.ndarray:
        bwidth = self.bwidth
        bins = np.arange(-np.pi, np.pi+bwidth, bwidth)
        return bins
        
        
    def get_degree_bins(self) -> np.ndarray:
        bwidth = np.rad2deg(self.bwidth)
        bins = np.arange(-180, 180+bwidth, bwidth)
        return bins
    

    
def main():
    print("passed")


if __name__ == "__main__":
    main()
