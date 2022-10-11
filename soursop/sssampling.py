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
from typing import List, Union
import pathlib
import os
import soursop.ssutils 
from .ssexceptions import SSException
from .sstrajectory import parallel_load_trjs
from glob import glob
from natsort import natsorted
import fnmatch
import seaborn as sns
import matplotlib.pyplot as plt

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
    if p.ndim == 3 and q.ndim == 3:
        hellingers = np.sqrt(np.sum((np.sqrt(p) - np.sqrt(q)) ** 2, axis=2)) / np.sqrt(2)
    elif p.ndim == 2 and q.ndim == 2:
        hellingers = np.sqrt(np.sum((np.sqrt(p) - np.sqrt(q)) ** 2, axis=1)) / np.sqrt(2)
    else:
        hellingers = np.sqrt(np.sum((np.sqrt(p) - np.sqrt(q)) ** 2)) / np.sqrt(2)
    return hellingers


def rel_entropy(p : np.ndarray, q : np.ndarray) -> np.ndarray:
    """Computes the relative entropy between two probability distributions p and q."""
    if p.ndim == 3 and q.ndim == 3:
        relative_entropy = np.sum(rel_entr(p,q), axis=2)
    elif p.ndim == 2 and q.ndim == 2:
        relative_entropy = np.sum(rel_entr(p,q), axis=1)
    else:
        relative_entropy = np.sum(rel_entr(p,q))
    return relative_entropy
    

def glob_traj_paths(root_dir : Union[str, pathlib.Path], num_reps : int, mode=None, traj_name="__traj.xtc",top_name="__START.pdb", exclude_dirs=None):
    """
    This function assembles the list of trajectory and topology paths for a set of simulations.

    Parameters
    ----------
    root_dir : Union[str, pathlib.path]
        Filepath or list of file paths
    num_reps : int
        number of replicates - will iterate over directories from [1,num_reps+1].
    mode : str, optional
        if "mega", the globbing will iterate over directories labeled both "coil_start" and "helical_start", by default None
    traj_name : str, optional
        trajectory filename, by default "__traj.xtc"
    top_name : str, optional
        topology filename, by default "__START.pdb"
    
    """
    top_paths, traj_paths = [], []
    if str(mode).lower() == "mega":
        cwd = pathlib.Path(f"{root_dir}").absolute().resolve()
        for directory in ["coil_start","helical_start"]:
            basepath = os.path.join(cwd,directory)
            for rep in range(1,num_reps+1):
                traj_paths.extend(glob(f"{basepath}/{rep}/{traj_name}"))
                top_paths.extend(glob(f"{basepath}/{rep}/{top_name}"))
    else:
        if not exclude_dirs:
            exclude_dirs = ["eq","FULL"]

        for root, dirs, files in os.walk(root_dir):
            if os.path.basename(root) in exclude_dirs:
                continue
            for filename in fnmatch.filter(files, traj_name):
                traj_paths.append(pathlib.Path(os.path.join(root, filename)).absolute().as_posix())
            for filename in fnmatch.filter(files, top_name):
                top_paths.append(pathlib.Path(os.path.join(root, filename)).absolute().as_posix())

    return natsorted(top_paths), natsorted(traj_paths)
class SamplingQuality:
    """Compare the sampling quality for a trajectory relative to some arbitrary referene model, usually a polymer limiting model.
    """
 
    def __init__(self, traj_list : List[str], 
                       polymer_model_traj_list : List[str],
                       top_file : str,
                       polymer_top : str,
                       method : str, 
                       bwidth : float = 0.2617993877991494,
                       proteinID : int = 0,
                       n_cpus : int = None,
                       truncate : bool = False,
                       **kwargs : dict,
                ):
        """_summary_

        Parameters
        ----------
        traj_list : List[str]
            a list of the trajectories associated with the simulated trajectories.
        polymer_model_traj_list : List[str]
            a list of the trajectories associated with the limiting polymer model.
        top_file : str
            path to the simulated trajectories topology file.
        polymer_top : str
            path to the polymer model's topology file.
        method : str
            The method used to compute the hellingers distance between the simulated trajectories and the polymer limiting model.
            options include: 'dihedral' and 'rmsd' or 'p_vects' [not currently implemented]
        bwidth : float, optional
            bin width parameter for segmenting histogrammed data into buckes\
            by default 0.2617993877991494 which corresponds to 15 degrees.
        proteinID : int, optional
            The ID of the protein where the ID is the proteins position
            in the ``self.proteinTrajectoryList`` list, by default 0.
        n_cpus : int, optional
            number of CPUs to use for parallel loading of SSTrajectory objects, by default None, which uses all available threads.
        truncate : bool, optional
            if True, will truncate trajectories arrays 

        Raises
        ------
        NotImplementedError
            _description_
        SSException
            _description_
        """

        super(SamplingQuality, self).__init__()
        self.traj_list = traj_list
        self.polymer_model_traj_list = polymer_model_traj_list
        self.top = top_file
        self.polymer_top = polymer_top
        self.proteinID = proteinID
        self.method = method
        self.bwidth = bwidth
        self.n_cpus = n_cpus
        self.truncate = truncate
        
        if not self.n_cpus:
            self.n_cpus = os.cpu_count()
        
        # Should probably add option to pass trajectories directly, and then also check for that optionality here 
        # best way to do this? idk alex halppp
        self.trajs = parallel_load_trjs(self.traj_list, top=self.top, n_procs=self.n_cpus,**kwargs)
        self.polymer_trajs = parallel_load_trjs(self.polymer_model_traj_list, top=self.polymer_top, n_procs=self.n_cpus, **kwargs)

        if str(self.method).lower() == "rmsd" or str(self.method).lower() == "p_vects":
            raise NotImplementedError("This functionality has not been implemented yet")

        soursop.ssutils.validate_keyword_option(method, ['dihedral', 'rmsd', 'p_vects'], 'method')
        if str(self.method).lower() == "dihedral":
            if self.bwidth > 2*np.pi or not (self.bwidth > 0):
                raise SSException(f'The bwidth parameter must be between 2*pi and greater than 0. Received {self.bwidth}')
            
            if truncate:
                lengths = []
                for trj, pol_trj in zip(self.trajs, self.polymer_trajs):
                    lengths.append([trj.n_frames,pol_trj.n_frames])
                self.min_length = np.min(lengths) - 1
                self.psi_angles, self.polymer_psi_angles, self.phi_angles, self.polymer_phi_angles = self.__compute_dihedrals(proteinID=self.proteinID, truncate=self.truncate)
            else:
                self.psi_angles, self.polymer_psi_angles, self.phi_angles, self.polymer_phi_angles = self.__compute_dihedrals(proteinID=self.proteinID)
        
    def __compute_dihedrals(self, proteinID : int = 0, truncate : bool = None) -> np.ndarray:
        """internal function to computes the phi/psi backbone dihedrals at a given index (proteinID) in the proteinTrajectoryList of an SSTrajectory.

        Parameters
        ----------
        proteinID : int, optional
            The ID of the protein where the ID is the proteins position
            in the ``self.proteinTrajectoryList`` list, by default 0.
        truncate  : bool, optional
            If True, will truncate each numpy array such that the longest array is the same length as the shortest trajectory.
        
        Returns
        -------
        np.ndarray
            Returns the psi and phi backbone dihedrals for the simulated trajectory and the limiting polyer model.
        """
        psi_angles = []
        phi_angles = []
        polymer_psi_angles = []
        polymer_phi_angles = []
        
        for trj, pol_trj in zip(self.trajs, self.polymer_trajs):
            psi_angles.append(trj.proteinTrajectoryList[proteinID].get_angles("psi")[1])
            phi_angles.append(trj.proteinTrajectoryList[proteinID].get_angles("phi")[1])
            polymer_psi_angles.append(pol_trj.proteinTrajectoryList[proteinID].get_angles("psi")[1])
            polymer_phi_angles.append(pol_trj.proteinTrajectoryList[proteinID].get_angles("phi")[1])
        
        if truncate:
            trunc_psi_angles, trunc_polymer_psi_angles, trunc_phi_angles,trunc_polymer_phi_angles = self.__truncate_arrays_by_min_length(psi_angles, polymer_psi_angles, phi_angles, polymer_phi_angles)
            return np.array([trunc_psi_angles, trunc_polymer_psi_angles, trunc_phi_angles, trunc_polymer_phi_angles])
        else:
            return np.array([psi_angles, polymer_psi_angles, phi_angles, polymer_phi_angles])

    def __truncate_arrays_by_min_length(self, 
                                        psi_angles : List[np.ndarray], 
                                        polymer_psi_angles : List[np.ndarray], 
                                        phi_angles : List[np.ndarray],
                                        polymer_phi_angles : List[np.ndarray]
                                        ) -> np.ndarray:
        """internal function to truncate arrays to ensure they're all the same shape. This is useful for when 
        performing analysis on trajectories that haven't completed, and have a varied number of conformations.

        Parameters
        ----------
        psi_angles : List[np.ndarray]
            _description_
        polymer_psi_angles : List[np.ndarray]
            _description_
        phi_angles : List[np.ndarray]
            _description_
        polymer_phi_angles : List[np.ndarray]
            _description_

        Returns
        -------
        np.ndarray
            _description_
        """
        trunc_psi_angles, trunc_phi_angles, trunc_polymer_psi_angles, trunc_polymer_phi_angles = [],[],[],[]
        for psi, pol_psi, phi, pol_phi in zip(psi_angles, polymer_psi_angles, phi_angles, polymer_phi_angles):
                trunc_psi_angles.append(psi[:, self.min_length])
                trunc_polymer_psi_angles.append(pol_psi[:, self.min_length])
                trunc_phi_angles.append(phi[:, self.min_length])
                trunc_polymer_phi_angles.append(pol_phi[:, self.min_length])
                
        return np.array([trunc_psi_angles,trunc_polymer_psi_angles,trunc_phi_angles,trunc_polymer_phi_angles])

    def compute_dihedral_hellingers(self) -> np.ndarray:
        """Compute the hellingers distance for both the phi and psi angles between a set of trajectories.

        Returns
        -------
        np.ndarray
            The hellinger distances between the probability density distributions for the phi and psi angles for a set of trajectories.
        """
        bins = self.get_degree_bins()

        psi_trj_pdfs = self.compute_pdf(self.psi_angles, bins=bins)
        psi_pol_trj_pdfs = self.compute_pdf(self.polymer_psi_angles, bins=bins)
        
        phi_trj_pdfs = self.compute_pdf(self.phi_angles, bins=bins)
        phi_pol_trj_pdfs = self.compute_pdf(self.polymer_phi_angles, bins=bins)

        psi_hellingers = hellinger_distance(psi_trj_pdfs, psi_pol_trj_pdfs)
        phi_hellingers = hellinger_distance(phi_trj_pdfs, phi_pol_trj_pdfs)

        return np.array([phi_hellingers, psi_hellingers])

    def compute_dihedral_rel_entropy(self) -> np.ndarray:
        """Compute the relative entropy for both the phi and psi angles between a set of trajectories.

        Returns
        -------
        np.ndarray
            The relative entropy between the probability density distributions for the phi and psi angles for a set of trajectories.
        """
        bins = self.get_degree_bins()

        psi_trj_pdfs = self.compute_pdf(self.psi_angles, bins=bins)
        psi_pol_trj_pdfs = self.compute_pdf(self.polymer_psi_angles, bins=bins)
        
        phi_trj_pdfs = self.compute_pdf(self.phi_angles, bins=bins)
        phi_pol_trj_pdfs = self.compute_pdf(self.polymer_phi_angles, bins=bins)

        psi_rel_entr = rel_entropy(psi_trj_pdfs, psi_pol_trj_pdfs)
        phi_rel_entr = rel_entropy(phi_trj_pdfs, phi_pol_trj_pdfs)

        return np.array([phi_rel_entr, psi_rel_entr])

    def compute_pdf(self, arr : np.ndarray, bins : np.ndarray) -> np.ndarray:
        """
        Computes a probability density by constructing a histogram of the data array with a specified set of bins. 

        Parameters
        ----------
        arr : np.ndarray
            A vector of shape (n_res x n_frames) or (traj x n_res x frames)
        bins : np.ndarray
            The set of bin edges that specify the range for the histogram buckets.

        Returns
        -------
        np.ndarray
            Returns a set of histograms of the probabilities densities for each residue in the amino acid sequence.
            Shape (n_res, len(bins) - 1) 
        """
        # Lambda function is used to ignore the bin edges returned by np.histogram at index 1
        if arr.ndim == 3:
            pdf = np.apply_along_axis(lambda col: np.histogram(col, bins=bins, density=True)[0], axis=2, arr=arr)
        else:
            pdf = np.apply_along_axis(lambda col: np.histogram(col, bins=bins, density=True)[0], axis=1, arr=arr)
        return pdf


    def get_radian_bins(self) -> np.ndarray:
        """Returns the edges of the bins in radians

        Returns
        -------
        np.ndarray
            an array of the bin edges in radians
        """
        bwidth = self.bwidth
        bins = np.arange(-np.pi, np.pi+bwidth, bwidth)
        return bins
        
    def get_degree_bins(self) -> np.ndarray:
        """Returns the edges of the bins in degrees

        Returns
        -------
        np.ndarray
            an array of the bin edges in degrees
        """
        bwidth = np.rad2deg(self.bwidth)
        bins = np.arange(-180, 180+bwidth, bwidth)
        return bins

    def plot_phi_psi_hellingers(self,figsize=(10,15),cmap=None, **kwargs):
        """plot heatmaps for phi and psi hellingers distances.\
        Optional keyword arguments are passed to 'plt.subplots'

        Parameters
        ----------
        figsize : tuple, optional
            dimensions of the figure to be rendered, by default (10,15)
        cmap : matplotlib colormap, optional
            The matplotlib colormap to be used for plotting the figure, by default None
        """
        if cmap == None:
            cmap = sns.color_palette("light:b", as_cmap=True)

        phi_hellingers, psi_hellingers = self.compute_dihedral_hellingers()
        fig, axes = plt.subplots(2,1,figsize=figsize,sharex=True,sharey=True, **kwargs)
        for i,ax in enumerate(axes):
            if i == 0:
                ax.set_xticks(np.arange(0,phi_hellingers[:,:].shape[1]+1))
                ax.set_yticks(np.arange(0,phi_hellingers[:,:].shape[0]+1))

                ax.set_xticklabels(np.arange(1,phi_hellingers[:,:].shape[1]+2),fontsize=16)
                ax.set_yticklabels(np.arange(1,phi_hellingers[:,:].shape[0]+2),fontsize=16)

                # ax.set_xlabel("Residue index",fontsize=24)
                ax.set_ylabel("Trajectory Index",fontsize=24)

                ax.set_title("Phi Hellinger's Distance",fontsize=24)
                im1 = ax.imshow(phi_hellingers[:,:],cmap=cmap)
            else:
                ax.set_xticks(np.arange(0,psi_hellingers[:,:].shape[1]+1))
                ax.set_yticks(np.arange(0,psi_hellingers[:,:].shape[0]+1))

                ax.set_xticklabels(np.arange(1,psi_hellingers[:,:].shape[1]+2),fontsize=16)
                ax.set_yticklabels(np.arange(1,psi_hellingers[:,:].shape[0]+2),fontsize=16)

                ax.set_xlabel("Residue index",fontsize=24)
                ax.set_ylabel("Trajectory Index",fontsize=24)

                ax.set_title("Psi Hellinger's Distance",fontsize=24)
                im2 = ax.imshow(psi_hellingers[:,:],cmap=cmap)
        plt.tight_layout()
        fig.colorbar(im1, ax=axes.ravel().tolist())

    @property
    def polymer_pdfs(self):
        """property for getting the pdfs computed from the phi/psi angles respectively

        Returns
        -------
        np.ndarray
            pdfs computed from the phi and psi angles with the specified bins.
        """
        bins = self.get_degree_bins()
        pol_phi_pdf = self.compute_pdf(self.phi_angles,bins=bins)
        pol_psi_pdf = self.compute_pdf(self.psi_angles,bins=bins)
        return np.array([pol_phi_pdf, pol_psi_pdf])

    @property
    def trj_pdfs(self):
        """property for getting the pdfs computed from the phi/psi angles respectively

        Returns
        -------
        np.ndarray
            pdfs computed from the phi and psi angles with the specified bins.
        """
        bins = self.get_degree_bins()
        trj_phi_pdf = self.compute_pdf(self.polymer_phi_angles, bins=bins)
        trj_psi_pdf = self.compute_pdf(self.polymer_psi_angles, bins=bins)
        return np.array([trj_phi_pdf, trj_psi_pdf])
    
    @property
    def hellingers_distances(self):
        """property for getting the pdfs computed from the phi/psi angles respectively

        Returns
        -------
        np.ndarray
            hellingers distance computed from the phi and psi angles with the specified bins.
        """
        return self.compute_dihedral_hellingers()
