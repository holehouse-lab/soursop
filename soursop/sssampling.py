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
from typing import List, Union, Tuple
import pathlib
import os
from soursop import ssutils
from .ssexceptions import SSException
from .sstrajectory import parallel_load_trjs, SSTrajectory
from glob import glob
from natsort import natsorted
import fnmatch
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from xhistogram.core import histogram
import itertools

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
    """Computes the relative entropy between two probability distributions p and q.
    Parameters
    ----------
    p : np.ndarray
        a probability distribution function, or series of probabiliy distribution functions, to compute the relative entropy distance on.
        p and q must be the same shape

    q : np.ndarray
        a probability distribution function, or series of probabiliy distribution functions, to compute the relative entropy on
        p and q must be the same shape
    
    Returns
    -------
    np.ndarray
        The relative entropy for each residue in the sequence
    """
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
        else: 
            np.unique(exclude_dirs.extend(["eq","FULL"]))

        for root, dirs, files in os.walk(root_dir):
            if os.path.basename(root) in exclude_dirs:
                continue
            for filename in fnmatch.filter(files, traj_name):
                traj_paths.append(pathlib.Path(os.path.join(root, filename)).absolute().resolve().as_posix())
            for filename in fnmatch.filter(files, top_name):
                top_paths.append(pathlib.Path(os.path.join(root, filename)).absolute().resolve().as_posix())

    return natsorted(top_paths), natsorted(traj_paths)

class SamplingQuality:
    """Compare the sampling quality for a trajectory relative to some arbitrary referene model, usually a polymer limiting model.
    """
 
    def __init__(self, traj_list : List[str], 
                       polymer_model_traj_list : List[str],
                       top_file : str,
                       polymer_top : str,
                       method : str, 
                       bwidth : float = np.deg2rad(15),
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
            bin width parameter for segmenting histogrammed data into buckets, by default 15 degrees.
        proteinID : int, optional
            The ID of the protein where the ID is the proteins position
            in the ``self.proteinTrajectoryList`` list, by default 0.
        n_cpus : int, optional
            number of CPUs to use for parallel loading of SSTrajectory objects, by default None, which uses all available threads.
        truncate : bool, optional
            if True, will slice all the trajectories such that they're all of the same minimum length. 

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
        
        if truncate:
            lengths = []
            for trj, pol_trj in zip(self.trajs, self.polymer_trajs):
                lengths.append([trj.n_frames, pol_trj.n_frames])

            # shift frames for np.array indexing purposes
            self.min_length = np.min(lengths) - 1
            
            print(lengths)
            f"The maximum trajectory length shared between all replicates is: {self.min_length}"
            self.trajs, self.polymer_trajs = self.__truncate_trajectories()

        if str(self.method).lower() == "rmsd" or str(self.method).lower() == "p_vects":
            raise NotImplementedError("This functionality has not been implemented yet")

        ssutils.validate_keyword_option(method, ['dihedral', 'rmsd', 'p_vects'], 'method')
        if str(self.method).lower() == "dihedral":
            if self.bwidth > 2*np.pi or not (self.bwidth > 0):
                raise SSException(f'The bwidth parameter must be between 2*pi and greater than 0. Received {self.bwidth}')
        
            self.psi_angles, self.polymer_psi_angles, self.phi_angles, self.polymer_phi_angles = self.__compute_dihedrals(proteinID=self.proteinID)

    def __truncate_trajectories(self) -> Tuple[List[SSTrajectory], List[SSTrajectory]]:
        """Internal function used to truncate the lengths of trajectories such that every trajectory has the same number of total frames.
        Useful for intermediary analysis of ongoing simulations.

        Returns
        -------
        Tuple[List[SSTrajectory], List[SSTrajectory]]
            A tuple containing two lists of SSTrajectory objects.\
                The first index corresponds to the empirical trajectories.\
                The second corresonds to the reference model - e.g., the polymer limiting model.
        """
        temp_trajs = []
        temp_pol_trjs = [] 
        for trj, pol_trj in zip(self.trajs,self.polymer_trajs):
            temp_trajs.append(
                    SSTrajectory(TRJ=trj.proteinTrajectoryList[self.proteinID].traj[0:self.min_length])
                )
            temp_pol_trjs.append(
                    SSTrajectory(TRJ=pol_trj.proteinTrajectoryList[self.proteinID].traj[0:self.min_length])
                )

        return (temp_trajs, temp_pol_trjs)


    def __compute_dihedrals(self, proteinID : int = 0) -> np.ndarray:
        """internal function to computes the phi/psi backbone dihedrals at a given index (proteinID) in the proteinTrajectoryList of an SSTrajectory.

        Parameters
        ----------
        proteinID : int, optional
            The ID of the protein where the ID is the proteins position
            in the ``self.proteinTrajectoryList`` list, by default 0.
        
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
        
        return np.array((psi_angles, polymer_psi_angles, phi_angles, polymer_phi_angles))

    def __compute_frac_helicity(self, proteinID : int = 0) -> np.ndarray:
        """internal function to computes the per residue fractional helicity at a given index (proteinID) in the proteinTrajectoryList of an SSTrajectory.

        Parameters
        ----------
        proteinID : int, optional
            The ID of the protein where the ID is the proteins position
            in the ``self.proteinTrajectoryList`` list, by default 0.
        
        Returns
        -------
        np.ndarray
            Returns the fractional helicity for the simulated trajectory and the reference model.
        """
        trj_helicity = []
        reference_helicity = []
                
        for trj, ref_traj in zip(self.trajs, self.polymer_trajs):
            trj_helicity.append(trj.proteinTrajectoryList[proteinID].get_secondary_structure_DSSP()[1])
            reference_helicity.append(ref_traj.proteinTrajectoryList[proteinID].get_secondary_structure_DSSP()[1])
            
        return np.array((trj_helicity, reference_helicity))


    def compute_dihedral_hellingers(self) -> np.ndarray:
        """Compute the hellingers distance for both the phi and psi angles between a set of trajectories.

        Returns
        -------
        np.ndarray
            The hellinger distances between the probability density distributions for the phi and psi angles for a set of trajectories.
        """
        bins = self.get_degree_bins()
        phi_trj_pdfs = self.compute_pdf(self.phi_angles, bins=bins)
        phi_pol_trj_pdfs = self.compute_pdf(self.polymer_phi_angles, bins=bins)

        psi_trj_pdfs = self.compute_pdf(self.psi_angles, bins=bins)
        psi_pol_trj_pdfs = self.compute_pdf(self.polymer_psi_angles, bins=bins)

        phi_hellingers = hellinger_distance(phi_trj_pdfs, phi_pol_trj_pdfs)
        psi_hellingers = hellinger_distance(psi_trj_pdfs, psi_pol_trj_pdfs)

        return np.array((phi_hellingers, psi_hellingers))

    def compute_dihedral_rel_entropy(self) -> np.ndarray:
        """Compute the relative entropy for both the phi and psi angles between a set of trajectories.

        Returns
        -------
        np.ndarray
            The relative entropy between the probability density distributions for the phi and psi angles for a set of trajectories.
        """
        bins = self.get_degree_bins()
        
        phi_trj_pdfs = self.compute_pdf(self.phi_angles, bins=bins)
        phi_pol_trj_pdfs = self.compute_pdf(self.polymer_phi_angles, bins=bins)

        psi_trj_pdfs = self.compute_pdf(self.psi_angles, bins=bins)
        psi_pol_trj_pdfs = self.compute_pdf(self.polymer_psi_angles, bins=bins)

        phi_rel_entr = rel_entropy(phi_trj_pdfs, phi_pol_trj_pdfs)
        psi_rel_entr = rel_entropy(psi_trj_pdfs, psi_pol_trj_pdfs)

        return np.array((phi_rel_entr, psi_rel_entr))

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
        # xhistogram is ~2x faster, but introduces depedency - keeping lambda function for legacy for now
        # tbh this isn't that slow with the lambda function after all 
        if arr.ndim == 3:
            # pdf = np.apply_along_axis(lambda col: np.histogram(col, bins=bins, density=True)[0], axis=2, arr=arr)*np.round(np.rad2deg(self.bwidth))
            
            # KEY POINT: multiplying by bin width to convert probability *density* to probabilty *mass*
            # implementation details may have to change here if supporting other methods.
            pdf = histogram(arr,bins=bins,axis=2,density=True)[0]*np.round(np.rad2deg(self.bwidth)) 
        else:
            # pdf = np.apply_along_axis(lambda col: np.histogram(col, bins=bins, density=True)[0], axis=1, arr=arr)*np.round(np.rad2deg(self.bwidth))
            pdf = histogram(arr,bins=bins,axis=1,density=True)[0]*np.round(np.rad2deg(self.bwidth))
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
        # have to round the conversion to handle floating point error so we get the right bins
        bwidth = np.round(np.rad2deg(self.bwidth))
        bins = np.arange(-180, 180+bwidth, bwidth)
        return bins

    def plot_phi_psi_heatmap(self, metric : str ="hellingers",
                                   figsize=(40,20), 
                                   annotate=True,
                                   cmap=None,
                                   vmin=0.0,
                                   vmax=1.0,
                                   filename : str="sampling_quality.png",
                                   save_dir=None,
                                   **kwargs,        
                            ):
        """Plot heatmaps for phi and psi metrics.\
        Optional keyword arguments are passed to 'plt.subplots'

        Parameters
        ----------
        metric : str, optional
            The distance metric to use - either "hellingers" or "relative entropy", by default "hellingers"
            Note: relative entropy is a divergence, and not a true distance metric. 
        figsize : tuple, optional
            dimensions of the figure to be rendered, by default (40,20)
        annotate : bool, optional
            Whether to display the data values from the metric in the plot, by default True
        cmap : str, optional
            The matplotlib colormap to be used for plotting the figure, by default None
        vmin : float, optional
            Minimum anchor point for colorbar, by default 0.0
        vmax : float, optional
            Maximum anchor point for colorbar, by default 1.0
        filename : str, optional
            _description_, by default "sampling_quality.png"
        save_dir : _type_, optional
            _description_, by default None

        Raises
        ------
        NotImplementedError
            _description_
        """
        if metric == "hellingers":
            phi_metric, psi_metric = self.compute_dihedral_hellingers()
        elif metric == "relative entropy":
            phi_metric, psi_metric = self.compute_dihedral_rel_entropy()
        else:
            raise NotImplementedError(f"The metric: {metric} is not implemented.")
        
        if cmap == None:
            cmap = sns.color_palette("light:b", as_cmap=True)

        fig, (ax1,ax2,axcb) = plt.subplots(1,3,figsize=figsize, gridspec_kw={'width_ratios':[1,1,0.05]},**kwargs)
        
        ax1.get_shared_y_axes().join(ax2)
        g1 = sns.heatmap(phi_metric,annot=annotate,annot_kws={"fontsize":24,"color":"k"},vmin=vmin,vmax=vmax,cmap=cmap,cbar=False,ax=ax1)
        g1.set_title(f'Phi {metric}',fontsize=36)
        g1.set_ylabel('Trajectory Index',fontsize=36)
        g1.set_xlabel('Residue Index',fontsize=36)
        g1.set_xticks(np.arange(0,phi_metric[:,:].shape[1])+0.5)
        g1.set_yticks(np.arange(0,phi_metric[:,:].shape[0])+0.5)

        g1.set_xticklabels(np.arange(1,phi_metric[:,:].shape[1]+1),fontsize=36)
        g1.set_yticklabels(np.arange(1,phi_metric[:,:].shape[0]+1),fontsize=36)

        g2 = sns.heatmap(psi_metric,annot=annotate,annot_kws={"fontsize":24,"color":"k"},vmin=vmin,vmax=vmax,cmap=cmap,cbar_ax=axcb,ax=ax2)
        g2.set_title(f'Psi {metric}',fontsize=36)
        g2.set_ylabel('Trajectory Index',fontsize=36)
        g2.set_xlabel('Residue Index',fontsize=36)
        g2.set_xticks(np.arange(0,psi_metric[:,:].shape[1])+0.5)
        g2.set_yticks(np.arange(0,psi_metric[:,:].shape[0])+0.5)

        g2.set_xticklabels(np.arange(1,psi_metric[:,:].shape[1]+1),fontsize=36)
        g2.set_yticklabels(np.arange(1,psi_metric[:,:].shape[0]+1),fontsize=36)
        axcb.tick_params(labelsize=36)
        plt.tight_layout()
        if save_dir is not None:
            os.makedirs(save_dir,exist_ok=True)
            outpath = os.path.join(save_dir,filename)
            fig.savefig(f"{outpath}",dpi=300)

    def get_all_to_all_trj_comparisons(self, metric : str ="hellingers") -> Tuple[pd.DataFrame,pd.DataFrame]:
        """function to aggregate an all-to-all comparison of pdfs 

        Returns
        -------
        Tuple[pd.DataFrame,pd.DataFrame]
            all-to-all trajectory comparisons for the hellingers distances in 'self.trj_pdfs'
        """
        phi_pdfs = self.trj_pdfs[0]
        psi_pdfs = self.trj_pdfs[1]
        phi_combinations = np.transpose(np.array(tuple(itertools.combinations(phi_pdfs,2))), axes=[1,0,2,3])
        psi_combinations = np.transpose(np.array(tuple(itertools.combinations(psi_pdfs,2))), axes=[1,0,2,3])
        if metric == "hellingers": 
            phi_metric = hellinger_distance(phi_combinations[0],phi_combinations[1])
            psi_metric = hellinger_distance(psi_combinations[0],psi_combinations[1])
        elif metric == "relative entropy":
            phi_metric = rel_entropy(phi_combinations[0],phi_combinations[1])
            psi_metric = rel_entropy(psi_combinations[0],psi_combinations[1])
        else: 
            raise NotImplementedError(f"The metric: {metric} is not implemented.")
            
        return pd.DataFrame(phi_metric), pd.DataFrame(psi_metric)

    def __get_all_to_all_ev_comparison(self) -> Tuple[pd.DataFrame,pd.DataFrame]:
        """internal function to aggregate an all-to-all comparison of pdfs 

        Returns
        -------
        Tuple[pd.DataFrame,pd.DataFrame]
            _description_
        """
        phi_pdfs = self.trj_pdfs[0]
        phi_ev_pdfs = self.polymer_pdfs[0]
        psi_pdfs = self.trj_pdfs[1]
        psi_ev_pdfs = self.polymer_pdfs[1]
        # if using this should revisit and do combinations like above
        phi_cartesian_prod = np.transpose(np.array(tuple(itertools.product(phi_pdfs, phi_ev_pdfs))),axes=[1,0,2,3])
        psi_cartesian_prod = np.transpose(np.array(tuple(itertools.product(psi_pdfs, psi_ev_pdfs))),axes=[1,0,2,3])
        phi_hellingers = hellinger_distance(phi_cartesian_prod[0],phi_cartesian_prod[1])
        psi_hellingers = hellinger_distance(psi_cartesian_prod[0],psi_cartesian_prod[1])
        return pd.DataFrame(phi_hellingers), pd.DataFrame(psi_hellingers)

    def plot_phi_psi_stripplot(self, metric : str ="hellingers",
                                    figsize : tuple =(20,30), 
                                    jitter : float = 0.4,
                                    cmap=None,
                                    vmin : float =0.0,
                                    vmax : float =1.0,
                                    filename : str="sampling_quality.png",
                                    save_dir : str =None,
                                    fontsize : int = 16,
                                    **kwargs,   
                            ):
        if metric == "hellingers":
            ax_label = "Hellinger's Distance"
            phi_metric, psi_metric = self.compute_dihedral_hellingers()
            phi_df, psi_df = pd.DataFrame(phi_metric), pd.DataFrame(psi_metric)
            all_phi_df, all_psi_df = self.__get_all_to_all_trj_comparisons()
        elif metric == "relative entropy":
            ax_label = "Relative Entropy"
            phi_metric, psi_metric = self.compute_dihedral_rel_entropy()
            phi_df, psi_df = pd.DataFrame(phi_metric), pd.DataFrame(psi_metric)
        else:
            raise NotImplementedError(f"The metric: {metric} is not implemented.")

        x_ticks = np.arange(0,phi_df.shape[1])
        x_ticklabels = np.arange(1,phi_df.shape[1]+1)
        y_ticks = np.arange(.2,1.2,0.2)
        y_ticklabels = np.round(np.arange(.2,1.2,0.2),1)

        # ref not being juplotted right now
        trj_frac_helicity, ref_traj_helicity = self.__compute_frac_helicity()
        helicity_df = pd.DataFrame(trj_frac_helicity)

        fig, axes = plt.subplots(5,1,facecolor="w",figsize=figsize)
        sns.stripplot(data=phi_df,jitter=jitter,ax=axes[0])
        axes[0].set_title("Phi Dihedrals",fontsize=fontsize)
        axes[0].set_xlabel("Residue",fontsize=fontsize)
        axes[0].set_ylabel(f"{ax_label}",fontsize=fontsize)
        axes[0].set_xticks(x_ticks,fontsize=fontsize)
        axes[0].set_xticklabels(x_ticklabels,fontsize=fontsize)
        axes[0].set_yticks(y_ticks,fontsize=fontsize)
        axes[0].set_yticklabels(y_ticklabels,fontsize=fontsize)

        sns.stripplot(data=psi_df,jitter=jitter,ax=axes[1])
        axes[1].set_title("Psi Dihedrals",fontsize=fontsize)
        axes[1].set_xlabel("Residue",fontsize=fontsize)
        axes[1].set_ylabel(f"{ax_label}",fontsize=fontsize)
        axes[1].set_xticks(x_ticks,fontsize=fontsize)
        axes[1].set_xticklabels(x_ticklabels,fontsize=fontsize)
        axes[1].set_yticks(y_ticks,fontsize=fontsize)
        axes[1].set_yticklabels(y_ticklabels,fontsize=fontsize)

        sns.stripplot(data=all_phi_df,jitter=jitter,ax=axes[2])        
        axes[2].set_title("All-to-All",fontsize=fontsize)
        axes[2].set_ylabel(f"{ax_label}",fontsize=fontsize)
        axes[2].set_xlabel("Residue",fontsize=fontsize)
        axes[2].set_xticks(x_ticks,fontsize=fontsize)
        axes[2].set_xticklabels(x_ticklabels,fontsize=fontsize)
        axes[2].set_yticks(y_ticks,fontsize=fontsize)
        axes[2].set_yticklabels(y_ticklabels,fontsize=fontsize)

        sns.stripplot(data=all_psi_df,jitter=jitter,ax=axes[3])        
        axes[3].set_title("All-to-All",fontsize=fontsize)
        axes[3].set_xlabel("Residue",fontsize=fontsize)
        axes[3].set_ylabel(f"{ax_label}",fontsize=fontsize)
        axes[3].set_xticks(x_ticks,fontsize=fontsize)
        axes[3].set_xticklabels(x_ticklabels,fontsize=fontsize)
        axes[3].set_yticks(y_ticks,fontsize=fontsize)
        axes[3].set_yticklabels(y_ticklabels,fontsize=fontsize)

        sns.stripplot(data=helicity_df,jitter=jitter,ax=axes[4])
        axes[4].set_xlabel("Residue",fontsize=fontsize)
        axes[4].set_ylabel(f"Fractional Helicity",fontsize=fontsize)
        axes[4].set_xticks(x_ticks,fontsize=fontsize)
        axes[4].set_xticklabels(x_ticklabels,fontsize=fontsize)
        axes[4].set_yticks(y_ticks,fontsize=fontsize)
        axes[4].set_yticklabels(y_ticklabels,fontsize=fontsize)

        plt.tight_layout()
        if save_dir is not None:
            os.makedirs(save_dir,exist_ok=True)
            outpath = os.path.join(save_dir,filename)
            fig.savefig(f"{outpath}",dpi=300)
        
        return (fig, axes)

    @property
    def trj_pdfs(self):
        """property for getting the pdfs computed from the phi/psi angles respectively

        Returns
        -------
        np.ndarray
            pdfs computed from the phi and psi angles with the specified bins.
        """
        bins = self.get_degree_bins()
        pol_phi_pdf = self.compute_pdf(self.phi_angles,bins=bins)
        pol_psi_pdf = self.compute_pdf(self.psi_angles,bins=bins)
        return np.array((pol_phi_pdf, pol_psi_pdf))

    @property
    def polymer_pdfs(self):
        """property for getting the pdfs computed from the phi/psi angles respectively

        Returns
        -------
        np.ndarray
            pdfs computed from the phi and psi angles with the specified bins.
        """
        bins = self.get_degree_bins()
        trj_phi_pdf = self.compute_pdf(self.polymer_phi_angles, bins=bins)
        trj_psi_pdf = self.compute_pdf(self.polymer_psi_angles, bins=bins)
        return np.array((trj_phi_pdf, trj_psi_pdf))
    
    @property
    def hellingers_distances(self):
        """property for getting the hellingers distances computed from the phi/psi angles respectively

        Returns
        -------
        np.ndarray
            hellingers distance computed from the phi and psi angles with the specified bins.
        """
        return self.compute_dihedral_hellingers()

    @property
    def fractional_helicity(self):
        """property for getting the per residue fractional helicity for all trajectories

        Returns
        -------
        np.ndarray
            The per residue fractional helicity for each trajectory in self.trajs and self.polymer_trajs.
        """
        return self.__compute_frac_helicity() 