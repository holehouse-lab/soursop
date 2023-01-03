#  _____  ____  _    _ _____   _____  ____  _____
# / ____|/ __ \| |  | |  __ \ / ____|/ __ \|  __ \
# | (___ | |  | | |  | | |__) | (___ | |  | | |__)|
# \___ \| |  | | |  | |  _  / \___ \| |  | |  ___/
# ____) | |__| | |__| | | \ \ ____) | |__| | |
#|_____/ \____/ \____/|_|  \_\_____/ \____/|_|
# Jeffrey M. Lotthammer (Holehouse Lab)
# Alex Holehouse (Pappu Lab and Holehouse Lab) and Jared Lalmansing (Pappu lab)
# Simulation analysis package
# Copyright 2014 - 2022
##
import fnmatch
import itertools
import os
import pathlib
from glob import glob
from typing import List, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import transforms
from natsort import natsorted
from scipy.special import rel_entr
# from xhistogram.core import histogram

from soursop import ssutils

from .ssexceptions import SSException
from .sstrajectory import SSTrajectory, parallel_load_trjs


def hellinger_distance(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    """
    Computes the hellinger distance between a set probability distributions p and q.
    The hellinger distances is defined by:
        H(P,Q) = r'\frac{1}{\sqrt{2}} \times \sqrt{\sum_{i=1}^{k}(\sqrt{p_i}-\sqrt{q_i})^2}'
    where k is the length of the probability vectors being compared.

    Parameters
    ----------
    p : np.ndarray
        a probability density function, or series of probabiliy density functions,
        to compute the hellingers distance on. p and q must be the same shape

    q : np.ndarray
        a probability density function, or series of probabiliy density functions,
        to compute the hellingers distance on. p and q must be the same shape

    Returns
    -------
    np.ndarray
        The hellingers distance for each residue in the sequence
    """
    if not isinstance(p, np.ndarray):
        p = np.array(p)
    if not isinstance(q, np.ndarray):
        q = np.array(q)

    if p.ndim == 3 and q.ndim == 3:
        hellingers = np.sqrt(np.sum((np.sqrt(p) - np.sqrt(q)) ** 2, axis=2)) / np.sqrt(2)
    elif p.ndim == 2 and q.ndim == 2:
        hellingers = np.sqrt(np.sum((np.sqrt(p) - np.sqrt(q)) ** 2, axis=1)) / np.sqrt(2)
    else:
        hellingers = np.sqrt(np.sum((np.sqrt(p) - np.sqrt(q)) ** 2)) / np.sqrt(2)

    return hellingers


def rel_entropy(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Computes the relative entropy between two probability distributions p and q.

    Parameters
    ----------
    p : np.ndarray
        a probability distribution function, or series of probabiliy distribution functions,
        to compute the relative entropy distance on. p and q must be the same shape

    q : np.ndarray
        a probability distribution function, or series of probabiliy distribution functions, to compute the relative entropy on
        p and q must be the same shape

    Returns
    -------
    np.ndarray
        The relative entropy for each residue in the sequence
    """
    if not isinstance(p, np.ndarray):
        p = np.array(p)
    if not isinstance(q, np.ndarray):
        q = np.array(q)

    if p.ndim == 3 and q.ndim == 3:
        relative_entropy = np.sum(rel_entr(p, q), axis=2)
    elif p.ndim == 2 and q.ndim == 2:
        relative_entropy = np.sum(rel_entr(p, q), axis=1)
    else:
        relative_entropy = np.sum(rel_entr(p, q))
    return relative_entropy


def glob_traj_paths(root_dir: Union[str, pathlib.Path],
                    num_reps=None, mode=None, traj_name="__traj.xtc",
                    top_name="__START.pdb", exclude_dirs=None
                    ):
    """
    This function assembles the list of trajectory and topology paths for a set of simulations.

    Parameters
    ----------
    root_dir : Union[str, pathlib.path]
        Filepath or list of file paths
    mode : str, optional
        if "mega", the globbing will iterate over directories labeled 
        both "coil_start" and "helical_start", by default None
    num_reps : int, optional
        if mode == 'mega' then this flag controls the number of replicates
        in each directory. Iterates over ["coil_start","helical_start"]
        gathering trajectories from child directories [1,num_reps+1], by default None
    traj_name : str, optional
        trajectory filename, by default "__traj.xtc"
    top_name : str, optional
        topology filename, by default "__START.pdb"
    """
    top_paths, traj_paths = [], []
    if str(mode).lower() == "mega":
        cwd = pathlib.Path(f"{root_dir}").absolute().resolve()
        for directory in ["coil_start", "helical_start"]:
            basepath = os.path.join(cwd, directory)
            for rep in range(1, num_reps+1):
                traj_paths.extend(glob(f"{basepath}/{rep}/{traj_name}"))
                top_paths.extend(glob(f"{basepath}/{rep}/{top_name}"))
    else:
        if not exclude_dirs:
            exclude_dirs = ["eq", "FULL"]
        else:
            np.unique(exclude_dirs.extend(["eq", "FULL"]))

        for root, _, files in os.walk(root_dir):
            if os.path.basename(root) in exclude_dirs:
                continue
            for filename in fnmatch.filter(files, traj_name):
                traj_paths.append(pathlib.Path(os.path.join(root, filename)).absolute().resolve().as_posix())
            for filename in fnmatch.filter(files, top_name):
                top_paths.append(pathlib.Path(os.path.join(root, filename)).absolute().resolve().as_posix())

    return natsorted(top_paths), natsorted(traj_paths)


class SamplingQuality:

    def __init__(self, 
                 traj_list: List[str],
                 reference_list: List[str],
                 top_file: str,
                 ref_top: str,
                 method: str,
                 bwidth: float = np.deg2rad(15),
                 proteinID: int = 0,
                 n_cpus: int = None,
                 truncate: bool = False,
                 **kwargs: dict,
                ):
        """SamplingQuality is a class to compare the sampling quality for a set of trajectories
            relative to some referene model. This is usually a polymer limiting model, but could
            be different proteins, mutations, or PTMs.

        Parameters
        ----------
        traj_list : List[str]
            a list of the trajectory paths associated with the simulated trajectories.
        reference_list : List[str]
            a list of the trajectory paths associated with the reference model.
            This is usually a limiting polymer model.
        top_file : str
            path to the simulated trajectories topology file.
        ref_top : str
            path to the reference topology - usually the same topology file
            if using a polymer model, but can differ if reference is e.g., a mutant.
        method : str
            The method used to compute the hellingers distance
            between the simulated trajectories and the polymer limiting model.
            options include: 'dihedral' and 'rmsd' or 'p_vects' [not currently implemented]
        bwidth : float, optional
            bin width parameter for segmenting histogrammed data into buckets,
            by default 15 degrees.
        proteinID : int, optional
            The ID of the protein where the ID is the proteins position
            in the ``self.proteinTrajectoryList`` list, by default 0.
        n_cpus : int, optional
            number of CPUs to use for parallel loading of SSTrajectory objects,
            by default None, which uses all available threads.
        truncate : bool, optional
            if True, will slice all the trajectories such that
            they're all of the same minimum length.
        kwargs : dict, optional
            if provided, key-value pairs will be passed to SSTrajectory for file
            loading. For example, can be useful for passing a stride. 
        Raises
        ------
        NotImplementedError
            Raised when the requested functionality has not yet been implemented.
        SSException
            Raised when the keyword methodology is not part of the validated options.
        """

        super(SamplingQuality, self).__init__()
        self.traj_list = traj_list
        self.reference_list = reference_list
        self.top = top_file
        self.ref_top = ref_top
        self.proteinID = proteinID
        self.method = method
        self.bwidth = bwidth
        self.n_cpus = n_cpus
        self.truncate = truncate

        if not self.n_cpus:
            self.n_cpus = os.cpu_count()

        # check to ensure trajectory lists are the same length and non-zero
        if len(self.traj_list) == 0 or len(self.reference_list) == 0\
            or len(self.traj_list) != len(self.reference_list):
            raise SSException(
                    f"Input trajectory and reference lists must\
                    be non-empty and equal in length.\
                    Received len(traj_list)={len(self.traj_list)},\
                    len(reference_list)={len(self.reference_list)}"
                )

        # weird thing I have to do to prevent issues
        # with multiprocessing parallel loading when there is only 1 trajectory to load
        # trajs/ref_trajs must be a list so they're iterables for __truncate_trajectories
        if len(self.traj_list) == 1:
            self.trajs = []
            self.trajs.append(SSTrajectory(self.traj_list,
                                            pdb_filename=self.top,
                                            **kwargs))
            self.ref_trajs = []
            self.ref_trajs.append(SSTrajectory(self.reference_list,
                                              pdb_filename=self.ref_top,
                                              **kwargs))
        else:
            # if many trajectories, load in parallel
            self.trajs = parallel_load_trjs(self.traj_list,
                                            top=self.top,
                                            n_procs=self.n_cpus,
                                            **kwargs)
            self.ref_trajs = parallel_load_trjs(self.reference_list, 
                                                    top=self.ref_top,
                                                    n_procs=self.n_cpus,
                                                    **kwargs)

        if truncate:
            lengths = []
            for trj, ref_trj in zip(self.trajs, self.ref_trajs):
                lengths.append([trj.n_frames, ref_trj.n_frames])

            # shift frames for np.array indexing purposes
            self.min_length = np.min(lengths)

            self.trajs, self.ref_trajs = self.__truncate_trajectories()
            print(
                f"Successfully truncated.\n\
                The shortest trajectory is: {self.min_length} frames.\
                All trajectories truncated to {self.min_length}"
                )

        if str(self.method).lower() == "rmsd" or str(self.method).lower() == "p_vects":
            raise NotImplementedError("This functionality has not been implemented yet")

        ssutils.validate_keyword_option(method, ['dihedral', 'rmsd', 'p_vects'], 'method')
        if str(self.method).lower() == "dihedral":
            if self.bwidth > 2*np.pi or not self.bwidth > 0:
                raise SSException(
                    f'The bwidth parameter must be between 0 and 2*pi.\
                        Received {self.bwidth}'
                )

            (self.psi_angles,
            self.ref_psi_angles,
            self.phi_angles,
            self.ref_phi_angles) = self.__compute_dihedrals(proteinID=self.proteinID)

    def __truncate_trajectories(self) -> Tuple[List[SSTrajectory], List[SSTrajectory]]:
        """Internal function used to truncate the lengths of trajectories
        such that every trajectory has the same number of total frames.
        Useful for intermediary analysis of ongoing simulations.

        Returns
        -------
        Tuple[List[SSTrajectory], List[SSTrajectory]]
            A tuple containing two lists of SSTrajectory objects.\
            The first index corresponds to the empirical trajectories.\
            The second corresonds to the reference model - e.g.,
            the polymer limiting model.
        """
        temp_trajs = []
        temp_ref_trjs = []
        for trj, ref_trj in zip(self.trajs, self.ref_trajs):
            temp_trajs.append(
                SSTrajectory(TRJ=trj.proteinTrajectoryList[self.proteinID].traj[0:self.min_length])
            )
            temp_ref_trjs.append(
                SSTrajectory(TRJ=ref_trj.proteinTrajectoryList[self.proteinID].traj[0:self.min_length])
            )

        return (temp_trajs, temp_ref_trjs)

    def __compute_dihedrals(self, proteinID: int = 0) -> np.ndarray:
        """internal function to computes the phi/psi backbone dihedrals
        at a given index proteinID in the proteinTrajectoryList of an SSTrajectory.

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
        ref_psi_angles = []
        ref_phi_angles = []

        for trj, ref_trj in zip(self.trajs, self.ref_trajs):
            psi_angles.append(trj.proteinTrajectoryList[proteinID].get_angles("psi")[1])
            phi_angles.append(trj.proteinTrajectoryList[proteinID].get_angles("phi")[1])
            ref_psi_angles.append(ref_trj.proteinTrajectoryList[proteinID].get_angles("psi")[1])
            ref_phi_angles.append(ref_trj.proteinTrajectoryList[proteinID].get_angles("phi")[1])

        return np.array((psi_angles, ref_psi_angles, phi_angles, ref_phi_angles))

    def __compute_frac_helicity(self, proteinID: int = 0) -> np.ndarray:
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

        for trj, ref_traj in zip(self.trajs, self.ref_trajs):
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
        phi_ref_trj_pdfs = self.compute_pdf(self.ref_phi_angles, bins=bins)

        psi_trj_pdfs = self.compute_pdf(self.psi_angles, bins=bins)
        psi_ref_trj_pdfs = self.compute_pdf(self.ref_psi_angles, bins=bins)

        phi_hellingers = hellinger_distance(phi_trj_pdfs, phi_ref_trj_pdfs)
        psi_hellingers = hellinger_distance(psi_trj_pdfs, psi_ref_trj_pdfs)

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
        phi_ref_trj_pdfs = self.compute_pdf(self.ref_phi_angles, bins=bins)

        psi_trj_pdfs = self.compute_pdf(self.psi_angles, bins=bins)
        psi_ref_trj_pdfs = self.compute_pdf(self.ref_psi_angles, bins=bins)

        phi_rel_entr = rel_entropy(phi_trj_pdfs, phi_ref_trj_pdfs)
        psi_rel_entr = rel_entropy(psi_trj_pdfs, psi_ref_trj_pdfs)

        return np.array((phi_rel_entr, psi_rel_entr))

    def compute_pdf(self, arr: np.ndarray, bins: np.ndarray) -> np.ndarray:
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
            Returns a set of histograms of the probabilities 
            densities for each residue in the amino acid sequence.
            Shape (n_res, len(bins) - 1)
        """
        # Lambda function is used to ignore the bin edges returned by np.histogram at index 1
        # xhistogram is ~2x faster, but introduces depedency - keeping lambda function for legacy for now
        # tbh this isn't that slow with the lambda function after all
        
        # if (traj x n_res x frames), histogram axis (2) associated all the frames 
        if arr.ndim == 3:
            pdf = np.apply_along_axis(lambda col: np.histogram(col, bins=bins, density=True)[0], axis=2, arr=arr)*np.round(np.rad2deg(self.bwidth))

            # KEY POINT: multiplying by bin width to convert probability *density* to probabilty *mass*
            # implementation details may have to change here if supporting other methods.
            # pdf = histogram(arr, bins=bins, axis=2, density=True)[0]*np.round(np.rad2deg(self.bwidth))
        # else (n_res x n_frames), histogram axis (1) associated with frames 
        else:
            pdf = np.apply_along_axis(lambda col: np.histogram(col, bins=bins, density=True)[0], axis=1, arr=arr)*np.round(np.rad2deg(self.bwidth))
            # pdf = histogram(arr, bins=bins, axis=1, density=True)[0]*np.round(np.rad2deg(self.bwidth))
        
        return pdf

    def get_all_to_all_trj_comparisons(self, metric: str = "hellingers") -> Tuple[pd.DataFrame, pd.DataFrame]:
        """function to aggregate an all-to-all comparison of pdfs

        Returns
        -------
        Tuple[pd.DataFrame, pd.DataFrame]
            all-to-all trajectory comparisons for the hellingers distances in 'self.trj_pdfs'
        """
        phi_pdfs = self.trj_pdfs[0]
        psi_pdfs = self.trj_pdfs[1]

        if phi_pdfs.shape[0] == 1 or psi_pdfs.shape[0] == 1:
            # if only 1 simulated traj and 1 ref traj all-to-all is just a 1:1 comparison.
            phi_combinations = np.transpose(np.array(tuple(itertools.combinations(phi_pdfs, 1))), axes=[1, 0, 2, 3])
            psi_combinations = np.transpose(np.array(tuple(itertools.combinations(psi_pdfs, 1))), axes=[1, 0, 2, 3])
        else:
            # returned array is (n_combinations, 2, num_resi, num_bins)
            # 2 because it's a pairwise head-to-head comparison of trajectories.
            # transposed for my sanity for indexing leaving final shape as:
            # (2, n_combinations, num_resi, num_bins)
            phi_combinations = np.transpose(np.array(tuple(itertools.combinations(phi_pdfs, 2))), 
                                            axes=[1, 0, 2, 3])
            psi_combinations = np.transpose(np.array(tuple(itertools.combinations(psi_pdfs, 2))),
                                            axes=[1, 0, 2, 3])

        if metric == "hellingers":
            # check if it's going to be a 1:1 comparison
            # note: i.e., the indexing changes in second variable if its a 1:1 comparison
            if phi_combinations.shape[0] == 1 and psi_combinations.shape[0] == 1:
                phi_metric = hellinger_distance(phi_combinations[0], phi_combinations[0])
                psi_metric = hellinger_distance(psi_combinations[0], psi_combinations[0])
            else:
                phi_metric = hellinger_distance(phi_combinations[0], phi_combinations[1])
                psi_metric = hellinger_distance(psi_combinations[0], psi_combinations[1])
        elif metric == "relative entropy":
            if phi_combinations.shape[0] == 1 and psi_combinations.shape[0] == 1:
                phi_metric = rel_entropy(phi_combinations[0], phi_combinations[0])
                psi_metric = rel_entropy(psi_combinations[0], psi_combinations[0])
            else:
                phi_metric = rel_entropy(phi_combinations[0], phi_combinations[1])
                psi_metric = rel_entropy(psi_combinations[0], psi_combinations[1])
        else:
            raise NotImplementedError(f"The metric: {metric} is not implemented.")

        return pd.DataFrame(phi_metric), pd.DataFrame(psi_metric)


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

    def quality_plot(self, 
                    dihedral : str = "phi", 
                    increment : int = 5,
                    figsize : Tuple[int, int] = (7, 5),
                    dpi : int = 400,
                    panel_labels : bool = False, 
                    fontsize : int = 10,
                    save_dir : str = None,
                    figname : str = "phi_hellingers.pdf"):
        """convenience plotting functionality for quick visual inspection of sampling quality

        Parameters
        ----------
        dihedral : str, optional
            the torsional angle to assess the quality of. Current options are phi and psi, by default "phi"
        increment : int, optional
             x axis stride, by default 5
        figsize : Tuple[int,int], optional
            the dimensions of the figure, by default (7,5)
        dpi : int, optional
            the dpi to save the figure, by default 400
        panel_labels : bool, optional
            whether or not to show manuscript figure panel labels, by default False
        fontsize : int, optional
            fontsize for all labels, by default 10
        save_dir : str, optional
            the directory to save the plotted graph in, by default None
        figname : str, optional
            name of the output figure, by default "phi_hellingers.pdf"

        Returns
        -------
        Tuple
            corresponding figure and axis of the subplot

        Raises
        ------
        NotImplementedError
            Raised when quality assessment is requested on a dihedral angle
            that is not currently supported
        """

        fig, axd = plt.subplot_mosaic(
            """AABB;CCDD""",
            sharex=True,
            figsize=figsize,
            dpi=dpi,
            facecolor="w",
            gridspec_kw={'height_ratios': [2, 2]}
        )
        if dihedral == "phi":
            metric = self.hellingers_distances[0]
        elif dihedral == "psi":
            metric = self.hellingers_distances[1]
        else:
            raise NotImplementedError

        phi_all_to_all, psi_all_to_all = self.get_all_to_all_trj_comparisons()

        trj_helicity, ref_helicity = self.fractional_helicity

        n_res = metric.shape[-1]
        idx = np.arange(1, n_res+1)
        xticks = np.arange(increment, idx[-1]+1, increment)
        xticklabels = np.arange(increment, idx[-1]+1, increment)

        for ax in axd:
            if ax == "A":
                axd[ax].set_yticks([0, 1], fontsize=fontsize)
                axd[ax].set_yticklabels([0, 1], fontsize=fontsize)
                axd[ax].set_ylim([0, 1])
                axd[ax].set_ylabel("Hellinger's Distance", fontsize=fontsize)
                axd[ax].set_title("Comparison to the Excluded Volume Limit", fontsize=fontsize)

                axd[ax].set_xticks(xticks, fontsize=fontsize)
                axd[ax].set_xticklabels(xticklabels, fontsize=fontsize)
                axd[ax].set_xlim([0, idx[-1]+1])

                # plot all red marks
                axd[ax].plot(idx, metric.transpose(), '.r', ms=4, alpha=0.3, mew=0)

                # plot mean
                axd[ax].plot(idx, np.mean(metric, axis=0), 'sk-', ms=2, alpha=1, mew=0, linewidth=0.5)

            elif ax == "B":
                axd[ax].set_yticks([0, 1], fontsize=fontsize)
                axd[ax].set_yticklabels([0, 1], fontsize=fontsize)
                axd[ax].set_ylim([0, 1])
                axd[ax].set_ylabel("Hellinger's Distance", fontsize=fontsize)

                axd[ax].set_title("All-to-All Trajectory Comparison", fontsize=fontsize)

                axd[ax].set_xticks(xticks, fontsize=fontsize)
                axd[ax].set_xticklabels(xticklabels, fontsize=fontsize)
                axd[ax].set_xlim([0, idx[-1]+1])

                axd[ax].plot(idx, phi_all_to_all.transpose(), '.r', ms=4, alpha=0.3, mew=0)

                # plot mean
                axd[ax].plot(idx, np.mean(phi_all_to_all, axis=0), 'sk-', ms=2, alpha=1, mew=0, linewidth=0.5)

            elif ax == "C":
                axd[ax].spines.right.set_visible(False)
                axd[ax].spines.top.set_visible(False)
                axd[ax].set_yticks([0, 1], fontsize=fontsize)
                axd[ax].set_yticklabels([0, 1], fontsize=fontsize)
                axd[ax].set_ylim([0, 1])
                axd[ax].set_ylabel("Hellinger's Distance\nmax - min", fontsize=fontsize)
                axd[ax].set_xlabel("Residue", fontsize=fontsize)

                max_minus_min = np.ptp(metric, axis=0)
                axd[ax].bar(idx, max_minus_min, width=0.8, color='k')

            elif ax == "D":
                axd[ax].set_yticks([0, 1], fontsize=fontsize)
                axd[ax].set_yticklabels([0, 1], fontsize=fontsize)
                axd[ax].set_ylim([0, 1])
                axd[ax].set_ylabel("Fractional Helicity", fontsize=fontsize)
                axd[ax].set_xlabel("Residue", fontsize=fontsize)

                axd[ax].set_xticks(xticks, fontsize=fontsize)
                axd[ax].set_xticklabels(xticklabels, fontsize=fontsize)
                axd[ax].set_xlim([0, idx[-1]+1])

                # plot red
                axd[ax].plot(idx, trj_helicity.transpose(), '.r', ms=4, alpha=0.3, mew=0)

                # plot line avg helicity
                axd[ax].plot(idx, np.mean(trj_helicity, axis=0), 'sk-', ms=2, alpha=1, mew=0, linewidth=0.5)

        if panel_labels:
            for ax in axd:
                trans = transforms.ScaledTranslation(-20/72, 7/72, fig.dpi_scale_trans)
                axd[ax].text(-0.0825, 1.10, ax, transform=axd[ax].transAxes+trans,
                             fontsize=fontsize, fontweight='bold', va='top', ha='right')
        plt.tight_layout()
        if save_dir is not None:
            os.makedirs(save_dir, exist_ok=True)
            outpath = os.path.join(save_dir, f"{figname}.pdf")
            fig.savefig(f"{outpath}", dpi=dpi)

        return fig, axd

    def plot_phi_psi_heatmap(self, 
                             metric: str = "hellingers",
                             figsize=(40, 20),
                             annotate=True,
                             cmap=None,
                             vmin=0.0,
                             vmax=1.0,
                             filename: str = "sampling_quality.png",
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
            filename, by default "sampling_quality.png"
        save_dir : _type_, optional
            directory path to save file, by default None

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

        if cmap is None:
            cmap = sns.color_palette("light:b", as_cmap=True)

        fig, (ax1, ax2, axcb) = plt.subplots(1, 3, figsize=figsize,
                                             gridspec_kw={'width_ratios': [1, 1, 0.05]}, **kwargs)

        ax1.get_shared_y_axes().join(ax2)
        g1 = sns.heatmap(phi_metric, annot=annotate, annot_kws={
                         "fontsize": 24, "color": "k"}, vmin=vmin, vmax=vmax, cmap=cmap, cbar=False, ax=ax1)
        g1.set_title(f'Phi {metric}', fontsize=36)
        g1.set_ylabel('Trajectory Index', fontsize=36)
        g1.set_xlabel('Residue Index', fontsize=36)
        g1.set_xticks(np.arange(0, phi_metric[:, :].shape[1])+0.5)
        g1.set_yticks(np.arange(0, phi_metric[:, :].shape[0])+0.5)

        g1.set_xticklabels(np.arange(1, phi_metric[:, :].shape[1]+1), fontsize=36)
        g1.set_yticklabels(np.arange(1, phi_metric[:, :].shape[0]+1), fontsize=36)

        g2 = sns.heatmap(psi_metric, annot=annotate, annot_kws={
                         "fontsize": 24, "color": "k"}, vmin=vmin, vmax=vmax, cmap=cmap, cbar_ax=axcb, ax=ax2)
        g2.set_title(f'Psi {metric}', fontsize=36)
        g2.set_ylabel('Trajectory Index', fontsize=36)
        g2.set_xlabel('Residue Index', fontsize=36)
        g2.set_xticks(np.arange(0, psi_metric[:, :].shape[1])+0.5)
        g2.set_yticks(np.arange(0, psi_metric[:, :].shape[0])+0.5)

        g2.set_xticklabels(np.arange(1, psi_metric[:, :].shape[1]+1), fontsize=36)
        g2.set_yticklabels(np.arange(1, psi_metric[:, :].shape[0]+1), fontsize=36)
        axcb.tick_params(labelsize=36)
        plt.tight_layout()
        if save_dir is not None:
            os.makedirs(save_dir, exist_ok=True)
            outpath = os.path.join(save_dir, filename)
            fig.savefig(f"{outpath}", dpi=300)


    @property
    def trj_pdfs(self):
        """property for getting the pdfs computed from the phi/psi angles respectively

        Returns
        -------
        np.ndarray
            pdfs computed from the phi and psi angles with the specified bins.
            returns (2, num_traj, n_res, n_bins)
        """
        bins = self.get_degree_bins()
        ref_phi_pdf = self.compute_pdf(self.phi_angles, bins=bins)
        ref_psi_pdf = self.compute_pdf(self.psi_angles, bins=bins)
        return np.array((ref_phi_pdf, ref_psi_pdf))

    @property
    def ref_pdfs(self):
        """property for getting the pdfs computed from the phi/psi angles respectively

        Returns
        -------
        np.ndarray
            pdfs computed from the phi and psi angles with the specified bins.
            returns (2, num_traj, n_res, n_bins) array
        """
        bins = self.get_degree_bins()
        trj_phi_pdf = self.compute_pdf(self.ref_phi_angles, bins=bins)
        trj_psi_pdf = self.compute_pdf(self.ref_psi_angles, bins=bins)
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
            The per residue fractional helicity for each trajectory in self.trajs and self.ref_trajs.
        """
        return self.__compute_frac_helicity()
