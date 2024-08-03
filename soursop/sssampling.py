#  _____  ____  _    _ _____   _____  ____  _____
# / ____|/ __ \| |  | |  __ \ / ____|/ __ \|  __ \
# | (___ | |  | | |  | | |__) | (___ | |  | | |__)|
# \___ \| |  | | |  | |  _  / \___ \| |  | |  ___/
# ____) | |__| | |__| | | \ \ ____) | |__| | |
# |_____/ \____/ \____/|_|  \_\_____/ \____/|_|
# Jeffrey M. Lotthammer (Holehouse Lab)
# Alex Holehouse (Pappu Lab and Holehouse Lab) and Jared Lalmansing (Pappu lab)
# Simulation analysis package
# Copyright 2014 - 2022
##
import itertools
import os
import pathlib
from typing import List, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import transforms
from natsort import natsorted
from scipy.special import rel_entr

from soursop import ssutils
from soursop.ssdata import (
    EV_RESIDUE_MAPPER,
    ONE_TO_THREE,
    PHI_EV_ANGLES_DICT,
    PSI_EV_ANGLES_DICT,
)

from .ssexceptions import SSException
from .sstrajectory import SSTrajectory, parallel_load_trjs


def compute_joint_hellinger_distance(p, q):
    # Compute the Bhattacharyya coefficient
    b_coefficient = np.sum(np.sqrt(p * q))

    # Compute the Hellinger's distance - note this doesn't need the normalization by sqrt(2)
    distance = np.sqrt(1 - b_coefficient)

    return distance

def hellinger_distance(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    """
    Computes the Hellinger distance between a set of probability distributions p and q.
    The Hellinger distances is defined as:
    
        :math:`H(P,Q) = \\frac{1}{\\sqrt{2}} \\times \\sqrt{\\sum_{i=1}^{k}(\\sqrt{p_i}-\\sqrt{q_i})^2}`
    
    where k is the length of the probability vectors being compared. 

    For sets of distributions, the datapoints should be in the last axis.
    
    Parameters
    ----------
    p : np.ndarray
        A probability distribution or set of probability distributions to compare.
    q : np.ndarray
        A probability distribution or set of probability distributions to compare.

    Returns
    -------
    np.ndarray
        The Hellinger distance(s) between p and q.
    """
    # Ensure that p and q are NumPy arrays
    p = np.asarray(p)
    q = np.asarray(q)

    # Compute the Hellinger distance
    numerator = np.sum(np.square(np.sqrt(p) - np.sqrt(q)), axis=-1)
    denominator = np.sqrt(2)
    return np.sqrt(numerator) / denominator


def rel_entropy(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Computes the relative entropy between two probability distributions p and q.
    For sets of distributions, the datapoints should be in the last axis.

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
    p = np.asarray(p)
    q = np.asarray(q)

    relative_entropy = np.sum(rel_entr(p, q), axis=-1)

    return relative_entropy


class SamplingQuality:
    def __init__(
        self,
        traj_list: List[str],
        reference_list: Union[List[str], None] = None,
        top_file: str = "__START.pdb",
        ref_top: Union[str, None] = None,
        method: str = "2D angle distributions",
        bwidth: float = np.deg2rad(15),
        proteinID: int = 0,
        n_cpus: int = None,
        truncate: bool = False,
        **kwargs: dict,
    ):
        """
        SamplingQuality is a class to compare the sampling quality for a set of trajectories
        relative to some referene model. This is usually a polymer limiting model, but could
        be different proteins, mutations, or PTMs.

        Parameters
        ----------
        traj_list : List[str]
            a list of the trajectory paths associated with the simulated trajectories.
        reference_list : Union[List[str], None]
            a list of the trajectory paths associated with the reference model.
            This is usually a limiting polymer model.
        top_file : str
            path to the simulated trajectories topology file.
        ref_top : Union[str, None]
            path to the reference topology - usually the same topology file
            if using a polymer model, but can differ if reference is e.g., a mutant.
        method : str
            The method used to compute the Hellinger distance
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
        self.kwargs = kwargs

        self.bins = self.get_degree_bins()
        self.__precomputed = {}

        self.__validate_arguments()

        self.__load_trajectories()

        # if reference trajectories have been provided
        # then self.ref_trajs should have been initialized.
        if self.reference_list:
            # if truncate is True,
            # then match the lengths of the trajectories before computing dihedrals
            if self.truncate:
                self.trajs, self.ref_trajs = self.__truncate_trajectories()

            # compute all dihedrals from trajectories and ref trajectories
            (
                self.psi_angles,
                self.ref_psi_angles,
                self.phi_angles,
                self.ref_phi_angles,
            ) = self.__compute_dihedrals(proteinID=self.proteinID)

        # if no reference trajectories have been provided
        else:
            if self.truncate:
                self.trajs, self.ref_trajs = self.__truncate_trajectories()

            (self.psi_angles, self.phi_angles) = self.__compute_dihedrals(
                proteinID=self.proteinID, precomputed=True
            )

            # if no reference list is provided, use precomputed reference dihedrals
            # for the limiting polymer model.

            ## NOTE this assumes that all trajectories will be the same sequence - this is implicit from the topology
            # anyway, so this is fine but just making it explicit.
            sequence = (
                self.trajs[0]
                .proteinTrajectoryList[self.proteinID]
                .get_amino_acid_sequence(oneletter=True)
            )

            # remove caps from sequence if present
            sequence = sequence.replace(">", "").replace("<", "")

            precomputed_interface = PrecomputedDihedralInterface(
                sequence,
                bins=self.bins,
                num_trajs=len(self.trajs),
                nsamples=len(self.trajs[0]),
            )

            self.ref_psi_angles = precomputed_interface.ref_psi_angles
            self.ref_phi_angles = precomputed_interface.ref_phi_angles

    def __validate_arguments(self):
        ssutils.validate_keyword_option(
            self.method, ["2D angle distributions", "1D angle distributions"], "method"
        )

        if self.bwidth > 2 * np.pi or not self.bwidth > 0:
            raise SSException(
                f"The bwidth parameter must be between 0 and 2*pi.\
                    Received {self.bwidth}"
            )

        if not self.n_cpus:
            self.n_cpus = os.cpu_count()

        if len(self.traj_list) == 0:
            raise SSException(
                f"Input trajectory list must be non-empty.\
                    Received len(traj_list)={len(self.traj_list)}"
            )

    def __load_trajectories(self):
        # weird thing I have to do to prevent issues with multiprocessing
        # parallel loading when there is only 1 trajectory to load
        # trajs/ref_trajs must be a list so they're iterables for __truncate_trajectories

        if len(self.traj_list) == 1:
            self.trajs = []
            self.trajs.append(
                SSTrajectory(self.traj_list, pdb_filename=self.top, **self.kwargs)
            )

            # if the reference list has been provided initialize the reference trajectories
            # else the reference dihedrals will be assigned from precomputed dihedrals later.
            if not self.reference_list:
                pass

            elif len(self.reference_list) == 1:
                self.ref_trajs = []
                self.ref_trajs.append(
                    SSTrajectory(
                        self.reference_list, pdb_filename=self.ref_top, **self.kwargs
                    )
                )

        else:
            # if many trajectories, load in parallel
            self.trajs = parallel_load_trjs(
                self.traj_list, self.top, n_procs=self.n_cpus, **self.kwargs
            )

            # if the reference list has been provided initialize the reference trajectories
            # else the reference dihedrals will be assigned from precomputed dihedrals later.
            if self.reference_list:
                self.ref_trajs = parallel_load_trjs(
                    self.reference_list,
                    self.ref_top,
                    n_procs=self.n_cpus,
                    **self.kwargs,
                )

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
        lengths = []
        # TODO: Make this work with Precomputed dihedrals
        if not self.reference_list:
            for trj in self.trajs:
                lengths.append(trj.n_frames)
            self.min_length = np.min(lengths)

            temp_trajs = []
            for trj in self.trajs:
                temp_trajs.append(
                    SSTrajectory(
                        TRJ=trj.proteinTrajectoryList[self.proteinID].traj[
                            0 : self.min_length
                        ]
                    )
                )
            print(
                f"Successfully truncated.\n\
                    The shortest trajectory is: {self.min_length} frames.\
                    All trajectories truncated to {self.min_length}"
            )
            return (temp_trajs, None)

        for trj, ref_trj in zip(self.trajs, self.ref_trajs):
            lengths.append([trj.n_frames, ref_trj.n_frames])

        # shift frames for np.array indexing purposes
        self.min_length = np.min(lengths)

        temp_trajs = []
        temp_ref_trjs = []
        for trj, ref_trj in zip(self.trajs, self.ref_trajs):
            temp_trajs.append(
                SSTrajectory(
                    TRJ=trj.proteinTrajectoryList[self.proteinID].traj[
                        0 : self.min_length
                    ]
                )
            )
            temp_ref_trjs.append(
                SSTrajectory(
                    TRJ=ref_trj.proteinTrajectoryList[self.proteinID].traj[
                        0 : self.min_length
                    ]
                )
            )

        print(
            f"Successfully truncated.\n\
                The shortest trajectory is: {self.min_length} frames.\
                All trajectories truncated to {self.min_length}"
        )

        return (temp_trajs, temp_ref_trjs)

    def __compute_dihedrals(
        self, proteinID: int = 0, precomputed: bool = False
    ) -> np.ndarray:
        """internal function to computes the phi/psi backbone dihedrals
        at a given index proteinID in the ``SSTrajectory.proteinTrajectoryList`` of an SSTrajectory.

        Parameters
        ----------
        proteinID : int, optional
            The ID of the protein where the ID is the proteins position
            in the ``SSTrajectory.proteinTrajectoryList`` list, by default 0.

        Returns
        -------
        np.ndarray
            Returns the psi and phi backbone dihedrals for the simulated trajectory and the limiting polyer model.
        """
        psi_angles = []
        phi_angles = []
        ref_psi_angles = []
        ref_phi_angles = []

        # if we're not using precomputed dihedrals, compute from the reference trajs
        if not precomputed:
            for trj, ref_trj in zip(self.trajs, self.ref_trajs):
                psi_angles.append(
                    trj.proteinTrajectoryList[proteinID].get_angles("psi")[1]
                )
                phi_angles.append(
                    trj.proteinTrajectoryList[proteinID].get_angles("phi")[1]
                )
                ref_psi_angles.append(
                    ref_trj.proteinTrajectoryList[proteinID].get_angles("psi")[1]
                )
                ref_phi_angles.append(
                    ref_trj.proteinTrajectoryList[proteinID].get_angles("phi")[1]
                )

            # return the angles for everything
            return np.array((psi_angles, ref_psi_angles, phi_angles, ref_phi_angles))

        # else only compute dihedrals from the simulated trajectories
        else:
            for trj in self.trajs:
                psi_angles.append(
                    trj.proteinTrajectoryList[proteinID].get_angles("psi")[1]
                )
                phi_angles.append(
                    trj.proteinTrajectoryList[proteinID].get_angles("phi")[1]
                )

            # return the angles for simulated trajectories only
            return np.array((psi_angles, phi_angles))

    def compute_frac_helicity(
        self, proteinID: int = 0, recompute: bool = False
    ) -> np.ndarray:
        """Function that computes the per residue fractional helicity at a given index (proteinID)
        in the proteinTrajectoryList of an SSTrajectory for all SSTrajectory objects provided
        in the ``self.trajs`` list.

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
        selectors = ("trj_helicity", "ref_helicity")
        if not recompute and all(
            selector in self.__precomputed for selector in selectors
        ):
            return self.__precomputed["trj_helicity"], self.__precomputed[
                "ref_helicity"
            ]

        trj_helicity = [
            trj.proteinTrajectoryList[proteinID].get_secondary_structure_DSSP()[1]
            for trj in self.trajs
        ]

        self.__precomputed["trj_helicity"] = np.array(trj_helicity)

        if self.reference_list:
            reference_helicity = [
                ref_trj.proteinTrajectoryList[proteinID].get_secondary_structure_DSSP()[
                    1
                ]
                for ref_trj in self.ref_trajs
            ]
        else:
            reference_helicity = np.zeros_like(self.__precomputed["trj_helicity"])

        self.__precomputed["ref_helicity"] = np.array(reference_helicity)

        return self.__precomputed["trj_helicity"], self.__precomputed["ref_helicity"]

    def compute_dihedral_hellingers(self) -> np.ndarray:
        """Compute the Hellinger distance for both the phi and psi angles between a set of trajectories.

        Returns
        -------
        np.ndarray
            The Hellinger distances between the probability density distributions for the phi and psi angles for a set of trajectories.
        """
        if self.method == "2D angle distributions":
            data = np.array([self.phi_angles, self.psi_angles])
            ref_data = np.array([self.ref_phi_angles, self.ref_psi_angles])

            pdfs = self.compute_series_of_histograms_along_axis(
                data, bins=self.bins, axis=2
            )
            ref_pdfs = self.compute_series_of_histograms_along_axis(
                ref_data, bins=self.bins, axis=2
            )

            joint_hellingers = self.__compute_2d_dihedral_hellingers(pdfs, ref_pdfs)

            return np.array(joint_hellingers)

        elif self.method == "1D angle distributions":
            phi_trj_pdfs = self.compute_pdf(self.phi_angles, bins=self.bins)
            phi_ref_trj_pdfs = self.compute_pdf(self.ref_phi_angles, bins=self.bins)

            psi_trj_pdfs = self.compute_pdf(self.psi_angles, bins=self.bins)
            psi_ref_trj_pdfs = self.compute_pdf(self.ref_psi_angles, bins=self.bins)

            phi_hellingers = hellinger_distance(phi_trj_pdfs, phi_ref_trj_pdfs)
            psi_hellingers = hellinger_distance(psi_trj_pdfs, psi_ref_trj_pdfs)

            return np.array((phi_hellingers, psi_hellingers))

        else:
            raise NotImplementedError(
                f"{self.method} is not defined!\
                                      Please use either 1D angle distributions\
                                      or 2D angle distributions"
            )

    def __compute_2d_dihedral_hellingers(self, trj_pdfs, ref_pdfs):
        """
        Helter function to Compute the Hellinger distances for
        2D dihedral angle probability density functions (PDFs).

        Parameters
        ----------
        trj_pdfs : ndarray
            Array of PDFs representing dihedral angle distributions for trajectory replicates.
        ref_pdfs : ndarray
            Array of PDFs representing reference dihedral angle distributions.

        Returns
        -------
        ndarray
            Array of Hellinger distances for each trajectory replica and dihedral angle.

        Notes
        -----
        - The input arrays trj_pdfs and ref_pdfs should have the same shape.
        - Each array has dimensions (num_replicates, num_angles, num_bins_phi, num_bins_psi),
        where num_replicates is the number of trajectory replicates,
        num_angles is the number of dihedral angles, and
        num_bins_phi and num_bins_psi are the number of bins in the phi and psi dimensions, respectively.
        - The function computes the Hellinger distances between the corresponding PDFs of each replicate and angle.
        - The Hellinger distance measures the similarity between two probability distributions.
        - 0 is returned if the two distributions are identical, and 1 is returned if the two distributions are completely different.
        - The computed distances are returned as an ndarray of shape (num_replicates, num_angles).
        """

        # Get the number of trajectory replicates
        num_replicates = trj_pdfs.shape[0]

        # Compute Hellinger's distances for each replicate
        hellinger_distances = []
        for replicate_idx in range(num_replicates):
            pdf1 = trj_pdfs[replicate_idx]
            pdf2 = ref_pdfs[replicate_idx]

            replicate_distances = []
            for angle_idx in range(pdf1.shape[0]):
                pdf1_angle = pdf1[angle_idx]
                pdf2_angle = pdf2[angle_idx]

                distance = compute_joint_hellinger_distance(pdf1_angle, pdf2_angle)
                replicate_distances.append(distance)

            hellinger_distances.append(replicate_distances)

        return hellinger_distances

    def compute_dihedral_rel_entropy(self) -> np.ndarray:
        """Compute the relative entropy for both the phi and psi angles between a set of trajectories.

        Returns
        -------
        np.ndarray
            The relative entropy between the probability density distributions for the phi and psi angles for a set of trajectories.
        """

        phi_trj_pdfs = self.compute_pdf(self.phi_angles, bins=self.bins)
        phi_ref_trj_pdfs = self.compute_pdf(self.ref_phi_angles, bins=self.bins)

        psi_trj_pdfs = self.compute_pdf(self.psi_angles, bins=self.bins)
        psi_ref_trj_pdfs = self.compute_pdf(self.ref_psi_angles, bins=self.bins)

        phi_rel_entr = rel_entropy(phi_trj_pdfs, phi_ref_trj_pdfs)
        psi_rel_entr = rel_entropy(psi_trj_pdfs, psi_ref_trj_pdfs)

        return np.array((phi_rel_entr, psi_rel_entr))

    def compute_series_of_histograms_along_axis(
        self, data: np.ndarray, bins: np.ndarray, axis: int = 0):
        """
        Compute a series of 2D histograms along an axis of a 4D array
        and convert them into probability density functions (PDFs).

        Parameters
        ----------
        data : ndarray
            4D array containing the joint phi/psi angle data.
    
        bins : ndarray
            1D array defining the bin edges for the histograms.
    
        axis : int, optional
            Axis along which to compute the histograms (default=0).

        Returns
        -------
        pdfs : ndarray
            Series of 2D probability density functions (PDFs) along the given axis.
        
        Notes
        -----
        - The input data array should have dimensions (2, n_trajs, n_residues, n_samples).
        - The bins array should contain the bin edges for the histograms.
        - The resulting PDFs will have dimensions (n_trajs, n_residues, num_bins_phi, num_bins_psi), where num_bins_phi and num_bins_psi are the number of bins in the phi and psi dimensions, respectively.        
        - The PDFs are computed by normalizing the histogram values and multiplying them by the corresponding bin widths.
        
        """
        # Get the shape of the input array
        shape = data.shape

        # Initialize an empty list to store the PDFs for each trajectory
        pdfs = []

        # Loop over the trajectories
        for traj_idx in range(shape[1]):
            traj_histograms = []

            # Loop over the residue indices
            for residue_idx in range(shape[2]):
                # Get the joint phi/psi angles for the current trajectory and residue
                angles = data[:, traj_idx, residue_idx, :]

                # Compute the 2D histogram for the joint phi/psi angles
                hist, x_edges, y_edges = np.histogram2d(
                    angles[0], angles[1], bins=bins, density=True
                )

                # Compute the bin widths along each dimension
                bin_width_phi = x_edges[1] - x_edges[0]
                bin_width_psi = y_edges[1] - y_edges[0]

                # Multiply the histogram values by the bin widths to obtain the PDF
                pdf = hist * (bin_width_phi * bin_width_psi)

                traj_histograms.append(pdf)

            pdfs.append(traj_histograms)

        return np.array(pdfs)

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
            pdf = np.apply_along_axis(
                lambda col: np.histogram(col, bins=bins, density=True)[0],
                axis=2,
                arr=arr,
            ) * np.round(np.rad2deg(self.bwidth))

            # KEY POINT: multiplying by bin width to convert probability *density* to probabilty *mass*
            # implementation details may have to change here if supporting other methods.
            # pdf = histogram(arr, bins=bins, axis=2, density=True)[0]*np.round(np.rad2deg(self.bwidth))
        # else (n_res x n_frames), histogram axis (1) associated with frames
        else:
            pdf = np.apply_along_axis(
                lambda col: np.histogram(col, bins=bins, density=True)[0],
                axis=1,
                arr=arr,
            ) * np.round(np.rad2deg(self.bwidth))
            # pdf = histogram(arr, bins=bins, axis=1, density=True)[0]*np.round(np.rad2deg(self.bwidth))

        return pdf

    def get_all_to_all_2d_trj_comparison(
        self, metric: str = "hellingers", recompute=False
    ) -> Tuple[pd.DataFrame]:
        # if self.method == "2D angle distributions":
        data = np.array([self.phi_angles, self.psi_angles])

        # shape = replicas, angles, phi_bins, psi_bins
        pdfs = self.compute_series_of_histograms_along_axis(
            data, bins=self.bins, axis=2
        )

        if pdfs.shape[0] == 1:
            # if only 1 simulated traj, an all-to-all is just a self:self comparison.
            # after transpose: [combinations, replicates, angles, phi_bins, psi_bins]
            pdf_combinations = np.transpose(
                np.array(tuple(itertools.combinations(pdfs, 1))), axes=[1, 0, 2, 3, 4]
            )
        else:
            # original shape is: [n_combinations, 2, angle, phi_bins, psi_bins]
            # 2 because it's a pairwise head-to-head comparison of trajectories.
            # transposed for my sanity for indexing leaving final shape as:
            # (2, n_combinations, num_resi, phi_bins, psi_bins)
            pdf_combinations = np.transpose(
                np.array(tuple(itertools.combinations(pdfs, 2))), axes=[1, 0, 2, 3, 4]
            )

        if metric == "hellingers":
            # check if it's going to be a 1:1 comparison
            # note: i.e., the indexing changes in second variable if its a 1:1 comparison
            if pdf_combinations.shape[0] == 1:
                dist_metric = []

                for replicate in range(pdf_combinations[0].shape[0]):
                    all_residue_replicate_distances = []
                    for angle in range(pdf_combinations[0][replicate].shape[0]):
                        # note the same index (0) for both pdfs because it's a self:self comparison
                        curr_residue_distance = compute_joint_hellinger_distance(
                            pdf_combinations[0][replicate][angle],
                            pdf_combinations[0][replicate][angle],
                        )

                        all_residue_replicate_distances.append(curr_residue_distance)

                    dist_metric.append(all_residue_replicate_distances)

                dist_metric = np.array(dist_metric)

            else:
                dist_metric = []

                for replicate in range(pdf_combinations[0].shape[0]):
                    all_residue_replicate_distances = []
                    for angle in range(pdf_combinations[0][replicate].shape[0]):
                        # note the different index (1) for both pdfs because it's a pairwise comparison
                        curr_residue_distance = compute_joint_hellinger_distance(
                            pdf_combinations[0][replicate][angle],
                            pdf_combinations[1][replicate][angle],
                        )

                        all_residue_replicate_distances.append(curr_residue_distance)

                    dist_metric.append(all_residue_replicate_distances)

        return np.array(dist_metric)

    def get_all_to_all_trj_comparisons(
        self, metric: str = "hellingers", recompute=False
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """function to aggregate an all-to-all comparison of pdfs

        Parameters
        ----------
        metric : str, optional
            The metric to use for the comparison, by default "hellingers"

        Returns
        -------
        Tuple[pd.DataFrame, pd.DataFrame]
            all-to-all trajectory comparisons for the Hellinger distances in 'self.trj_pdfs'
        """
        phi_pdfs = self.trj_pdfs(recompute=recompute, dihedral="trj_phi_pdfs")
        psi_pdfs = self.trj_pdfs(recompute=recompute, dihedral="trj_psi_pdfs")

        if phi_pdfs.shape[0] == 1 or psi_pdfs.shape[0] == 1:
            # if only 1 simulated traj and 1 ref traj all-to-all is just a 1:1 comparison.
            phi_combinations = np.transpose(
                np.array(tuple(itertools.combinations(phi_pdfs, 1))), axes=[1, 0, 2, 3]
            )
            psi_combinations = np.transpose(
                np.array(tuple(itertools.combinations(psi_pdfs, 1))), axes=[1, 0, 2, 3]
            )
        else:
            # returned array is (n_combinations, 2, num_resi, num_bins)
            # 2 because it's a pairwise head-to-head comparison of trajectories.
            # transposed for my sanity for indexing leaving final shape as:
            # (2, n_combinations, num_resi, num_bins)
            phi_combinations = np.transpose(
                np.array(tuple(itertools.combinations(phi_pdfs, 2))), axes=[1, 0, 2, 3]
            )
            psi_combinations = np.transpose(
                np.array(tuple(itertools.combinations(psi_pdfs, 2))), axes=[1, 0, 2, 3]
            )

        if metric == "hellingers":
            # check if it's going to be a 1:1 comparison
            # note: i.e., the indexing changes in second variable if its a 1:1 comparison
            if phi_combinations.shape[0] == 1 and psi_combinations.shape[0] == 1:
                phi_metric = hellinger_distance(
                    phi_combinations[0], phi_combinations[0]
                )
                psi_metric = hellinger_distance(
                    psi_combinations[0], psi_combinations[0]
                )
            else:
                phi_metric = hellinger_distance(
                    phi_combinations[0], phi_combinations[1]
                )
                psi_metric = hellinger_distance(
                    psi_combinations[0], psi_combinations[1]
                )
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
        bins = np.arange(-180, 180 + bwidth, bwidth)
        return bins

    def quality_plot(
        self,
        increment: int = 5,
        figsize: Tuple[int, int] = (7, 5),
        dpi: int = 400,
        panel_labels: bool = False,
        fontsize: int = 10,
        save_dir: str = None,
        dihedral: Union[None, str] = "2D",
        figname: str = "hellingers.pdf",
    ):
        """
        This function enables a convenient plotting functionality for quick visual 
        inspection of the sampling quality.

        Parameters
        ----------
        dihedral : str, optional
            the torsional angle to assess the quality of. Current options are 2D, phi, and psi, by default "2D"
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
            gridspec_kw={"height_ratios": [2, 2]},
        )
        if self.method == "1D angle distributions" and dihedral == "2D":
            raise ValueError(
                f"Cannot plot 1D angle distributions with dihedral = {dihedral} selector.\
                             Please set dihedral to phi or psi"
            )

        selector = {
            "2D": self.compute_dihedral_hellingers(),
            "phi": self.compute_dihedral_hellingers()[0],
            "psi": self.compute_dihedral_hellingers()[1],
        }

        all_to_all_selector = {
            "2D": self.get_all_to_all_2d_trj_comparison(),
            "phi": self.get_all_to_all_trj_comparisons()[0],
            "psi": self.get_all_to_all_trj_comparisons()[1],
        }

        metric = selector[dihedral]
        print("metric shape: ", metric.shape)
        all_to_all = all_to_all_selector[dihedral]

        trj_helicity, ref_helicity = self.fractional_helicity()

        # if self.method == "2D angle distributions" and dihedral == "2D":
        #     metric = selector["2D"]
        #     joint_all_to_all = self.get_all_to_all_2d_trj_comparison()
        # elif self.method == "1D angle distributions" and dihedral == "phi":
        #     metric = selector["phi"]
        #     phi_all_to_all, psi_all_to_all = self.get_all_to_all_trj_comparisons()
        # elif self.method == "1D angle distributions" and dihedral == "psi":
        #     metric = selector["psi"]
        #     phi_all_to_all, psi_all_to_all = self.get_all_to_all_trj_comparisons()
        # else:
        #     raise NotImplementedError(f"{self.method} cannot be used with {dihedral}." +
        #                               f"Currently supported options are:\
        #                               1D angle distributions and phi/psi or 2D angle distributions and 2D")

        n_res = metric.shape[-1]
        idx = np.arange(1, n_res + 1)
        xticks = np.arange(increment, idx[-1] + 1, increment)
        xticklabels = np.arange(increment, idx[-1] + 1, increment)

        yticks = [0, 0.2, 0.4, 0.6, 0.8, 1]
        ytick_labels = [0, 0.2, 0.4, 0.6, 0.8, 1]
        for ax in axd:
            if ax == "A":
                axd[ax].set_yticks(yticks)
                axd[ax].set_yticklabels(ytick_labels, fontsize=fontsize)
                axd[ax].set_ylim([0, 1])
                axd[ax].set_ylabel("Hellinger's Distance", fontsize=fontsize)
                axd[ax].set_title(
                    "Comparison to the Excluded Volume Limit", fontsize=fontsize
                )

                axd[ax].set_xticks(
                    xticks,
                )
                axd[ax].set_xticklabels(xticklabels, fontsize=fontsize)
                axd[ax].set_xlim([0, idx[-1] + 1])

                # plot all red marks
                axd[ax].plot(idx, metric.transpose(), ".r", ms=4, alpha=0.3, mew=0)

                # plot mean
                axd[ax].plot(
                    idx,
                    np.mean(metric, axis=0),
                    "sk-",
                    ms=2,
                    alpha=1,
                    mew=0,
                    linewidth=0.5,
                )

            elif ax == "B":
                axd[ax].set_yticks(yticks)
                axd[ax].set_yticklabels(ytick_labels, fontsize=fontsize)
                axd[ax].set_ylim([0, 1])
                axd[ax].set_ylabel("Hellinger's Distance", fontsize=fontsize)

                axd[ax].set_title("All-to-All Trajectory Comparison", fontsize=fontsize)

                axd[ax].set_xticks(xticks)
                axd[ax].set_xticklabels(xticklabels, fontsize=fontsize)
                axd[ax].set_xlim([0, idx[-1] + 1])

                axd[ax].plot(idx, all_to_all.transpose(), ".r", ms=4, alpha=0.3, mew=0)

                # plot mean
                axd[ax].plot(
                    idx,
                    np.mean(all_to_all, axis=0),
                    "sk-",
                    ms=2,
                    alpha=1,
                    mew=0,
                    linewidth=0.5,
                )

            elif ax == "C":
                # axd[ax].spines.right.set_visible(False)
                # axd[ax].spines.top.set_visible(False)
                axd[ax].set_yticks(yticks)
                axd[ax].set_yticklabels(ytick_labels, fontsize=fontsize)
                axd[ax].set_ylim([0, 1])
                axd[ax].set_ylabel("Hellinger's Distance\nmax - min", fontsize=fontsize)
                axd[ax].set_xlabel("Residue", fontsize=fontsize)

                max_minus_min = np.ptp(metric, axis=0)
                axd[ax].bar(idx, max_minus_min, width=0.8, color="k")

            elif ax == "D":
                axd[ax].set_yticks(yticks)
                axd[ax].set_yticklabels(ytick_labels, fontsize=fontsize)
                axd[ax].set_ylim([0, 1])
                axd[ax].set_ylabel("Fractional Helicity", fontsize=fontsize)
                axd[ax].set_xlabel("Residue", fontsize=fontsize)

                axd[ax].set_xticks(
                    xticks,
                )
                axd[ax].set_xticklabels(xticklabels, fontsize=fontsize)
                axd[ax].set_xlim([0, idx[-1] + 1])

                # plot red
                axd[ax].plot(
                    idx, trj_helicity.transpose(), ".r", ms=4, alpha=0.3, mew=0
                )

                # plot line avg helicity
                axd[ax].plot(
                    idx,
                    np.mean(trj_helicity, axis=0),
                    "sk-",
                    ms=2,
                    alpha=1,
                    mew=0,
                    linewidth=0.5,
                )

        if panel_labels:
            for ax in axd:
                trans = transforms.ScaledTranslation(
                    -20 / 72, 7 / 72, fig.dpi_scale_trans
                )
                axd[ax].text(
                    -0.0825,
                    1.10,
                    ax,
                    transform=axd[ax].transAxes + trans,
                    fontsize=fontsize,
                    fontweight="bold",
                    va="top",
                    ha="right",
                )
        plt.tight_layout()
        if save_dir is not None:
            os.makedirs(save_dir, exist_ok=True)
            outpath = os.path.join(save_dir, f"{dihedral}_{figname}.pdf")
            fig.savefig(f"{outpath}", dpi=dpi)

        return fig, axd

    def trj_pdfs(self, dihedral: str = "joint", recompute: bool = False):
        """Function to return the pdfs computed from the phi/psi angles, respectively

        Parameters
        ----------
        dihedral : str, optional
              method to use to return specific PDFs. Options are: "trj_phi_pdfs", "trj_psi_pdfs",
              "joint". By default this is "joint".
        
        recompute : bool, optional
            Whether or not to recompute the PDFs, by default False.

        Returns
        -------
        np.ndarray
            PDFs computed from the phi and psi angles with the specified bins.
            returns (2, num_traj, n_res, n_bins)

        Raises
        ------
        NotImplementedError
            Raised if the selector is not one of the implemented options.
        """
        selectors = ["trj_phi_pdfs", "trj_psi_pdfs", "joint"]
        if dihedral not in selectors:
            raise NotImplementedError(
                f"Should not arrive here: {selectors} is not implemented."
                + "Please try one of trj_phi_pdfs, trj_psi_pdfs, joint instead."
            )

        for selector in selectors:
            if selector not in self.__precomputed or recompute is True:
                if selector == "trj_phi_pdfs":
                    self.__precomputed[selector] = self.compute_pdf(
                        self.phi_angles, bins=self.bins
                    )
                elif selector == "trj_psi_pdfs":
                    self.__precomputed[selector] = self.compute_pdf(
                        self.psi_angles, bins=self.bins
                    )
                elif selector == "joint":
                    data = np.array([self.phi_angles, self.psi_angles])
                    pdfs = self.compute_series_of_histograms_along_axis(
                        data, bins=self.bins, axis=2
                    )
                    self.__precomputed[selector] = pdfs

        return self.__precomputed[dihedral]


    def ref_pdfs(self, dihedral="joint", recompute=False):
        """Function to return the pdfs computed from the phi/psi angles, respectively

        Parameters
        ----------
        dihedral : str, optional
              method to use to return specific PDFs. Options are: "trj_phi_pdfs", "trj_psi_pdfs",
              and "joint". By default this is "joint".
        
        recompute : bool, optional
            Whether or not to recompute the PDFs, by default False.

        Returns
        -------
        np.ndarray
            PDFs computed from the phi and psi angles with the specified bins.
            returns (2, num_traj, n_res, n_bins)

        Raises
        ------
        NotImplementedError
            Raised if the selector is not one of the implemented options.
        """
        selectors = ["ref_phi_pdfs", "ref_psi_pdfs", "joint"]

        if dihedral not in selectors:
            raise NotImplementedError(
                f"Should not arrive here: {dihedral} is not implemented."
                + "Please try one of ref_phi_pdfs, ref_psi_pdfs, joint instead."
            )

        for selector in selectors:
            if selector not in self.__precomputed or recompute is True:
                if selector == "ref_phi_pdfs":
                    self.__precomputed[selector] = self.compute_pdf(
                        self.ref_phi_angles, bins=self.bins
                    )
                elif selector == "ref_psi_pdfs":
                    self.__precomputed[selector] = self.compute_pdf(
                        self.ref_psi_angles, bins=self.bins
                    )
                elif selector == "joint":
                    data = np.array([self.ref_phi_angles, self.ref_psi_angles])
                    pdfs = self.compute_series_of_histograms_along_axis(
                        data, bins=self.bins, axis=2
                    )
                    self.__precomputed[selector] = pdfs

        return self.__precomputed[dihedral]


    def hellingers_distances(self, recompute=False):
        """property for getting the Hellinger distances computed from the phi/psi angles, respectively

        Returns
        -------
        np.ndarray
            Hellinger distances computed from the phi and psi angles with the specified bins.
            2 x n_reps x n_angles
            where 0 is phi and 1 is psi
        """
        selector = "hellingers"

        if selector not in self.__precomputed or recompute is True:
            self.__precomputed[selector] = self.compute_dihedral_hellingers()

        return self.__precomputed[selector]

    def fractional_helicity(self, recompute=False):
        """
        Property for getting the per residue fractional helicity for all trajectories.

        Returns:
        ----------
        np.ndarray
            The per residue fractional helicity for each trajectory in self.trajs and self.ref_trajs.
        
        """
        selectors = ("trj_helicity", "ref_helicity")
        if not recompute and all(
            selector in self.__precomputed for selector in selectors
        ):
            return self.__precomputed["trj_helicity"], self.__precomputed[
                "ref_helicity"
            ]

        trj_helicity, ref_helicity = self.compute_frac_helicity()

        return trj_helicity, ref_helicity


# Interface to separate computation of dihedrals from SamplingQuality class
# will serve to return precomputed excluded volume dihedral angle distributions
# if no EV trajectories are provided.


class PrecomputedDihedralInterface:
    """docstring for PrecomputedDihedralInterface."""

    def __init__(self, sequence, bins, num_trajs, nsamples):
        self.sequence = sequence
        self.num_trajs = num_trajs
        self.nsamples = nsamples
        self.bins = bins

        self.tmp_phi_angles = self.gather_phi_reference_dihedrals(self.sequence)
        self.tmp_psi_angles = self.gather_psi_reference_dihedrals(self.sequence)

        # ensure len ref angles is equal to number of angles found in traj arrays.

        # test case used to ensure we match exactly when not sampling
        # self.ref_phi_angles = np.tile(self.gather_phi_reference_dihedrals(sequence), (self.num_trajs, 1, 1))
        # self.ref_psi_angles = np.tile(self.gather_psi_reference_dihedrals(sequence), (self.num_trajs, 1, 1))

        # sampling introduces a small amount of error from sampling, but this error is inconsequential
        # and will asymtotically decrease with larger trajectories
        # and is easier than me refactoring...
        self.ref_phi_angles = np.tile(self.sample_angles("phi"), (self.num_trajs, 1, 1))
        self.ref_psi_angles = np.tile(self.sample_angles("psi"), (self.num_trajs, 1, 1))

    def sample_angles(self, angle):
        dist_selector = {
            "phi": self.gather_phi_reference_dihedrals(self.sequence),
            "psi": self.gather_psi_reference_dihedrals(self.sequence),
        }

        dihedral_hist = []

        for dihedral in range(dist_selector[angle].shape[0]):
            dihedral_angles = dist_selector[angle][dihedral, :]

            # GOAL: Generate samples that adhere to the underlying distribution
            # Step 1: Compute the distribution & bin centers

            hist, bin_edges = np.histogram(
                dihedral_angles, bins=self.bins, density=True
            )

            bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

            # Step 2: Calculate the cumulative distribution function (CDF)
            cdf = np.cumsum(hist * np.diff(bin_edges))
            cdf /= cdf[-1]  # Normalize the CDF

            # Step 3: Generate random values between 0 and 1 and interpolate to get corresponding bin values
            rand_values = np.random.random(size=self.nsamples)
            sampled_dihedrals = np.interp(rand_values, cdf, bin_centers)

            dihedral_hist.append(sampled_dihedrals)

        return np.array(dihedral_hist)

    def gather_phi_reference_dihedrals(self, sequence: str) -> np.ndarray:
        """Gather the reference phi dihedral angles for a given sequence.

        Parameters
        ----------
        sequence : str
            The amino acid sequence of the current trajectory.

        Returns
        -------
        np.ndarray
            The reference phi dihedral angles for the given sequence.
        """

        phi_angles = []
        for i, residue in enumerate(sequence):
            if i == 0:
                phi_preceeding_context = "A"
            else:
                phi_preceeding_context = sequence[i - 1]

            three_letter_residue = ONE_TO_THREE[phi_preceeding_context]
            approximate_residue = EV_RESIDUE_MAPPER[three_letter_residue]

            phi_angles.append(
                PHI_EV_ANGLES_DICT[three_letter_residue][approximate_residue]
            )

        return np.array(phi_angles)

    def gather_psi_reference_dihedrals(self, sequence: str) -> np.ndarray:
        """Gather the reference psi dihedral angles for a given sequence.

        Parameters
        ----------
        sequence : str
            The amino acid sequence of the current trajectory.

        Returns
        -------
        np.ndarray
            The reference psi dihedral angles for the given sequence.
        """

        psi_angles = []
        for i, residue in enumerate(sequence):
            if i == len(sequence) - 1:
                psi_subsequent_context = "A"
            else:
                psi_subsequent_context = sequence[i + 1]

            three_letter_residue = ONE_TO_THREE[psi_subsequent_context]
            approximate_residue = EV_RESIDUE_MAPPER[three_letter_residue]
            psi_angles.append(
                PSI_EV_ANGLES_DICT[three_letter_residue][approximate_residue]
            )

        return np.array(psi_angles)
