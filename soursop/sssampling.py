#  _____  ____  _    _ _____   _____  ____  _____
# / ____|/ __ \| |  | |  __ \ / ____|/ __ \|  __ \
# | (___ | |  | | |  | | |__) | (___ | |  | | |__)|
# \___ \| |  | | |  | |  _  / \___ \| |  | |  ___/
# ____) | |__| | |__| | | \ \ ____) | |__| | |
# |_____/ \____/ \____/|_|  \_\_____/ \____/|_|
# Jeffrey M. Lotthammer (Holehouse Lab)
# Alex Holehouse (Pappu Lab and Holehouse Lab) and Jared Lalmansing (Pappu lab)
# Simulation analysis package
## Copyright 2014 - 2026
##

import itertools
import os
from typing import List, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import transforms
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
    """Hellinger distance between two joint probability distributions.

    Defined via the Bhattacharyya coefficient
    :math:`BC = \\sum_i \\sqrt{p_i q_i}` as
    :math:`H = \\sqrt{1 - BC}`. The normalisation factor of
    :math:`1/\\sqrt{2}` used in :func:`hellinger_distance` is omitted here
    because the inputs are already joint (2D) probability surfaces that
    sum to 1.

    Parameters
    ----------
    p, q : array_like
        Joint probability distributions of the same shape. Each must sum
        to 1 (within numerical tolerance) for the result to be a valid
        Hellinger distance.

    Returns
    -------
    float
        Hellinger distance in ``[0, 1]``; 0 means identical distributions
        and 1 means disjoint supports.

    Example
    -------
    >>> import numpy as np
    >>> from soursop.sssampling import compute_joint_hellinger_distance
    >>> p = np.eye(2) / 2
    >>> compute_joint_hellinger_distance(p, p)
    0.0
    """
    # Compute the Bhattacharyya coefficient
    b_coefficient = np.sum(np.sqrt(p * q))

    # Compute the Hellinger's distance - note this doesn't need the normalization by sqrt(2)
    distance = np.sqrt(1 - b_coefficient)

    return distance


def hellinger_distance(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Hellinger distance(s) between pairs of 1D probability distributions.

    For each pair the distance is

    .. math::
        H(P, Q) = \\frac{1}{\\sqrt{2}} \\sqrt{\\sum_{i=1}^{k} (\\sqrt{p_i} - \\sqrt{q_i})^2}

    which lies in ``[0, 1]``. The reduction is taken along the last axis
    of ``p`` / ``q``, so passing higher-rank arrays computes one distance
    per leading-axis slice.

    Parameters
    ----------
    p, q : np.ndarray
        Probability distributions of identical shape. The last axis is
        treated as the distribution axis; any leading axes are broadcast.

    Returns
    -------
    np.ndarray
        Hellinger distance(s) with shape ``p.shape[:-1]``.

    Example
    -------
    >>> import numpy as np
    >>> from soursop.sssampling import hellinger_distance
    >>> hellinger_distance(np.array([0.5, 0.5]), np.array([0.5, 0.5]))
    0.0
    >>> # per-residue distances for an (n_residues, n_bins) PDF stack
    >>> pdf_a = np.full((10, 20), 1/20)
    >>> pdf_b = np.full((10, 20), 1/20)
    >>> hellinger_distance(pdf_a, pdf_b).shape
    (10,)
    """
    # Ensure that p and q are NumPy arrays
    p = np.asarray(p)
    q = np.asarray(q)

    # Compute the Hellinger distance
    numerator = np.sum(np.square(np.sqrt(p) - np.sqrt(q)), axis=-1)
    denominator = np.sqrt(2)
    return np.sqrt(numerator) / denominator


def rel_entropy(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Kullback-Leibler relative entropy :math:`D_{KL}(P || Q)`.

    Computed via ``scipy.special.rel_entr`` (which handles ``p == 0`` and
    ``q == 0`` correctly), summed along the last axis. Asymmetric in
    ``p`` and ``q``; the result is always non-negative and is 0 only when
    ``p == q`` almost everywhere.

    Parameters
    ----------
    p, q : np.ndarray
        Probability distributions of identical shape. The last axis is
        the distribution axis; leading axes are broadcast.

    Returns
    -------
    np.ndarray
        Relative entropy values with shape ``p.shape[:-1]``, in nats.

    Example
    -------
    >>> import numpy as np
    >>> from soursop.sssampling import rel_entropy
    >>> rel_entropy(np.array([0.5, 0.5]), np.array([0.5, 0.5]))
    0.0
    >>> rel_entropy(np.array([0.9, 0.1]), np.array([0.5, 0.5]))
    0.368
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
        force_sequential: bool = False,
        **kwargs: dict,
    ):
        """Compare sampling quality of one or more trajectories against a reference.

        The reference can be a limiting-polymer-model ensemble, a wild-type
        simulation, or any other set of trajectories. If a ``reference_list``
        is not supplied, SOURSOP falls back to the precomputed excluded-
        volume (EV) limiting polymer angles tabulated in ``ssdata``.

        On construction, the class loads (or truncates, if requested) every
        trajectory, computes phi/psi dihedrals for the chosen ``proteinID``,
        and stores them for downstream methods such as
        :meth:`compute_dihedral_hellingers`, :meth:`compute_frac_helicity`,
        and :meth:`quality_plot`.

        Parameters
        ----------
        traj_list : list of str
            Trajectory file paths (xtc / dcd) for the simulated ensembles.
        reference_list : list of str or None, optional
            Trajectory file paths for the reference ensembles. If ``None``,
            the precomputed EV limiting-polymer dihedrals are used as the
            reference.
        top_file : str, optional
            Topology PDB for the simulated trajectories. Default
            ``"__START.pdb"``.
        ref_top : str or None, optional
            Topology PDB for the reference trajectories. Only required when
            ``reference_list`` is supplied.
        method : {'2D angle distributions', '1D angle distributions'}, optional
            Histogram strategy used when computing Hellinger distances and
            relative entropies. Default ``'2D angle distributions'``.
        bwidth : float, optional
            Histogram bin width in radians. Default ``deg2rad(15)``.
        proteinID : int, optional
            Index into each trajectory's ``proteinTrajectoryList`` that
            picks the chain to analyse. Default 0.
        n_cpus : int or None, optional
            Number of worker processes for parallel trajectory loading.
            None (default) uses all CPUs reported by ``os.cpu_count()``.
        truncate : bool, optional
            If True, slice every trajectory to the minimum length across
            the input set before computing dihedrals. Useful for mid-run
            analysis. Default False.
        force_sequential : bool, optional
            If True, load trajectories one-by-one rather than in parallel.
            Default False.
        **kwargs : dict
            Extra keyword arguments forwarded to :class:`SSTrajectory` (e.g.
            ``stride``).

        Raises
        ------
        SSException
            If ``method`` is not one of the allowed options, ``bwidth`` is
            out of range, or ``traj_list`` is empty.

        Example
        -------
        >>> from soursop.sssampling import SamplingQuality
        >>> sq = SamplingQuality(
        ...     traj_list=['rep0/traj.xtc', 'rep1/traj.xtc'],
        ...     top_file='topology.pdb',
        ... )
        >>> hellingers = sq.compute_dihedral_hellingers()
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
        self.force_sequential = force_sequential
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
                # __truncate_trajectories rebuilds each trajectory from the
                # single selected protein, so after truncation that protein is
                # at index 0 of the new proteinTrajectoryList. Reset proteinID
                # so downstream indexing (dihedrals, sequence) is not out of
                # range for the original proteinID > 0.
                self.proteinID = 0

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
                # see note above: after truncation the selected protein is at
                # index 0 of every rebuilt trajectory.
                self.proteinID = 0

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

            # Align the EV reference count with the trajectory dihedral count.
            # The reference is gathered once per residue, but get_angles('phi')
            # / ('psi') return one fewer angle than there are residues for an
            # UNCAPPED chain (phi is undefined at the N-terminal residue, psi at
            # the C-terminal residue). Without this the EV path crashed on
            # uncapped inputs with a broadcast error (n_res vs n_res-1). Trim
            # the reference from the phi-start / psi-end so the arrays align.
            self.ref_phi_angles = self.__align_reference_angles(
                precomputed_interface.ref_phi_angles, self.phi_angles, drop="first"
            )
            self.ref_psi_angles = self.__align_reference_angles(
                precomputed_interface.ref_psi_angles, self.psi_angles, drop="last"
            )

    @staticmethod
    def __align_reference_angles(reference, trajectory, drop):
        """Trim a per-residue reference angle array to the trajectory's count.

        ``reference`` and ``trajectory`` are ``(n_trajs, n_residues, n_samples)``
        arrays. When the reference has exactly one extra residue (the uncapped
        case), drop it from the terminus where the corresponding backbone
        dihedral is undefined (``'first'`` for phi, ``'last'`` for psi).

        Parameters
        ----------
        reference : np.ndarray
            Reference dihedral angles, ``(n_trajs, n_ref_res, n_samples)``.
        trajectory : np.ndarray
            Trajectory dihedral angles, ``(n_trajs, n_traj_res, n_samples)``.
        drop : {'first', 'last'}
            Terminus to trim from when the reference has one extra residue.

        Returns
        -------
        np.ndarray
            The reference array trimmed to match ``trajectory``'s residue axis.
        """
        n_ref = reference.shape[1]
        n_traj = trajectory.shape[1]
        if n_ref == n_traj + 1:
            return reference[:, 1:, :] if drop == "first" else reference[:, :-1, :]
        return reference

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
            if self.force_sequential:
                # Load trajectories sequentially
                self.trajs = []
                for traj in self.traj_list:
                    self.trajs.append(
                        SSTrajectory([traj], pdb_filename=self.top, **self.kwargs)
                    )

                if self.reference_list:
                    self.ref_trajs = []
                    for ref_traj in self.reference_list:
                        self.ref_trajs.append(
                            SSTrajectory(
                                [ref_traj], pdb_filename=self.ref_top, **self.kwargs
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
        """Per-residue fractional helicity for every loaded trajectory and reference.

        Helicity is taken directly from
        :meth:`SSProtein.get_secondary_structure_DSSP` (the helix column of
        the DSSP summary). If no reference trajectories were supplied,
        the reference helicity is zeros — the precomputed EV polymer
        reference is dihedral-only and has no DSSP equivalent.

        Results are cached on the instance; pass ``recompute=True`` to
        bypass the cache.

        Parameters
        ----------
        proteinID : int, optional
            Index of the chain in each trajectory's ``proteinTrajectoryList``.
            Default 0.
        recompute : bool, optional
            If True, ignore any cached result and recompute. Default False.

        Returns
        -------
        tuple of (np.ndarray, np.ndarray)
            ``(trj_helicity, ref_helicity)`` each of shape
            ``(n_trajectories, n_residues)``. When no reference was
            supplied, ``ref_helicity`` is all zeros.

        Example
        -------
        >>> trj_h, ref_h = sq.compute_frac_helicity()
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
        """Per-residue Hellinger distance between simulated and reference dihedrals.

        Behaviour depends on the ``method`` set at construction:

        * ``'2D angle distributions'``: histograms the joint
          ``(phi, psi)`` distribution per residue, then computes one
          joint-Hellinger distance per (trajectory, residue) pair via
          :func:`compute_joint_hellinger_distance`. Returns a shape
          ``(n_trajectories, n_residues)`` array.
        * ``'1D angle distributions'``: histograms phi and psi
          separately and computes a Hellinger distance for each. Returns
          a shape ``(2, n_trajectories, n_residues)`` array stacked as
          ``[phi_hellingers, psi_hellingers]``.

        Returns
        -------
        np.ndarray
            Hellinger distances, shape depending on ``method`` (see above).

        Raises
        ------
        NotImplementedError
            If ``method`` is not one of the two supported strings.

        Example
        -------
        >>> H = sq.compute_dihedral_hellingers()
        >>> H.shape         # for 2D angle distributions
        (3, 56)
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
        """Per-residue Kullback-Leibler relative entropy between simulated and reference dihedrals.

        Histograms phi and psi 1D distributions independently, then
        computes :math:`D_{KL}(P || Q)` per residue via :func:`rel_entropy`.

        Returns
        -------
        np.ndarray
            Array of shape ``(2, n_trajectories, n_residues)`` stacked as
            ``[phi_rel_entropy, psi_rel_entropy]``. Values are in nats.

        Example
        -------
        >>> rel_e = sq.compute_dihedral_rel_entropy()
        >>> rel_e.shape
        (2, 3, 56)
        """

        phi_trj_pdfs = self.compute_pdf(self.phi_angles, bins=self.bins)
        phi_ref_trj_pdfs = self.compute_pdf(self.ref_phi_angles, bins=self.bins)

        psi_trj_pdfs = self.compute_pdf(self.psi_angles, bins=self.bins)
        psi_ref_trj_pdfs = self.compute_pdf(self.ref_psi_angles, bins=self.bins)

        phi_rel_entr = rel_entropy(phi_trj_pdfs, phi_ref_trj_pdfs)
        psi_rel_entr = rel_entropy(psi_trj_pdfs, psi_ref_trj_pdfs)

        return np.array((phi_rel_entr, psi_rel_entr))

    def compute_series_of_histograms_along_axis(
        self, data: np.ndarray, bins: np.ndarray, axis: int = 0
    ):
        """2D ``(phi, psi)`` PDFs for every (trajectory, residue) pair.

        Builds an ``n_trajectories x n_residues`` grid of 2D joint
        histograms (one per residue) and normalises them so each is a
        probability density. The result is the per-pair PDF stack
        consumed by :meth:`compute_dihedral_hellingers` in 2D mode.

        Parameters
        ----------
        data : np.ndarray
            4D array of shape ``(2, n_trajectories, n_residues, n_frames)``
            where the leading axis stacks ``[phi_angles, psi_angles]``.
        bins : np.ndarray
            1D bin edges shared by both phi and psi axes.
        axis : int, optional
            Retained for API compatibility; the function always reduces
            over the frame axis internally. Default 0.

        Returns
        -------
        np.ndarray
            PDFs of shape
            ``(n_trajectories, n_residues, len(bins)-1, len(bins)-1)``.
            Each ``[i, j]`` slice is a normalised joint phi/psi
            distribution that sums to 1.

        Example
        -------
        >>> pdfs = sq.compute_series_of_histograms_along_axis(
        ...     np.array([sq.phi_angles, sq.psi_angles]),
        ...     bins=sq.bins,
        ... )
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
        """Per-residue 1D probability density histograms.

        Operates on either a 2D ``(n_residues, n_frames)`` array or a 3D
        ``(n_trajectories, n_residues, n_frames)`` stack. Each
        residue-level histogram is normalised by ``np.histogram(...,
        density=True)`` and rescaled by the bin width (in degrees), so
        each row sums to ~1.

        Parameters
        ----------
        arr : np.ndarray
            2D or 3D angle array. The last axis is the frame axis.
        bins : np.ndarray
            1D array of bin edges.

        Returns
        -------
        np.ndarray
            * Input 2D -> output ``(n_residues, len(bins) - 1)``.
            * Input 3D -> output
              ``(n_trajectories, n_residues, len(bins) - 1)``.

        Example
        -------
        >>> phi_pdfs = sq.compute_pdf(sq.phi_angles, bins=sq.bins)
        """
        # Lambda function is used to ignore the bin edges returned by np.histogram at index 1
        # xhistogram is ~2x faster, but introduces depedency - keeping lambda function for legacy for now

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
        """All-vs-all 2D joint-dihedral Hellinger distances across trajectories.

        Histograms the joint ``(phi, psi)`` distribution per residue for every
        loaded trajectory, then forms every pairwise trajectory combination
        (using ``itertools.combinations``) and computes a per-residue
        Hellinger distance for each pair. With a single trajectory this
        degenerates to a 1:1 self-comparison (which is always zero) and is
        useful only as a sanity check.

        Parameters
        ----------
        metric : str, optional
            Currently only ``'hellingers'`` is implemented. Default
            ``'hellingers'``.
        recompute : bool, optional
            Currently unused (accepted for API symmetry with
            :meth:`get_all_to_all_trj_comparisons`). Default False.

        Returns
        -------
        np.ndarray
            Shape ``(n_combinations, n_residues)`` of pairwise per-residue
            Hellinger distances in ``[0, 1]``.

        Example
        -------
        >>> mat = sq.get_all_to_all_2d_trj_comparison()
        >>> mat.shape   # 3 trajs -> C(3,2) == 3 pairs
        (3, 56)
        """
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
        """All-vs-all per-residue dihedral comparisons (separate phi and psi).

        Builds the per-residue 1D phi and psi PDFs for every trajectory,
        enumerates every pairwise trajectory combination, and computes a
        per-residue Hellinger distance or relative entropy for each pair.
        The two dihedrals are kept separate (unlike
        :meth:`get_all_to_all_2d_trj_comparison`).

        Parameters
        ----------
        metric : {'hellingers', 'relative entropy'}, optional
            Which divergence to compute. Default ``'hellingers'``.
        recompute : bool, optional
            If True, ignore cached PDFs from ``self.trj_pdfs`` and rebuild
            them. Default False.

        Returns
        -------
        tuple of (pd.DataFrame, pd.DataFrame)
            ``(phi_df, psi_df)`` each of shape
            ``(n_combinations, n_residues)`` containing the chosen metric
            for every pairwise comparison.

        Raises
        ------
        NotImplementedError
            If ``metric`` is not one of the two supported strings.

        Example
        -------
        >>> phi_df, psi_df = sq.get_all_to_all_trj_comparisons()
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
        """Histogram bin edges spanning ``[-180, 180]`` degrees.

        Constructs the bin edges used by every histogram-based method on
        this class. Uses ``self.bwidth`` (in radians) converted to
        degrees and rounded to handle floating-point error so the final
        edge lands cleanly on 180.

        Returns
        -------
        np.ndarray
            1D array of bin edges in degrees, monotonically increasing,
            starting at -180 and ending at 180.

        Example
        -------
        >>> sq.get_degree_bins()       # bwidth = 15 degrees
        array([-180., -165., ...,  165.,  180.])
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
        """Plot a four-panel sampling-quality summary figure.

        The four panels are:

        * **A** - per-residue Hellinger distance vs. the chosen reference
          (e.g. excluded-volume limit) with per-trajectory points and
          across-trajectory mean.
        * **B** - per-residue all-vs-all trajectory Hellinger distances.
        * **C** - fractional helicity (simulated trajectories + reference).
        * **D** - paired comparison panel (configurable).

        Layout is mosaic ``"AABB;CCDD"``. The chosen ``dihedral`` selector
        controls which of phi / psi / joint 2D is shown.

        Parameters
        ----------
        increment : int, optional
            X-axis tick stride (residues). Default 5.
        figsize : tuple of (int, int), optional
            Figure dimensions in inches. Default ``(7, 5)``.
        dpi : int, optional
            Output DPI for ``savefig``. Default 400.
        panel_labels : bool, optional
            If True, add A/B/C/D panel labels for manuscript figures.
            Default False.
        fontsize : int, optional
            Font size used for tick labels, titles, and axis labels.
            Default 10.
        save_dir : str or None, optional
            If given, write the figure to ``<save_dir>/<figname>``.
            Default None (no file written; figure is returned only).
        dihedral : {'2D', 'phi', 'psi'} or None, optional
            Which dihedral comparison to plot. ``'2D'`` requires
            ``method='2D angle distributions'``. Default ``'2D'``.
        figname : str, optional
            File name (joined with ``save_dir``). Default
            ``'hellingers.pdf'``.

        Returns
        -------
        tuple
            ``(fig, axd)`` — the matplotlib figure and the mosaic Axes
            dictionary keyed by ``'A'``, ``'B'``, ``'C'``, ``'D'``.

        Raises
        ------
        ValueError
            If ``method='1D angle distributions'`` is paired with
            ``dihedral='2D'``.
        NotImplementedError
            If a requested combination of method and dihedral isn't yet
            supported.

        Example
        -------
        >>> fig, axd = sq.quality_plot(dihedral='phi', save_dir='./figs')
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
        """Per-residue PDFs from the simulated trajectories' dihedral angles.

        Builds (and memoises) the three PDF stacks used by Hellinger /
        relative-entropy calculations against the trajectories:

        * ``'trj_phi_pdfs'`` — 1D phi histogram per residue.
        * ``'trj_psi_pdfs'`` — 1D psi histogram per residue.
        * ``'joint'`` (default) — 2D ``(phi, psi)`` histogram per residue.

        On every call all three are populated on the cache; the selector
        determines which is returned. ``recompute=True`` bypasses the
        cache.

        Parameters
        ----------
        dihedral : {'joint', 'trj_phi_pdfs', 'trj_psi_pdfs'}, optional
            Which PDF stack to return. Default ``'joint'``.
        recompute : bool, optional
            If True, ignore any cached PDFs and rebuild. Default False.

        Returns
        -------
        np.ndarray
            * For ``'trj_phi_pdfs'`` / ``'trj_psi_pdfs'``:
              ``(n_trajectories, n_residues, n_bins)``.
            * For ``'joint'``:
              ``(n_trajectories, n_residues, n_bins, n_bins)``.

        Raises
        ------
        NotImplementedError
            If ``dihedral`` is not one of the three allowed strings.

        Example
        -------
        >>> phi_pdfs = sq.trj_pdfs(dihedral='trj_phi_pdfs')
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
        """Per-residue PDFs from the reference trajectories' dihedral angles.

        The reference analogue of :meth:`trj_pdfs`. The three accepted
        selectors here are ``'ref_phi_pdfs'``, ``'ref_psi_pdfs'``, and
        ``'joint'`` (default). When no reference trajectories were
        supplied to the constructor, the reference angles came from the
        precomputed excluded-volume polymer model.

        Parameters
        ----------
        dihedral : {'joint', 'ref_phi_pdfs', 'ref_psi_pdfs'}, optional
            Which PDF stack to return. Default ``'joint'``.
        recompute : bool, optional
            If True, ignore any cached PDFs and rebuild. Default False.

        Returns
        -------
        np.ndarray
            * For ``'ref_phi_pdfs'`` / ``'ref_psi_pdfs'``:
              ``(n_trajectories, n_residues, n_bins)``.
            * For ``'joint'``:
              ``(n_trajectories, n_residues, n_bins, n_bins)``.

        Raises
        ------
        NotImplementedError
            If ``dihedral`` is not one of the three allowed strings.

        Example
        -------
        >>> ref_phi = sq.ref_pdfs(dihedral='ref_phi_pdfs')
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
        """Cached accessor for per-residue Hellinger distances.

        Thin wrapper around :meth:`compute_dihedral_hellingers` that
        memoises the result on first call. Pass ``recompute=True`` to
        invalidate the cache and force a fresh computation.

        Parameters
        ----------
        recompute : bool, optional
            If True, rebuild the Hellinger distances from scratch.
            Default False.

        Returns
        -------
        np.ndarray
            Shape and meaning depend on the SamplingQuality method:

            * ``'2D angle distributions'``:
              ``(n_trajectories, n_residues)`` joint Hellinger distances.
            * ``'1D angle distributions'``:
              ``(2, n_trajectories, n_residues)`` stacked as
              ``[phi_hellingers, psi_hellingers]``.

        Example
        -------
        >>> H = sq.hellingers_distances()
        """
        selector = "hellingers"

        if selector not in self.__precomputed or recompute is True:
            self.__precomputed[selector] = self.compute_dihedral_hellingers()

        return self.__precomputed[selector]

    def fractional_helicity(self, recompute=False):
        """Cached accessor for per-residue fractional helicity.

        Thin wrapper around :meth:`compute_frac_helicity` that returns
        results from the instance cache if available. The
        ``recompute=True`` flag is forwarded so the underlying
        computation re-runs.

        Parameters
        ----------
        recompute : bool, optional
            If True, bypass the cache and recompute. Default False.

        Returns
        -------
        tuple of (np.ndarray, np.ndarray)
            ``(trj_helicity, ref_helicity)`` each of shape
            ``(n_trajectories, n_residues)``.

        Example
        -------
        >>> trj_h, ref_h = sq.fractional_helicity()
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
    """Reference dihedral provider backed by the precomputed EV polymer model.

    Used by :class:`SamplingQuality` when no explicit ``reference_list``
    of reference trajectories is supplied. For each residue the
    appropriate phi / psi distribution is looked up from the
    excluded-volume polymer reference tables in ``ssdata``, then
    inverse-CDF sampled to produce a synthetic per-trajectory angle
    array of the same shape as the simulated trajectories' dihedral
    arrays — so the downstream Hellinger and relative-entropy code
    treats it identically.

    Parameters
    ----------
    sequence : str
        Single-letter amino-acid sequence of the trajectory chain
        (caps already stripped).
    bins : np.ndarray
        Histogram bin edges (in degrees) used for the inverse-CDF
        sampling. Should match the SamplingQuality instance's bins.
    num_trajs : int
        Number of simulated-trajectory replicas to mimic. The
        precomputed angles are tiled across this dimension.
    nsamples : int
        Number of synthetic frames to generate per replica.

    Attributes
    ----------
    ref_phi_angles, ref_psi_angles : np.ndarray
        Inverse-CDF-sampled reference angles of shape
        ``(num_trajs, n_residues, nsamples)``.

    Example
    -------
    >>> from soursop.sssampling import PrecomputedDihedralInterface
    >>> ev = PrecomputedDihedralInterface(
    ...     sequence='AAAAAAAA', bins=np.arange(-180, 181, 15),
    ...     num_trajs=3, nsamples=1000,
    ... )
    >>> ev.ref_phi_angles.shape
    (3, 8, 1000)
    """

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
        """Inverse-CDF sample reference dihedrals to match a target sample count.

        Builds a per-residue histogram from the precomputed reference
        angles, normalises it to a CDF, and inverse-CDF samples
        ``self.nsamples`` synthetic angles per residue using
        ``numpy.random``. Used at construction time to populate
        :attr:`ref_phi_angles` and :attr:`ref_psi_angles`.

        Parameters
        ----------
        angle : {'phi', 'psi'}
            Which backbone dihedral to sample.

        Returns
        -------
        np.ndarray
            Array of shape ``(n_residues, self.nsamples)`` of synthetic
            dihedral angles (in degrees) drawn from the precomputed
            reference distributions.

        Example
        -------
        >>> ev.sample_angles('phi').shape
        (8, 1000)
        """
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
        """Look up the excluded-volume reference phi distribution for each residue.

        For each position ``i``, the relevant phi distribution depends on
        the chemical context of residue ``i-1`` (the residue preceding
        the rotatable phi bond). The lookup maps the preceding residue
        to its EV-table key via :data:`EV_RESIDUE_MAPPER` and pulls the
        distribution from :data:`PHI_EV_ANGLES_DICT`. For position 0 we
        substitute alanine as the preceding context.

        Parameters
        ----------
        sequence : str
            One-letter amino-acid sequence (no caps).

        Returns
        -------
        np.ndarray
            Array of shape ``(n_residues, n_reference_samples)`` where
            each row is the EV reference phi distribution for that
            residue.

        Example
        -------
        >>> ev.gather_phi_reference_dihedrals('AAAA').shape
        (4, 50000)
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
        """Look up the excluded-volume reference psi distribution for each residue.

        Mirror of :meth:`gather_phi_reference_dihedrals` for psi: the
        distribution at position ``i`` depends on the *following*
        residue ``i+1`` (because psi rotates the i-(i+1) bond). For the
        last position we substitute alanine as the following context.
        Reference distributions come from :data:`PSI_EV_ANGLES_DICT`.

        Parameters
        ----------
        sequence : str
            One-letter amino-acid sequence (no caps).

        Returns
        -------
        np.ndarray
            Array of shape ``(n_residues, n_reference_samples)`` where
            each row is the EV reference psi distribution for that
            residue.

        Example
        -------
        >>> ev.gather_psi_reference_dihedrals('AAAA').shape
        (4, 50000)
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
