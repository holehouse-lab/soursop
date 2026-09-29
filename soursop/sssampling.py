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
    :math:`H = \\sqrt{1 - BC}`. For any pair of normalised distributions
    this is algebraically identical to the
    :math:`\\frac{1}{\\sqrt{2}}\\lVert\\sqrt{p} - \\sqrt{q}\\rVert_2` form used
    in :func:`hellinger_distance` (nothing is omitted); the 2D inputs are
    simply summed over both axes.

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
    # For identical distributions floating point can push BC a hair above 1,
    # which would give sqrt of a negative number (nan); clamp at 0.
    distance = np.sqrt(max(1.0 - b_coefficient, 0.0))

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


def smooth_pdf(pdf: np.ndarray, pseudocount: float = 1e-6) -> np.ndarray:
    """Add a pseudocount to every bin of a PDF and renormalise.

    Histogrammed dihedral distributions routinely contain empty bins, and
    the Kullback-Leibler divergence is undefined (infinite) wherever the
    reference bin is empty but the trajectory bin is not. Adding a small
    constant mass to every bin and renormalising along the last axis
    gives every distribution full support so the divergence is finite.

    Parameters
    ----------
    pdf : np.ndarray
        Probability mass array. The last axis is the distribution axis;
        leading axes are treated independently.
    pseudocount : float, optional
        Mass added to every bin before renormalising. ``0`` returns the
        input unchanged. Default ``1e-6``.

    Returns
    -------
    np.ndarray
        Smoothed array of the same shape, each row summing to 1.

    Example
    -------
    >>> import numpy as np
    >>> from soursop.sssampling import smooth_pdf
    >>> smooth_pdf(np.array([1.0, 0.0]), pseudocount=0.5)
    array([0.75, 0.25])
    """
    pdf = np.asarray(pdf, dtype=float)
    if not np.isfinite(pseudocount) or pseudocount < 0:
        raise SSException(
            f"pseudocount must be a finite, non-negative number. Received {pseudocount}"
        )
    if pseudocount == 0:
        return pdf
    smoothed = pdf + pseudocount
    return smoothed / np.sum(smoothed, axis=-1, keepdims=True)


def rel_entropy(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Kullback-Leibler relative entropy :math:`D_{KL}(P || Q)`.

    Computed via ``scipy.special.rel_entr`` summed along the last axis.
    Asymmetric in ``p`` and ``q``; the result is always non-negative and
    is 0 only when ``p == q`` almost everywhere.

    Note that KL divergence is undefined on disjoint support: any bin
    with ``p > 0`` and ``q == 0`` contributes ``inf``. This function
    returns that raw value and does no smoothing. The
    :class:`SamplingQuality` methods that call it
    (:meth:`~SamplingQuality.compute_dihedral_rel_entropy` and
    :meth:`~SamplingQuality.get_all_to_all_trj_comparisons`) apply a
    pseudocount via :func:`smooth_pdf` by default, because histogrammed
    dihedral distributions almost always contain empty bins.

    Parameters
    ----------
    p, q : np.ndarray
        Probability distributions of identical shape. The last axis is
        the distribution axis; leading axes are broadcast.

    Returns
    -------
    np.ndarray
        Relative entropy values with shape ``p.shape[:-1]``, in nats.
        ``inf`` wherever ``q`` is zero on part of ``p``'s support.

    Example
    -------
    >>> import numpy as np
    >>> from soursop.sssampling import rel_entropy
    >>> rel_entropy(np.array([0.5, 0.5]), np.array([0.5, 0.5]))
    0.0
    >>> round(float(rel_entropy(np.array([0.9, 0.1]), np.array([0.5, 0.5]))), 4)
    0.3681
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
        seed: Union[int, None] = None,
        verbose: bool = True,
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
        top_file : str or list of str, optional
            Topology PDB for the simulated trajectories. Either a single
            path shared by every trajectory, or a list with one entry per
            element of ``traj_list``. Default ``"__START.pdb"``.
        ref_top : str or list of str or None, optional
            Topology PDB for the reference trajectories; a single path or
            one per element of ``reference_list``. Only required when
            ``reference_list`` is supplied.
        method : {'2D angle distributions', '1D angle distributions'}, optional
            Histogram strategy used when computing Hellinger distances and
            relative entropies. Default ``'2D angle distributions'``.
        bwidth : float, optional
            Histogram bin width in radians. Must lie in ``(0, pi]`` (at
            least two bins) and divide 360 degrees into an integer number
            of bins (so 15, 10, 7.5 or 5 degrees are fine, 7 degrees is
            not). Default ``deg2rad(15)``.
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
        seed : int or None, optional
            Seed for the random number generator used to resample the
            precomputed EV reference dihedrals when no ``reference_list``
            is supplied. Passing an integer makes the EV reference (and
            hence every downstream Hellinger distance) reproducible. None
            (default) draws a fresh, unseeded generator.
        verbose : bool, optional
            If True, print a short status message when trajectories are
            truncated. Default True.
        **kwargs : dict
            Extra keyword arguments forwarded verbatim to
            :class:`SSTrajectory`, e.g. ``extra_valid_residue_names``,
            ``explicit_residue_checking``, ``protein_grouping`` or
            ``pdblead``. Note that there is no stride or frame-weighting
            support here: every frame of every trajectory is used.

        Attributes
        ----------
        residue_indices : list of int
            0-based residue indices (in the analysed chain's own numbering) of
            the columns of every dihedral array. Every residue that carries
            both a phi and a psi: all non-cap residues on a capped chain, the
            interior residues on an uncapped one.

        Raises
        ------
        SSException
            If ``method`` is not one of the allowed options; if ``bwidth``
            is outside ``(0, pi]`` or does not divide 360 degrees into a
            whole number of bins; if ``traj_list`` is empty; if
            ``reference_list`` is given without ``ref_top``; if
            ``reference_list`` (or a per-trajectory ``top_file`` /
            ``ref_top`` list) does not match the length of ``traj_list``;
            if ``seed`` is not None or a non-negative integer; if
            ``proteinID`` does not index a protein in every loaded
            trajectory; if the trajectories (and references) are not all
            the same chain with the same caps; or if the loaded trajectories
            do not all have the same number of frames and ``truncate`` is
            False.

        Example
        -------
        >>> from soursop.sssampling import SamplingQuality
        >>> # the default loader uses multiprocessing; in a script, guard the
        >>> # call with `if __name__ == "__main__":` (needed on macOS and
        >>> # Windows), or pass force_sequential=True
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
        self.seed = seed
        self.verbose = verbose
        self.kwargs = kwargs

        self.__precomputed = {}

        # validate before building the bin grid: get_degree_bins relies on
        # bwidth being sane, and a tiny bwidth would otherwise try to
        # allocate an enormous array before we ever got to check it.
        self.__validate_arguments()

        self.bins = self.get_degree_bins()

        self.__load_trajectories()

        # proteinID can only be checked once we know how many proteins each
        # trajectory holds; an out-of-range value used to raise a bare
        # IndexError deep inside the dihedral code
        self.__check_protein_id(self.proteinID)

        self.__check_frame_counts()

        # if reference trajectories have been provided
        # then self.ref_trajs should have been initialized.
        if self.reference_list:
            # if truncate is True,
            # then match the lengths of the trajectories before computing dihedrals
            if self.truncate:
                # rebuilds every trajectory with the same proteins, so
                # proteinID still refers to the same chain
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

            # remove caps from sequence if present: ACE '<', NME '>', and the
            # FOR '(' / NH2 ')' caps, which used to be left in
            sequence = "".join(c for c in sequence if c not in "<>()")

            precomputed_interface = PrecomputedDihedralInterface(
                sequence,
                bins=self.bins,
                num_trajs=len(self.trajs),
                nsamples=len(self.trajs[0]),
                seed=self.seed,
            )

            # The reference tables hold one row per (cap-stripped) sequence
            # position. Select the rows for exactly the residues that carry
            # both a phi and a psi in the trajectories (self.residue_indices),
            # so trajectory column k and reference column k are the same
            # residue. On a capped chain that is every non-cap residue; on an
            # uncapped chain the two terminal residues are excluded.
            protein = self.trajs[0].proteinTrajectoryList[self.proteinID]
            offset = 1 if protein.ncap else 0
            rows = [r - offset for r in self.residue_indices]
            self.ref_phi_angles = precomputed_interface.ref_phi_angles[:, rows, :]
            self.ref_psi_angles = precomputed_interface.ref_psi_angles[:, rows, :]

    def __check_protein_id(self, proteinID):
        """Raise unless ``proteinID`` indexes a protein in every trajectory."""
        if isinstance(proteinID, bool) or not isinstance(proteinID, (int, np.integer)):
            raise SSException(f"proteinID must be an integer. Received {proteinID!r}")
        all_trajs = list(self.trajs)
        if self.reference_list:
            all_trajs.extend(self.ref_trajs)
        for trj in all_trajs:
            n_proteins = len(trj.proteinTrajectoryList)
            if not 0 <= proteinID < n_proteins:
                raise SSException(
                    f"proteinID={proteinID} is out of range: a loaded trajectory has "
                    f"{n_proteins} protein(s), so valid values are 0..{n_proteins - 1}."
                )

    def __validate_arguments(self):
        ssutils.validate_keyword_option(
            self.method, ["2D angle distributions", "1D angle distributions"], "method"
        )

        # a single path is a one-element list; a bare string used to be
        # iterated character by character
        if isinstance(self.traj_list, (str, os.PathLike)):
            self.traj_list = [self.traj_list]
        if isinstance(self.reference_list, (str, os.PathLike)):
            self.reference_list = [self.reference_list]
        self.traj_list = [os.fspath(t) for t in self.traj_list]
        if self.reference_list:
            self.reference_list = [os.fspath(t) for t in self.reference_list]

        # numpy's generator needs a non-negative integer seed; checking here
        # means a bad seed fails before any trajectory has been loaded
        if self.seed is not None and (
            isinstance(self.seed, bool)
            or not isinstance(self.seed, (int, np.integer))
            or self.seed < 0
        ):
            raise SSException(
                f"seed must be None or a non-negative integer. Received {self.seed!r}"
            )

        if self.reference_list and self.ref_top is None:
            raise SSException(
                "reference_list was given without ref_top: a topology file is needed "
                "to read the reference trajectories."
            )

        if (
            not np.isfinite(self.bwidth)
            or self.bwidth > 2 * np.pi
            or not self.bwidth > 0
        ):
            raise SSException(
                f"The bwidth parameter must be between 0 and 2*pi. Received {self.bwidth}"
            )

        # The bin grid must tile [-180, 180] exactly, so 360 degrees has to be
        # an integer multiple of the bin width. Anything else either leaves a
        # ragged final bin or (previously) got silently rounded to whole
        # degrees.
        n_bins = 2 * np.pi / self.bwidth
        # a single bin puts every angle in the same place, so every distance
        # is trivially zero
        if np.round(n_bins) < 2:
            raise SSException(
                "The bwidth parameter must give at least two bins (bwidth <= pi). "
                f"Received {np.rad2deg(self.bwidth):.6g} degrees."
            )
        if abs(n_bins - np.round(n_bins)) > 1e-6:
            raise SSException(
                "The bwidth parameter must divide 360 degrees into a whole "
                f"number of bins. Received {np.rad2deg(self.bwidth):.6g} degrees, "
                f"which gives {n_bins:.6g} bins. Use e.g. 15, 10, 7.5 or 5 degrees."
            )

        if not self.n_cpus:
            self.n_cpus = os.cpu_count()

        if len(self.traj_list) == 0:
            raise SSException(
                f"Input trajectory list must be non-empty.\
                    Received len(traj_list)={len(self.traj_list)}"
            )

        # Trajectories and references are paired one-to-one (zip), so the
        # two lists must have the same length. Previously a mismatch either
        # silently dropped the extra entries or, for a single trajectory with
        # several references, never set self.ref_trajs and failed later with
        # an AttributeError.
        if self.reference_list and len(self.reference_list) != len(self.traj_list):
            raise SSException(
                "reference_list must contain one reference trajectory per "
                "input trajectory (they are paired one-to-one). Received "
                f"len(traj_list)={len(self.traj_list)} and "
                f"len(reference_list)={len(self.reference_list)}."
            )

        # per-trajectory topology lists must line up with the trajectories. A
        # single str or pathlib.Path is shared by every trajectory (prior to
        # 2.0.6 a Path was mistaken for a list and raised a TypeError)
        if isinstance(self.top, os.PathLike):
            self.top = os.fspath(self.top)
        if isinstance(self.ref_top, os.PathLike):
            self.ref_top = os.fspath(self.ref_top)
        if not isinstance(self.top, str) and len(self.top) != len(self.traj_list):
            raise SSException(
                "top_file must be a single path or a list with one entry per "
                f"trajectory. Received len(top_file)={len(self.top)} and "
                f"len(traj_list)={len(self.traj_list)}."
            )
        if (
            self.reference_list
            and self.ref_top is not None
            and not isinstance(self.ref_top, str)
            and len(self.ref_top) != len(self.reference_list)
        ):
            raise SSException(
                "ref_top must be a single path or a list with one entry per "
                f"reference trajectory. Received len(ref_top)={len(self.ref_top)} "
                f"and len(reference_list)={len(self.reference_list)}."
            )

    def __load_trajectories(self):
        # weird thing I have to do to prevent issues with multiprocessing
        # parallel loading when there is only 1 trajectory to load
        # trajs/ref_trajs must be a list so they're iterables for __truncate_trajectories

        if len(self.traj_list) == 1:
            self.trajs = []
            self.trajs.append(
                SSTrajectory(
                    self.traj_list,
                    pdb_filename=self.__topology_for(self.top, 0),
                    **self.kwargs,
                )
            )

            # if the reference list has been provided initialize the reference trajectories
            # else the reference dihedrals will be assigned from precomputed dihedrals later.
            if not self.reference_list:
                pass

            elif len(self.reference_list) == 1:
                self.ref_trajs = []
                self.ref_trajs.append(
                    SSTrajectory(
                        self.reference_list,
                        pdb_filename=self.__topology_for(self.ref_top, 0),
                        **self.kwargs,
                    )
                )

        else:
            if self.force_sequential:
                # Load trajectories sequentially
                self.trajs = []
                for i, traj in enumerate(self.traj_list):
                    self.trajs.append(
                        SSTrajectory(
                            [traj],
                            pdb_filename=self.__topology_for(self.top, i),
                            **self.kwargs,
                        )
                    )

                if self.reference_list:
                    self.ref_trajs = []
                    for i, ref_traj in enumerate(self.reference_list):
                        self.ref_trajs.append(
                            SSTrajectory(
                                [ref_traj],
                                pdb_filename=self.__topology_for(self.ref_top, i),
                                **self.kwargs,
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

    @staticmethod
    def __topology_for(top, index):
        """Resolve the topology path for the trajectory at ``index``.

        ``top`` may be a single path (shared by every trajectory), a
        per-trajectory list, or None. Length agreement is checked in
        ``__validate_arguments``.
        """
        if top is None or isinstance(top, str):
            return top
        return top[index]

    def __check_frame_counts(self):
        """Raise if the loaded trajectories cannot be stacked into one array.

        Dihedral angles from every trajectory (and its paired reference,
        or the resampled EV reference) are stacked into a single array, so
        every trajectory must contain the same number of frames unless
        ``truncate=True`` was requested. Previously a mismatch surfaced as
        an opaque inhomogeneous-shape ValueError from numpy.
        """
        if self.truncate:
            return

        lengths = [trj.n_frames for trj in self.trajs]
        if self.reference_list:
            lengths.extend(ref_trj.n_frames for ref_trj in self.ref_trajs)

        if len(set(lengths)) > 1:
            raise SSException(
                "All trajectories (and reference trajectories) must contain the "
                f"same number of frames. Received frame counts {lengths}. Pass "
                "truncate=True to slice every trajectory to the shortest length."
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
        # Rebuild each trajectory from the first min_length frames of the
        # *whole* system, with the same protein-definition options it was
        # loaded with. Prior to 2.0.6 only the selected protein was kept and
        # the SSTrajectory keywords (protein_grouping,
        # extra_valid_residue_names, ...) were dropped, so chain detection
        # re-ran on the sliced protein: a protein_grouping protein spanning
        # two chains lost half its residues, and proteinID had to be reset
        # to 0, making self.trajs differ in structure depending on truncate.
        rebuild_keys = (
            "protein_grouping",
            "debug",
            "extra_valid_residue_names",
            "explicit_residue_checking",
            "swan_trajectory",
            "check_whole_molecules",
        )
        rebuild_kwargs = {k: v for k, v in self.kwargs.items() if k in rebuild_keys}
        # the full trajectories were already checked for split molecules on
        # load, so don't warn a second time for their first min_length frames
        rebuild_kwargs.setdefault("check_whole_molecules", False)

        def _truncated(trj):
            return SSTrajectory(TRJ=trj.traj[0 : self.min_length], **rebuild_kwargs)

        lengths = [trj.n_frames for trj in self.trajs]
        if self.reference_list:
            lengths.extend(ref_trj.n_frames for ref_trj in self.ref_trajs)
        self.min_length = int(np.min(lengths))

        temp_trajs = [_truncated(trj) for trj in self.trajs]
        temp_ref_trjs = (
            [_truncated(ref_trj) for ref_trj in self.ref_trajs]
            if self.reference_list
            else None
        )

        if self.verbose:
            print(
                f"Successfully truncated.\n\
                The shortest trajectory is: {self.min_length} frames.\
                All trajectories truncated to {self.min_length}"
            )

        return (temp_trajs, temp_ref_trjs)

    @staticmethod
    def __aligned_dihedrals(protein):
        """Phi and psi of one protein, restricted to residues that have both.

        ``SSProtein.get_angles('phi')`` is undefined for the first residue of
        an uncapped chain and ``get_angles('psi')`` for the last, so the two
        arrays are offset by one residue relative to each other. PENGUIN
        compares the *joint* (phi, psi) distribution per residue, so both
        arrays are cut down to the residues that carry both angles - every
        non-cap residue on a capped chain, the interior residues on an
        uncapped one.

        Parameters
        ----------
        protein : SSProtein

        Returns
        -------
        (phi, psi, residue_indices)
            ``phi`` and ``psi`` are ``(n_residues, n_frames)`` arrays in
            degrees; ``residue_indices`` is the list of 0-based residue
            indices (in the protein's own numbering) of the columns.
        """
        phi_atoms, phi = protein.get_angles("phi")
        psi_atoms, psi = protein.get_angles("psi")

        # the central residue of phi = [C(i-1), N(i), CA(i), C(i)] and of
        # psi = [N(i), CA(i), C(i), N(i+1)] is the residue of atom 1
        phi_res = [quad[1].residue.index for quad in phi_atoms]
        psi_res = [quad[1].residue.index for quad in psi_atoms]
        common = sorted(set(phi_res) & set(psi_res))
        if len(common) == 0:
            raise SSException(
                "No residue carries both a phi and a psi backbone dihedral; "
                "PENGUIN needs at least one such residue."
            )

        phi_lookup = {r: k for k, r in enumerate(phi_res)}
        psi_lookup = {r: k for k, r in enumerate(psi_res)}
        phi = np.asarray(phi)[[phi_lookup[r] for r in common]]
        psi = np.asarray(psi)[[psi_lookup[r] for r in common]]
        return phi, psi, common

    def __compute_dihedrals(
        self, proteinID: int = 0, precomputed: bool = False
    ) -> np.ndarray:
        """Compute residue-aligned phi/psi dihedrals for every loaded trajectory.

        Angles are taken from chain ``proteinID`` of each trajectory (and each
        reference trajectory when ``precomputed`` is False) and restricted to
        the residues that carry both a phi and a psi (see
        ``__aligned_dihedrals``), so that column ``k`` of every returned array
        is the same residue. The residue indices of the columns are stored on
        ``self.residue_indices``.

        Parameters
        ----------
        proteinID : int, optional
            Position of the chain in each ``SSTrajectory.proteinTrajectoryList``.
        precomputed : bool, optional
            If True only the simulated trajectories are processed (the
            reference comes from the precomputed EV tables).

        Returns
        -------
        np.ndarray
            ``(psi, ref_psi, phi, ref_phi)`` when ``precomputed`` is False,
            else ``(psi, phi)``; each element is
            ``(n_trajectories, n_residues, n_frames)``.
        """
        psi_angles, phi_angles = [], []
        ref_psi_angles, ref_phi_angles = [], []

        # Column k must be the same residue in every array. Matching residue
        # index lists is not enough on its own: e.g. an ACE-capped,
        # C-terminally uncapped chain and an uncapped, NME-capped one both
        # give [1, ..., n-2] but are offset by one residue, so the cap state
        # has to match as well. Residue *names* are deliberately not compared
        # (a wild-type reference for a mutant is a documented use).
        self.residue_indices = None
        caps = None
        for k, trj in enumerate(self.trajs):
            protein = trj.proteinTrajectoryList[proteinID]
            phi, psi, resids = self.__aligned_dihedrals(protein)
            if self.residue_indices is None:
                self.residue_indices = resids
                caps = (protein.ncap, protein.ccap)
            elif resids != self.residue_indices or (protein.ncap, protein.ccap) != caps:
                raise SSException(
                    f"Trajectory {k} ({self.traj_list[k]}) has phi/psi residues "
                    f"{resids} (N-cap {protein.ncap}, C-cap {protein.ccap}), but "
                    f"trajectory 0 has {self.residue_indices} (N-cap {caps[0]}, "
                    f"C-cap {caps[1]}). Every trajectory must be of the same chain so "
                    "that column k of the dihedral arrays is the same residue."
                )
            phi_angles.append(phi)
            psi_angles.append(psi)

        if precomputed:
            return np.array((psi_angles, phi_angles))

        for k, ref_trj in enumerate(self.ref_trajs):
            ref_protein = ref_trj.proteinTrajectoryList[proteinID]
            phi, psi, resids = self.__aligned_dihedrals(ref_protein)
            # the reference is compared residue by residue, so it must carry
            # the same residues (and the same caps, see above); previously a
            # length mismatch surfaced as an inhomogeneous-shape ValueError
            # from numpy
            if (
                resids != self.residue_indices
                or (ref_protein.ncap, ref_protein.ccap) != caps
            ):
                raise SSException(
                    f"Reference trajectory {k} ({self.reference_list[k]}) has phi/psi "
                    f"residues {resids} (N-cap {ref_protein.ncap}, C-cap "
                    f"{ref_protein.ccap}), but the trajectories have "
                    f"{self.residue_indices} (N-cap {caps[0]}, C-cap {caps[1]}). "
                    "References must be of the same chain as the trajectories they "
                    "are paired with."
                )
            ref_phi_angles.append(phi)
            ref_psi_angles.append(psi)

        return np.array((psi_angles, ref_psi_angles, phi_angles, ref_phi_angles))

    def compute_frac_helicity(
        self, proteinID: Union[int, None] = None, recompute: bool = False
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
        proteinID : int or None, optional
            Index of the chain in each trajectory's ``proteinTrajectoryList``.
            If None (default) the chain the :class:`SamplingQuality` instance
            was constructed with (``self.proteinID``) is used.
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
        # Default to the chain this SamplingQuality instance was built on.
        if proteinID is None:
            proteinID = self.proteinID
        self.__check_protein_id(proteinID)

        # the cache is keyed on the chain; prior to 2.0.6 a single entry was
        # shared, so asking for a second chain returned the first chain's
        # helicity (and vice versa for quality_plot)
        key = ("helicity", int(proteinID))
        if not recompute and key in self.__precomputed:
            return self.__precomputed[key]

        trj_helicity = np.array(
            [
                trj.proteinTrajectoryList[proteinID].get_secondary_structure_DSSP()[1]
                for trj in self.trajs
            ]
        )

        if self.reference_list:
            reference_helicity = np.array(
                [
                    ref_trj.proteinTrajectoryList[
                        proteinID
                    ].get_secondary_structure_DSSP()[1]
                    for ref_trj in self.ref_trajs
                ]
            )
        else:
            reference_helicity = np.zeros_like(trj_helicity)

        self.__precomputed[key] = (trj_helicity, reference_helicity)

        return self.__precomputed[key]

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
        SSException
            If ``method`` is not one of the two supported strings.

        Example
        -------
        >>> H = sq.compute_dihedral_hellingers()
        >>> H.shape         # 2D method, 3 trajectories of NTL9 (54 phi/psi residues)
        (3, 54)
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
            raise SSException(
                f"method '{self.method}' is not defined. Please use either "
                "'1D angle distributions' or '2D angle distributions'."
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

    def compute_dihedral_rel_entropy(self, pseudocount: float = 1e-6) -> np.ndarray:
        """Per-residue Kullback-Leibler relative entropy between simulated and reference dihedrals.

        Histograms phi and psi 1D distributions independently, then
        computes :math:`D_{KL}(P || Q)` per residue via :func:`rel_entropy`,
        with ``P`` the trajectory distribution and ``Q`` the reference.

        KL divergence is infinite whenever a reference bin is empty but
        the corresponding trajectory bin is populated, which with 15
        degree bins on finite ensembles is the norm rather than the
        exception. To keep the values usable, both distributions are
        smoothed with :func:`smooth_pdf` (a small pseudocount added to
        every bin, then renormalised) before the divergence is taken.

        Parameters
        ----------
        pseudocount : float, optional
            Mass added to every histogram bin before renormalising. Set to
            ``0`` to recover the raw, unsmoothed divergence (which may be
            ``inf``). Default ``1e-6``.

        Returns
        -------
        np.ndarray
            Array of shape ``(2, n_trajectories, n_residues)`` stacked as
            ``[phi_rel_entropy, psi_rel_entropy]``. Values are in nats.

        Example
        -------
        >>> rel_e = sq.compute_dihedral_rel_entropy()
        >>> rel_e.shape     # 3 trajectories of NTL9 (54 phi/psi residues)
        (2, 3, 54)
        """

        phi_trj_pdfs = self.compute_pdf(self.phi_angles, bins=self.bins)
        phi_ref_trj_pdfs = self.compute_pdf(self.ref_phi_angles, bins=self.bins)

        psi_trj_pdfs = self.compute_pdf(self.psi_angles, bins=self.bins)
        psi_ref_trj_pdfs = self.compute_pdf(self.ref_psi_angles, bins=self.bins)

        phi_rel_entr = rel_entropy(
            smooth_pdf(phi_trj_pdfs, pseudocount),
            smooth_pdf(phi_ref_trj_pdfs, pseudocount),
        )
        psi_rel_entr = rel_entropy(
            smooth_pdf(psi_trj_pdfs, pseudocount),
            smooth_pdf(psi_ref_trj_pdfs, pseudocount),
        )

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

                # Compute the 2D histogram for the joint phi/psi angles. The
                # edges are passed once per axis: a bare length-2 array (one
                # bin) would otherwise be read as [nx, ny] bin *counts*
                hist, x_edges, y_edges = np.histogram2d(
                    angles[0], angles[1], bins=[bins, bins], density=True
                )

                # Multiply each density by the area of its own bin to obtain
                # the per-bin probability (prior to 2.0.6 the first bin's area
                # was used for every bin, so non-uniform bins did not sum to 1)
                pdf = hist * np.outer(np.diff(x_edges), np.diff(y_edges))

                traj_histograms.append(pdf)

            pdfs.append(traj_histograms)

        return np.array(pdfs)

    def compute_pdf(self, arr: np.ndarray, bins: np.ndarray) -> np.ndarray:
        """Per-residue 1D probability density histograms.

        Operates on either a 2D ``(n_residues, n_frames)`` array or a 3D
        ``(n_trajectories, n_residues, n_frames)`` stack. Each
        residue-level histogram is normalised by ``np.histogram(...,
        density=True)`` and then multiplied by the width of each bin in
        ``bins`` to convert density to probability mass, so each row sums
        to 1 whatever bin edges are passed.

        Parameters
        ----------
        arr : np.ndarray
            2D or 3D angle array. The last axis is the frame axis.
        bins : np.ndarray
            1D array of bin edges (in the same units as ``arr``). Need not
            match ``self.bins``, but note that the precomputed EV reference
            angles are drawn uniformly *within* ``self.bins`` bins, so a
            finer histogram of ``ref_phi_angles`` / ``ref_psi_angles`` shows
            a staircase rather than extra detail.

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

        # KEY POINT: multiplying by the per-bin width converts probability
        # *density* to probability *mass*. This must use the bins actually
        # passed in, not self.bwidth, otherwise rows only sum to 1 when the
        # two happen to agree.
        bins = np.asarray(bins, dtype=float)
        bin_widths = np.diff(bins)

        # if (traj x n_res x frames), histogram axis (2) associated all the frames
        if arr.ndim == 3:
            pdf = (
                np.apply_along_axis(
                    lambda col: np.histogram(col, bins=bins, density=True)[0],
                    axis=2,
                    arr=arr,
                )
                * bin_widths
            )
        # else (n_res x n_frames), histogram axis (1) associated with frames
        else:
            pdf = (
                np.apply_along_axis(
                    lambda col: np.histogram(col, bins=bins, density=True)[0],
                    axis=1,
                    arr=arr,
                )
                * bin_widths
            )

        return pdf

    def get_all_to_all_2d_trj_comparison(
        self, metric: str = "hellingers", recompute: bool = False
    ) -> np.ndarray:
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
            If True, ignore the cached joint PDFs from :meth:`trj_pdfs`
            and rebuild them. Default False.

        Returns
        -------
        np.ndarray
            Shape ``(n_combinations, n_residues)`` of pairwise per-residue
            Hellinger distances in ``[0, 1]``. Rows follow
            ``itertools.combinations`` order, i.e. ``(0, 1), (0, 2), ...,
            (1, 2), ...``; column ``k`` is residue
            ``self.residue_indices[k]``.

        Raises
        ------
        SSException
            If ``metric`` is anything other than ``'hellingers'``.

        Example
        -------
        >>> mat = sq.get_all_to_all_2d_trj_comparison()
        >>> mat.shape   # 3 trajs -> C(3,2) == 3 pairs; NTL9 has 54 phi/psi residues
        (3, 54)
        """
        if metric != "hellingers":
            raise SSException(
                f"metric '{metric}' is not implemented for the 2D comparison; "
                "only 'hellingers' is available."
            )

        # shape = replicas, angles, phi_bins, psi_bins
        pdfs = self.trj_pdfs(dihedral="joint", recompute=recompute)

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
        self,
        metric: str = "hellingers",
        recompute: bool = False,
        pseudocount: float = 1e-6,
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """All-vs-all per-residue dihedral comparisons (separate phi and psi).

        Builds the per-residue 1D phi and psi PDFs for every trajectory,
        enumerates every pairwise trajectory combination, and computes a
        per-residue Hellinger distance or relative entropy for each pair.
        The two dihedrals are kept separate (unlike
        :meth:`get_all_to_all_2d_trj_comparison`).

        For ``metric='relative entropy'`` the PDFs are smoothed with
        :func:`smooth_pdf` first, since KL divergence is infinite wherever
        the second distribution has an empty bin that the first does not.

        Parameters
        ----------
        metric : {'hellingers', 'relative entropy'}, optional
            Which divergence to compute. Default ``'hellingers'``.
        recompute : bool, optional
            If True, ignore cached PDFs from ``self.trj_pdfs`` and rebuild
            them. Default False.
        pseudocount : float, optional
            Mass added to every bin before renormalising, used only for
            ``'relative entropy'``. ``0`` gives the raw (possibly infinite)
            divergence. Default ``1e-6``.

        Returns
        -------
        tuple of (pd.DataFrame, pd.DataFrame)
            ``(phi_df, psi_df)`` each of shape
            ``(n_combinations, n_residues)`` containing the chosen metric
            for every pairwise comparison. Rows follow
            ``itertools.combinations`` order, i.e. ``(0, 1), (0, 2), ...,
            (1, 2), ...``, and for ``'relative entropy'`` row ``(i, j)`` is
            :math:`D_{KL}(P_i \\| P_j)` with ``i < j``. The DataFrame columns
            are 0-based column positions, not residue indices: column ``k``
            is residue ``self.residue_indices[k]``.

        Raises
        ------
        SSException
            If ``metric`` is not one of the two supported strings.

        Example
        -------
        >>> phi_df, psi_df = sq.get_all_to_all_trj_comparisons()
        """
        if metric not in ("hellingers", "relative entropy"):
            raise SSException(
                f"The metric '{metric}' is not implemented. Use 'hellingers' or "
                "'relative entropy'."
            )

        phi_pdfs = self.trj_pdfs(recompute=recompute, dihedral="trj_phi_pdfs")
        psi_pdfs = self.trj_pdfs(recompute=recompute, dihedral="trj_psi_pdfs")

        if metric == "relative entropy":
            phi_pdfs = smooth_pdf(phi_pdfs, pseudocount)
            psi_pdfs = smooth_pdf(psi_pdfs, pseudocount)

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
        else:
            if phi_combinations.shape[0] == 1 and psi_combinations.shape[0] == 1:
                phi_metric = rel_entropy(phi_combinations[0], phi_combinations[0])
                psi_metric = rel_entropy(psi_combinations[0], psi_combinations[0])
            else:
                phi_metric = rel_entropy(phi_combinations[0], phi_combinations[1])
                psi_metric = rel_entropy(psi_combinations[0], psi_combinations[1])

        return pd.DataFrame(phi_metric), pd.DataFrame(psi_metric)

    def get_degree_bins(self) -> np.ndarray:
        """Histogram bin edges spanning ``[-180, 180]`` degrees.

        Constructs the bin edges used by every histogram-based method on
        this class. The number of bins is ``2*pi / self.bwidth`` (which
        ``__validate_arguments`` guarantees is a whole number) and the
        edges are laid out with ``np.linspace`` so the first and last
        edges land exactly on -180 and 180 with no rounding of the width.

        Returns
        -------
        np.ndarray
            1D array of ``n_bins + 1`` bin edges in degrees, monotonically
            increasing, starting at -180 and ending at 180.

        Example
        -------
        >>> sq.get_degree_bins()       # bwidth = 15 degrees
        array([-180., -165., ...,  165.,  180.])
        """
        n_bins = int(np.round(2 * np.pi / self.bwidth))
        return np.linspace(-180.0, 180.0, n_bins + 1)

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

        * **A** - per-residue Hellinger distance vs. the reference (the
          supplied ``reference_list``, or the precomputed excluded-volume
          limit) with per-trajectory points and across-trajectory mean.
        * **B** - per-residue all-vs-all trajectory Hellinger distances.
        * **C** - per-residue spread (max minus min) of the panel A
          Hellinger distances across trajectories.
        * **D** - per-residue fractional helicity for each trajectory with
          the across-trajectory mean; when a ``reference_list`` was
          supplied the mean reference helicity is drawn as a dashed line.

        Layout is mosaic ``"AABB;CCDD"``. The chosen ``dihedral`` selector
        controls which of phi / psi / joint 2D is shown. All four panels
        share one x-axis: the 1-based position of each residue among the
        CA-bearing residues of the chain (caps are not counted). The
        dihedral panels cover only the residues that carry both a phi and a
        psi (``self.residue_indices``), so on an uncapped chain they run
        from position 2 to ``n - 1`` while the helicity panel runs from 1 to
        ``n``; on a capped chain all four cover positions 1 to ``n``.

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
            If given, write the figure to
            ``<save_dir>/<dihedral>_<figname>`` (e.g.
            ``figs/2D_hellingers.pdf``). Default None (no file written;
            figure is returned only).
        dihedral : {'2D', 'phi', 'psi'}, optional
            Which dihedral comparison to plot. ``'2D'`` requires
            ``method='2D angle distributions'``. Default ``'2D'``.
        figname : str, optional
            File name including extension, prefixed with ``dihedral`` and
            joined with ``save_dir``. Default ``'hellingers.pdf'``.

        Returns
        -------
        tuple
            ``(fig, axd)`` — the matplotlib figure and the mosaic Axes
            dictionary keyed by ``'A'``, ``'B'``, ``'C'``, ``'D'``.

        Raises
        ------
        SSException
            If the requested ``dihedral`` is not compatible with the
            instance's ``method`` (``'2D'`` requires
            ``'2D angle distributions'``; ``'phi'`` / ``'psi'`` require
            ``'1D angle distributions'``).

        Example
        -------
        >>> fig, axd = sq.quality_plot(save_dir='./figs')   # default 2D method
        >>> # 'phi' / 'psi' panels need method='1D angle distributions'
        >>> fig, axd = sq_1d.quality_plot(dihedral='phi')
        """

        fig, axd = plt.subplot_mosaic(
            """AABB;CCDD""",
            sharex=True,
            figsize=figsize,
            dpi=dpi,
            facecolor="w",
            gridspec_kw={"height_ratios": [2, 2]},
        )
        # Only the combinations below are meaningful: the 2D method yields a
        # single joint (phi, psi) Hellinger distance per residue, while the 1D
        # method yields separate phi and psi distances. Resolve the branch
        # before computing anything so an unsupported combination fails with
        # a clear message rather than a downstream shape error.
        if self.method == "2D angle distributions":
            if dihedral != "2D":
                plt.close(fig)
                raise SSException(
                    f"dihedral='{dihedral}' is not available with "
                    "method='2D angle distributions'; use dihedral='2D', or "
                    "construct SamplingQuality with "
                    "method='1D angle distributions' for phi/psi."
                )
            metric = self.hellingers_distances()
            all_to_all = self.get_all_to_all_2d_trj_comparison()
        elif self.method == "1D angle distributions":
            if dihedral not in ("phi", "psi"):
                plt.close(fig)
                raise SSException(
                    f"dihedral='{dihedral}' is not available with "
                    "method='1D angle distributions'; use dihedral='phi' "
                    "or dihedral='psi'."
                )
            index = 0 if dihedral == "phi" else 1
            metric = self.hellingers_distances()[index]
            all_to_all = self.get_all_to_all_trj_comparisons()[index]
        else:
            plt.close(fig)
            raise SSException(f"Unsupported method '{self.method}'")

        metric = np.asarray(metric)

        trj_helicity, ref_helicity = self.fractional_helicity()
        trj_helicity = np.asarray(trj_helicity)
        ref_helicity = np.asarray(ref_helicity)

        # Every panel shares one residue axis: the 1-based position of a
        # residue among the CA-bearing residues of the chain (caps are not
        # counted). The dihedral panels only cover residues that carry both a
        # phi and a psi (self.residue_indices), so on an uncapped chain they
        # start at position 2 and stop at n-1, while the helicity panel covers
        # every CA-bearing residue. Previously both were plotted at 1..len(),
        # which put residue k of the dihedral panels one position to the left
        # of residue k in the helicity panel on uncapped chains.
        protein = self.trajs[0].proteinTrajectoryList[self.proteinID]
        ca_resids = list(protein.resid_with_CA)
        n_pos = len(ca_resids)
        idx = np.array([ca_resids.index(r) + 1 for r in self.residue_indices])
        helix_idx = np.arange(1, trj_helicity.shape[-1] + 1)
        xticks = np.arange(increment, n_pos + 1, increment)
        xticklabels = np.arange(increment, n_pos + 1, increment)

        if self.reference_list:
            panel_a_title = "Comparison to reference ensemble"
        else:
            panel_a_title = "Comparison to the Excluded Volume Limit"

        yticks = [0, 0.2, 0.4, 0.6, 0.8, 1]
        ytick_labels = [0, 0.2, 0.4, 0.6, 0.8, 1]
        for ax in axd:
            if ax == "A":
                axd[ax].set_yticks(yticks)
                axd[ax].set_yticklabels(ytick_labels, fontsize=fontsize)
                axd[ax].set_ylim([0, 1])
                axd[ax].set_ylabel("Hellinger's Distance", fontsize=fontsize)
                axd[ax].set_title(panel_a_title, fontsize=fontsize)

                axd[ax].set_xticks(
                    xticks,
                )
                axd[ax].set_xticklabels(xticklabels, fontsize=fontsize)
                axd[ax].set_xlim([0, n_pos + 1])

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
                axd[ax].set_xlim([0, n_pos + 1])

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
                axd[ax].set_xlim([0, n_pos + 1])

                # plot red
                axd[ax].plot(
                    helix_idx, trj_helicity.transpose(), ".r", ms=4, alpha=0.3, mew=0
                )

                # plot line avg helicity
                axd[ax].plot(
                    helix_idx,
                    np.mean(trj_helicity, axis=0),
                    "sk-",
                    ms=2,
                    alpha=1,
                    mew=0,
                    linewidth=0.5,
                    label="trajectories",
                )

                # reference helicity is only meaningful when real reference
                # trajectories were supplied (the EV tables carry no DSSP).
                if self.reference_list:
                    axd[ax].plot(
                        helix_idx,
                        np.mean(ref_helicity, axis=0),
                        "--",
                        color="tab:blue",
                        linewidth=0.8,
                        label="reference",
                    )
                    axd[ax].legend(fontsize=fontsize - 2, frameon=False)

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
            outpath = os.path.join(save_dir, f"{dihedral}_{figname}")
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
        SSException
            If ``dihedral`` is not one of the three allowed strings.

        Example
        -------
        >>> phi_pdfs = sq.trj_pdfs(dihedral='trj_phi_pdfs')
        """
        selectors = ["trj_phi_pdfs", "trj_psi_pdfs", "joint"]
        if dihedral not in selectors:
            raise SSException(
                f"dihedral '{dihedral}' is not a valid selector. "
                "Use one of 'trj_phi_pdfs', 'trj_psi_pdfs' or 'joint'."
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
        SSException
            If ``dihedral`` is not one of the three allowed strings.

        Example
        -------
        >>> ref_phi = sq.ref_pdfs(dihedral='ref_phi_pdfs')
        """
        selectors = ["ref_phi_pdfs", "ref_psi_pdfs", "joint"]

        if dihedral not in selectors:
            raise SSException(
                f"dihedral '{dihedral}' is not a valid selector. "
                "Use one of 'ref_phi_pdfs', 'ref_psi_pdfs' or 'joint'."
            )

        # The public selector is 'joint' (to mirror trj_pdfs) but the cache
        # key is 'ref_joint' so the reference and trajectory joint PDFs never
        # overwrite one another.
        cache_keys = {
            "ref_phi_pdfs": "ref_phi_pdfs",
            "ref_psi_pdfs": "ref_psi_pdfs",
            "joint": "ref_joint",
        }

        for selector in selectors:
            key = cache_keys[selector]
            if key not in self.__precomputed or recompute is True:
                if selector == "ref_phi_pdfs":
                    self.__precomputed[key] = self.compute_pdf(
                        self.ref_phi_angles, bins=self.bins
                    )
                elif selector == "ref_psi_pdfs":
                    self.__precomputed[key] = self.compute_pdf(
                        self.ref_psi_angles, bins=self.bins
                    )
                elif selector == "joint":
                    data = np.array([self.ref_phi_angles, self.ref_psi_angles])
                    pdfs = self.compute_series_of_histograms_along_axis(
                        data, bins=self.bins, axis=2
                    )
                    self.__precomputed[key] = pdfs

        return self.__precomputed[cache_keys[dihedral]]

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
        # keyed on method so changing self.method after construction does not
        # hand back an array of the other method's shape
        selector = ("hellingers", self.method)

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
        return self.compute_frac_helicity(proteinID=self.proteinID, recompute=recompute)


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

    Note that the tables hold the phi and psi *marginals* and the two are
    resampled independently, so the synthetic reference has no phi/psi
    correlation: its joint distribution is the product of its marginals.
    Under the 2D method any phi/psi correlation present in the simulated
    ensemble therefore adds to the distance from this reference. Supply
    excluded-volume trajectories via ``reference_list`` if the joint
    reference distribution matters.

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
    seed : int or None, optional
        Seed for the ``numpy.random.default_rng`` generator used for the
        resampling. None (default) gives a fresh unseeded generator, so
        two instances built with ``seed=None`` will differ.

    Attributes
    ----------
    ref_phi_angles, ref_psi_angles : np.ndarray
        Inverse-CDF-sampled reference angles of shape
        ``(num_trajs, n_residues, nsamples)``.
    rng : np.random.Generator
        The generator used for sampling.

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

    def __init__(
        self,
        sequence: str,
        bins: np.ndarray,
        num_trajs: int,
        nsamples: int,
        seed: Union[int, None] = None,
    ):
        self.sequence = sequence
        self.num_trajs = num_trajs
        self.nsamples = nsamples
        self.bins = np.asarray(bins, dtype=float)
        self.seed = seed
        self.rng = np.random.default_rng(seed)

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
        angles on ``self.bins``, converts it to a CDF, and inverse-CDF
        samples ``self.nsamples`` synthetic angles per residue: each draw
        picks a bin in proportion to its mass and then a uniform position
        within that bin's edges, so the resampled histogram reproduces the
        reference histogram (including the bins at -180 and 180) up to
        sampling noise. Draws come from ``self.rng``. Used at construction
        time to populate :attr:`ref_phi_angles` and :attr:`ref_psi_angles`.

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
        # only phi and psi have EV tables (anything else used to raise a
        # bare KeyError); only the requested angle is gathered
        ssutils.validate_keyword_option(angle, ["phi", "psi"], "angle")
        if angle == "phi":
            reference = self.gather_phi_reference_dihedrals(self.sequence)
        else:
            reference = self.gather_psi_reference_dihedrals(self.sequence)

        dihedral_hist = []

        for dihedral in range(reference.shape[0]):
            dihedral_angles = reference[dihedral, :]

            # GOAL: Generate samples that adhere to the underlying distribution
            # Step 1: Compute the histogram and its CDF on self.bins
            hist, bin_edges = np.histogram(
                dihedral_angles, bins=self.bins, density=True
            )

            cdf = np.cumsum(hist * np.diff(bin_edges))
            cdf /= cdf[-1]  # Normalize the CDF

            # Step 2: Inverse-CDF sample. Pick the bin whose cumulative mass
            # first exceeds the uniform draw, then place the sample uniformly
            # inside that bin. Interpolating between bin centres (the old
            # approach) biased the tails and could never reach the outermost
            # half-bins at -180 and 180.
            rand_values = self.rng.random(size=self.nsamples)
            bin_idx = np.searchsorted(cdf, rand_values, side="right")
            bin_idx = np.minimum(bin_idx, len(hist) - 1)

            lower = bin_edges[bin_idx]
            upper = bin_edges[bin_idx + 1]
            sampled_dihedrals = lower + self.rng.random(size=self.nsamples) * (
                upper - lower
            )

            dihedral_hist.append(sampled_dihedrals)

        return np.array(dihedral_hist)

    @staticmethod
    def __three_letter(residue: str, position: int) -> str:
        """Map a one-letter code to three letters, raising a clear error.

        The EV tables only cover the 20 standard amino acids, so a
        non-standard letter in the sequence (e.g. a modified residue that
        survived caps stripping) is reported by name and position rather
        than surfacing as a bare KeyError.
        """
        # ONE_TO_THREE also maps the cap symbols (<, >, (, )), but the EV
        # tables only cover the 20 standard residues, so check against those
        # (prior to 2.0.6 an NH2/FOR cap slipped through here and the
        # defaultdict tables raised a bare KeyError after gaining an empty
        # entry for the cap)
        three = ONE_TO_THREE.get(residue)
        if three is not None and three in EV_RESIDUE_MAPPER:
            return three
        else:
            raise SSException(
                f"Residue '{residue}' at position {position} is not one of the "
                "20 standard amino acids, so no precomputed excluded-volume "
                "reference is available for it. Supply your own reference "
                "trajectories via reference_list instead."
            ) from None

    def gather_phi_reference_dihedrals(self, sequence: str) -> np.ndarray:
        """Look up the excluded-volume reference phi distribution for each residue.

        For each position ``i`` the phi distribution is that of residue
        ``i`` itself, in the context of residue ``i-1`` (the residue
        preceding the rotatable phi bond). The tables in
        :data:`PHI_EV_ANGLES_DICT` are keyed first by the three-letter
        code of residue ``i`` and second by the coarse-grained identity
        of the preceding residue (via :data:`EV_RESIDUE_MAPPER`, which
        collapses the 20 amino acids to ALA/LEU/PRO-like contexts). For
        position 0 we substitute alanine as the preceding context.

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
        (4, 4000)
        """

        phi_angles = []
        for i, residue in enumerate(sequence):
            if i == 0:
                phi_preceeding_context = "A"
            else:
                phi_preceeding_context = sequence[i - 1]

            # The EV tables are keyed first by the residue whose phi we want
            # and second by the (coarse-grained) identity of the preceding
            # residue, i.e. PHI_EV_ANGLES_DICT[res_i][mapper(res_{i-1})].
            three_letter_residue = self.__three_letter(residue, i)
            approximate_context = EV_RESIDUE_MAPPER[
                self.__three_letter(phi_preceeding_context, i - 1)
            ]

            phi_angles.append(
                PHI_EV_ANGLES_DICT[three_letter_residue][approximate_context]
            )

        return np.array(phi_angles)

    def gather_psi_reference_dihedrals(self, sequence: str) -> np.ndarray:
        """Look up the excluded-volume reference psi distribution for each residue.

        Mirror of :meth:`gather_phi_reference_dihedrals` for psi: the
        distribution is that of residue ``i`` in the context of the
        *following* residue ``i+1`` (because psi rotates the i-(i+1)
        bond), looked up as
        ``PSI_EV_ANGLES_DICT[res_i][EV_RESIDUE_MAPPER[res_{i+1}]]``. For
        the last position we substitute alanine as the following context.

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
        (4, 4000)
        """

        psi_angles = []
        for i, residue in enumerate(sequence):
            if i == len(sequence) - 1:
                psi_subsequent_context = "A"
            else:
                psi_subsequent_context = sequence[i + 1]

            # As for phi: PSI_EV_ANGLES_DICT[res_i][mapper(res_{i+1})].
            three_letter_residue = self.__three_letter(residue, i)
            approximate_context = EV_RESIDUE_MAPPER[
                self.__three_letter(psi_subsequent_context, i + 1)
            ]
            psi_angles.append(
                PSI_EV_ANGLES_DICT[three_letter_residue][approximate_context]
            )

        return np.array(psi_angles)
