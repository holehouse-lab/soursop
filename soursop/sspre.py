##     _____  ____  _    _ _____   _____  ____  _____
##   / ____|/ __ \| |  | |  __ \ / ____|/ __ \|  __ \
##  | (___ | |  | | |  | | |__) | (___ | |  | | |__) |
##   \___ \| |  | | |  | |  _  / \___ \| |  | |  ___/
##   ____) | |__| | |__| | | \ \ ____) | |__| | |
##  |_____/ \____/ \____/|_|  \_\_____/ \____/|_|

## Alex Holehouse (Pappu Lab and Holehouse Lab) and Jared Lalmansing (Pappu lab)
## Simulation analysis package
## Copyright 2014 - 2026
##

import numpy as np
from .ssexceptions import SSWarning, SSException
from .ssprotein import SSProtein

### SSPRE contains all the functionality associated with calculating
### PRE profiles.
###
###
###
###


original_K = 1.2300e-32  # K constant in cm6*s-2
K_IN_NM6 = original_K * 1e42  # K constant in nm6 s-2
# W_H        = 267530000        # Proton Larmor frequency
# W_H_SQUARED = W_H * W_H

# Default coarse-grained spin-label cloud geometry.
#
# These defaults were measured directly from the MTSSL 175K X1X2 rotamer
# library shipped with DEER-PREdict (the library used in the HyRes paper,
# bioRxiv 2026.06.23.734133). Over the 46 rotamers, the (Boltzmann-weighted)
# distance from CB to the nitroxide electron (N-O midpoint) is 7.0 +/- 1.0 A
# and the angle between the CA->CB direction and the CB->electron vector is
# 84 +/- 25 deg (range 34-134 deg). A single CG bead placed ~7 A from CB in a
# wide cone that excludes the backward (toward-backbone) direction therefore
# reproduces the geometry of the paramagnetic centre without the cost or the
# all-atom requirement of an explicit rotamer library.
DEFAULT_LABEL_DISTANCE = 7.0  # Angstroms, CB -> effective electron position
DEFAULT_LABEL_CONE_ANGLE = 130.0  # degrees, half-angle of the cloud about CA->CB
DEFAULT_N_LABEL_CONFORMERS = 100
# Effective steric radius of the coarse-grained label bead, in Angstroms.
# Steric exclusion is ALWAYS applied when use_label is True: for every frame
# any bead whose centre lies within this radius of a non-label CA is treated
# as clashing with the chain and dropped. This mimics the steric pruning
# DEER-PREdict applies to its rotamer library and stops the cloud from being
# placed through the protein interior (critical for folded/compact chains).
# The radius is not a CA-CA contact distance: it represents the exclusion
# between the bulky nitroxide bead and the volume occupied by each residue
# (whose side-chain centroid sits beyond CA). The 6.0 A default is the value
# that minimised the mean RMSD of the CG cloud against DEER-PREdict's explicit
# MTSL rotamer PRE profiles in a joint distance x radius grid search over
# folded (NTL9) and disordered (LAF-1) all-atom benchmarks with multiple
# label sites each. Folded chains favour a slightly larger radius and IDPs a
# slightly smaller one; 6.0 A is the balanced optimum, paired with the 7.0 A
# label_distance (itself the measured mean CB->electron distance). The 5.5 A
# default below is the value that, with the default soft steric wall, minimised
# the dataset-balanced RMSD against DEER-PREdict across 85 disordered
# (trajectory) and ~2700 folded (AlphaFold) label sites (see the soft-vs-hard
# calibration). The soft wall is flat enough in this region that 5.5-6.0 A are
# essentially equivalent.
DEFAULT_LABEL_BEAD_RADIUS = 5.5

# How the label bead's steric exclusion against the chain is applied:
#   'hard' - a bead is kept (weight 1) only if it clears every non-label CA by
#            label_bead_radius, otherwise dropped (weight 0).
#   'soft' - every bead gets a smooth Boltzmann-style weight from a WCA-like
#            repulsive wall against the CA atoms,
#                w_i = exp( - sum_j (label_bead_radius / d_ij)^wall_stiffness ),
#            so beads grazing the chain are down-weighted rather than discarded.
#            Reduces to 'hard' as wall_stiffness -> infinity. In a calibration
#            against DEER-PREdict the soft wall is less sensitive to the radius
#            and gives a single value that works across folded and disordered
#            chains, whereas the hard cutoff needs a regime-dependent radius.
# A full-scale calibration against DEER-PREdict (85 disordered + ~2700 folded
# label sites) found the soft wall to be marginally more accurate than the hard
# cutoff on a dataset-balanced basis and, more importantly, about three times
# less sensitive to label_bead_radius - a single value works across folded and
# disordered chains. 'soft' with a stiffness of 6.0 is therefore the default.
DEFAULT_LABEL_STERIC = "soft"
DEFAULT_LABEL_WALL_STIFFNESS = 6.0


def _fibonacci_cap_directions(n, half_angle_deg):
    """Deterministic, near-uniform unit vectors on a spherical cap about +z.

    Generates ``n`` unit vectors whose polar angle (from the +z axis) is
    distributed uniformly in *solid angle* between 0 and ``half_angle_deg``,
    using a Fibonacci lattice for even azimuthal coverage. With
    ``half_angle_deg = 180`` this returns points spread over the full sphere.

    Parameters
    ----------
    n : int
        Number of directions (conformers in the label cloud).

    half_angle_deg : float
        Half-angle of the cap, in degrees. The cloud is restricted to
        directions within this angle of +z.

    Returns
    -------
    numpy.ndarray
        Array of shape ``(n, 3)`` of unit vectors.
    """

    n = int(n)
    cos_max = np.cos(np.radians(half_angle_deg))
    golden = np.pi * (3.0 - np.sqrt(5.0))  # golden angle
    i = np.arange(n) + 0.5
    # z uniform in [cos(half_angle), 1] -> uniform in solid angle over the cap
    z = 1.0 - (i / n) * (1.0 - cos_max)
    r = np.sqrt(np.clip(1.0 - z * z, 0.0, None))
    phi = i * golden
    return np.column_stack((r * np.cos(phi), r * np.sin(phi), z))


def _rotation_from_z(axis):
    """Rotation matrix mapping the +z axis onto a target unit vector.

    Parameters
    ----------
    axis : numpy.ndarray
        Length-3 unit vector (the desired image of +z).

    Returns
    -------
    numpy.ndarray
        A ``(3, 3)`` rotation matrix ``R`` such that ``R @ [0, 0, 1] == axis``.
    """

    z = np.array([0.0, 0.0, 1.0])
    axis = axis / np.linalg.norm(axis)
    v = np.cross(z, axis)
    c = float(np.dot(z, axis))
    s = np.linalg.norm(v)
    if s < 1e-8:
        # axis is parallel (identity) or antiparallel (flip about x) to +z
        return np.eye(3) if c > 0 else np.diag([1.0, -1.0, -1.0])
    vx = np.array([[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]])
    return np.eye(3) + vx + vx @ vx * ((1.0 - c) / (s * s))


class SSPRE:
    """Synthetic paramagnetic relaxation enhancement (PRE) calculations.

    Wraps an ``SSProtein`` ensemble and computes intramolecular PRE
    intensity ratios and the associated transverse relaxation
    enhancement (gamma) profiles for a nitroxide spin label placed at an
    arbitrary sequence position. The model follows the approach used by
    Meng & Lyle (PNAS 2013), Das et al. (PNAS 2016) and Peran &
    Holehouse (PNAS 2019), and is fast enough to apply to large
    (multi-thousand frame) disordered-protein ensembles.

    A typical workflow is to extract an ``SSProtein`` from an
    ``SSTrajectory``, construct an ``SSPRE`` object with the relevant
    experimental parameters (effective correlation time, INEPT delay,
    diamagnetic transverse relaxation rate and proton Larmor frequency)
    and then call ``generate_PRE_profile`` for one or more label
    positions. The bound ``SSProtein`` is treated as read-only and is not
    modified.

    Example
    -------
    >>> P = traj.proteinTrajectoryList[0]
    >>> spre = SSPRE(P, tau_c=5.0, t_delay=10.0, R_2D=10.0, W_H=600000000)
    >>> profile, gamma = spre.generate_PRE_profile(label_position=42)
    """

    # ........................................................................
    #
    def __init__(self, SSProteinObject, tau_c, t_delay, R_2D, W_H):
        """Build an SSPRE object bound to a single-chain protein ensemble.

        Stores the bound ``SSProtein`` and the experimental parameters,
        runs physiological sanity checks on the inputs (emitting an
        ``SSWarning`` for values far outside their normal range without
        blocking analysis) and precomputes the spectral-density prefactor
        used by ``generate_PRE_profile``. The resulting object can then be
        used to calculate PRE profiles from the underlying ensemble; the
        calculation is extremely fast.

        Parameters
        ----------
        SSProteinObject : soursop.ssprotein.SSProtein
            SSProtein object extracted from an ``SSTrajectory``. This is
            the main object over which protein-centric analysis is
            performed in SOURSOP. It is treated as read-only.

        tau_c : float
            Effective rotational correlation time in nanoseconds,
            typically between 1 and 30.

        t_delay : float
            Total duration of the INEPT delays in the PRE experiment, in
            milliseconds. Depends on the pulse sequence but is typically
            ~1-30 ms for an HSQC.

        R_2D : float
            Transverse relaxation rate of the backbone amide protons in
            the diamagnetic form of the protein, in Hz (per second). A
            value of around 10 is typical.

        W_H : float
            Proton Larmor frequency of the magnet, in Hz - i.e. the
            nominal "MHz" rating of the magnet expressed in Hz (a 600 MHz
            magnet uses ``600000000``). The proton Larmor frequency at
            1 Tesla is 267530000 per second per Tesla.

        Raises
        ------
        soursop.ssexceptions.SSException
            If ``SSProteinObject`` is not an ``SSProtein`` instance.

        Example
        -------
        >>> P = traj.proteinTrajectoryList[0]
        >>> spre = SSPRE(P, tau_c=5.0, t_delay=10.0, R_2D=10.0, W_H=600000000)
        """

        # NOTE: NO CHANGES SHOULD BE MADE TO THE SSPO by the object - this should be treated
        # as a read only object (no such explicit control in Python)
        self.SSPO = SSProteinObject

        if type(self.SSPO) is not SSProtein:
            raise SSException(
                f"SSPRE requires an SSProtein object to be passed, but instead the first argument is of type {type(self.SSPO)}"
            )

        # set the INEPT delay value and the backbone amide transverse relaxation rate which
        # is used explicitly later
        self.t_delay = float(t_delay)  # in ms - INDEPT delay
        self.R_2D = float(R_2D)  # in Hz - backbone amide transverse relaxation rate
        self.tau_c = tau_c  # in ns - effecive correlation time
        self.W_H = W_H  # in Hz - Proton Larmor frequency in the magnet

        # ------------------------------------------
        # sanity checks to warn if any of the input values seem dratsically wrong. NOTE that these won't block the analysis but will
        # throw up errors
        if self.R_2D < 0.01 or self.R_2D > 100:
            SSWarning(
                "WARNING: The value of R_2D (bacbone amide transverese relaxation rate) is far from the normal expected value of ~10 (R_2D = %4.4e) - recal this value is in units of Herz"
                % (self.R_2D)
            )

        if self.t_delay < 0.01 or self.t_delay > 100:
            SSWarning(
                "WARNING: The value of t_delay (INEPT delay) is far from the normal expected value of ~15 (t_delay = %4.4e) - recal this is in units of ms"
                % (self.t_delay)
            )

        if tau_c < 0.01 or tau_c > 100:
            SSWarning(
                "WARNING: The value of tau_c (effective correlation time) is far from the normal expected value of ~5 (t_delay = %4.4e) - recal this is in units of ns"
                % (self.tau_c)
            )

        # if Larmor frequency less than 100 MhZ or above 2 GHz assume something is wrong
        if W_H < 50000000 or W_H > 2000000000:
            SSWarning(
                f"WARNING: The value of W_h {self.W_H} (proton Larmor frequency) is far from the normal expected value of ~600 000 000 - recal this value should be provided in Herz"
            )

        # # convert tau_c to seconds and calculate tau_c squared
        tau_c = float(tau_c) / 1000000000  # tau c in seconds
        tau_c_squared = tau_c * tau_c  #

        # compute the prefactor term which will be used when computing the PRE dependent relaxation profile
        # by the generate_PRE_profile function
        W_H_SQUARED = W_H * W_H
        PREFACTOR = (3 * tau_c) / (1 + W_H_SQUARED * tau_c_squared)
        PREFACTOR = 4 * tau_c + PREFACTOR
        self.PREFACTOR = PREFACTOR * K_IN_NM6

    # ........................................................................
    #
    def __repr__(self):
        """Return a concise string summary of the SSPRE object.

        The representation includes the object's memory address and the
        four experimental parameters (diamagnetic transverse relaxation
        rate, INEPT delay, correlation time and proton Larmor frequency).

        Returns
        -------
        str
            Human-readable single-line description of this SSPRE object.
        """
        return (
            "["
            + hex(id(self))
            + "]: SSPRE OBJ - (R_2D = %3.2f Hz, t_delay = %3.2f ms, tau_c = %3.2f ns, H1 Larmor = %3.3e Hz)"
            % (self.R_2D, self.t_delay, self.tau_c, self.W_H)
        )

    # ........................................................................
    #
    def generate_PRE_profile(
        self,
        label_position,
        spin_label_atom="CB",
        target_relaxation_atom="N",
        use_label=True,
        label_distance=DEFAULT_LABEL_DISTANCE,
        label_cone_angle=DEFAULT_LABEL_CONE_ANGLE,
        n_label_conformers=DEFAULT_N_LABEL_CONFORMERS,
        label_bead_radius=DEFAULT_LABEL_BEAD_RADIUS,
        label_steric=DEFAULT_LABEL_STERIC,
        label_wall_stiffness=DEFAULT_LABEL_WALL_STIFFNESS,
    ):
        """Compute the PRE intensity ratio and gamma profile for a spin label.

        Places a nitroxide spin label on the ``spin_label_atom`` of
        residue ``label_position`` and, for every residue that has a
        cached CA atom, computes the spin-label-induced transverse
        relaxation enhancement (the gamma profile) and the corresponding
        paramagnetic/diamagnetic intensity ratio (the PRE profile). By
        default the label is on CB and relaxation is assessed from the
        CB-to-backbone-N distances, matching the parameterisation used by
        Meng, Lyle, Luan, Raleigh & Pappu (PNAS 2013, 110:2123-2128),
        Das, Huang, Phillips, Kriwacki & Pappu (PNAS 2016, 113:5616-5621)
        and Peran, Holehouse, Carrico, Pappu, Bilsel & Raleigh (PNAS
        2019, 116:12301-12310).

        The PRE profile is the intensity ratio
        (I\\ :sub:`paramagnetic` / I\\ :sub:`diamagnetic`) and typically
        lies between 0 and 1: a ratio near 0 means the spin label
        dominates relaxation (the residue is, on average, close to the
        label) while a ratio near 1 means relaxation is dominated by
        non-spin-label mechanisms (the residue is far from the label).
        The gamma profile is the spin-label-induced amide-proton
        relaxation rate in per-second units, a second observable directly
        comparable with experiment.

        The r^-6 distance dependence is averaged *after* the per-frame
        relaxation is computed (not by averaging distances first), because
        the distance-to-relaxation mapping is non-linear. This routine is
        extremely fast (sub-10 s on a ~6000-frame ensemble).

        .. note::

            **Breaking change (SOURSOP 2.0.2).** ``use_label`` now defaults
            to ``True``, so by default the paramagnetic centre is modelled as
            a coarse-grained spin-label cloud calibrated against DEER-PREdict
            (see ``use_label`` below) rather than sitting directly on the
            ``spin_label_atom``. This produces different (and more accurate)
            profiles than SOURSOP <= 2.0.1. To recover the previous
            point-at-CB behaviour exactly, pass ``use_label=False``.

        Parameters
        ----------
        label_position : int
            Sequence position at which the spin label is located. Should
            contain a CB atom (i.e. not be glycine) unless
            ``spin_label_atom`` is changed to ``'CA'`` (see below).

        spin_label_atom : str, optional
            Name of the atom on which the spin label sits. Should be
            ``'CB'`` but may be changed if the labelled residue lacks a CB
            (e.g. a glycine in place of the cysteine nitroxide site).
            Default ``'CB'``.

        target_relaxation_atom : str, optional
            Name of the atom at which relaxation is evaluated. This should
            be ``'N'`` (backbone amide), as that is how the approach is
            parameterised; an atom of this name is searched for in every
            residue. Changing it is strongly discouraged. Default
            ``'N'``.

        use_label : bool, optional
            If ``True`` (the default as of SOURSOP 2.0.2), the paramagnetic
            centre is not taken to sit on the ``spin_label_atom`` itself but
            is instead modelled as a coarse-grained cloud of
            ``n_label_conformers`` beads placed a fixed ``label_distance``
            away from that atom, mimicking the geometry of an MTSL nitroxide
            side chain without requiring an all-atom rotamer library. For
            every frame the relaxation is averaged over the whole cloud (and,
            as always, over frames), so the r^-6 non-linearity is respected
            across both the conformer cloud and the ensemble. This works on
            both all-atom and coarse-grained (CA-only) trajectories, and the
            default ``label_*`` parameters below are calibrated against
            DEER-PREdict. Set ``use_label=False`` to recover the classic
            point-at-``spin_label_atom`` behaviour used by SOURSOP <= 2.0.1.

        label_distance : float, optional
            Distance, in Angstroms, from ``spin_label_atom`` to the cloud
            of label beads. The default of 7.0 A is the mean CB-to-electron
            distance of the MTSSL 175K X1X2 rotamer library. Only used when
            ``use_label`` is ``True``.

        label_cone_angle : float, optional
            Half-angle, in degrees, of the cone (about the CA->CB direction)
            within which the label cloud is distributed. A wide cone
            (default 130 deg) reproduces the broad lateral spread of the
            MTSL side chain while excluding directions that point back into
            the backbone. If the labelled residue has no CB distinct from
            the anchor atom (e.g. a CA-only coarse-grained bead, or the
            anchor *is* CA), the cloud falls back to a full isotropic sphere
            and this argument is ignored. Only used when ``use_label`` is
            ``True``.

        n_label_conformers : int, optional
            Number of beads in the label cloud. The cloud is generated
            deterministically (a Fibonacci lattice), so results are
            reproducible. Default 100. Only used when ``use_label`` is
            ``True``.

        label_bead_radius : float, optional
            Effective steric radius of the coarse-grained label bead, in
            Angstroms. Steric exclusion is *always* applied in the label
            model. For ``label_steric='hard'`` any bead whose centre lies
            within this radius of a CA atom of a *different* residue is
            dropped from that frame's average; for ``'soft'`` it is the
            distance scale of the repulsive wall (see ``label_steric``). This
            mimics the steric pruning that DEER-PREdict applies to its
            rotamer library and stops the cloud being placed through the
            protein interior (critical for folded/compact structures). The
            default was calibrated to most closely reproduce DEER-PREdict's
            explicit MTSL rotamer PRE profiles; it is not a bare CA-CA
            contact distance. Only used when ``use_label`` is ``True``.

        label_steric : {'hard', 'soft'}, optional
            How the bead cloud is sterically excluded from the chain.
            ``'hard'`` (default) keeps a bead only if it clears every
            non-label CA by ``label_bead_radius`` and averages uniformly over
            the survivors. ``'soft'`` instead weights every bead by a smooth
            WCA-like repulsive wall against the CA atoms,
            ``w_i = exp(-sum_j (label_bead_radius / d_ij)**label_wall_stiffness)``,
            and takes the weighted cloud average. The soft wall reduces to the
            hard cutoff as ``label_wall_stiffness`` grows, but is less
            sensitive to ``label_bead_radius`` and transfers better between
            folded and disordered chains. Only used when ``use_label`` is
            ``True``.

        label_wall_stiffness : float, optional
            Exponent of the soft repulsive wall (only used when
            ``label_steric='soft'``). Larger values give a sharper wall that
            approaches the hard cutoff; smaller values give a softer, longer
            ranged penalty. Default 6.0 (calibrated against DEER-PREdict).
            Only used when ``use_label`` is ``True``.

        Returns
        -------
        tuple of (list of float, list of float)
            A 2-tuple ``(profile, gamma)`` where ``profile`` is the PRE
            intensity ratio and ``gamma`` is the spin-label-induced amide
            proton relaxation rate (per second), each with one entry per
            CA-containing residue (in ascending residue order).

        Example
        -------
        >>> profile, gamma = spre.generate_PRE_profile(label_position=42)
        >>> len(profile) == len(gamma)
        True
        """

        # get index value of all residues
        # residue_list = self.SSPO.get_residue_index_list()

        # Residues carrying a single CA, in ascending order. Using the public
        # resid_with_CA property (rather than the private CA cache) means this
        # also works for coarse-grained one-bead-per-residue chains, where the
        # CA cache is populated by a different code path.
        residue_list = sorted(self.SSPO.resid_with_CA)

        # first calculate mean rij distance for pair-residue distances
        gamma = []

        if use_label:
            if label_steric not in ("hard", "soft"):
                raise SSException(
                    f"label_steric must be 'hard' or 'soft', got '{label_steric}'"
                )
            # Coarse-grained label-cloud model: instead of the point-at-CB
            # distance, the paramagnetic centre is a cloud of beads placed
            # label_distance away from the anchor atom, and relaxation is
            # averaged over the whole cloud as well as over frames.
            gamma = self.__generate_label_cloud_gamma(
                label_position,
                residue_list,
                spin_label_atom,
                target_relaxation_atom,
                label_distance,
                label_cone_angle,
                n_label_conformers,
                label_bead_radius,
                label_steric,
                label_wall_stiffness,
            )
        else:
            # finally for each residue calculate the r^6 distances associated with each frame and for EACH FRAME calculate
            # the PREFACTOR / r^6 value and then take the mean. THIS gives a different answer to if you take the mean distance
            # and calculate the PREFACTOR/<R^6> value because there is a non-linear mapping between relaxation and distance so
            # it's important the former method is used (i.e. only average at the end). This calculates the gamma coefficient for
            # each residue, which measures relaxation
            for idx in residue_list:
                r_6_nm = np.power(
                    0.1
                    * self.SSPO.get_inter_residue_atomic_distance(
                        label_position,
                        idx,
                        A1=spin_label_atom,
                        A2=target_relaxation_atom,
                    ),
                    6,
                )
                gamma.append(np.mean(self.PREFACTOR / r_6_nm))

        # convert the t_delay from ms to seconds
        t_delay_in_seconds = self.t_delay / 1000

        # for each gamma compute the intensity ration
        profile = []
        for g in gamma:
            profile.append(
                (self.R_2D * np.exp(-g * t_delay_in_seconds)) / (self.R_2D + g)
            )

        return (profile, gamma)

    # ........................................................................
    #
    def __get_atom_xyz_nm(self, resid, atom_name):
        """Return per-frame coordinates (nm) of a named atom in a residue.

        Parameters
        ----------
        resid : int
            Residue index to select from.

        atom_name : str
            Atom name to select (e.g. ``'CB'`` or ``'N'``).

        Returns
        -------
        numpy.ndarray or None
            Array of shape ``(n_frames, 3)`` in nanometres, or ``None`` if
            the residue has no atom of that name.
        """

        # Use SOURSOP's memoised O(1) residue/atom lookup rather than a fresh
        # topology.select() per call - the latter is very slow and, when the
        # label cloud fetches CA/target atoms for every residue and radius,
        # dominates the runtime on large proteins.
        sel = self.SSPO._SSProtein__residue_atom_lookup(resid, atom_name)
        if sel is None or len(sel) == 0:
            return None
        # mdtraj xyz is already in nm; take the first matching atom
        return self.SSPO.traj.xyz[:, sel[0], :]

    # ........................................................................
    #
    def __generate_label_cloud_gamma(
        self,
        label_position,
        residue_list,
        spin_label_atom,
        target_relaxation_atom,
        label_distance,
        label_cone_angle,
        n_label_conformers,
        label_bead_radius,
        label_steric,
        label_wall_stiffness,
    ):
        """Compute the per-residue gamma profile using the CG label cloud.

        For each frame a cloud of ``n_label_conformers`` beads is placed a
        distance ``label_distance`` from the anchor atom (``spin_label_atom``
        of ``label_position``). When a CB distinct from a CA is available the
        cloud is a cone of half-angle ``label_cone_angle`` about the CA->CB
        direction; otherwise it is a full isotropic sphere (so coarse-grained
        CA-only chains are supported). The spin-label-induced relaxation
        ``PREFACTOR / r^6`` is averaged over every (frame, cloud-bead, target)
        combination, preserving the r^-6 non-linearity across both the cloud
        and the ensemble.

        Parameters
        ----------
        label_position : int
            Residue index carrying the spin label.

        residue_list : list of int
            Sorted residue indices at which to evaluate relaxation.

        spin_label_atom : str
            Anchor atom name the cloud emanates from (``'CB'`` for all-atom,
            ``'CA'`` for a coarse-grained bead).

        target_relaxation_atom : str
            Atom name relaxation is evaluated at (``'N'``).

        label_distance : float
            Anchor-to-bead distance in Angstroms.

        label_cone_angle : float
            Half-angle in degrees of the cloud about the CA->CB axis.

        n_label_conformers : int
            Number of beads in the cloud.

        label_bead_radius : float
            Effective steric radius / wall distance scale of the label bead,
            in Angstroms.

        label_steric : {'hard', 'soft'}
            Hard cutoff (drop clashing beads) or soft WCA-like weighting of
            all beads. Steric exclusion is always applied either way.

        label_wall_stiffness : float
            Exponent of the soft repulsive wall (used only for ``'soft'``).

        Returns
        -------
        list of float
            The gamma coefficient (per second) for each residue in
            ``residue_list``.

        Raises
        ------
        soursop.ssexceptions.SSException
            If the anchor atom cannot be found on the labelled residue.
        """

        label_distance_nm = label_distance / 10.0

        # anchor atom the label cloud hangs off (CB for AA, CA for CG)
        anchor_xyz = self.__get_atom_xyz_nm(label_position, spin_label_atom)
        if anchor_xyz is None:
            raise SSException(
                f"Unable to find spin-label anchor atom [{spin_label_atom}] "
                f"in residue {label_position} for the label-cloud model."
            )

        # CA is used to orient the cloud away from the backbone. If the anchor
        # already is CA (coarse-grained bead) or no separate CB exists, we fall
        # back to an isotropic sphere.
        ca_xyz = self.__get_atom_xyz_nm(label_position, "CA")
        oriented = (
            spin_label_atom != "CA"
            and ca_xyz is not None
            and not np.allclose(anchor_xyz, ca_xyz)
        )

        if oriented:
            base_dirs = _fibonacci_cap_directions(n_label_conformers, label_cone_angle)
        else:
            base_dirs = _fibonacci_cap_directions(n_label_conformers, 180.0)

        # pre-fetch target atom coordinates (nm) for every residue, once
        target_xyz = {}
        for idx in residue_list:
            target_xyz[idx] = self.__get_atom_xyz_nm(idx, target_relaxation_atom)

        # Steric exclusion is always applied: gather CA coordinates of all
        # residues except the labelled one (stacked into (n_frames, n_ca, 3))
        # so beads that clash with the chain can be dropped frame by frame.
        bead_radius_nm = label_bead_radius / 10.0
        ca_stack = None
        ca_list = []
        for idx in residue_list:
            if idx == label_position:
                continue
            ca = self.__get_atom_xyz_nm(idx, "CA")
            if ca is not None:
                ca_list.append(ca)
        if len(ca_list) > 0:
            ca_stack = np.stack(ca_list, axis=1)  # (n_frames, n_ca, 3)

        n_frames = anchor_xyz.shape[0]

        # Accumulate a per-frame, within-cloud weighted average of
        # PREFACTOR / r^6 for each residue, then average uniformly over frames.
        # This matches the frame averaging used by SOURSOP's point model and by
        # DEER-PREdict (which normalises its rotamer weights within each frame
        # and then averages gamma over frames). For 'hard' steric the bead
        # weights are 0/1 (kept or dropped); for 'soft' they are the smooth WCA
        # wall weights. Frames in which the label is fully excluded contribute
        # nothing and are dropped from the frame count.
        relax_sum = {idx: 0.0 for idx in residue_list}
        frame_count = 0

        for f in range(n_frames):
            anchor = anchor_xyz[f]  # (3,)
            if oriented:
                axis = anchor - ca_xyz[f]
                norm = np.linalg.norm(axis)
                if norm < 1e-8:
                    dirs = base_dirs
                else:
                    dirs = base_dirs @ _rotation_from_z(axis / norm).T
            else:
                dirs = base_dirs

            # cloud bead positions for this frame: (n_conf, 3), in nm
            beads = anchor + label_distance_nm * dirs

            # per-bead steric weight against the non-label CA atoms
            if ca_stack is not None:
                # (n_conf, n_ca) distances from each bead to each CA this frame
                dca = np.linalg.norm(
                    beads[:, None, :] - ca_stack[f][None, :, :], axis=2
                )
                if label_steric == "soft":
                    # WCA-like repulsive wall; clip the exponent to avoid overflow
                    clash = np.sum(
                        np.power(bead_radius_nm / dca, label_wall_stiffness), axis=1
                    )
                    weights = np.exp(-np.clip(clash, 0.0, 700.0))
                else:  # 'hard'
                    weights = (np.min(dca, axis=1) >= bead_radius_nm).astype(float)
            else:
                # no other CA to clash against (e.g. a very short chain)
                weights = np.ones(beads.shape[0])

            frame_weight = float(np.sum(weights))
            if frame_weight <= 1e-12:
                # every bead is excluded this frame; nothing to contribute
                continue
            frame_count += 1

            for idx in residue_list:
                t = target_xyz[idx]
                if t is None:
                    continue
                # distances from every bead to this target atom
                d = np.linalg.norm(beads - t[f], axis=1)  # (n_conf,)
                # within-cloud weighted average of PREFACTOR / r^6 this frame
                relax_sum[idx] += (
                    np.sum(weights * (self.PREFACTOR / np.power(d, 6))) / frame_weight
                )

        gamma = []
        for idx in residue_list:
            if frame_count == 0 or target_xyz[idx] is None:
                gamma.append(np.nan)
            else:
                gamma.append(relax_sum[idx] / frame_count)
        return gamma
