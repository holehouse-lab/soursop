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

import mdtraj as md
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
        self, label_position, spin_label_atom="CB", target_relaxation_atom="N"
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
            r_6_nm = np.power(
                0.1
                * self.SSPO.get_inter_residue_atomic_distance(
                    label_position, idx, A1=spin_label_atom, A2=target_relaxation_atom
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
