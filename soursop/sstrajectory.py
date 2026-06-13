##     _____  ____  _    _ _____   _____  ____  _____
##   / ____|/ __ \| |  | |  __ \ / ____|/ __ \|  __ \
##  | (___ | |  | | |  | | |__) | (___ | |  | | |__) |
##   \___ \| |  | | |  | |  _  / \___ \| |  | |  ___/
##   ____) | |__| | |__| | | \ \ ____) | |__| | |
##  |_____/ \____/ \____/|_|  \_\_____/ \____/|_|

## Alex Holehouse (Pappu Lab and Holehouse Lab) and Jared Lalmansingh (Pappu lab)
## Simulation analysis package
## Copyright 2014 - 2026
##

import mdtraj as md
import numpy as np
from .configs import *
from .ssdata import ALL_VALID_RESIDUE_NAMES

from .ssprotein import SSProtein
from .ssexceptions import SSException
from . import ssutils
from . import ssio
from . import sstools

## Order of standard args:
## 1. stride
## 2. weights
## 3. verbose

from functools import partial, wraps
from multiprocessing import Pool, cpu_count


def lazy_loading_single_protein_trajectory(func):
    """
    Deocrator function that means we lazyily load the full
    trajectory (once) only if we needed it. This is a stand-
    alone decorator function that can be used on any function
    in the SSTrajectory class that needs the hidden
    single_protein_traj object.

    Parameters
    -----------
    func : function
        Function to be decorated

    Returns
    --------
    function
        Returns a function that will load the full protein trajectory
        if it's not already loaded.
    """

    @wraps(func)
    def wrapper(self, *args, **kwargs):
        if self._SSTrajectory__single_protein_traj is None:
            self._SSTrajectory__single_protein_traj = (
                self._SSTrajectory__get_all_proteins(
                    self.traj, self._SSTrajectory__explicit_residue_checking
                )
            )
        return func(self, *args, **kwargs)

    return wrapper


class SSTrajectory:
    # oxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxo
    #
    #
    def __init__(
        self,
        trajectory_filename=None,
        pdb_filename=None,
        TRJ=None,
        protein_grouping=None,
        pdblead=False,
        debug=False,
        extra_valid_residue_names=None,
        explicit_residue_checking=False,
        print_warnings=False,
        swan_trajectory=None,
    ):
        """
        SSTrajectory is the class used to read in a work with simulation
        trajectories in SOURSOP.

        There are two ways new SSTrajectory objects can be generated;

           1. By reading from disk (i.e. passing in a trajectory file and a topology file.

           2. By passing in an existing trajectory object that has already been read in.

        Both of these are described below.

        SOURSOP will, by default, extract out the protein component from
        your trajectory automatically, which lets you ask questions about the
        protein only (i.e. without salt ions getting in the way).

        You can also explicitly define which protein groups should be considered
        independently.

        Note that by default the mechanism by which individual proteins are
        identified is by cycling through the unique chains and determining
        if they are protein or not. You can also provide manual grouping
        via the protein_grouping option, which lets you define which
        residues should make up an individual protein. This can be useful
        if you have multiple proteins associated with the same chain, which
        happens in CAMPARI if you have more than 26 separate chains (i.e.
        every protein after the 26th is the 'Z' chain).

        **Class Variables**

        .proteinTrajectoryList : list
            List of SSProtein objects which themselves contain a large set of associated functions

        .traj : mdtraj.trajectory
            The underlying mdtraj.trajectory object enables any/all mdtraj analyses to be performed
            as one might want normally.

        .swan_trajectory : bool
            True if the trajectory was detected (or forced) to be a SWAN 2-bead
            (CA backbone / CB sidechain) coarse-grained model. This is auto-detected
            on load and propagated to every SSProtein in proteinTrajectoryList, which
            switches sidechain-vector and secondary-structure analyses to SWAN-aware
            CA/CB definitions.


        Parameters
        -------------
        trajectory_filename : str
            Filename which contains the trajectory file of interest. Normally
            this is `__traj.xtc` or `__traj.dcd`.

        pdb_filename : str
            Filename which contains the pdb file associated with the trajectory
            of interest. Normally this is `__START.pdb`.

        TRJ : mdtraj.Trajectory
            It is sometimes useful to re-defined a trajectory and create a new
            SSTrajectory  object from that trajectory. This could be done by
            writing that new trajectory  to file, but this is extremely slow
            due to the I/O impact of reading/writing  from disk. If an mdtraj
            trajectory objected is passed, this is used as the new trajectory
            from which the SSTrajectory object is constructed. Default = None

        protein_grouping : list of lists of ints
            Lets you manually define protein groups to be considered independently.
            Default = None

        pdblead : bool
            Lets you set the PDB file (which is normally ONLY used as a topology
            file) to be the first frame of the trajectory. This is useful when
            the first PDB file holds some specific reference information which
            you want to use (e.g. RMSD or Q). Default = False

        debug : book
            Prints warning/help information to help debug weird stuff during
            initial  trajectory read-in. Default = False.

        extra_valid_residue_names : list
            By default, SOURSOP identifies chains as proteins based on the a set
            of normally seen protein residue names. These are defined in
            soursop/ssdata, and are listed below::

                'ALA', 'CYS', 'ASP', 'ASH', 'GLU', 'GLU', 'PHE', 'GLY',
                'HIE', 'HIS', 'HID', 'HIP', 'ILE', 'LEU', 'LYS', 'LYD',
                'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL',
                'TRP', 'TYR', 'AIB', 'ABA', 'NVA', 'NLE', 'ORN', 'DAB',
                'PTR', 'TPO', 'SEP', 'KAC', 'KM1', 'KM2'  'KM3', 'ACE',
                'NME', 'FOR', 'NH2'

            This keyword allows the user to pass a list of ADDITIONAL residues
            that we want SOURSOP to recognize as valid residues to extract a
            chain as a protein molecule. This can be especially useful if you
            want to trick SOURSOP into analyzing polymer simulations where your
            PDB file may have non-standard residue names (e.g., XXX).

        explicit_residue_checking : bool
            In early versions of SOURSOP we operate under the assumption that
            every residue in a chain will be protein if the first residue is
            a protein. In general this is a good assumption, especially because
            it means we can implicitly deal with non-standard residue names that
            are internal to a chain. However, if you're reading in a .gro file, or
            any kind of topology file that lacks explicit chains, then any
            non-protein residues end up being brought in and counted which if you
            have 10000s of solvent molecules makes parsing impossible. If
            explicit_residue_checking is set to True, then each residue in each
            chain is explicitly checked, meaning solvent molecules/atoms in
            a single chain are discarded.

        print_warnings : bool
            Print warnings if the unit cell lengths are zero or not set. This
            is a common issue with old CAMPARI trajectories or trajectories
            generated without CRYSTAL records, but is generally not an issue.
            Default = False

        swan_trajectory : bool or None
            Controls SWAN 2-bead (CA/CB) coarse-grained handling. If None
            (default) SOURSOP auto-detects whether the topology is a SWAN model
            (a single CA per residue plus a single CB per non-glycine residue,
            and nothing else) on load. Pass True or False to force the behaviour
            and skip auto-detection. The resulting value is stored as the
            ``.swan_trajectory`` attribute and propagated to every SSProtein.
            Default = None

        Example
        -----------
        Example of reading in an XTC trajectory file::

            from soursop.sstrajectory import SSTrajectory

            TrajOb = SSTrajector('traj.xtc','start.pdb')


        """

        self.valid_residue_names = []
        self.valid_residue_names.extend(ALL_VALID_RESIDUE_NAMES)

        # save this so we can lazy load the full protein trajectory as needed
        self.__explicit_residue_checking = explicit_residue_checking

        if extra_valid_residue_names is not None:
            try:
                self.valid_residue_names.extend(extra_valid_residue_names)
            except Exception:
                print(
                    "Unable to use the extra_valid_residue_names - this must be a list of strings"
                )

        # first we decide if we're reading from file or from an existing trajectory
        if (trajectory_filename is None) and (pdb_filename is None):
            if TRJ is None:
                raise SSException(
                    "No input provided! Please provide ether a pdb and trajectory file OR a pre-formed traj object"
                )

            # note the [:] means this is a COPY!
            self.traj = TRJ[:]
        else:
            if trajectory_filename is None:
                raise SSException("No trajectory file provided!")
            if pdb_filename is None:
                raise SSException("No PDB file provided!")

            # read in the raw trajectory
            self.traj = self.__readTrajectory(
                trajectory_filename, pdb_filename, pdblead, print_warnings
            )

        # determine whether this is a SWAN 2-bead (CA/CB) coarse-grained model.
        # This is auto-detected from the topology unless explicitly forced via
        # the swan_trajectory keyword, and must be resolved BEFORE the per-protein
        # SSProtein objects are built so the flag can be threaded into them.
        if swan_trajectory is None:
            self.swan_trajectory = ssutils.is_swan_topology(self.traj.topology)
        else:
            self.swan_trajectory = bool(swan_trajectory)

        # Next, having read in the trajectory we parse out into proteins
        # extract a list of protein trajectories where each protein is assumed
        # to be in its own chain
        if protein_grouping == None:
            self.proteinTrajectoryList = self.__get_proteins(
                self.traj, debug, explicit_residue_checking=explicit_residue_checking
            )
        else:
            self.proteinTrajectoryList = self.__get_proteins_by_residue(
                self.traj, protein_grouping, debug
            )

        # this is initialized to None and then gets defined using the lazy_loading_single_protein_trajectory()
        # decorator when it's first needed
        self.__single_protein_traj = None

    def __repr__(self):
        return "SSTrajectory (%s): %i proteins and %i frames" % (
            hex(id(self)),
            self.n_proteins,
            self.n_frames,
        )

    def __len__(self):
        # Edited to mimic the behavior of `mdtraj` trajectory objects.
        # Originally: `return (self.n_proteins, self.n_frames)`
        return self.n_frames

    @property
    def n_frames(self):
        """Number of frames in the underlying mdtraj trajectory.

        Returns
        -------
        int

        Example
        -------
        >>> traj.n_frames
        1000
        """
        return len(self.traj)

    @property
    def n_proteins(self):
        """Number of distinct protein chains identified during loading.

        For a single-chain PDB this is 1; for multi-chain inputs (or when a
        manual ``protein_grouping`` list was provided to ``SSTrajectory``)
        the value matches the number of chains.

        Returns
        -------
        int

        Example
        -------
        >>> two_chain_traj.n_proteins
        2
        """
        return len(self.proteinTrajectoryList)

    @property
    def length(self):
        """Convenience ``(n_proteins, n_frames)`` summary tuple.

        Preserves the historic behaviour of ``__len__`` from before
        SOURSOP switched its ``len(traj)`` semantics to match mdtraj
        (frame count only).

        Returns
        -------
        tuple of (int, int)
            ``(n_proteins, n_frames)``.

        Example
        -------
        >>> traj.length
        (2, 1000)
        """
        return (self.n_proteins, self.n_frames)

    @property
    def unitcell(self):
        """Unit-cell box lengths of the first frame, in Angstroms.

        SOURSOP reports distances in Angstroms, so the values stored in
        nanometres by mdtraj are multiplied by 10 here. Raises
        ``TypeError`` if the trajectory has no periodic box (some
        coarse-grained or vacuum simulations).

        Returns
        -------
        np.ndarray
            Array of shape ``(3,)`` with box lengths in Angstroms.

        Example
        -------
        >>> traj.unitcell
        array([60.0, 60.0, 60.0])
        """
        return self.traj.unitcell_lengths[0] * 10

    # oxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxo
    #
    #
    def __readTrajectory(
        self, trajectory_filename, pdb_filename, pdblead, print_warnings
    ):
        """
        Internal function which parses and reads in a CAMPARI trajectory

        Read a trajectory file. This was separated out into its own
        function in case we want to add additional sanity checks during
        the file loading.

        Notably older versions of CAMPARI mess up the unitcell length
        vectors, so will cause problems, but you can get around this by
        having GROMACS rebuild the trajectory. If this happens an error
        pops up but instructions on how to use GROMACS to fix it are
        presented.

        Parameters
        -----------
        trajectory_filename : str
            Filename which contains the trajectory file of interest. File type
            is automatically detected and dealt with mdtraj' 'load' command
            (i.e. md.load(filename, top=pdb_filename))

        pdb_filename : str
            Filename which contains the pdb file associated with the trajectory
            of interest. This defines the topology of the system and must match
            the trajectory in terms of number of atomas


        pdblead : bool
            Also extract the coordinates from the PDB file and append it to
            the front of the trajectory. This is useful if you are starting
            an analysis where that first structure should be a reference frame
            but it's not actually included in the trajectory file.

        print_warnings : bool
            Print warnings if the unit cell lengths are zero or not set. This
            is a common issue with old CAMPARI trajectories or trajectories
            generated without CRYSTAL records, but is generally not an issue.

        Returns
        --------
        mdtraj.traj
            Returns an mdtraj trajectory object

        """

        # straight up read the trajectory first using mdtraj's awesome
        # trajectory reading facility
        traj = md.load(trajectory_filename, top=pdb_filename)

        # check unit cell lengths
        try:
            uc_lengths = traj.unitcell_lengths[0]

            # this is s custom warning for a specific edge-case we encounter a lot
            if uc_lengths[0] == 0 or uc_lengths[1] == 0 or uc_lengths[2] == 0:
                if print_warnings:
                    ssio.warning_message(
                        "Trajectory file unit cell lengths are zero for at least one dimension. This is a probably a bug with an FRC generated __START.pdb file, because in the old version of CAMPARI used to do grid based FRC calculations the unit cell dimensions are not written correctly. This may cause issues but we're going to assume everything is OK for now. Check the validity of any analysis output. If you're worried, you can use the following workaround.\n\n:::: WORK AROUNDS ::::\nSimply run\n\ntrjconv -f __traj.xtc -s __START.pdb -box a b c -o frc.xtc \n\nAn then \n\ntrjconv -f frc.xtc -s __START.pdb -box a b c -o start.pdb -dump 0\n\n\nHere\n-f n__traj.xtc   : defines the trajectory file\n-s __start.pdb   : defines the pdb file used to parse the topology\n-box a b c       : defines the box lengths **in nanometers**\n-o frc.xtc       : is the name of the new trajectory file with updated box lengths\nSelect 0 (system) when asked to 'Select group for output'.The second step creates the equivalent PDB file with the header-line correctly defining the box unit cell lengths and angles. These two new files should then be used for analysis.\n\nAs an example, if my FRC simulation had a sphere radius of 100 angstroms then my correction command would look something like \n\ntrjconv -f __traj.xtc -s __START.pdb -box 20 20 20 -o frc.xtc\ntrjconv -f frc.xtc -s __START.pdb -box 20 20 20 -o start.pdb -dump 0"
                    )

        except TypeError:
            if print_warnings:
                ssio.warning_message(
                    "Warning: UnitCell lengths were not provided... This may cause issues but we're going to assume everything is OK for now..."
                )

        # if pdbLead is true then load the pdb_filename as a trajectory
        # and then add it to the front (the PDB file is its own topology
        # file so no need to specificy the top= file here!
        if pdblead:
            pdbtraj = md.load(pdb_filename)
            traj = pdbtraj + traj

            # having added the PDB file now check all the unit-cells match up!
            try:
                uc_lengths_1 = traj.unitcell_lengths[0]
                uc_lengths_2 = traj.unitcell_lengths[1]

                if (
                    (uc_lengths_1[0] != uc_lengths_2[0])
                    or (uc_lengths_1[1] != uc_lengths_2[1])
                    or (uc_lengths_1[2] != uc_lengths_2[2])
                ):
                    ssio.warning_message(
                        f"........................\nWARNING:\nThe unit cell dimensions of the PDB file and trajectory file did not match, specifically\nPDB file = [ {uc_lengths_1} ]\nXTC file = [ {uc_lengths_2} ]\nThis may cause issues if native MDTraj utilities are used (and potentially for SOURSOP utilities that are based on these). It is not necessarily an issue, but PLEASE sanity check your outcome. To be safe we recommend editing the PDB-file unitcell dimensions to match."
                    )

            except TypeError:
                ssio.warning_message(
                    "Warning: UnitCell lengths were not provided... This may cause issues but we're going to assume everything is OK for now..."
                )

        return traj

    # oxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxo
    #
    #
    def __get_all_proteins(self, trajectory, explicit_residue_checking=False):
        """
        Internal function that builds a single trajectory which contains all protein
        residues.

        Parameters
        -----------
        trajectory : mdtraj.Trajectory
            An already parsed trajectory object (i.e. checked for CAMPARI-
            relevant defects such as unitcell issues etc)

        Returns
        ----------
        ssprotein.SSProtein
            Returns a single SSProtein object

        """

        # extract full system topology
        topology = trajectory.topology

        protein_atoms = []

        # for each chain in this toplogy determine if the
        # first residue is protein or not. If it's protein we parse it if
        # not it gets skipped

        for chain in topology.chains:
            # if explicit residue checking is False we assume the first residue of the
            # chain is representative of the whole chain. This is generally
            # fair assumption
            if explicit_residue_checking is False:
                if chain.residue(0).name in self.valid_residue_names:
                    # get every atom in this chain
                    protein_atoms.extend(atom.index for atom in chain.atoms)
            else:
                # initialize an empty list of atoms
                local_atoms = []

                for res in chain.residues:
                    # for each residue in that chain ask if that residue is valid
                    if res.name in self.valid_residue_names:
                        # if yes
                        local_atoms.extend([a.index for a in res.atoms])

                protein_atoms.extend(local_atoms)

        return SSProtein(
            trajectory.atom_slice(protein_atoms), swan=self.swan_trajectory
        )

    # oxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxo
    #
    #
    def __get_proteins(self, trajectory, debug, explicit_residue_checking=False):
        """
        Internal function that takes an MDTraj trajectory and returns a list
        of mdtraj trajectory objects corresponding to each protein in the
        system, ASSUMING that each protein is in its own chain.

        The way this works is to cycle through each chain, identify if that
        chain contains protein or not, and if it does grab all the atoms in
        that chain and perform an atomslice using those atoms on the main
        trajectory.

        Parameters
        -----------
        trajectory : mdtraj.Trajectory
            An already parsed trajectory object (i.e. checked for CAMPARI-
            relevant defects such as unitcell issues etc)

        Returns
        ---------
        list
            Returns a proteinTrajectoryList - a list with one or more
            SSProtein objects in it

        """

        # extract full system topology
        topology = trajectory.topology

        chainAtoms = []

        # for each chain in this toplogy determine if the
        # first residue is protein or not. If it's protein we parse it if
        # not it gets skipped
        for chain in topology.chains:
            # if hybrid chains is False we assume the first residue of the
            # chain is representative of the whole chain. This is generally
            # fair assumption
            if explicit_residue_checking is False:
                # if the first residue in the chain is protein
                # so we include an edgecase here for that
                if chain.residue(0).name in self.valid_residue_names:
                    # intialize an empty list of atoms
                    local_atoms = []

                    # get every atom in this chain
                    for atom in chain.atoms:
                        local_atoms.append(atom.index)

                    chainAtoms.append(local_atoms)

                else:
                    if debug:
                        ssio.debug_message(
                            "Skipping residue %s from %s"
                            % (chain.residue(0).name, chain)
                        )
            else:
                # initialize an empty list of atoms
                local_atoms = []

                for res in chain.residues:
                    # for each residue in that chain ask if that residue is valid
                    if res.name in self.valid_residue_names:
                        # if yes
                        local_atoms.extend([a.index for a in res.atoms])

                chainAtoms.append(local_atoms)

        # for each protein chain that we have atomic indices
        # for (hopefully all of them!) cycle through and create
        # sub-trajectories
        proteinTrajectoryList = []

        for local_chain_atoms in chainAtoms:
            # generate a trajectory composed of *JUST* the
            # $chain atoms. The PT object created now is self
            # consistent and contains an associated and fully
            # correct .topology object (NOTE this fixes a
            # previous bug in CAMPARITraj 0.1.4)

            PT = trajectory.atom_slice(local_chain_atoms)

            # WA
            # gets the resid offset in a way that is ensures internal
            # consistency for the SSProtein object
            first_resid = PT.topology.chain(0).residue(0).index

            if first_resid != 0:
                raise SSException(
                    "After extracting a protein subtrajectory, the first resid is not 0. This may reflect a bug, or you may not be using MDTraj 1.9.5"
                )

            # add that SSProtein to the ever-growing proteinTrajectory list. NOTE that
            # THIS is the line that is the source of most of the slowness for trajectory
            # loading....
            proteinTrajectoryList.append(SSProtein(PT, swan=self.swan_trajectory))

        if len(proteinTrajectoryList) == 0:
            ssio.warning_message("No protein chains found in the trajectory")

        return proteinTrajectoryList

    # oxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxo
    #
    #
    def __get_proteins_by_residue(self, trajectory, residue_grouping, debug):
        """
        Internal function which returns a list of mdtraj trajectory objects
        corresponding to each protein where we *explicitly* define the residues
        in each protein.


        Unlike the `__get_proteins()` function, which doesn't require any
        manual input in identifying the proteins, here we provide a list of
        groups, where each group is the set of residues associated with a
        protein.

        The way this works is to cycle through each group, and for each
        residue  in each group grabs all the atoms and uses these to carry
        out an  atomslice on the full trajectory.

        Parameters
        -----------
        trajectory : mdtraj.Trajectory
            An already parsed trajectory object (i.e. checked for CAMPARI-
            relevant defects such as unitcell issues etc)

        residue_grouping : list of list of integers
            Must be a list containing one or more lists, where each internal
            list contains a set of monotonically increasing residues (which
            correspond to the full protein trajectory). In other words, each
            sub-list  defines a single protein. The integer indexing here -
            importantly - uses the  CAMPARITraj internal residue indexing,
            meaning that  indexing begins at 0 from the first residue in the
            PDB file.

        Returns
        ---------
        list
            Returns a proteinTrajectoryList - a list with one or more SSProtein
            objects in it.

        """

        group_atoms = []

        # extract full system topology
        topology = trajectory.topology

        # for each chain in this toplogy determine if the
        # first residue is protein or not
        for group in residue_grouping:
            # build a string of the resids
            res_string = ""
            for r in group:
                res_string = res_string + " %i" % (int(r))

            # select atoms based on the resid string
            local_atoms = topology.select("resid %s" % (res_string))

            if len(local_atoms) == 0:
                ssio.warning_message(
                    "In residue group [%s ...] no residues in the trajectory were found..."
                    % (str(group)[0:8])
                )

            else:
                group_atoms.append(local_atoms)

        # for each protein group that we have atomic indices
        # for cycle through and create sub-trajectories
        proteinTrajectoryList = []

        # cycle through each group of atoms as was defined by the residue_grouping
        # indices

        for local_group_atoms in group_atoms:
            # generate a trajectory composed of *JUST* the
            # $chain atoms. The PT object created now is self
            # consistent and contains an associated and fully
            # correct .topology object (NOTE this fixes a
            # previous bug in CAMPARITraj 0.1.4)
            PT = trajectory.atom_slice(local_group_atoms)

            # gets the resid offset in a way that is ensures internal
            # consistency for the SSProtein object
            resid_offset = PT.topology.chain(0).residue(0).index

            if resid_offset != 0:
                raise SSException(
                    "After extracting a protein subtrajectory, the first resid is not 0. This may reflect a bug, or you may not be using MDTraj 1.9.5"
                )

            proteinTrajectoryList.append(SSProtein(PT, swan=self.swan_trajectory))

        if len(proteinTrajectoryList) == 0:
            ssio.warning_message("No protein chains found in the trajectory")

        return proteinTrajectoryList

    # oxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxo
    #
    #
    @lazy_loading_single_protein_trajectory
    def get_overall_radius_of_gyration(self, weights=False, etol=0.0000001):
        """Per-frame radius of gyration computed across every protein chain.

        For multi-chain systems the chains are combined into a single
        "pseudo-chain" before computing :math:`R_g`. For single-chain
        systems the result is identical to
        ``self.proteinTrajectoryList[0].get_radius_of_gyration()``.

        .. warning::
            This does NOT apply periodic-boundary corrections. If your
            chains are split across PBC images, centre the molecule first.

        Parameters
        ----------
        weights : array_like or False, optional
            Per-frame re-weighting vector forwarded to
            :meth:`SSProtein.get_radius_of_gyration`. ``False`` (default)
            returns the per-frame array; if supplied the scalar
            deterministic weighted-mean :math:`R_g` is returned.
        etol : float, optional
            Tolerance on ``|sum(weights) - 1|``. Default ``1e-7``.

        Returns
        -------
        np.ndarray or float
            Per-frame :math:`R_g` (length ``n_frames``), or the scalar
            weighted mean if ``weights`` is supplied. Angstroms.

        Example
        -------
        >>> rg_all = traj.get_overall_radius_of_gyration()
        """

        return self.__single_protein_traj.get_radius_of_gyration(
            weights=weights, etol=etol
        )

    # oxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxo
    #
    #
    @lazy_loading_single_protein_trajectory
    def get_overall_asphericity(self, weights=False, etol=0.0000001):
        """Per-frame asphericity computed across every protein chain.

        For multi-chain systems the chains are combined into a single
        "pseudo-chain" before computing the asphericity. For single-chain
        systems the result equals
        ``self.proteinTrajectoryList[0].get_asphericity()``.

        .. warning::
            This does NOT apply periodic-boundary corrections.

        Parameters
        ----------
        weights : array_like or False, optional
            Per-frame re-weighting vector forwarded to
            :meth:`SSProtein.get_asphericity`. ``False`` (default) returns
            the per-frame array; if supplied the scalar deterministic
            weighted-mean asphericity is returned.
        etol : float, optional
            Tolerance on ``|sum(weights) - 1|``. Default ``1e-7``.

        Returns
        -------
        np.ndarray or float
            Per-frame asphericity (length ``n_frames``), or the scalar
            weighted mean if ``weights`` is supplied.

        Example
        -------
        >>> asph_all = traj.get_overall_asphericity()
        """

        return self.__single_protein_traj.get_asphericity(weights=weights, etol=etol)

    # oxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxo
    #
    #
    @lazy_loading_single_protein_trajectory
    def get_overall_hydrodynamic_radius(self, weights=False, etol=0.0000001):
        """Per-frame hydrodynamic radius computed across every protein chain.

        For multi-chain systems the chains are combined into a single
        "pseudo-chain" before computing :math:`R_h`. Uses the default
        Nygaard-et-al estimator. For single-chain systems the result equals
        ``self.proteinTrajectoryList[0].get_hydrodynamic_radius()``.

        Parameters
        ----------
        weights : array_like or False, optional
            Per-frame re-weighting vector forwarded to
            :meth:`SSProtein.get_hydrodynamic_radius`. ``False`` (default)
            returns the per-frame array; if supplied the scalar
            deterministic weighted-mean :math:`R_h` is returned.
        etol : float, optional
            Tolerance on ``|sum(weights) - 1|``. Default ``1e-7``.

        Returns
        -------
        np.ndarray or float
            Per-frame :math:`R_h` (length ``n_frames``), or the scalar
            weighted mean if ``weights`` is supplied. Angstroms.

        Example
        -------
        >>> rh_all = traj.get_overall_hydrodynamic_radius()
        """

        return self.__single_protein_traj.get_hydrodynamic_radius(
            weights=weights, etol=etol
        )

    # oxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxo
    #
    #
    def get_interchain_distance_map(
        self,
        proteinID1,
        proteinID2,
        mode="CA",
        periodic=False,
        weights=False,
        etol=0.0000001,
    ):
        """Mean and std inter-chain distance map between two protein chains.

        Returns an ``(n_res_P1, n_res_P2)`` matrix of per-pair *mean*
        inter-residue distances along with the matching per-pair standard
        deviations. Distances are in Angstroms.

        Calling ``get_interchain_distance_map(i, i)`` produces the
        intra-chain distance map of protein ``i``, which is the same as
        ``self.proteinTrajectoryList[i].get_distance_map()`` and makes a
        useful sanity check.

        Parameters
        ----------
        proteinID1 : int
            Index of the first protein in ``self.proteinTrajectoryList``.
        proteinID2 : int
            Index of the second protein in ``self.proteinTrajectoryList``.
        mode : {'CA', 'COM'}, optional
            * ``'CA'`` (default) - distances between alpha-carbon atoms.
            * ``'COM'`` - distances between residue centres of mass.
        periodic : bool, optional
            If True, apply the minimum-image convention. Requires a
            recorded periodic box and a cubic cell. Generally it is better
            to centre the molecule first and leave this False. Default
            False.
        weights : array_like or False, optional
            Per-frame re-weighting vector (validated against the shared
            trajectory frame count). ``False`` (default) gives the
            ordinary per-pair mean/std (byte-identical to before); if
            supplied each pair's mean/std is the deterministic weighted
            mean / weighted (population) std over frames.
        etol : float, optional
            Tolerance on ``|sum(weights) - 1|``. Default ``1e-7``.

        Returns
        -------
        tuple of (np.ndarray, np.ndarray)
            ``(distance_map, std_map)`` each of shape
            ``(n_res_P1, n_res_P2)``, in Angstroms.

        Example
        -------
        >>> dmap, smap = two_chain_traj.get_interchain_distance_map(0, 1)
        >>> dmap.shape
        (58, 58)
        """

        ssutils.validate_keyword_option(mode, ["CA", "COM"], "mode")

        # get SSProtein objects for the two IDs passed (could be the same)
        P1 = self.proteinTrajectoryList[proteinID1]
        P2 = self.proteinTrajectoryList[proteinID2]

        # optional deterministic per-frame re-weighting of every pair's
        # mean/std (validated against the shared trajectory frame count).
        wv = ssutils.validate_weights(weights, P1.n_frames, 1, etol)

        # create the empty distance maps
        p1_residues = P1.resid_with_CA
        p2_residues = P2.resid_with_CA
        map_shape = (len(p1_residues), len(p2_residues))
        distanceMap = np.zeros(map_shape)
        stdMap = np.zeros(map_shape)

        # Hoist the per-residue center-of-mass computation out of the
        # nested loop. A residue's COM is independent of the residue it
        # is paired with, so computing it once per residue (O(n1 + n2))
        # rather than once per pair (O(n1 * n2)) is numerically identical
        # and removes the dominant cost. The atom_name selection per mode
        # is exactly that of the original per-pair calls.
        if mode == "COM":
            com1 = [P1.get_residue_COM(r1) for r1 in p1_residues]
            com2 = [P2.get_residue_COM(r2) for r2 in p2_residues]
        else:
            com1 = [P1.get_residue_COM(r1, atom_name="CA") for r1 in p1_residues]
            com2 = [P2.get_residue_COM(r2, atom_name="CA") for r2 in p2_residues]

        # the (unchanged) ncap-shifted output indices
        p1_indices = [(r1 - 1) if P1.ncap else r1 for r1 in p1_residues]
        p2_indices = [(r2 - 1) if P2.ncap else r2 for r2 in p2_residues]

        if periodic:
            # preserve the exact original per-pair minimum-image call
            for i, COM_1 in enumerate(com1):
                p1_index = p1_indices[i]
                for j, COM_2 in enumerate(com2):
                    d = sstools.get_distance_periodic(
                        COM_1, COM_2, self.unitcell[0], "cube"
                    )
                    if wv is False:
                        distanceMap[p1_index, p2_indices[j]] = np.mean(d, 0)
                        stdMap[p1_index, p2_indices[j]] = np.std(d, 0)
                    else:
                        d = np.asarray(d)
                        distanceMap[p1_index, p2_indices[j]] = ssutils.weighted_mean(
                            d, wv
                        )
                        stdMap[p1_index, p2_indices[j]] = ssutils.weighted_std(d, wv)
        else:
            # com2 stacked once -> (n2, F, 3); broadcasting COM_1 (F, 3)
            # against it reproduces the per-pair np.linalg.norm exactly.
            com2_stack = np.stack(com2, axis=0)
            for i, COM_1 in enumerate(com1):
                d = np.linalg.norm(COM_1 - com2_stack, axis=-1)  # (n2, F)
                if wv is False:
                    row_mean = np.mean(d, axis=1)
                    row_std = np.std(d, axis=1)
                else:
                    row_mean = ssutils.weighted_mean(d, wv, axis=1)
                    row_std = ssutils.weighted_std(d, wv, axis=1)
                p1_index = p1_indices[i]
                for j, p2_index in enumerate(p2_indices):
                    distanceMap[p1_index, p2_index] = row_mean[j]
                    stdMap[p1_index, p2_index] = row_std[j]

        return (distanceMap, stdMap)

    # oxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxo
    #
    #
    def get_interchain_contact_map(
        self,
        proteinID1,
        proteinID2,
        threshold=5.0,
        mode="atom",
        A1="CA",
        A2="CA",
        periodic=False,
        stride=1,
        verbose=False,
    ):
        """Inter-chain contact-fraction map between two protein chains.

        Returns an ``(n_res_P1, n_res_P2)`` matrix where entry ``[i, j]`` is
        the fraction of (strided) frames in which residue ``i`` of chain
        ``proteinID1`` is within ``threshold`` Angstroms of residue ``j``
        of chain ``proteinID2`` under the chosen distance mode. Calling
        with ``proteinID1 == proteinID2`` produces the chain's intra-chain
        contact map (caps and glycines contribute zero rows/columns for
        modes that can't be evaluated for them).

        Parameters
        ----------
        proteinID1 : int
            Index of the first protein in ``self.proteinTrajectoryList``.
        proteinID2 : int
            Index of the second protein in ``self.proteinTrajectoryList``.
        threshold : float, optional
            Contact distance cutoff in Angstroms. Default 5.0.
        mode : {'atom', 'ca', 'closest', 'closest-heavy', 'sidechain', 'sidechain-heavy'}, optional
            How inter-residue distance is defined (see
            :meth:`SSProtein.get_inter_residue_atomic_distance` for full
            descriptions). For ``mode='atom'`` the ``A1`` / ``A2`` atom
            names are used; for every other mode they are ignored.
            Default ``'atom'``.
        A1 : str, optional
            Atom name in residue R1 when ``mode='atom'``. Default ``'CA'``.
        A2 : str, optional
            Atom name in residue R2 when ``mode='atom'``. Default ``'CA'``.
        periodic : bool, optional
            If True, apply the minimum-image convention to the
            mdtraj-driven modes (closest, closest-heavy, sidechain,
            sidechain-heavy). Default False.
        stride : int, optional
            Use every ``stride``-th frame when computing the per-pair
            distances; the contact fraction is then
            ``sum(distances < threshold) / n_strided_frames``. Default 1.
        verbose : bool, optional
            If True, print one status line per residue in the outer loop.
            Default False.

        Returns
        -------
        np.ndarray
            Array of shape ``(n_res_P1, n_res_P2)`` of contact fractions in
            ``[0, 1]``.

        Raises
        ------
        SSException
            If either ``proteinID`` is out of range, ``mode`` is invalid,
            ``stride`` is not a positive integer, or every per-residue
            inner call failed (e.g. an invalid atom name was supplied).

        Example
        -------
        >>> cmap = two_chain_traj.get_interchain_contact_map(0, 1, mode='closest-heavy')
        >>> cmap.shape
        (58, 58)
        """

        # Eagerly validate proteinIDs so an out-of-range value raises
        # SSException up front instead of an IndexError partway through.
        n_proteins = len(self.proteinTrajectoryList)
        for label, pid in (("proteinID1", proteinID1), ("proteinID2", proteinID2)):
            if not (0 <= pid < n_proteins):
                raise SSException(
                    f"In get_interchain_contact_map(): {label}={pid} is out of "
                    f"range; valid indices are 0..{n_proteins - 1}"
                )

        # Validate the mode keyword up front. Doing this here avoids the
        # silent all-zero matrix that would otherwise be produced by the
        # per-residue try/except below when every call raises with the
        # same "bad mode" error.
        allowed_modes = [
            "atom",
            "ca",
            "closest",
            "closest-heavy",
            "sidechain",
            "sidechain-heavy",
        ]
        if mode not in allowed_modes:
            raise SSException(
                f"In get_interchain_contact_map(): mode keyword must be one of "
                f"{allowed_modes}. Provided keyword was [{mode}]"
            )

        if not (isinstance(stride, (int, np.integer)) and stride >= 1):
            raise SSException(
                f"In get_interchain_contact_map(): stride must be a positive "
                f"integer, got {stride!r}"
            )

        # get number of residues/bases for the two proteins
        n_res_P1 = self.proteinTrajectoryList[proteinID1].n_residues
        n_res_P2 = self.proteinTrajectoryList[proteinID2].n_residues

        P1 = self.proteinTrajectoryList[proteinID1]
        P2 = self.proteinTrajectoryList[proteinID2]

        # Fast path: for mode='atom' (non-periodic) the named-atom
        # position of a residue is independent of the residue it is
        # paired with. Compute each residue's atom position once
        # (O(n1 + n2)) instead of re-deriving it inside an O(n1 * n2)
        # loop of get_interchain_distance() calls. This reproduces the
        # per-pair atom math exactly: same residue/atom selection, same
        # frame stride, the same 10x compute_center_of_mass, Euclidean
        # norm, threshold, and the same zero-fill for residues whose
        # A1/A2 atom is absent (caps, missing atom names). The periodic
        # branch is intentionally left to the original per-pair path so
        # its behaviour is byte-for-byte unchanged.
        if mode == "atom" and not periodic:

            def _residue_atom_positions(P, n_res, atom_name):
                # positions[r] = 10x-COM (F,3) of the named atom of
                # residue r over the strided frames; ok[r] mirrors the
                # exact success/failure of the original
                # get_interchain_distance() 'atom' selection for that
                # residue (the conditions the caller's except catches).
                positions = [None] * n_res
                ok = [False] * n_res
                for r in range(n_res):
                    try:
                        sel = P.topology.select("resid %i" % r)
                        if len(sel) == 0:
                            continue
                        sub = P.traj.atom_slice(sel)
                        if stride > 1:
                            sub = sub[::stride]
                        a = sub.topology.select('resid 0 and name "%s"' % atom_name)
                        if len(a) != 1:
                            continue
                        positions[r] = 10 * md.compute_center_of_mass(sub.atom_slice(a))
                        ok[r] = True
                    except (SSException, ValueError, IndexError):
                        continue
                return positions, ok

            pos1, ok1 = _residue_atom_positions(P1, n_res_P1, A1)
            pos2, ok2 = _residue_atom_positions(P2, n_res_P2, A2)

            # The original per-pair loop succeeds for a pair iff both
            # residues' atom selections succeed (independently), so the
            # number of successful pairs is the product of the per-chain
            # success counts. Preserve the original "all failed -> raise"
            # behaviour and exception message.
            success_count = int(np.sum(ok1)) * int(np.sum(ok2))
            if success_count == 0:
                raise SSException(
                    f"In get_interchain_contact_map(): no residue pair could be "
                    f"computed for mode={mode!r} (A1={A1!r}, A2={A2!r}). Check that "
                    f"the mode/atom selections are valid for the residues in "
                    f"protein {proteinID1} and protein {proteinID2}."
                )

            cmap = np.zeros((n_res_P1, n_res_P2))
            ok2_idx = [j for j in range(n_res_P2) if ok2[j]]
            # success_count > 0 guarantees ok2_idx is non-empty here
            pos2_stack = np.stack([pos2[j] for j in ok2_idx], axis=0)  # (m, F, 3)
            for i in range(n_res_P1):
                if verbose:
                    print(f"On {i} of {n_res_P1}")
                if not ok1[i]:
                    continue
                d = np.linalg.norm(pos1[i] - pos2_stack, axis=-1)  # (m, F)
                frac = np.sum(d < threshold, axis=1) / d.shape[1]
                for k, j in enumerate(ok2_idx):
                    cmap[i, j] = frac[k]

            return cmap

        all_contact_fractions = []
        # Track successful (R1, R2) computations. If zero pairs succeed (e.g.,
        # mode='atom' with a non-existent atom name, or mode='sidechain' on a
        # poly-glycine chain) we raise rather than return a silent zero matrix.
        success_count = 0

        # cycle over each residue in protein 1
        for p1_res_idx in range(0, n_res_P1):
            if verbose:
                print(f"On {p1_res_idx} of {n_res_P1}")

            tmp = []

            # cycle over each residue in protein 2
            for p2_res_idx in range(0, n_res_P2):
                # Per-residue computations can legitimately fail in two cases:
                #   1. Cap residues (ACE / NME) lack a CA atom, so mode='atom'
                #      with A1='CA' (default) and mode='ca' both raise.
                #   2. Glycines lack a sidechain, so mode='sidechain' and
                #      mode='sidechain-heavy' raise (and mdtraj's empty
                #      sidechain selector emits a ValueError further down
                #      the stack).
                # In both cases the correct contact fraction is 0 — there is
                # nothing to be in contact via the requested mode. Catch and
                # zero-fill so the matrix shape (n_res_P1, n_res_P2) is
                # preserved and the rest of the matrix still works.
                try:
                    distances = self.get_interchain_distance(
                        proteinID1,
                        proteinID2,
                        p1_res_idx,
                        p2_res_idx,
                        A1=A1,
                        A2=A2,
                        mode=mode,
                        periodic=periodic,
                        stride=stride,
                    )
                except (SSException, ValueError):
                    tmp.append(0.0)
                    continue

                # distances is already strided by get_interchain_distance, so
                # len(distances) is the number of frames actually evaluated.
                contact_fraction = np.sum(distances < threshold) / len(distances)
                tmp.append(contact_fraction)
                success_count += 1

            all_contact_fractions.append(tmp)

        if success_count == 0:
            raise SSException(
                f"In get_interchain_contact_map(): no residue pair could be "
                f"computed for mode={mode!r} (A1={A1!r}, A2={A2!r}). Check that "
                f"the mode/atom selections are valid for the residues in "
                f"protein {proteinID1} and protein {proteinID2}."
            )

        return np.array(all_contact_fractions)

    # oxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxo
    #
    #
    def get_interchain_distance(
        self,
        proteinID1,
        proteinID2,
        R1,
        R2,
        A1="CA",
        A2="CA",
        mode="atom",
        periodic=False,
        stride=1,
    ):
        """Per-frame distance between two residues, one in each chain.

        Resids are interpreted as they would be inside the corresponding
        :class:`SSProtein` (so the same R1=5 means "residue 5 of chain
        ``proteinID1``" rather than a global topology index).

        For ``mode='atom'`` the distance is between the named atoms
        ``A1`` of R1 in chain 1 and ``A2`` of R2 in chain 2. For every
        other mode the ``A1`` / ``A2`` arguments are ignored and the
        distance comes from ``mdtraj.compute_contacts``.

        Distance is returned in Angstroms.

        Parameters
        ----------
        proteinID1 : int
            Index of the first protein in ``self.proteinTrajectoryList``.
        proteinID2 : int
            Index of the second protein.
        R1, R2 : int
            Resids to measure between (R1 in chain 1, R2 in chain 2).
        A1, A2 : str, optional
            For ``mode='atom'``, the atom names. Defaults ``'CA'``.
        mode : {'atom', 'ca', 'closest', 'closest-heavy', 'sidechain', 'sidechain-heavy'}, optional
            How the distance is defined. See
            :meth:`SSProtein.get_inter_residue_atomic_distance` for full
            descriptions. Default ``'atom'``.
        periodic : bool, optional
            If True, apply the minimum-image convention for the
            mdtraj-driven modes. Default False.
        stride : int, optional
            Use every ``stride``-th frame. Returned array has length
            ``ceil(n_frames / stride)``. Must be a positive integer.
            Default 1.

        Returns
        -------
        np.ndarray
            1D array of length ``ceil(n_frames / stride)`` with the
            inter-chain distance in each (strided) frame, in Angstroms.

        Raises
        ------
        SSException
            If ``mode`` is invalid, a protein ID is out of range, a resid
            has no atoms, an atom name cannot be found, or ``stride`` is
            not a positive integer.

        Example
        -------
        >>> d = two_chain_traj.get_interchain_distance(0, 1, R1=10, R2=15)
        >>> d.shape
        (1000,)
        """

        # check mode keyword is valid
        allowed_modes = [
            "atom",
            "ca",
            "closest",
            "closest-heavy",
            "sidechain",
            "sidechain-heavy",
        ]
        if mode not in allowed_modes:
            raise SSException(
                "Provided mode keyword must be one of 'atom', 'ca', 'closest', 'closest-heavy', 'sidechain', or 'sidechain-heavy'. Provided keyword was [%s]"
                % (mode)
            )

        # get SSProtein objects for the two IDs passed (could be the same)
        try:
            P1 = self.proteinTrajectoryList[proteinID1]
            P2 = self.proteinTrajectoryList[proteinID2]

        except IndexError as e:
            raise SSException(
                "In get_interchain_distance(): When selecting protein indices %i and %i at least one of these was out of range (indices are from 0...%i)"
                % (proteinID1, proteinID2, len(self.proteinTrajectoryList) - 1)
            )

        # next build a new trajectory that contains ONLY the two residues selected
        local_atoms1 = P1.topology.select("resid %i" % (R1))
        local_atoms2 = P2.topology.select("resid %i" % (R2))

        if len(local_atoms1) == 0:
            raise SSException(
                "In get_interchain_distance(): When selecting resid %i from proteinID1 found no atoms"
                % (R1)
            )

        if len(local_atoms2) == 0:
            raise SSException(
                "In get_interchain_distance(): When selecting resid %i from proteinID2 found no atoms"
                % (R2)
            )

        subtraj_p1 = P1.traj.atom_slice(local_atoms1)
        subtraj_p2 = P2.traj.atom_slice(local_atoms2)

        # Apply frame stride before any expensive mdtraj compute so the
        # downstream COM / compute_contacts / compute_distances calls operate
        # on n_frames // stride frames rather than the full trajectory.
        if not (isinstance(stride, (int, np.integer)) and stride >= 1):
            raise SSException(
                f"In get_interchain_distance(): stride must be a positive "
                f"integer, got {stride!r}"
            )
        if stride > 1:
            subtraj_p1 = subtraj_p1[::stride]
            subtraj_p2 = subtraj_p2[::stride]

        # this is now a subtrajectory which in principle contains just two residues. We can check
        # this to ensure that the trajectory has exactly 2 residues
        full_subtraj = subtraj_p1.stack(subtraj_p2)
        full_subtraj_residues = [i for i in full_subtraj.topology.residues]
        if len(full_subtraj_residues) != 2:
            raise SSException(
                "In get_interchain_distance(): When passed in two residues (R1=%i, R2=%i) in proteins %i and %i found multiple residues (%i)...these resids could not be found "
                % (R1, R2, proteinID1, proteinID2, len(full_subtraj_residues))
            )

        # if we're looking at a specific pair of atoms (note we use resid 0 and 1 because we KNOW this trajectory only has 2 residues and we know R1 is 0 and R2 is 1
        if mode == "atom":
            atom1 = full_subtraj.topology.select('resid 0 and name "%s"' % A1)
            if len(atom1) != 1:
                raise SSException(
                    "In get_interchain_distance() when selecting atom %s from residue %i in protein %i no atoms were found "
                    % (A1, R1, proteinID1)
                )

            COM_1 = 10 * md.compute_center_of_mass(full_subtraj.atom_slice(atom1))

            atom2 = full_subtraj.topology.select('resid 1 and name "%s"' % A2)
            if len(atom2) != 1:
                raise SSException(
                    "In get_interchain_distance() when selecting atom %s from residue %i in protein %i no atoms were found "
                    % (A2, R2, proteinID2)
                )

            COM_2 = 10 * md.compute_center_of_mass(full_subtraj.atom_slice(atom2))

            # finally compute distances. Use minimum image convention if the periodic keyword is passed
            if periodic:
                distances = sstools.get_distance_periodic(
                    COM_1, COM_2, self.unitcell[0], "cube"
                )

            else:
                # revised way
                distances = np.linalg.norm(COM_1 - COM_2, axis=1)

                # old way
                # distances = np.sqrt(np.square(np.transpose(COM_1)[0] - np.transpose(COM_2)[0]) + np.square(np.transpose(COM_1)[1] - np.transpose(COM_2)[1])+np.square(np.transpose(COM_1)[2] - np.transpose(COM_2)[2]))

        else:
            # TODO: Documentation missing!
            # use the compute_contacts() function from mdtraj, multiplying by 10 because this will
            # by default give you numbers that...
            distances = (
                10
                * md.compute_contacts(
                    full_subtraj, [[0, 1]], scheme=mode, periodic=periodic
                )[0].ravel()
            )

        return distances


def __load_trajectory(trj_filename, top_filename, **kwargs):
    """Private helper function for loading trajectories in parallel

    Parameters
    ----------
    trj_filename : str
        Filename of trajectory to be loaded by SSTrajectory
    top_filenames : str
        Topology filename to be used for loading the trajectory
    **kwargs: dict, optional
        Key value pairs to be passed directly to SSTrajectory.

    Returns
    -------
    SSTrajectory
        Returns an SSTrajectory object for the given trajectory and topology.
    """

    try:
        return SSTrajectory(trj_filename, pdb_filename=top_filename, **kwargs)
    except ValueError as e:
        raise RuntimeError(
            f"Failed to load SSTrajectory for trajectory '{trj_filename}' with topology '{top_filename}'. "
            f"Perhaps you supplied the wrong topology?: {e}"
        )


def parallel_load_trjs(trj_filenames, top_filenames, n_procs=None, **kwargs):
    """
    Parallel loading of trajectories with optional kwargs.

    Parameters
    ----------
    trj_filenames : list of str
        A list of strings containing the trajectory file paths to be loaded.
    top_filenames : list of str
        A list of strings containing the topology file paths corresponding to the trajectories.
    n_procs : int, optional
        Number of separate processors to use for loading, by default None.
        If None, it will use the number of available CPU cores.
    **kwargs: dict, optional
        Key value pairs to be passed directly to SSTrajectory.

    Returns
    -------
    list
        Returns a list of SSTrajectory objects.
    """
    if n_procs is None:
        n_procs = cpu_count()

    # Normalize the topology argument. Accept either a single shared
    # topology (str) or a per-trajectory list. A single string -- or a
    # one-element list -- is broadcast across all trajectories.
    if isinstance(top_filenames, str):
        top_filenames = [top_filenames] * len(trj_filenames)
    elif len(top_filenames) == 1 and len(trj_filenames) > 1:
        top_filenames = list(top_filenames) * len(trj_filenames)

    if len(top_filenames) != len(trj_filenames):
        raise SSException(
            f"parallel_load_trjs(): number of topology files "
            f"({len(top_filenames)}) does not match number of "
            f"trajectory files ({len(trj_filenames)})"
        )

    # Partially apply load_trajectory with topology file path and kwargs
    partial_load = partial(__load_trajectory, **kwargs)

    # Parallelize load with **kwargs
    with Pool(processes=n_procs) as pool:
        trjs = pool.starmap(partial_load, zip(trj_filenames, top_filenames))

    return trjs
