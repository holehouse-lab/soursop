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

import os

import mdtraj as md
import numpy as np
from .ssdata import ALL_VALID_RESIDUE_NAMES

from .ssprotein import SSProtein
from .ssexceptions import SSException, SSWarning
from . import ssutils
from . import ssio

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
        check_whole_molecules=True,
    ):
        """
        SSTrajectory is the class used to read in a work with simulation
        trajectories in SOURSOP.

        **Input contract: whole molecules, no periodic-boundary corrections.**
        SOURSOP computes every distance, size and contact directly from the
        coordinates it is given. It never applies the minimum-image
        convention (the only exception is the opt-in ``periodic=True``
        flag on the inter-chain methods of this class) and it never
        re-images or unwraps coordinates. A chain that has been wrapped
        back into the primary cell by the simulation engine, so that part
        of it sits on the far side of the box, will therefore give wrong
        radii of gyration, distance maps, contact maps and so on with no
        error. Make molecules whole before loading (e.g. ``gmx trjconv -pbc
        mol -center``, or in Python ``traj = traj.make_molecules_whole()``
        / ``traj = traj.image_molecules()`` on the mdtraj trajectory; both
        return a new trajectory unless ``inplace=True`` is passed). On
        load, SOURSOP checks every protein chain for the tell-tale
        signature of a wrapped molecule - two bonded atoms (or, without
        bond information, two consecutive CA atoms) further apart than half
        the shortest box vector - and warns (or raises, see
        ``check_whole_molecules``) if it finds one; the same check is
        available at any time as :meth:`check_molecules_whole`.

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
            True if any protein was detected (or forced) to be a two-bead
            (CA backbone / CB sidechain) coarse-grained model. Detection is
            done per protein, on that protein's own atoms, so ions, ligands
            or solvent elsewhere in the system do not affect it; check
            ``SSProtein.is_swan`` for an individual chain. Two-bead chains
            switch sidechain-vector and secondary-structure analyses to the
            two-bead CA/CB definitions.


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
            Lets you manually define protein groups to be considered
            independently: each inner list holds the (0-indexed) residue
            indices of one protein. The groups are validated on load - each
            must be a non-empty iterable of integer residue indices in
            strictly increasing order (no duplicates), and no residue may
            appear in more than one group; violations raise an
            ``SSException``. Note that prior to 2.0.4 malformed groups were
            accepted and silently reordered (or duplicated across proteins)
            by the underlying topology selection. Default = None

        pdblead : bool
            Lets you set the PDB file (which is normally ONLY used as a topology
            file) to be the first frame of the trajectory. This is useful when
            the first PDB file holds some specific reference information which
            you want to use (e.g. RMSD or Q). Default = False

        debug : book
            Prints warning/help information to help debug weird stuff during
            initial  trajectory read-in. Default = False.

        extra_valid_residue_names : str or list of str
            By default, SOURSOP identifies chains as proteins based on the a set
            of normally seen protein residue names. These are defined in
            soursop/ssdata, and are listed below::

                'ALA', 'CYS', 'ASP', 'ASH', 'GLU', 'GLH', 'PHE', 'GLY',
                'HIE', 'HIS', 'HID', 'HIP', 'HSD', 'HSE', 'HSP', 'CYX',
                'CYM', 'ILE', 'LEU', 'LYS', 'LYD', 'LYN', 'MET', 'ASN',
                'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR',
                'AIB', 'ABA', 'NVA', 'NLE', 'ORN', 'DAB', 'PTR', 'TPO',
                'SEP', 'KAC', 'KM1', 'KM2', 'KM3', 'ACE', 'NME', 'FOR',
                'NH2'

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
            Controls two-bead (CA/CB) coarse-grained handling. If None
            (default) SOURSOP auto-detects, separately for each protein,
            whether its topology is a two-bead model (a single CA per residue
            plus a single CB per non-glycine residue, and nothing else). Pass
            True or False to force the behaviour for every protein and skip
            auto-detection. ``.swan_trajectory`` is then True if any protein
            is treated as two-bead.
            Default = None

        check_whole_molecules : bool or str
            Whether to test each protein chain for being split across the
            periodic boundary once it has been loaded (see
            :meth:`check_molecules_whole`). ``True`` (default) runs the
            check and emits an ``SSWarning`` naming the chain, the residue
            pair and the frame count if a chain is split; ``'raise'`` runs
            the check and raises an ``SSException`` instead; ``False``
            skips it. The check needs a unit cell; trajectories without one
            (or with zero-length box vectors) are not tested.

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
            # a bare string is almost certainly a single residue name; wrap it
            # rather than letting list.extend() split it into characters
            if isinstance(extra_valid_residue_names, str):
                extra_valid_residue_names = [extra_valid_residue_names]

            try:
                extra_names = list(extra_valid_residue_names)
            except TypeError:
                raise SSException(
                    "extra_valid_residue_names must be a string or an iterable of strings (got %s)"
                    % type(extra_valid_residue_names).__name__
                )

            for name in extra_names:
                if not isinstance(name, str):
                    raise SSException(
                        "extra_valid_residue_names must contain only strings (found %s)"
                        % type(name).__name__
                    )

            self.valid_residue_names.extend(extra_names)

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

        # determine whether this is a two-bead (CA/CB) coarse-grained model.
        # If swan_trajectory is forced (True/False) that value is threaded into
        # every SSProtein. If it is None each SSProtein auto-detects from its
        # own (protein-only) topology. Prior to 2.0.6 detection ran once on
        # the whole system, so a single ion, ligand or solvent bead switched
        # the two-bead paths off for every chain and DSSP then reported the
        # backbone-less chains as 100% coil.
        if swan_trajectory is None:
            self.__swan_request = None
        else:
            self.__swan_request = bool(swan_trajectory)

        # atom indices (in self.traj) of each protein, filled in when the
        # proteins are extracted below
        self.__protein_atom_indices = None

        # Next, having read in the trajectory we parse out into proteins
        # extract a list of protein trajectories where each protein is assumed
        # to be in its own chain
        if protein_grouping is None:
            self.proteinTrajectoryList = self.__get_proteins(
                self.traj, debug, explicit_residue_checking=explicit_residue_checking
            )
        else:
            self.proteinTrajectoryList = self.__get_proteins_by_residue(
                self.traj, protein_grouping, debug
            )

        # the public flag is True if any protein was forced or detected to be a
        # two-bead model
        if self.__swan_request is None:
            self.swan_trajectory = any(P.is_swan for P in self.proteinTrajectoryList)
        else:
            self.swan_trajectory = self.__swan_request

        # this is initialized to None and then gets defined using the lazy_loading_single_protein_trajectory()
        # decorator when it's first needed
        self.__single_protein_traj = None

        # SOURSOP applies no periodic-boundary corrections anywhere, so a
        # chain that the simulation engine wrapped back into the box gives
        # silently wrong results. Look for the signature of that up front.
        if check_whole_molecules is not False:
            if check_whole_molecules not in (True, "raise"):
                raise SSException(
                    "check_whole_molecules must be True, False or 'raise'; "
                    f"received {check_whole_molecules!r}"
                )
            self.__report_split_molecules(
                self.check_molecules_whole(),
                raise_on_split=(check_whole_molecules == "raise"),
            )

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
    @staticmethod
    def __same_pdb_file(trajectory_filename, pdb_filename):
        """Is the trajectory file the same PDB file as the topology file?

        Returns True only when both names resolve to the same file on disk AND
        that file is a PDB. This is the (common) case of loading a multi-model
        PDB - e.g. the output of idpconformergenerator, PED downloads, or any
        NMR-style multi-MODEL entry - where the trajectory is its own topology.

        Parameters
        ----------
        trajectory_filename : str
            Filename of the trajectory file.

        pdb_filename : str
            Filename of the PDB (topology) file.

        Returns
        -------
        bool
            True if both arguments name the same on-disk PDB file.
        """

        if trajectory_filename is None or pdb_filename is None:
            return False

        if not str(trajectory_filename).lower().endswith(".pdb"):
            return False

        # compare resolved paths so './x.pdb' and 'x.pdb' (or a symlink) match.
        # If either path cannot be resolved just fall back to the safe default
        # of False, which preserves the original md.load() behaviour.
        try:
            traj_path = os.path.realpath(trajectory_filename)
            pdb_path = os.path.realpath(pdb_filename)
        except (OSError, ValueError, TypeError):
            return False

        return traj_path == pdb_path

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

        # Is the trajectory file its own topology (i.e. a multi-model PDB passed
        # as both arguments)? If so we take a faster path - see below.
        same_pdb = self.__same_pdb_file(trajectory_filename, pdb_filename)

        # straight up read the trajectory first using mdtraj's awesome
        # trajectory reading facility.
        #
        # NOTE: when the trajectory file IS the topology file and both are a
        # PDB, we call md.load_pdb() directly rather than md.load(..., top=...).
        # md.load() pre-parses the topology via _parse_topology(top), which for
        # a PDB calls load_frame(top, 0); mdtraj's PDB reader has no random
        # access, so it reads EVERY model and then discards all but the first.
        # The file is then parsed a second time for the coordinates. For a large
        # multi-model PDB this doubles an already-expensive pure-Python text
        # parse (~2x wall-clock; e.g. 59 s -> 31 s for a 520 MB, 1000-model,
        # 6.4M-atom file). load_pdb() parses once and returns an identical
        # trajectory. Every other combination of inputs is unaffected.
        if same_pdb:
            traj = md.load_pdb(trajectory_filename)
        else:
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
                    "UnitCell lengths were not provided... This may cause issues but we're going to assume everything is OK for now..."
                )

        # if pdbLead is true then load the pdb_filename as a trajectory
        # and then add it to the front (the PDB file is its own topology
        # file so no need to specificy the top= file here!
        if pdblead:
            # If the PDB *is* the trajectory file we have already parsed it, so
            # reuse that result rather than paying for a third full parse of the
            # same file. md.load(pdb_filename) on that file returns exactly the
            # trajectory we just built, so this is behaviour-preserving.
            #
            # Only the first model is prepended: the docs promise the PDB
            # becomes "the first frame", but prior to 2.0.6 every model of a
            # multi-model PDB was prepended (duplicating the whole trajectory
            # when the PDB was also the trajectory file)
            if same_pdb:
                pdbtraj = traj[0]
            else:
                pdbtraj = md.load_frame(pdb_filename, 0)
            try:
                traj = pdbtraj + traj
            except ValueError as e:
                raise SSException(
                    f"pdblead=True: could not prepend the PDB file to the trajectory ({e}). This usually means one of the two files has periodic box information and the other does not; add a CRYST1 record to the PDB file (or remove the box from the trajectory), or load with pdblead=False."
                ) from None

            # having added the PDB file now check all the unit-cells match up!
            try:
                uc_lengths_1 = traj.unitcell_lengths[0]
                uc_lengths_2 = traj.unitcell_lengths[1]

                # PDB CRYST1 records store box lengths to 0.001 A while
                # xtc/dcd files store float32 nm, so compare with a tolerance
                # rather than exactly (an exact compare warned on every file)
                if print_warnings and not np.allclose(
                    uc_lengths_1, uc_lengths_2, rtol=0, atol=1e-3
                ):
                    ssio.warning_message(
                        f"The unit cell dimensions of the PDB file and trajectory file did not match, specifically\nPDB file = [ {uc_lengths_1} ]\nXTC file = [ {uc_lengths_2} ]\nThis may cause issues if native MDTraj utilities are used (and potentially for SOURSOP utilities that are based on these). It is not necessarily an issue, but PLEASE sanity check your outcome. To be safe we recommend editing the PDB-file unitcell dimensions to match."
                    )

            except TypeError:
                if print_warnings:
                    ssio.warning_message(
                        "UnitCell lengths were not provided... This may cause issues but we're going to assume everything is OK for now..."
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

        # Use exactly the atoms of the proteins in proteinTrajectoryList.
        # Prior to 2.0.6 this re-ran chain detection, so with
        # protein_grouping the overall_* observables silently included
        # residues outside every group (or failed outright)
        if self.__protein_atom_indices is not None:
            protein_atoms = sorted(
                int(a) for group in self.__protein_atom_indices for a in group
            )
            return SSProtein(
                trajectory.atom_slice(protein_atoms), swan=self.__swan_request
            )

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

        return SSProtein(trajectory.atom_slice(protein_atoms), swan=self.__swan_request)

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

                # only keep chains that actually contained protein residues -
                # a chain with none (e.g. a solvent chain, the very case this
                # flag exists for) would otherwise contribute an empty atom
                # list, and atom_slice([]) builds an empty topology that
                # crashes downstream residue lookups
                if len(local_atoms) > 0:
                    chainAtoms.append(local_atoms)
                elif debug:
                    ssio.debug_message(
                        f"Skipping {chain}: no valid protein residues found"
                    )

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
            proteinTrajectoryList.append(SSProtein(PT, swan=self.__swan_request))

        # remember which atoms of the full system make up the proteins so the
        # overall_* observables use exactly the same set
        self.__protein_atom_indices = [list(a) for a in chainAtoms]

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
            Must be a list containing one or more non-empty lists, where each
            internal list contains a set of strictly increasing (no
            duplicates) integer residue indices, and no residue index may
            appear in more than one list. In other words, each sub-list
            defines a single protein. Violations of any of these rules raise
            an ``SSException`` (previously such input was silently reordered
            or duplicated by the topology selector). The integer indexing
            here - importantly - uses the CAMPARITraj internal residue
            indexing, meaning that indexing begins at 0 from the first
            residue in the PDB file.

        Returns
        ---------
        list
            Returns a proteinTrajectoryList - a list with one or more SSProtein
            objects in it.

        """

        # extract full system topology
        topology = trajectory.topology

        # A topology selection returns atoms in topology order, regardless of
        # the order in the selection string.  Previously this meant malformed
        # groups such as [2, 1, 0] appeared to work while being silently
        # reordered; overlapping groups also duplicated residues into multiple
        # SSProtein objects.  Validate the documented contract before selecting
        # any atoms so invalid input cannot produce plausible-but-wrong output.
        try:
            normalized_groups = [list(group) for group in residue_grouping]
        except TypeError as exc:
            raise SSException(
                "protein_grouping must be an iterable of residue-index iterables"
            ) from exc

        seen_residues = set()
        for group_index, group in enumerate(normalized_groups):
            if len(group) == 0:
                raise SSException(f"protein_grouping group {group_index} is empty")

            normalized = []
            for resid in group:
                if isinstance(resid, bool) or not isinstance(resid, (int, np.integer)):
                    raise SSException(
                        f"protein_grouping group {group_index} contains a non-integer "
                        f"residue index: {resid!r}"
                    )
                # out-of-range indices used to be dropped silently (or, if a
                # whole group was out of range, the group was skipped and
                # every later protein index shifted); negative ones reached
                # mdtraj's selection parser
                if not (0 <= int(resid) < topology.n_residues):
                    raise SSException(
                        f"protein_grouping group {group_index} contains residue index "
                        f"{resid}, but valid indices are 0..{topology.n_residues - 1}"
                    )
                normalized.append(int(resid))

            if any(a >= b for a, b in zip(normalized, normalized[1:])):
                raise SSException(
                    f"protein_grouping group {group_index} must contain strictly "
                    f"increasing, non-duplicate residue indices; received {group!r}"
                )

            overlap = seen_residues.intersection(normalized)
            if overlap:
                raise SSException(
                    "protein_grouping groups must not overlap; residue indices "
                    f"{sorted(overlap)} occur in more than one group"
                )

            seen_residues.update(normalized)
            normalized_groups[group_index] = normalized

        group_atoms = []

        # for each chain in this toplogy determine if the
        # first residue is protein or not
        for group in normalized_groups:
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

            proteinTrajectoryList.append(SSProtein(PT, swan=self.__swan_request))

        # remember which atoms of the full system make up the proteins so the
        # overall_* observables use exactly the same set
        self.__protein_atom_indices = [list(a) for a in group_atoms]

        if len(proteinTrajectoryList) == 0:
            ssio.warning_message("No protein chains found in the trajectory")

        return proteinTrajectoryList

    # oxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxo
    #
    #
    def check_molecules_whole(self, chunk_size=1000):
        """Test every protein chain for being split across the periodic boundary.

        SOURSOP never applies periodic-boundary corrections, so it relies on
        the molecules it is given being whole. A chain (or part of one) that
        a simulation engine has wrapped back into the primary cell leaves a
        clear signature: two covalently bonded atoms sit on opposite sides
        of the box, so they are separated by roughly a box vector rather
        than the 1-2 A of a bond. This method looks for exactly that - a
        bonded-atom distance greater than half the shortest box vector of
        that frame - in every chain and every frame.

        The bonds tested are those in the topology that join two atoms of
        the same residue, or of two sequence-adjacent residues of the same
        chain (consecutive residue indices whose PDB residue numbers also
        differ by one). This catches a split sidechain or cap as well as a
        split backbone, and never pairs residues that were never bonded
        (e.g. the end of one chain and the start of the next when
        ``protein_grouping`` or a ``.gro`` file puts them side by side). If
        the topology has no bonds (e.g. one-bead coarse-grained models read
        without CONECT records), consecutive CA atoms that pass the same
        adjacency test are used instead.

        A whole chain that is simply larger than the box is *not* flagged
        (its bonded atoms remain bonded distances apart), and neither is a
        genuine chain break shorter than half the box, which is not a
        periodic-boundary problem.

        The check needs a unit cell. Trajectories without one, or whose box
        vectors are zero or non-finite (old CAMPARI files), cannot be tested
        and are reported as such rather than as whole.

        Parameters
        ----------
        chunk_size : int, optional
            Number of frames processed at a time, to bound memory on very
            long trajectories. Must be a positive integer. Default 1000.

        Returns
        -------
        list of dict
            One entry per protein in ``proteinTrajectoryList`` with the keys
            ``'protein'`` (its index), ``'tested'`` (False when there is no
            usable unit cell or nothing to test), ``'split'`` (True if any
            tested pair exceeds half the box in any frame),
            ``'n_frames_split'``, ``'worst_frame'``, ``'worst_pair'`` (the
            resids, in the chain's own 0-based numbering, of the two atoms
            in the most stretched pair; both entries are the same resid for
            an intra-residue bond), ``'worst_atoms'`` (the two atom names),
            ``'max_distance'`` (the largest tested distance, Angstroms) and
            ``'half_box'`` (the smallest half box vector over the trajectory,
            Angstroms, or ``None``).

        Raises
        ------
        SSException
            If ``chunk_size`` is not a positive integer.

        Example
        -------
        >>> report = traj.check_molecules_whole()
        >>> any(r["split"] for r in report)
        False
        """

        if (
            isinstance(chunk_size, bool)
            or not isinstance(chunk_size, (int, np.integer))
            or chunk_size < 1
        ):
            raise SSException(
                f"chunk_size must be a positive integer; received {chunk_size!r}"
            )
        chunk_size = int(chunk_size)

        lengths = self.traj.unitcell_lengths
        usable_box = (
            lengths is not None and np.all(np.isfinite(lengths)) and np.all(lengths > 0)
        )
        # half the shortest box vector in every frame, in Angstroms
        half_box = 0.5 * 10.0 * lengths.min(axis=1) if usable_box else None

        def _bonded_neighbours(res_a, res_b):
            # same residue, or sequence neighbours in the same chain. The
            # resSeq test stops a .gro file (one chain for everything) or a
            # protein_grouping group that spans chains from pairing the end
            # of one chain with the start of the next
            if res_a.index == res_b.index:
                return True
            if res_a.chain.index != res_b.chain.index:
                return False
            return (
                abs(res_a.index - res_b.index) == 1
                and abs(res_a.resSeq - res_b.resSeq) == 1
            )

        report = []
        for k, protein in enumerate(self.proteinTrajectoryList):
            entry = {
                "protein": k,
                "tested": False,
                "split": False,
                "n_frames_split": 0,
                "worst_frame": None,
                "worst_pair": None,
                "worst_atoms": None,
                "max_distance": None,
                "half_box": None if half_box is None else float(np.min(half_box)),
            }

            top = protein.topology
            if top.n_bonds > 0:
                # every intra-residue and sequence-adjacent bond, so a split
                # sidechain or cap is caught as well as a split backbone
                atom_pairs = [
                    [a.index, b.index]
                    for a, b in top.bonds
                    if _bonded_neighbours(a.residue, b.residue)
                ]
            else:
                # no bonds recorded: fall back to consecutive CA atoms
                resids = protein.resid_with_CA
                atom_pairs = [
                    [protein.get_CA_index(a), protein.get_CA_index(b)]
                    for a, b in zip(resids[:-1], resids[1:])
                    if _bonded_neighbours(top.residue(a), top.residue(b))
                ]

            if not usable_box or len(atom_pairs) == 0:
                report.append(entry)
                continue

            atom_pairs = np.array(atom_pairs)

            entry["tested"] = True
            n_frames = protein.n_frames
            frames_split = np.zeros(n_frames, dtype=bool)
            max_d, worst = -1.0, (None, None)
            for start in range(0, n_frames, chunk_size):
                stop = min(start + chunk_size, n_frames)
                d = 10.0 * md.compute_distances(
                    protein.traj[start:stop], atom_pairs, periodic=False
                )
                over = d > half_box[start:stop, None]
                frames_split[start:stop] = over.any(axis=1)
                local_max = d.max()
                if local_max > max_d:
                    max_d = float(local_max)
                    f, p = np.unravel_index(np.argmax(d), d.shape)
                    worst = (start + int(f), int(p))

            atom_a = top.atom(int(atom_pairs[worst[1], 0]))
            atom_b = top.atom(int(atom_pairs[worst[1], 1]))
            # report the pair in sequence order
            if atom_b.residue.index < atom_a.residue.index:
                atom_a, atom_b = atom_b, atom_a

            entry["split"] = bool(frames_split.any())
            entry["n_frames_split"] = int(frames_split.sum())
            entry["max_distance"] = max_d
            entry["worst_frame"] = worst[0]
            entry["worst_pair"] = (atom_a.residue.index, atom_b.residue.index)
            entry["worst_atoms"] = (atom_a.name, atom_b.name)
            report.append(entry)

        return report

    # oxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxo
    #
    #
    def __report_split_molecules(self, report, raise_on_split=False):
        """Warn (or raise) for every chain that check_molecules_whole flagged."""
        split = [r for r in report if r["split"]]
        if len(split) == 0:
            return

        lines = []
        for r in split:
            a, b = r["worst_pair"]
            name_a, name_b = r["worst_atoms"]
            where = f"residue {a}" if a == b else f"residues {a}-{b}"
            lines.append(
                f"protein {r['protein']}: bonded atoms {name_a} and {name_b} ({where}) "
                f"are {r['max_distance']:.1f} A apart in frame {r['worst_frame']} (half "
                f"the shortest box vector is {r['half_box']:.1f} A); "
                f"{r['n_frames_split']} of {self.n_frames} frames affected"
            )
        message = (
            "Protein chain(s) appear to be split across the periodic boundary - "
            + "; ".join(lines)
            + ". SOURSOP applies no periodic-boundary corrections and expects whole "
            "molecules, so sizes, distances and contacts computed from this trajectory "
            "will be wrong. Make the molecules whole before loading (e.g. "
            "'gmx trjconv -pbc mol -center', or traj = traj.make_molecules_whole() / "
            "traj = traj.image_molecules() in mdtraj - both return a new trajectory "
            "rather than editing in place), or pass check_whole_molecules=False "
            "to skip this check."
        )
        if raise_on_split:
            raise SSException(message)
        SSWarning(message)

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

        return self.__single_protein_traj.get_asphericity(
            verbose=False, weights=weights, etol=etol
        )

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
        intra-chain distance map of protein ``i``. The per-pair values
        match ``self.proteinTrajectoryList[i].get_distance_map()``, but
        note that this function returns a full symmetric matrix whereas
        :meth:`SSProtein.get_distance_map` returns an upper-triangular
        one (zeros below the diagonal), so compare against the upper
        triangle when using this as a sanity check.

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
            If True, apply the minimum-image convention using each frame's
            own box. Requires a recorded rectangular (orthorhombic) periodic
            box; triclinic boxes raise. Generally it is better to make the
            system whole and centred first and leave this False. Default
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
        P1 = self.__validated_protein(
            proteinID1, "proteinID1", "get_interchain_distance_map"
        )
        P2 = self.__validated_protein(
            proteinID2, "proteinID2", "get_interchain_distance_map"
        )

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

        # row/column k of the map is the k-th CA-bearing residue of each
        # chain (the same convention as SSProtein.get_distance_map). This
        # used to be derived as ``resid - 1 if ncap else resid``, which is
        # only equivalent when every residue after the cap has a CA; a chain
        # with any other CA-less residue indexed past the end of the map.
        p1_indices = list(range(len(p1_residues)))
        p2_indices = list(range(len(p2_residues)))

        # minimum image in a rectangular box, using each frame's own box
        # lengths. Prior to 2.0.6 this assumed a cube with frame 0's x length
        # and only handled separations below 1.5 box lengths, so slab boxes,
        # NPT boxes and chains that had drifted out of the primary cell all
        # gave silently wrong distances
        if periodic:
            box = self.__orthorhombic_box_lengths("get_interchain_distance_map")

        # com2 stacked once -> (n2, F, 3); broadcasting COM_1 (F, 3)
        # against it reproduces the per-pair np.linalg.norm exactly.
        com2_stack = np.stack(com2, axis=0)
        for i, COM_1 in enumerate(com1):
            delta = COM_1 - com2_stack  # (n2, F, 3)
            if periodic:
                delta = delta - box * np.round(delta / box)
            d = np.linalg.norm(delta, axis=-1)  # (n2, F)
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
    def __validated_protein(self, proteinID, label, function_name):
        """Return ``proteinTrajectoryList[proteinID]`` after validating the ID.

        Negative indices used to wrap silently to the end of the list, and
        an out-of-range index raised a bare ``IndexError`` in some methods.

        Parameters
        ----------
        proteinID : int
            Index into ``self.proteinTrajectoryList``.
        label : str
            Name of the argument, used in the error message.
        function_name : str
            Name of the calling method, used in the error message.

        Returns
        -------
        SSProtein
            The selected protein.

        Raises
        ------
        SSException
            If ``proteinID`` is not an integer in
            ``[0, len(proteinTrajectoryList) - 1]``.
        """
        n_proteins = len(self.proteinTrajectoryList)
        if isinstance(proteinID, bool) or not isinstance(proteinID, (int, np.integer)):
            raise SSException(
                f"In {function_name}(): {label} must be an integer; received {proteinID!r}"
            )
        if not (0 <= proteinID < n_proteins):
            raise SSException(
                f"In {function_name}(): {label}={proteinID} is out of range; valid indices are 0..{n_proteins - 1}"
            )
        return self.proteinTrajectoryList[proteinID]

    # oxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxo
    #
    #
    def __orthorhombic_box_lengths(self, function_name):
        """Per-frame box lengths (Angstroms) for a rectangular periodic box.

        Parameters
        ----------
        function_name : str
            Name of the calling method, used in error messages.

        Returns
        -------
        np.ndarray
            Array of shape ``(n_frames, 3)`` holding each frame's
            ``(Lx, Ly, Lz)`` in Angstroms.

        Raises
        ------
        SSException
            If the trajectory has no unit cell, if any box length is not
            positive, or if the box is not rectangular (any angle differs
            from 90 degrees by more than 0.001 degrees).
        """
        lengths = self.traj.unitcell_lengths
        angles = self.traj.unitcell_angles
        if lengths is None or angles is None:
            raise SSException(
                f"{function_name}(): periodic=True needs a periodic box, but this trajectory has no unit cell information."
            )
        if not np.allclose(angles, 90.0, rtol=0, atol=1e-3):
            raise SSException(
                f"{function_name}(): periodic=True only supports rectangular (orthorhombic) boxes, but this trajectory has box angles {angles[0]}. Make the molecules whole and use periodic=False instead."
            )
        lengths = 10.0 * np.asarray(lengths, dtype=np.float64)
        if np.any(lengths <= 0):
            raise SSException(
                f"{function_name}(): periodic=True needs positive box lengths, but at least one frame has a zero or negative box length."
            )
        return lengths

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

        Returns an ``(n_CA_P1, n_CA_P2)`` matrix where entry ``[i, j]`` is
        the fraction of (strided) frames in which the ``i``-th CA-bearing
        residue of chain ``proteinID1`` is within ``threshold`` Angstroms of
        the ``j``-th CA-bearing residue of chain ``proteinID2`` under the
        chosen distance mode. Calling with ``proteinID1 == proteinID2``
        produces the chain's intra-chain contact map.

        Rows and columns follow ``resid_with_CA`` of each chain, so caps
        (ACE, NME, FOR, NH2, ...) and any other residue without a CA atom are
        excluded - the same convention as :meth:`get_interchain_distance_map`
        and :meth:`SSProtein.get_contact_map`. Row ``k`` is therefore residue
        ``P.resid_with_CA[k]``, not resid ``k``. **Breaking change in
        2.0.6:** earlier versions indexed every residue by resid, caps
        included, returning an ``(n_residues, n_residues)`` matrix whose rows
        were offset by one from the other per-residue maps on a chain with an
        N-terminal cap.

        A residue for which the chosen mode cannot be evaluated (e.g. one
        lacking the atom named by ``A1``/``A2`` in ``mode='atom'``, or a
        coarse-grained bead with a sidechain mode) contributes a row/column
        of zeros. Note that for the ``'sidechain'`` modes mdtraj uses a
        glycine's HA atoms as its sidechain, so glycine rows are generally
        not zero.

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
            If True, apply the minimum-image convention (using each frame's
            own box) for every mode. Requires a recorded periodic box.
            Default False.
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
            Array of shape ``(n_CA_P1, n_CA_P2)`` of contact fractions in
            ``[0, 1]``, where ``n_CA`` is ``len(P.resid_with_CA)`` for each
            chain.

        Raises
        ------
        SSException
            If either ``proteinID`` is not a valid index, ``mode`` is
            invalid, ``stride`` is invalid, ``periodic=True`` is requested
            without a periodic box, either chain has no CA-bearing residues,
            or every per-residue inner call failed (e.g. an invalid atom name
            was supplied).

        Example
        -------
        >>> cmap = two_chain_traj.get_interchain_contact_map(0, 1, mode='closest-heavy')
        >>> cmap.shape
        (58, 58)
        """

        # Eagerly validate proteinIDs so an out-of-range value raises
        # SSException up front instead of an IndexError partway through.
        self.__validated_protein(proteinID1, "proteinID1", "get_interchain_contact_map")
        self.__validated_protein(proteinID2, "proteinID2", "get_interchain_contact_map")

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

        # the shared stride check (also rejects bools and strides longer than
        # the trajectory, which used to be accepted here)
        stride = ssutils.validate_stride(stride, self.n_frames)

        # check this up front; otherwise every per-pair call fails, is
        # zero-filled, and the error blames the mode/atom selection
        if periodic and self.traj.unitcell_vectors is None:
            raise SSException(
                "In get_interchain_contact_map(): periodic=True needs a periodic box, but this trajectory has no unit cell information."
            )

        P1 = self.proteinTrajectoryList[proteinID1]
        P2 = self.proteinTrajectoryList[proteinID2]

        # Rows/columns are the CA-bearing residues of each chain (caps and any
        # other CA-less residue excluded), matching get_interchain_distance_map
        # and SSProtein.get_contact_map. Prior to 2.0.6 every residue was
        # included and the map was indexed by resid, caps and all.
        res_P1 = list(P1.resid_with_CA)
        res_P2 = list(P2.resid_with_CA)
        for label, pid, res in (
            ("proteinID1", proteinID1, res_P1),
            ("proteinID2", proteinID2, res_P2),
        ):
            if len(res) == 0:
                raise SSException(
                    f"In get_interchain_contact_map(): protein {pid} ({label}) has no "
                    "CA-bearing residues, so the contact map has no rows/columns."
                )
        n_res_P1 = len(res_P1)
        n_res_P2 = len(res_P2)

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

            def _residue_atom_positions(P, residues, atom_name):
                # positions[k] = 10x-COM (F,3) of the named atom of the k-th
                # residue in `residues` over the strided frames; ok[k]
                # mirrors the exact success/failure of the original
                # get_interchain_distance() 'atom' selection for that
                # residue (the conditions the caller's except catches).
                positions = [None] * len(residues)
                ok = [False] * len(residues)
                for k, r in enumerate(residues):
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
                        positions[k] = 10 * md.compute_center_of_mass(sub.atom_slice(a))
                        ok[k] = True
                    except (SSException, ValueError, IndexError):
                        continue
                return positions, ok

            pos1, ok1 = _residue_atom_positions(P1, res_P1, A1)
            pos2, ok2 = _residue_atom_positions(P2, res_P2, A2)

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

        # cycle over each CA-bearing residue in protein 1
        for k1, p1_res_idx in enumerate(res_P1):
            if verbose:
                print(f"On {k1} of {n_res_P1}")

            tmp = []

            # cycle over each CA-bearing residue in protein 2
            for p2_res_idx in res_P2:
                # Per-residue computations can legitimately fail when a
                # residue lacks the atoms the mode needs: e.g. a residue
                # without a sidechain atom mdtraj recognises (such as a
                # coarse-grained bead) makes the sidechain modes raise. Note
                # mdtraj treats a glycine's HA atoms as its "sidechain", so
                # glycine rows are generally NOT zero. In the failing cases
                # the contact fraction is set to 0 - there is nothing to be in
                # contact via the requested mode - so the matrix shape
                # (n_CA_P1, n_CA_P2) is preserved.
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
            If True, apply the minimum-image convention using each frame's
            own box (triclinic boxes included), for every mode. Requires a
            recorded periodic box. Default False.
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
            If ``mode`` is invalid, a protein ID or resid is not an integer
            or is out of range, a resid has no atoms, an atom name cannot be
            found, mdtraj cannot compute the requested mode for the two
            residues (e.g. a sidechain mode on glycine or a coarse-grained
            bead), ``periodic=True`` is requested without a periodic box, or
            ``stride`` is invalid.

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
        P1 = self.__validated_protein(
            proteinID1, "proteinID1", "get_interchain_distance"
        )
        P2 = self.__validated_protein(
            proteinID2, "proteinID2", "get_interchain_distance"
        )

        # validate the residue indices; a float used to be truncated by the
        # "%i" selection string and a negative one reached mdtraj's parser
        for label, R, P in (("R1", R1, P1), ("R2", R2, P2)):
            if isinstance(R, bool) or not isinstance(R, (int, np.integer)):
                raise SSException(
                    f"In get_interchain_distance(): {label} must be an integer residue index; received {R!r}"
                )
            if not (0 <= R < P.n_residues):
                raise SSException(
                    f"In get_interchain_distance(): {label}={R} is out of range; valid indices are 0..{P.n_residues - 1}"
                )

        # the shared stride check (also rejects bools and strides longer than
        # the trajectory, which used to be accepted here)
        stride = ssutils.validate_stride(stride, self.n_frames)

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
                # let mdtraj apply the minimum image: it uses each frame's
                # own box and handles triclinic cells. Prior to 2.0.6 this
                # assumed a cube with frame 0's x length and only handled
                # separations below 1.5 box lengths
                if full_subtraj.unitcell_vectors is None:
                    raise SSException(
                        "In get_interchain_distance(): periodic=True needs a periodic box, but this trajectory has no unit cell information."
                    )
                distances = 10 * md.compute_distances(
                    full_subtraj, [[atom1[0], atom2[0]]], periodic=True
                )[:, 0].astype(float)

            else:
                # revised way
                distances = np.linalg.norm(COM_1 - COM_2, axis=1)

                # old way
                # distances = np.sqrt(np.square(np.transpose(COM_1)[0] - np.transpose(COM_2)[0]) + np.square(np.transpose(COM_1)[1] - np.transpose(COM_2)[1])+np.square(np.transpose(COM_1)[2] - np.transpose(COM_2)[2]))

        else:
            # use the compute_contacts() function from mdtraj (which works in
            # nm, hence the factor of 10). mdtraj raises a bare ValueError when
            # a residue lacks the atoms a scheme needs (e.g. sidechain modes on
            # glycine, caps or coarse-grained beads), so translate that
            if periodic and full_subtraj.unitcell_vectors is None:
                raise SSException(
                    "In get_interchain_distance(): periodic=True needs a periodic box, but this trajectory has no unit cell information."
                )
            try:
                distances = (
                    10
                    * md.compute_contacts(
                        full_subtraj, [[0, 1]], scheme=mode, periodic=periodic
                    )[0].ravel()
                )
            except ValueError as e:
                raise SSException(
                    f"In get_interchain_distance(): mdtraj could not compute a '{mode}' distance between residue {R1} of protein {proteinID1} ({full_subtraj_residues[0].name}) and residue {R2} of protein {proteinID2} ({full_subtraj_residues[1].name}); one of them probably lacks the atoms this mode needs (e.g. a sidechain). mdtraj said: {e}"
                ) from None

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
    trj_filenames : list of str or pathlib.Path
        The trajectory file paths to be loaded. A single path (str or
        ``pathlib.Path``) is treated as a one-element list.
    top_filenames : str, pathlib.Path, or list of these
        The topology file paths corresponding to the trajectories. A single
        path (or a one-element list) is used for every trajectory.
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

    # A bare path used to be iterated character by character (str) or
    # raise a TypeError (pathlib.Path), so wrap single paths in a list
    if isinstance(trj_filenames, (str, os.PathLike)):
        trj_filenames = [trj_filenames]
    trj_filenames = [os.fspath(f) for f in trj_filenames]

    # Normalize the topology argument. Accept either a single shared
    # topology (str / Path) or a per-trajectory list. A single path -- or a
    # one-element list -- is broadcast across all trajectories.
    if isinstance(top_filenames, (str, os.PathLike)):
        top_filenames = [os.fspath(top_filenames)] * len(trj_filenames)
    elif len(top_filenames) == 1 and len(trj_filenames) > 1:
        top_filenames = list(top_filenames) * len(trj_filenames)

    top_filenames = [os.fspath(f) for f in top_filenames]

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
