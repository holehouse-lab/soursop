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

import time
import mdtraj as md
import numpy as np
from numpy import linalg as LA
from itertools import combinations
from scipy import stats
from scipy.special import expit
import scipy.optimize as SPO
from scipy.spatial import ConvexHull
from numpy.random import choice
from .configs import DEBUGGING
from .ssdata import THREE_TO_ONE, DEFAULT_SIDECHAIN_VECTOR_ATOMS, ALL_VALID_RESIDUE_NAMES
from .ssexceptions import SSException
from . import ssmutualinformation, ssio, sstools, sspolymer, ssutils
#from .sstrajectory import SSTrajectory

from . _internal_data import BBSEG2

import scipy.cluster.hierarchy


## Order of standard args:
## 1. stride
## 2. weights
#  3. verbose
##
##


class SSProtein:

    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ##
    ## A note for the code:
    ## As of version 0.1.3 we assume that for each protein the first residue resid
    ## will ALWAYS index from 0 onwards. This is explicitly checked in SSTrajectory where
    ## SSProtein objects are built. This is a change from prior versions were an offset
    ## backend was built to give the illusion of resid = 0 indexinging, while on the
    ## level of the underlying mdtraj topology object this was not a given. In MDTraj
    ## 1.9.5 this has been resolved, allowing the codebase to become substantially
    ## simpler.
    ##
    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


    # ........................................................................
    #
    def __init__(self, traj, debug=DEBUGGING, check_one_bead_per_residue=True):
        """SSProtein objects are initialized with a trajectory subset that
        contains only the atoms a specific, single protein. This means that a
        SSProtein object allows operations to performed on a single protein. A
        single trajectory may have multiple proteins in it.

        Indexing with a protein object assumes that the protein is indexed from
        0 to `n`, `n` is the number of residues.

        **This is an important idea to emphasize - it means that you (the user)
        will need to determine the correct residue index for a region of interest
        being examined. The region will NOT (necessarily) correspond to the
        residue index in the PDB file used.**

        To make this easier the function ``SSProtein.print_residues()`` will
        print the mapping of residue index to residue name and residue number.

        To re-iterate:

        **Residue number is the number of the residue in the PDB file.**

        Residue index is the index value associated with a residue in a specific
        protein and will always begin from 0 - note this will include the peptide
        caps (ACE/NME) if present.

        Initialize a SSProtein object instance using trajectory information.

        Parameters
        ----------
        traj: `sstrajectory.SSTrajectory`
            An instance of a system's trajectory populated via
            `sstrajectory.SSTrajectory`.

        debug: bool {False}
            If set to `True` the code will print out debug messages to the
            screen. This is useful for debugging purposes.

        check_one_bead_per_residue: bool {True}
            If set to `True` the code will check to see if this passed trajectory
            makes sense as a 1 bead per residue (number of atoms = number of amino
            acids) and if yes initialize accordingly.

        """

        
        # This is necessary to support sstrajectory.Trajectory as well as the default `mdtraj`.
        # note that this is not super elegant - we'd rather use isinstance(), but this would
        # neceisstate a circular import which is not great so
        if str(type(traj)) == "<class 'soursop.sstrajectory.SSTrajectory'>":
            self.traj       = traj.traj
            self.topology   = traj.traj.topology

        elif isinstance(traj, md.core.trajectory.Trajectory):
            # set the trajectory object for easy access
            self.traj     = traj
            self.topology = traj.topology
        else:
            raise RuntimeError('The argument passed as `traj` is not a supported Trajectory object. Please use an mdtraj or SSTrajectory object.')

        if debug:
            ssio.debug_message("Creating protein")
            residues_strings = list()
            for r in self.topology.chain(0).residues:
                residues_strings.append(str(r))
            r_string = '-'.join(residues_strings)
            ssio.debug_message(f"Residue string from residues in self.topology.chain(0).residues: {r_string}")

            # delete the vaiable to avoid any possible introduction of this var into the namespace
            del r_string

        # initialze various protein-centric data
        self.__num_residues       = sum( 1 for _ in self.topology.residues)

        # remember the constructor's CG-detection preference so that
        # reset_cache() can reproduce an identical initialization later
        self.__check_one_bead_per_residue = check_one_bead_per_residue

        ## DEVELOPMENT NOTES
        # Everything set up by __initialize_memoization_and_topology() is the
        # state that reset_cache() discards and rebuilds.
        self.__initialize_memoization_and_topology()

        
        

    # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    #
    # IMPORTANT FUNCTION - enables the user to manually override stored
    def reset_cache(self):
        """Clear every value the SSProtein has memoised so far.

        Many SSProtein methods cache their first-call result (sequences, atom
        index tables, per-residue centres of mass, SASA results, dihedral
        angles). This method discards those caches and re-initialises the
        ``resid_with_CA`` / cap flags from the underlying topology. It is
        rarely necessary in normal use; reach for it if you have mutated the
        trajectory in-place and need the SSProtein to re-derive its tables,
        or if you are profiling cold-vs-warm performance.

        Returns
        -------
        None

        Example
        -------
        >>> protein.get_radius_of_gyration()    # first call: populates caches
        >>> protein.reset_cache()               # discard all cached values
        >>> protein.get_radius_of_gyration()    # cold-cache call again
        """

        self.__initialize_memoization_and_topology()


    # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    #
    def __initialize_memoization_and_topology(self):
        """Initialize (or re-initialize) every memoised lookup table and the
        topology-derived CA / cap state.

        Called once from ``__init__`` and again by ``reset_cache()``; routing
        both through this single method guarantees a freshly-constructed
        object and a reset object are in an identical state (previously
        ``reset_cache`` forced ``__cg_onechain`` to False, silently changing
        behaviour for coarse-grained chains).

        Returns
        -------
        None
        """

        # empty values populated on demand by functions that drive local memoization
        self.__amino_acids_3LTR   = None
        self.__amino_acids_1LTR   = None
        self.__residue_index_list = None
        self.__CA_residue_atom    = {}
        self.__residue_atom_table = {}
        self.__residue_COM        = {}
        self.__residue_atom_COM   = {}
        self.__SASA_saved         = {}
        self.__all_angles         = {} # different dihedral (backbone and sidechain) angles

        # determine whether __cg_onechain can be set True. This ONLY serves to
        # dramatically speed up initialization but has no effect beyond that.
        if self.__check_one_bead_per_residue is True:
            self.__cg_onechain = self.__check_cg_onebead()
        else:
            self.__cg_onechain = False

        # build the list of resids that have a CA atom
        self.__resid_with_CA = self.__get_resid_with_CA()

        # caps are inferred from resid 0 / the last resid lacking a CA
        self.__ncap = 0 not in self.resid_with_CA
        self.__ccap = (self.n_residues - 1) not in self.resid_with_CA



    # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    #
    # Properties


    @property
    def resid_with_CA(self):
        """Resids of residues that have a CA (alpha-carbon) atom.

        ACE and NME capping residues, water, ions, and any non-standard
        residue without a CA are excluded. Use this whenever you need to
        iterate over "real" residues and not the cap pseudo-residues.

        Returns
        -------
        list of int
            Resids (0-indexed) of residues with a CA atom, in ascending order.

        Example
        -------
        >>> protein.resid_with_CA
        [1, 2, 3, ..., 92]   # ctl9_AA has ACE at 0 and NME at 93
        """
        return self.__resid_with_CA


    @property
    def ncap(self):
        """True if an N-terminal capping residue (e.g. ACE) is present.

        Inferred from whether resid 0 has a CA atom: if it does not, a cap is
        assumed to occupy resid 0.

        Returns
        -------
        bool
            True if an N-terminal cap is present, False otherwise.

        Example
        -------
        >>> protein.ncap
        True
        """
        return self.__ncap

    @property
    def ccap(self):
        """True if a C-terminal capping residue (e.g. NME) is present.

        Inferred from whether the last resid has a CA atom: if it does not, a
        cap is assumed to occupy that resid.

        Returns
        -------
        bool
            True if a C-terminal cap is present, False otherwise.

        Example
        -------
        >>> protein.ccap
        True
        """

        return self.__ccap

    @property
    def n_frames(self):
        """Number of frames in the underlying mdtraj trajectory.

        Returns
        -------
        int
            Total frame count of the simulation trajectory.

        Example
        -------
        >>> protein.n_frames
        1000
        """

        return self.traj.n_frames

    @property
    def n_residues(self):
        """Number of residues in the protein, including any cap residues.

        Returns
        -------
        int
            Total residue count, counting ACE / NME caps if present.

        Example
        -------
        >>> protein.n_residues
        94
        """

        return self.__num_residues

    @property
    def residue_index_list(self):
        """All resids associated with this protein, extracted from the topology.

        As of v0.1.3 the list always starts at 0 and increments by one up to
        ``n_residues - 1``. It is extracted live from the underlying
        ``mdtraj.topology`` (and then cached) so it is useful as a
        ground-truth reference when debugging residue-index issues.

        Returns
        -------
        list of int
            Resids in ascending order; length equals ``n_residues``.

        Example
        -------
        >>> protein.residue_index_list[:5]
        [0, 1, 2, 3, 4]
        """

        if self.__residue_index_list is None:
            reslist = []
            for res in self.topology.residues:
                reslist.append(res.index)

            self.__residue_index_list = reslist
        return self.__residue_index_list


    @property
    def unitcell(self):
        """Unit cell box lengths of the first frame, in Angstroms.

        SOURSOP follows the convention of reporting distances in Angstroms, so
        the values stored by mdtraj (which are in nanometres) are multiplied
        by 10 here. Raises ``TypeError`` if the trajectory has no recorded
        periodic box (some coarse-grained runs).

        Returns
        -------
        np.ndarray
            Array of shape (3,) giving the box lengths along x, y, z in
            Angstroms.

        Example
        -------
        >>> protein.unitcell
        array([60.0, 60.0, 60.0])
        """
        return self.traj.unitcell_lengths[0]*10
    



    def  __repr__(self):
        return f"SSProtein ({hex(id(self))}): {self.n_residues} res and {self.n_frames} frames"


    def __len__(self):
        # Edited to mimic the behavior of `mdtraj` trajectory objects.
        # Originally: `return (self.n_residues, self.n_frames)`
        return self.n_frames


    def length(self):
        """Return a ``(n_residues, n_frames)`` tuple summarising the protein.

        Convenience method that preserves the historic behaviour of
        ``__len__`` from before SOURSOP switched its ``len(protein)``
        semantics to match mdtraj (frame count only).

        Returns
        -------
        tuple of (int, int)
            ``(n_residues, n_frames)`` — number of residues including caps,
            and number of trajectory frames.

        Example
        -------
        >>> protein.length()
        (94, 1000)
        """
        return (self.n_residues, self.n_frames)




    # ........................................................................
    #
    def __check_weights(self, weights, stride=1, etol=0.0000001):
        """Validate a frame-weights vector before it is used in averaging.

        Ensures the weights are a numeric numpy array, that there is one
        weight per frame, and that they sum to 1 within ``etol`` (after
        any stride subsampling). Type-casts to ``np.float64`` so the
        result can be used as numpy indexing or weighting input
        downstream.

        Parameters
        ----------
        weights : array_like or False
            Per-frame weights. The literal value ``False`` is treated as
            "no weights" and returned unchanged.
        stride : int, optional
            Stepsize used while iterating across frames. If ``> 1`` a
            warning is emitted and the strided slice of the weights is
            returned (since strided weighting is rarely what you want).
            Default 1.
        etol : float, optional
            Tolerance for the ``|sum(weights) - 1| < etol`` check.
            Default ``1e-7``.

        Returns
        -------
        np.ndarray or False
            * ``False`` if the input was ``False`` (no weights).
            * Otherwise a ``np.float64`` array of length
              ``ceil(n_frames / stride)``.
        """

        if weights is not False:
            # convert to an array of doubles
            try:
                weights = np.array(weights, dtype=np.float64)
            except ValueError as e:
                ssio.exception_message(f"Unable to convert passed weights to a np.array(). Likely means the passed value is not numerical (printed below):\n\n{weights}", e, with_frills=True, raise_exception=True)


            if len(weights) != self.n_frames:
                raise SSException(f'Passed frame weights array is {len(weights)} in length, while there are actually {self.n_frames} frames - these must match')


            if stride > 1:
                ssio.warning_message("WARNING: Using stride with weights is ALMOST certainly not a good idea unless the weights are\ncalculated for every stride-th residue", with_frills=True)
                return weights[list(range(0,self.n_frames,stride))]
            else:
                weights_sum = np.sum(weights)
                weights_difference = weights_sum - 1.0
                abs_weights_difference = abs(weights_difference)
                if abs(np.sum(weights) - 1.0) < etol:
                    return weights
                else:
                    ssio.exception_message(f"The passed weights are not within the specified floating point epsilon tolerance (etol={etol:f}). | sum of weights - 1 | = {abs_weights_difference:f}\n\n", with_frills=True, raise_exception=True)


        return False


    # ........................................................................
    #
    def __get_first_and_last(self, R1, R2, withCA=False):
        """Resolve a possibly-``None`` ``(R1, R2)`` pair into a valid residue range.

        Handles three things in one place:

        * Defaults: ``R1=None`` becomes the first residue (or first CA
          residue if ``withCA``); ``R2=None`` becomes the last.
        * Cap-awareness: when ``withCA=True`` the range is shrunk by one
          on each cap-bearing end.
        * Validation: out-of-range or boolean inputs raise SSException;
          ``R1 > R2`` is silently swapped.

        Parameters
        ----------
        R1, R2 : int or None
            Inclusive endpoints of the residue range. ``None`` means
            "use the chain endpoint".
        withCA : bool, optional
            If True, restrict default endpoints to CA-bearing residues
            (i.e. exclude ACE / NME caps). Default False.

        Returns
        -------
        tuple of (int, int, str)
            ``(R1, R2, selection_string)`` where the third element is an
            mdtraj-compatible ``"resid R1 to R2"`` string that can be
            handed directly to ``topology.select``.

        Raises
        ------
        SSException
            If ``R1`` / ``R2`` are booleans (deprecated), if either is
            negative, or if either exceeds ``n_residues - 1``.
        """

        # this is a defensive sanity check to revert a potentially bug-causing
        # flexibility in a previous version of the code
        if isinstance(R1, bool) or isinstance(R2, bool):
            raise SSException(f'Deprecation error: Prior to 0.1.9x versions __get_first_and_last() could take boolean values for R1 and R2. Starting with 0.2.0 it can only take integers and None. This message reflects a bug in SOURSOP, please contact Alex directly or raise an issue on GitHub')

        # first define as if we're starting from first and last residue with/without caps
        if R1 == None:
            R1 = 0
            if withCA:
                if self.ncap:
                    R1 = 1

        if R2 == None:
            R2 = self.n_residues - 1
            if withCA:
                if self.ccap:
                   R2 = self.n_residues - 2

        # finally flip around if R1 is larger than R2
        if R1 > R2:
            tmp = R2
            R2 = R1
            R1 = tmp

        if R1 < 0:
            raise SSException('ERROR: Requested a residue index (resid) less than 0')

        if R2 >= self.n_residues:
            raise SSException(f'ERROR: This protein only has {self.n_residues:d} residues, SO valid indices for selection are between 0 and {self.n_residues-1:d}, yet [{R2}] was passed to function.')

        return (R1, R2, f"resid {R1} to {R2}")


    # ........................................................................
    #
    def __check_stride(self, stride):
        """Checks that a passed stride value doesn't break everything. Returns
        `None` or raises a `SSException`.

        Parameters
        ----------

        stride: int
            The non-zero number of steps to perform while iterating across a
            trajectory.

        Raises
        ------
        SSException
            When the stride is larger than the number of available frames in the
            trajectory, or less than 1.
        """
        if stride > self.n_frames:
            raise SSException(f'stride ({stride}) is larger than the number of frames ({self.n_frames})')

        if stride < 1:
            raise SSException(f'stride ({stride}) is less than 1')



    # ........................................................................
    #
    def __check_single_residue(self, R1):
        """Internal function that checks that a single residue provided makes
        sense in the context of this protein.

        Returns `None` or raises a `SSException`.

        Parameters
        ----------
        R1: int
            The resid (`int`) of the residue to be  checked and validated.

        Raises
        ------
        SSException
            When the residue ID is greater than the chain length, or when
            the distances explored are greater than the chain size.
        """

        if R1 < 0:
            raise SSException(f"Trying to use a negative residue index [residue index = {R1}]")

        if R1  >= self.n_residues:
            raise SSException(f"Trying to use a residue ID greater than the chain length [residue index = {R1}, chain length = {self.n_residues}] ")

    # ........................................................................
    #
    def __check_cg_onebead(self):
        """Internal function that checks if the chain is a single chain and
        coarse grained (1 bead per residue).

        Returns
        --------------
        bool
            Returns True if the chain is a single chain and coarse grained
            (1 bead per residue), False otherwise.
        """

        # the only condition in which we have a single chain and coarse grained
        # is if the number of residues is equal to the number of atoms.
        # Use the already-computed residue count rather than building the
        # one-letter sequence: the latter raises KeyError on non-standard
        # residue names and would crash construction of CG/non-standard chains.
        return self.traj.n_atoms == self.__num_residues




        
        

    # ........................................................................
    #
    def __check_contains_CA(self, R1):
        """Function which checks if residue R1 contains an alpha carbon (CA)
        atom.

        Returns `None` or raises a `SSException`.

        Parameters
        ----------
        R1: int
            The resid to be checked (recall that resids always start from 0).

        Raises
        ------
        SSException
            If the residue is not found in the resid_with_CA list an exception
            is raised
        """
        if R1 not in self.resid_with_CA:
            raise SSException(f"Resid {R1} lacks an alpha carbon atom")

        return None

    # ........................................................................
    #
    def __get_subtrajectory(self, traj, stride):
        """Internal function which returns a subtrajectory. Expects `traj` to
        be an `mdtraj` trajectory object and `stride` to be an `int`.

        Parameters
        ----------
        traj: mdtraj.Trajectory
            An instance of an `mdtraj.Trajectory` which is non-empty - i.e.
            contains at least 1 frame.

        stride: int
            The non-zero number of steps to perform while iterating across
            the input trajectory, `traj`.

        Returns
        -----------
        mdtraj.Trajectory
            A sliced trajectory which contains the frames selected every
            `stride` step from the input trajectory, `traj`.
        """

        stride = int(stride)
        self.__check_stride(stride)

        if stride == 1:
            return traj
        else:
            return traj[::stride]


    # ........................................................................
    #
    def __get_resid_with_CA(self):
        """Internal function which should only be needed during initialization.
        Defines the list of residues where CA atoms are present, and the list
        of zero-indexed residue indices which contain CA atoms.

        In the case that the first resid in the self.topology object is 0 then
        these two lists are the same (which for MDTraj 1.9.5 or higher should
        always be the same) but for systems with multiple protein chains in
        MDTraj 1.9.4 or lower the residuesWithCA and idxWithCA willn differe
        for the second chain and onwards.

        This list is then assigned to the property variable ``self.resid_with_CA`` .

        Note this is quite slow for large trajectories

        Returns
        -----------
        list of int
            Ascending list of resids that have exactly one C-alpha atom
            (cap / non-standard residues without a unique CA are excluded).

        See also
        ------------
        resid_with_CA
        """

        # if this is a one-bead-per-residue coarse-grained chain we initialize
        # both the CA and 'all atoms' tables for every residue here. Every
        # residue is then its own CA, so resid_with_CA is simply 0..n-1.
        if self.__cg_onechain:
            self.__initialize_cg_atoms()
            return list(range(self.__num_residues))

        # All-atom path: a single pass over the topology atoms is dramatically
        # faster than one topology.select() per residue (the historical
        # initialization bottleneck). We then reproduce *exactly* the cache
        # side effects the old per-residue get_CA_index()/__residue_atom_lookup()
        # path produced, so downstream behaviour is byte-for-byte identical:
        #   * every residue gets a __residue_atom_table[ridx] dict whose 'CA'
        #     entry is the same array topology.select(...) would have returned
        #     (length-1 for a single CA, empty for caps, length-N otherwise);
        #   * __CA_residue_atom[ridx] is set only for residues with exactly
        #     one CA, to that array's [0] element (a numpy integer scalar,
        #     not a Python int - several callers/tests depend on the dtype);
        #   * resid_with_CA keeps residues with exactly one CA, in ascending
        #     topology order.
        ca_atoms_for_resid = {}
        for atom in self.topology.atoms:
            if atom.name == 'CA':
                ca_atoms_for_resid.setdefault(atom.residue.index, []).append(atom.index)

        resid_with_CA = []
        for res in self.topology.residues:
            ridx = res.index
            ca_atoms = ca_atoms_for_resid.get(ridx, [])

            if ridx not in self.__residue_atom_table:
                self.__residue_atom_table[ridx] = {}
            ca_array = np.array(ca_atoms, dtype=np.int64)
            self.__residue_atom_table[ridx]['CA'] = ca_array

            if len(ca_array) == 1:
                resid_with_CA.append(ridx)
                self.__CA_residue_atom[ridx] = ca_array[0]

        return resid_with_CA
    
    def __initialize_cg_atoms(self):
        """Internal function that initializes the CA and 'all atoms' for each residue
        in the topology. This is only called if the chain is a single chain and
        coarse grained (1 bead per residue).

        Returns
        -----------
        None
        """

        for res in self.topology.residues:
            self.__residue_atom_table[res.index] = {}
            self.__residue_atom_table[res.index]['all_atoms'] = np.array([res.index])
            self.__residue_atom_table[res.index]['CA'] = np.array([res.index])

        

    # ........................................................................
    #
    def __residue_atom_lookup(self, resid, atom_name=None):
        """Memoisation function to lookup the atomic index of a specific
        residues atom. Originally I'd assumed the underlying MDTraj
        ``topology.select()`` operation was basically a lookup, BUT it turns
        out it's actually *really* expensive, so this method converts
        atom/residue lookup information into a dynamic O(1) operation, greatly
        improving the performance of a number of different methods in the
        processes.
        
        As of Nov 2024, this function also explicitly will short circuit for 
        coarse grained chains (1 bead per residue) and raise an exception if 
        the atom name is not 'CA' (since that's the only atom we can select 
        for a coarse grained chain).

        Parameters
        ----------
        resid: int
            The residue index to lookup. If the residue has not been cached, it will be added to the
            lookup table for later reuse.

        atom_name: str or None {None}
            The name of the atom to lookup which will return the corresponding residue ID. Like the previous parameter,
            if that residue does not exist in the lookup table it will be added for later reuse.

        Returns
        -------
        list
            A list containing all the atoms corresponding to a given residue id that match the input residue id (`resid`)
            or, the residue corresponding to the atom name (`atom_name`).

        Raises
        ------
        SSException
            -  If the atom name is not 'CA' and the chain is coarse grained (1 bead per residue).
        """

        # if resid is not yet in table, create an empty dicitionary
        if resid not in self.__residue_atom_table:
            self.__residue_atom_table[resid] = {}
            

        # if we haven't specified an atom look up ALL the atoms associated with this residue
        if atom_name is None:

            # if all_atoms not yet associated with this residue
            if 'all_atoms' not in self.__residue_atom_table[resid]:
                self.__residue_atom_table[resid]['all_atoms'] = self.topology.select(f'resid {resid}')

            # return set of all atoms
            return self.__residue_atom_table[resid]['all_atoms']

        # if atom_ame not yet associated with this resid lookup
        # the atom_name from the underlying topology
        if atom_name not in self.__residue_atom_table[resid]:

            # if our chain is a single chain and coarse grained (1 bead per residue)
            # we actually raise an exception if that atom name is not CA
            if self.__cg_onechain:
                if atom_name != 'CA':
                    raise SSException("Trying to select a single atom for a coarse grained chain, but the atom name is not 'CA'. SOURSOP requires all bead atoms to be defined as 'CA'")
            
            self.__residue_atom_table[resid][atom_name] = self.topology.select(f'resid {resid} and name "{atom_name}"')

        # at this point we know the resid-atom_name pair is in the table
        # so goahead and look it up!
        return self.__residue_atom_table[resid][atom_name]


    # ........................................................................
    #
    def __get_selection_atoms(self, region=None, backbone=True, heavy=False):
        """Function which returns a list of atoms associated with the residues
        defined by the region keyword. If no region is supplied this returns
        the entire region (NME/ACE caps included).

        Parameters
        ----------

        region : `np.array`, `list`, or `tuple` {None}
            An array_like object of size 2 which defines the first and last
            residue (INCLUSIVE) for a region to be examined.

        backbone: bool {True}
            Boolean flag to determine if only the backbone atoms should be
            returned, or if all the full chain's atoms should be included
            (i.e. including sidechain).

        heavy: bool {False}
            Boolean flag to determine if we should only select heavy atoms
            or not (i.e. not H).


        Returns
        -------
        selectionatoms
            A `numpy.array` comprised of atom indices corresponding to the
            residues in a given region.

        Raises
        ------
        SSException
            When the input region is larger than 2.
        """



        # if a valid regino was passed...
        if region is None:
            R1 = 0
            R2 = self.n_residues - 1

            if backbone:
                if heavy:
                    selectionatoms = self.topology.select(f'backbone and resid {R1} to {R2} and not type "H"')
                else:
                    selectionatoms = self.topology.select(f'backbone and resid {R1} to {R2}')

            else:
                if heavy:
                    selectionatoms = self.topology.select(f'resid {R1} to {R2} and not type "H"')
                else:
                    selectionatoms = self.topology.select(f'resid {R1} to {R2}')

        elif len(region) == 2:
            if backbone:
                if heavy:
                    selectionatoms = self.topology.select(f'backbone and resid {region[0]} to {region[1]} and not type "H"')
                else:
                    selectionatoms = self.topology.select(f'backbone and resid {region[0]} to {region[1]}')
            else:
                if heavy:
                    selectionatoms = self.topology.select(f'resid {region[0]} to {region[1]} and not type "H"')
                else:
                    selectionatoms = self.topology.select(f'resid {region[0]} to {region[1]}')

        else:
            raise SSException(f"Trying to select a subsection of atoms, but the provided 'region' tuple/list is not of exactly length two [region={region}].\nCould indicate a problem, so be safe raising an exception")

        return selectionatoms





    # ........................................................................
    #
    def print_residues(self, verbose=True):
        """Map each zero-indexed residue position to its PDB resname-resid label.

        Useful when debugging cap/non-cap offsets or aligning SOURSOP residue
        indices against a published PDB numbering.

        Parameters
        ----------
        verbose : bool, optional
            If True (default), also prints each ``"<idx> --> <resname-resid>"``
            line to stdout. If False, the list is returned silently.

        Returns
        -------
        list of [int, str]
            One ``[index, resname-resid]`` pair per residue, in topology order.

        Example
        -------
        >>> mapping = protein.print_residues(verbose=False)
        >>> mapping[:3]
        [[0, 'ACE-0'], [1, 'MET-1'], [2, 'ALA-2']]
        """


        AA = self.get_amino_acid_sequence()
        return_list = []
        for i in range(0, len(AA)):
            if verbose is True:
                print(f"{i} --> {AA[i]}")
            return_list.append([i,AA[i]])
        return return_list



    # ........................................................................
    #
    def get_residue_atom_indices(self, resid, atom_name=None):
        """Look up the atom indices that belong to a specific residue.

        Results are memoised on first lookup so repeated calls for the same
        ``(resid, atom_name)`` pair are essentially free. The previous name of
        this method was ``residue_atom_com()``.

        Parameters
        ----------
        resid : int
            Resid of the residue to look up.
        atom_name : str, optional
            If given, restrict the result to atoms whose name matches this
            string (e.g. ``"CA"``, ``"CB"``, ``"N"``). If None (default), all
            atoms of the residue are returned.

        Returns
        -------
        list of int
            Atom indices into the underlying mdtraj topology that satisfy the
            ``resid`` / ``atom_name`` filter.

        Example
        -------
        >>> ca_atoms = protein.get_residue_atom_indices(5, atom_name='CA')
        >>> all_atoms = protein.get_residue_atom_indices(5)
        """

        return self.__residue_atom_lookup(resid, atom_name)



    # ........................................................................
    #
    def get_residue_COM(self, resid, atom_name=None):
        """Per-frame centre of mass (COM) of a residue or a specific atom of it.

        The COM is computed once and cached, so repeated calls for the same
        ``(resid, atom_name)`` pair are essentially free. All positions are
        returned in Angstroms.

        Parameters
        ----------
        resid : int
            Resid of the residue to extract.
        atom_name : str, optional
            If given, only the named atom contributes to the COM (a single-atom
            COM is just the atom's position). If None (default), every atom of
            the residue contributes.

        Returns
        -------
        np.ndarray
            Array of shape ``(n_frames, 3)`` giving the COM ``(x, y, z)`` of
            the residue in each frame.

        Example
        -------
        >>> com = protein.get_residue_COM(5)              # full-residue COM
        >>> ca  = protein.get_residue_COM(5, 'CA')        # CA atom position
        """

        # check its ok...
        self.__check_single_residue(resid)

        # for whole residues (atom_name not specified)
        if atom_name is None:
            if resid not in self.__residue_COM:
                atoms = self.__residue_atom_lookup(resid)
                TRJ_1 = self.traj.atom_slice(atoms)
                self.__residue_COM[resid] = 10*md.compute_center_of_mass(TRJ_1)
                del TRJ_1

            return self.__residue_COM[resid]

        else:

            # if resid is not yet in table, create an empty dicitionary
            if resid not in self.__residue_atom_COM:
                self.__residue_atom_COM[resid] = {}

            if atom_name not in self.__residue_atom_COM[resid]:
                atoms = self.__residue_atom_lookup(resid, atom_name)
                TRJ_1 = self.traj.atom_slice(atoms)
                self.__residue_atom_COM[resid][atom_name] = 10*md.compute_center_of_mass(TRJ_1)
                del TRJ_1

            return self.__residue_atom_COM[resid][atom_name]





    # ........................................................................
    #
    def get_amino_acid_sequence(self, oneletter=False, numbered=True):
        """Return the amino-acid sequence of the protein in several formats.

        Caps (ACE / NME) and other non-standard residues are included in the
        sequence using their three-letter or one-letter codes as recorded by
        the topology.

        Parameters
        ----------
        oneletter : bool, optional
            If True, use single-letter amino-acid codes. If False (default),
            use three-letter codes.
        numbered : bool, optional
            If True (default), include the resid in each entry as a
            ``"RESNAME-RESID"`` string. If False, just the resname.

        Returns
        -------
        list of str OR str
            * ``oneletter=False, numbered=True`` -> ``['MET-0', 'ALA-1', ...]``
            * ``oneletter=False, numbered=False`` -> ``['MET', 'ALA', ...]``
            * ``oneletter=True, numbered=True`` -> the full one-letter string,
              e.g. ``"MA..."``
            * ``oneletter=True, numbered=False`` -> ``['M', 'A', ...]``

        Example
        -------
        >>> protein.get_amino_acid_sequence(oneletter=True, numbered=True)
        'MAAEELANAK...'
        >>> protein.get_amino_acid_sequence()[:3]
        ['MET-0', 'ALA-1', 'ALA-2']
        """

        if oneletter:

            if self.__amino_acids_1LTR is None:
                res = []
                for i in self.topology.residues:
                    res.append(THREE_TO_ONE[str(i)[0:3].upper()])
                self.__amino_acids_1LTR = "".join(res)

        else:

            if self.__amino_acids_3LTR is None:
                res = []
                for i in self.topology.residues:
                    res.append(str(i)[0:3]+"-"+str(i)[3:])
                self.__amino_acids_3LTR = res

        # if numbered requsted
        if numbered:
            if oneletter:
                return self.__amino_acids_1LTR
            else:
                return self.__amino_acids_3LTR

        # else strip out the numbering with list comprehension
        else:
            if oneletter:
                return [x.split('-')[0] for x in self.__amino_acids_1LTR]
            else:
                return [x.split('-')[0] for x in self.__amino_acids_3LTR]




    # ........................................................................
    #
    def get_CA_index(self, resid):
        """Return the atom index of the alpha-carbon (CA) for a residue.

        Results are memoised on first lookup so repeated calls are essentially
        free. Raises ``SSException`` if the residue does not have exactly one
        CA atom (e.g. cap residues like ACE / NME have none).

        Parameters
        ----------
        resid : int
            Resid whose CA atom should be returned.

        Returns
        -------
        int
            The mdtraj atom index of the CA atom of residue ``resid``.

        Raises
        ------
        SSException
            If the residue has no CA atom (caps, non-standard residues, or
            coarse-grained beads named differently).

        Example
        -------
        >>> protein.get_CA_index(5)
        37
        """

        if resid not in self.__CA_residue_atom:
            return_val = self.__residue_atom_lookup(resid, 'CA')

            if len(return_val) == 1:
                self.__CA_residue_atom[resid] = return_val[0]
            else:
                raise SSException(f"get_CA_index - unable to find residue {resid}")

        return self.__CA_residue_atom[resid]


    # ........................................................................
    #
    def get_all_atomic_indices(self, resid):
        """Return every atom index that belongs to a residue.

        Results are memoised on first lookup so repeated calls for the same
        ``resid`` are essentially free.

        Parameters
        ----------
        resid : int
            Resid whose atoms should be returned.

        Returns
        -------
        list of int
            Atom indices (into the underlying mdtraj topology) for every atom
            of residue ``resid``.

        Example
        -------
        >>> protein.get_all_atomic_indices(5)
        [35, 36, 37, ..., 49]
        """

        return self.__residue_atom_lookup(resid)


    # ........................................................................
    #
    def get_multiple_CA_index(self, resID_list=None):
        """Return CA atom indices for many residues at once.

        Convenience wrapper around :meth:`get_CA_index` that handles three
        input forms: a single integer, a list of integers, or ``None``
        (meaning "every residue with a CA"). Residues without a CA are
        silently skipped rather than raising.

        Parameters
        ----------
        resID_list : int, list of int, or None, optional
            * ``None`` (default) -> use :attr:`resid_with_CA` (every residue
              that has a CA).
            * A single ``int`` -> wrap in a list, return one-element list.
            * A ``list[int]`` -> look up each in turn, skipping any that
              raise.

        Returns
        -------
        list of int
            CA atom indices, in the same order as the input residues.

        Example
        -------
        >>> protein.get_multiple_CA_index([5, 10, 15])
        [37, 86, 134]
        >>> all_ca = protein.get_multiple_CA_index()   # every CA
        """

        # if we've just passed a single unlisted integer
        # then just return the single residue associated
        # with
        if type(resID_list) == int:
            return ([self.get_CA_index(resID_list)])

        # if no value passed grab all the residues
        if resID_list == None:
            resID_list = self.resid_with_CA

        CAlist = []
        for res in resID_list:
            try:
                CAlist.append(self.get_CA_index(res))
            except SSException as e:
                print(e)
                continue

        return CAlist


    # ........................................................................
    #
    def calculate_all_CA_distances(self, residueIndex, mode='CA', only_C_terminal_residues=True, stride=1):
        """Distances between a reference residue and every other CA residue.

        Computes either CA-CA or COM-COM distances from ``residueIndex`` to
        all other residues in the chain that have a CA atom. By default only
        residues C-terminal of the reference are considered, which is the
        right choice when building an all-vs-all upper-triangular matrix.

        Distance is returned in Angstroms.

        Parameters
        ----------
        residueIndex : int
            Reference residue. Must have a CA atom; if not, the function
            returns ``-1`` rather than raising.
        mode : {'CA', 'COM'}, optional
            * ``'CA'`` (default) - distance between alpha-carbon atoms.
            * ``'COM'`` - distance between residue centres of mass.
        only_C_terminal_residues : bool, optional
            If True (default) only consider residues with index strictly
            greater than ``residueIndex``. If False, also include residues
            N-terminal of the reference.
        stride : int, optional
            Use every ``stride``-th frame. Default is 1.

        Returns
        -------
        np.ndarray or int
            * On success: array of shape ``(n_frames_after_stride, M)`` where
              ``M`` is the number of other residues considered.
            * If ``residueIndex`` has no CA: the integer ``-1``.

        Raises
        ------
        SSException
            If ``mode`` is not one of ``'CA'`` or ``'COM'``.

        Example
        -------
        >>> d = protein.calculate_all_CA_distances(5)
        >>> d.shape
        (1000, 86)
        """

        # validate input
        ssutils.validate_keyword_option(mode, ['CA', 'COM'], 'mode')

        # determine atomic index of CA atom for the residue you passed in
        try:
            CA_base = self.get_CA_index(residueIndex)
        except SSException:

            # if we couldln't find a C-alpha for this residue then nothing
            # makes sense so return -1
            return -1

        # list of atomic indices for C-alpha atoms we care about
        CAlist = []

        ##
        ## Block to use if using CA as base for computing distances
        ##
        if mode == 'CA':

            for residue in self.resid_with_CA:

                # only compute distances bigger than the residueIndex
                # such that by default full iteration calculates the
                # non-redundant map (i.e. the half diagonal)
                if residue <= residueIndex and only_C_terminal_residues:
                    continue

                # build a list of c-alpah atomic indices
                CAlist.append(self.get_CA_index(residue))

            # now we construct a nested list of lists of pairs to compute distances
            # between
            pairs=[]
            for CA in CAlist:
                pairs.append([CA_base, CA])

            if len(pairs) == 0:
                return np.array([])


            local_traj = self.__get_subtrajectory(self.traj, stride)

            # 10* for angstroms
            return 10*md.compute_distances(local_traj, np.array(pairs), periodic=False)

        if mode == 'COM':
            ##
            ## Block to use if using COM as base for computing distances
            ##
            return_distances = []

            # cycle over each resid for a residue with a CA atom
            for residue in self.resid_with_CA:

                # only compute distances bigger than the residueIndex
                # such that by default full iteration calculates the
                # non-redundant map (i.e. the half diagonal)
                if residue <= residueIndex and only_C_terminal_residues:
                    continue

                return_distances.append(self.get_inter_residue_COM_distance(residueIndex, residue))

            # finally convert list to numpy array and flip so returns in same format as CA mode
            return np.array(return_distances).transpose()



    # ........................................................................
    #
    def get_distance_map(self, mode='CA', RMS=False, stride=1, return_instantaneous_maps=False, weights=False, verbose=True):
        """Inter-residue distance map (and its standard deviation).

        Builds the full N x N matrix of inter-residue distances where N is
        the number of residues with a CA atom (ACE / NME caps are therefore
        excluded). Distances are returned in Angstroms. The matrix is upper
        triangular (entries below the diagonal are zero).

        Parameters
        ----------
        mode : {'CA', 'COM'}, optional
            * ``'CA'`` (default) - distances between alpha-carbon atoms.
            * ``'COM'`` - distances between residue centres of mass.
        RMS : bool, optional
            If True, report root-mean-square distances
            :math:`\\sqrt{\\langle r_{ij}^2 \\rangle}`, the polymer-physics
            order parameter. If False (default), report the ensemble mean
            :math:`\\langle r_{ij} \\rangle`.
        stride : int, optional
            Use every ``stride``-th frame. Default is 1 (every frame).
        return_instantaneous_maps : bool, optional
            If True, the first element of the returned tuple is a 3D array
            ``(n_frames, N, N)`` of per-frame distance maps rather than the
            2D ensemble-average map. Default is False.
        weights : list or np.ndarray, optional
            Per-frame weights for re-weighted averaging (e.g. T-WHAM output).
            Default ``False`` means uniform weighting; in that case the
            std-map is also returned, otherwise it is None.
        verbose : bool, optional
            If True (default), print one status line per row of the matrix.
            This calculation can be slow on long trajectories.

        Returns
        -------
        tuple of (np.ndarray, np.ndarray)
            ``(distance_map, std_distance_map)`` where:

            * ``distance_map`` is ``(N, N)`` ensemble-mean distances (or
              ``(n_frames, N, N)`` per-frame distances if
              ``return_instantaneous_maps=True``).
            * ``std_distance_map`` is ``(N, N)`` per-pair standard
              deviations across frames; ``None`` if ``weights`` was given.

        Example
        -------
        >>> mean_map, std_map = protein.get_distance_map(verbose=False)
        >>> mean_map.shape
        (92, 92)
        >>> rms_map, _ = protein.get_distance_map(mode='COM', RMS=True, verbose=False)
        """

        ssutils.validate_keyword_option(mode, ['CA', 'COM'], 'mode')

        weights = self.__check_weights(weights, stride)

        # use the previously identified residues with CA
        residuesWithCA = self.resid_with_CA

        # initialize empty matrices that we're gonna fill up
        n_res = len(self.resid_with_CA)

        # if we want to return ALL the distance maps...
        if return_instantaneous_maps is True:
            # use this to empircally work out dimensions of return dime
            test_data = self.calculate_all_CA_distances(self.resid_with_CA[0], mode=mode, stride=stride)
            distance_map = np.zeros([n_res, test_data.shape[0], n_res])

        # if we want to return JUST the average...
        else:
            distance_map = np.zeros([n_res, n_res])

        # return same standard deviation either way...
        std_distance_map = np.zeros([n_res, n_res])

        # cycle over CA-containing residues - note we have to define this SM_index
        # because the matrix itself indices from 0 onwards (not from the first
        # residue index)
        SM_index = 0
        for resIndex in self.resid_with_CA[0:-1]:

            ssio.status_message(f"On protein residue {SM_index} (overall residue index = {resIndex}) of {int(len(residuesWithCA))} [distance calculations]", verbose)

            # get all CA-CA distances between the residue of index resIndex and every other residue.
            # Note this gives the non-redudant upper triangle.
            full_data = self.calculate_all_CA_distances(resIndex, mode=mode, stride=stride)

            # if we want root mean square then NOW square each distances
            if RMS:
                full_data = np.power(full_data,2)

            # calculate mean and standard deviation
            if weights is not False:
                
                mean_data = np.average(full_data,0,weights=weights)

                # if we want RMS then NOW take square root of <rij^2>
                if RMS:
                    mean_data = np.sqrt(mean_data)
                std_data  = None
            else:

                mean_data = np.mean(full_data,0)

                # if we want RMS then NOW take square root of <rij^2>
                if RMS:
                    mean_data = np.sqrt(mean_data)

                std_data = np.std(full_data,0)
                
            # update the maps appropriately and increment the counter
            if return_instantaneous_maps is True:
                distance_map[SM_index].transpose()[1+SM_index:] = full_data.transpose()

            else:
                distance_map[SM_index][1+SM_index:len(residuesWithCA)] = mean_data

            # updated std map
            std_distance_map[SM_index][1+SM_index:len(residuesWithCA)] = std_data

            SM_index = SM_index + 1

        if return_instantaneous_maps is True:

            # note we have to transpose the distance map so the 1st index is the
            # frame index
            return (np.transpose(distance_map, axes=[1,0,2]), std_distance_map)
        else:
            return (distance_map, std_distance_map)



    # ........................................................................
    #
    def get_polymer_scaled_distance_map(self, nu=None, A0=None, min_separation=10, mode='fractional-change', stride=1, weights=False, etol=0.0000001, verbose=True):
        """Quantify how each inter-residue distance deviates from a homopolymer model.

        For a standard polymer the equilibrium distance scales as
        :math:`\\langle r_{ij} \\rangle = A_0 |i-j|^{\\nu}`. This method
        builds an N x N matrix where each entry quantifies how far the
        observed mean distance is from that homopolymer prediction. The
        deviation can be expressed in several ways (see ``mode``).

        Note that :math:`A_0` here is the prefactor of the inter-residue
        distance scaling and is NOT the same numerical value as the
        :math:`R_0` prefactor that defines :math:`R_g = R_0 N^{\\nu}`. The
        scaling exponent :math:`\\nu` should typically lie in
        ``[0.33, 0.598]`` for real polymers.

        If ``nu`` and ``A0`` are not supplied, the homopolymer model is fit
        automatically by :meth:`get_scaling_exponent` (with default settings)
        and that fit is used as the reference model.

        Parameters
        ----------
        nu : float, optional
            Scaling exponent used as the homopolymer reference. If supplied,
            ``A0`` must also be supplied. If ``None`` (default), both are
            obtained by calling :meth:`get_scaling_exponent`.
        A0 : float, optional
            Scaling prefactor for the homopolymer reference; paired with
            ``nu``.
        min_separation : int, optional
            Minimum sequence separation ``|i-j|`` for which a deviation is
            computed. Entries with separations below this threshold are filled
            with the mode-dependent default value (0 for the change modes,
            1 for ``'scaled'``). Default is 10.
        mode : {'fractional-change', 'signed-fractional-change', 'signed-absolute-change', 'scaled'}, optional
            How to quantify the deviation between the observed mean distance
            :math:`r_{ij}` and the homopolymer prediction :math:`p_{ij}`:

            * ``'fractional-change'`` (default): :math:`|r_{ij} - p_{ij}| / p_{ij}`.
              Unsigned magnitude.
            * ``'signed-fractional-change'``: :math:`(r_{ij} - p_{ij}) / p_{ij}`.
              Positive => expansion vs. the homopolymer, negative => contraction.
            * ``'signed-absolute-change'``: :math:`r_{ij} - p_{ij}` (Angstroms).
            * ``'scaled'``: :math:`r_{ij} / p_{ij}` (default value 1 on near-diagonal).
        stride : int, optional
            Use every ``stride``-th frame for the distance map. Default is 1.
        weights : list or np.ndarray, optional
            Per-frame weights for re-weighted analysis (e.g. T-WHAM output).
            Default ``False`` means uniform weighting.
        etol : float, optional
            Tolerance for the weights-sum-to-one check. Default ``1e-7``.
        verbose : bool, optional
            If True (default), print status updates during the distance map
            and (when needed) the scaling-exponent fit.

        Returns
        -------
        tuple of (np.ndarray, float, float, float)
            ``(map, nu, A0, redchi)`` where:

            * ``map`` is an ``(N, N)`` numpy matrix of deviations (units
              depend on ``mode``).
            * ``nu`` is the homopolymer scaling exponent used.
            * ``A0`` is the homopolymer scaling prefactor used.
            * ``redchi`` is the reduced chi-squared of the homopolymer fit
              (``-1`` if ``nu`` / ``A0`` were supplied rather than fit).

        Raises
        ------
        SSException
            If ``mode`` is not one of the four allowed values, if
            ``min_separation < 1``, if exactly one of ``nu``/``A0`` is given,
            or if the supplied ``nu``/``A0`` are out of physical range.

        Example
        -------
        >>> m, nu, A0, chi2 = protein.get_polymer_scaled_distance_map(verbose=False)
        >>> m.shape, round(nu, 2)
        ((92, 92), 0.6)
        """

        # First validate keyword
        ssutils.validate_keyword_option(mode, ['fractional-change','scaled', 'signed-fractional-change', 'signed-absolute-change'], 'mode')

        # next check that the minimum separation requested makes sense... (this is only a partial check)
        if min_separation < 1:
            raise SSException("Minimum separation to be used must be greater than 0")

        ## next see if nu or A0 have been provided...

        # If NEITHER provided then do fitting here and now!
        if nu == None and A0 == None:
            ssio.status_message("Fitting data to homopolymer mode...", verbose)

            # remind that this is the old get_scaling_exponent_v2()
            # ALSO note this uses COM which is also how the distance map is computed
            SE = self.get_scaling_exponent(verbose=False, weights = weights, stride = stride, etol=etol)
            nu = SE[0]
            A0 = SE[1]
            REDCHI = SE[7]

        elif nu is None:
            raise SSException(f"A0 parameter provided [{A0:.5f}] but nu was not. Must provide BOTH or neither (in which case fitting is done)")

        elif A0 is None:
            raise SSException(f"nu parameter provided [{nu:.5f}] but A0 was not. Must provide BOTH or neither (in which case fitting is done)")

        # else both were provided so we double check they're valid..
        else:
            # sanity check nu
            if nu <= 0 or nu > 1:
                raise SSException("Nu parameter must be in interval 0 < nu <= 1 (and probably should be between 0.33 and 1.0...)")

            # sanity check A0
            if A0 <= 0:
                raise SSException("A0 paameter must be greater than 0")

            # not computing a reduced chi
            REDCHI = -1

        # We now define a function which will evaluate how we assess the deviation (or lack thereof) from the
        # traditional polymer scaling behaviour. This just saves us having if/then statements that get evaluated
        # on each loop of the analysis script below, and also makes it easy to add in additional 'modes'

        if mode  == 'fractional-change':
            def d_funct(dMap_val, p_val):
                return abs(dMap_val - p_val)/p_val

            default_val = 0

        elif mode == 'signed-fractional-change':
            def d_funct(dMap_val, p_val):
                return (dMap_val - p_val)/p_val

            default_val = 0

        elif mode == 'signed-absolute-change':
            def d_funct(dMap_val, p_val):
                return (dMap_val - p_val)

            default_val = 0

        elif mode == 'scaled':
            def d_funct(dMap_val, p_val):
                return dMap_val/p_val

            default_val = 1

        ## Now we've set everything up we can actually compute some numbers

        # first compute and get the distance map (note [0] means we get first element
        # which is mean distance ([1] is STDEV)
        distance_map = self.get_distance_map(mode='COM', RMS=True, stride=stride, weights=weights, verbose=False)[0]

        # get distance map dimensions (will be a square so just take X-dim)
        dimensions = distance_map.shape[0]

        if dimensions <= min_separation:
            raise SSException('The minimum separation is shorter than the chain length')

        # compute expected distance given the standard polymer scaling model
        expected_distances = sstools.powermodel(list(range(0,dimensions)), nu, A0)

        # initialize the return matrix and then populate for distances that are
        # above the minimum threshold. We only populate upper right triangle
        return_matrix = np.zeros((dimensions, dimensions))
        for i in range(0,dimensions):
            for j in range(0,dimensions):
                if j-i < min_separation:
                    return_matrix[i,j] = default_val
                else:
                    return_matrix[i,j] = d_funct(distance_map[i,j], expected_distances[(j-i)])

        return (return_matrix, nu, A0, REDCHI)



    # ........................................................................
    #
    def get_local_heterogeneity(self, fragment_size=10, bins=None, stride=20, verbose=True):
        """Sliding-window heterogeneity: per-position distribution of intra-window RMSDs.

        At each starting residue ``i`` (from 0 to ``n_residues - fragment_size``)
        the method computes the RMSD of the local fragment ``[i, i+fragment_size]``
        across every (strided) frame, then summarises that distribution by
        its mean, standard deviation, and histogram. The Phi parameter of
        Lyle, Das & Pappu (2013) is built from this same family of D-values.

        Computational cost scales with ``n_residues * (n_frames / stride)``;
        the default ``stride=20`` keeps things tractable on long trajectories.

        Parameters
        ----------
        fragment_size : int, optional
            Size of the sliding window (residues). Must be ``>= 2`` and
            ``<= n_residues``. Default is 10.
        bins : np.ndarray or list, optional
            Histogram bin edges. If None (default), uses
            ``np.arange(0, 10, 0.01)``. Must contain at least two evenly-
            spaced floats.
        stride : int, optional
            Use every ``stride``-th frame in the inner RMSD pass. Default 20.
        verbose : bool, optional
            If True (default), print one status line per starting residue.

        Returns
        -------
        tuple of (list, list, list, np.ndarray)
            ``(mean, std, histo, bins)`` where:

            * ``mean`` (list of float, length ``n_residues - fragment_size``):
              mean intra-window RMSD per starting residue.
            * ``std`` (list of float, same length): standard deviation of the
              same RMSD distribution.
            * ``histo`` (list of np.ndarray, same length): histogram counts
              per starting residue.
            * ``bins`` (np.ndarray): the bin edges used (echoed back).

        Raises
        ------
        SSException
            If ``fragment_size`` is out of range, or ``bins`` is not a valid
            evenly-spaced 1D vector.

        Example
        -------
        >>> mean, std, hist, bins = protein.get_local_heterogeneity(fragment_size=8)

        References
        ----------
        Lyle N, Das RK, Pappu RV. A quantitative measure for protein
        conformational heterogeneity. J Chem Phys. 2013;139(12):121907.
        """

        # validate bins
        if bins is None:
            bins = np.arange(0,10,0.01)
        else:
            try:
                if len(bins)  < 2:
                    raise SSException('Bins should be a numpy defined vector of values - e.g. np.arange(0,1,0.01)')
            except TypeError:
                raise SSException('Bins should be a list, vector, or numpy array of evenly spaced values')

            try:
                bins = np.array(bins, dtype=float)
            except ValueError:
                raise SSException('Passed bins could not be converted to a numpy array of floats')

        # check stride is ok
        self.__check_stride(stride)

        # get the residue IDXs were going to use
        #res_idx_list = self.residue_index_list
        n_frames = self.n_frames


        # check the fragment_size is appropriate
        if fragment_size > self.n_residues:
            raise SSException('fragment_size is larger than the number of residues')
        if fragment_size < 2:
            raise SSException('fragment_size must be 2 or larger')


        meanData = []
        stdData  = []
        histo    = []

        # cycle over each sub-region in the sequence

        # untested - for loop used to be frag_idx in res_idx_list[0:-fragment_size]:
        for frag_idx in range(0, self.n_residues - fragment_size):
            tmp = []
            ssio.status_message(f"On range {frag_idx}", verbose)

            # for each frame in ensemble, calculate RMSD for that sub-region compared to
            # all other sub-regions (i.e. we're doing a 1-vs-all RMSD calculation for EACH
            # frame (after adjusting for stride) for a subregion of the protein
            for j in range(0, n_frames, stride):
                tmp.extend(self.get_RMSD(j ,-1, region=[frag_idx, frag_idx+fragment_size]))

            # compute a histogram for this large dataset
            (b,c) = np.histogram(tmp,bins)
            histo.append(b)

            meanData.append(np.mean(tmp))
            stdData.append(np.std(tmp))

        return (meanData, stdData, histo, bins)



    # ........................................................................
    #
    def get_D_vector(self, stride=20, verbose=True):
        """Pairwise frame-vs-frame conformational dissimilarity vector :math:`D`.

        For every pair of (strided) frames :math:`(A, B)` the method
        computes a single dissimilarity scalar
        :math:`D_{AB} = 1 - V_A \\cdot V_B / (||V_A|| ||V_B||)`
        where :math:`V` is the flattened upper-triangular CA-CA distance
        vector for that frame. This is the input to the Phi parameter of
        Lyle, Das, Pappu (2013).

        Using CA-CA distances (rather than full inter-atomic distances)
        makes the measure backbone-specific and *much* faster than the
        original formulation.

        Parameters
        ----------
        stride : int, optional
            Use every ``stride``-th frame. Default is 20.
        verbose : bool, optional
            If True (default), print one status line per outer-loop frame.

        Returns
        -------
        np.ndarray
            1D array of length :math:`\\binom{n_{frames}/\\text{stride}}{2}`
            containing the pairwise dissimilarity values.

        Example
        -------
        >>> D = protein.get_D_vector(stride=20)
        >>> D.mean()   # mean pairwise dissimilarity
        0.27

        References
        ----------
        Lyle N, Das RK, Pappu RV. A quantitative measure for protein
        conformational heterogeneity. J Chem Phys. 2013;139(12):121907.
        """

        # get the list of residues which have CA (typically this means we exlcude
        # ACE and NME)
        residuesWithCA = self.resid_with_CA

        # compute number of frames exactly (this is empirical, but ensures we're always consistent with the projected
        # dimensions in the first for-loop)
        tmp = self.calculate_all_CA_distances(residuesWithCA[0], stride=stride, only_C_terminal_residues=True)
        n_frames =  np.shape(tmp)[0]

        all_distances = np.zeros([len(residuesWithCA), len(residuesWithCA), n_frames])

        # first compute upper triangle only (lower traingle is identical and doesn't change the answer
        # so we stick with the upper traingle only)
        SM_index=0
        for resIndex in residuesWithCA[0:-1]:
            ssio.status_message(f"Calculating non redundant distance for res. {resIndex} ", verbose)

            vals = self.calculate_all_CA_distances(resIndex, stride=stride, only_C_terminal_residues=True)

            # have to include a -1 here because we don't have a self:self distance
            all_distances[SM_index][0:(len(residuesWithCA)-1)-SM_index] = vals.transpose()
            SM_index=SM_index+1

        # number of residues we're calculating distances between
        n_res = np.shape(all_distances)[0]

        # reshape to n_frame x n_rij 2D array
        all_distance_tmp = all_distances.transpose().reshape(n_frames,n_res*n_res)

        # find the idx of non-zero elements in the first frame (will be true for all frames - zeros
        # originate because we only computed the upper triangle, so 1/2 of elements in each row of
        # all_distance_tmp are zero
        non_zero_idx = np.nonzero(all_distance_tmp[0])[0]

        # finally extract out the positions that are nonzero, setting us up for the phi analysis
        non_zero_distance = all_distance_tmp[:,non_zero_idx]

        # calculate the D-vector of all frames
        D_vector = []
        for A in range(0, n_frames):
            ssio.status_message(f"Running PHI calculation on frame {A} of {n_frames}", verbose)


            for B in range(A+1, n_frames):

                """# This is the old less efficient implementation of the
                algorithm, kept in case,

                # for some reason, the new implementation has issues...

                # get the vector of distances for frame A and frame B - this extracts all the
                # inter-residue distances for frame A and frame B
                VA = all_distances.transpose()[A].flatten()
                VB = all_distances.transpose()[B].flatten()

                # remove zero entries (all real distances MUST be greater
                # than zero because there's a hardsphere potential stopping things actually
                # reaching zero distance (distance is center-of-mass calculated for Ca)

                # we need to calculate the allowed set
                NZ_A = np.nonzero(VA)[0]
                NZ_B = np.nonzero(VB)[0]

                # this finds the positions that in both vectors were non-zero
                all_zero = np.intersect1d(NZ_A, NZ_B)

                # now for those positions extract a new length-match vector for frame
                # A and B that contains only non-zero positions
                VA = VA[all_zero]
                VB = VB[all_zero]
                """

                VA = non_zero_distance[A]
                VB = non_zero_distance[B]

                # and compute the D value for comparing these two frames
                D_vector.append(1 - np.dot(VA,VB)/(np.linalg.norm(VA)*np.linalg.norm(VB)))

        return np.array(D_vector)



    # ........................................................................
    #
    def get_RMSD(self, frame1, frame2=-1, region=None, backbone=True, stride=1):
        """Aligned root-mean-square deviation (RMSD) of one frame vs. another (or all).

        Both modes superpose the target onto the reference before computing
        the RMSD, so the result is rotation/translation invariant. Units are
        Angstroms.

        Parameters
        ----------
        frame1 : int
            Index of the reference frame.
        frame2 : int, optional
            Index of the comparison frame. ``-1`` (default) means "every
            frame" — i.e. compute an RMSD trace of ``frame1`` against the
            entire trajectory (subsampled by ``stride``).
        region : list or tuple of length 2, optional
            ``[first_resid, last_resid]`` (inclusive) restricting the atoms
            used to compute and align the RMSD. ``None`` (default) uses the
            full chain.
        backbone : bool, optional
            If True (default), use only the backbone heavy atoms (N, CA, C, O).
            If False, use every atom in the selection.
        stride : int, optional
            Used only when ``frame2 == -1``. Compare ``frame1`` against every
            ``stride``-th frame of the trajectory. Default is 1.

        Returns
        -------
        np.ndarray
            * 1-element array if ``frame2`` is a specific frame index.
            * Length-``ceil(n_frames/stride)`` array if ``frame2 == -1``.

            Values are in Angstroms.

        Example
        -------
        >>> rmsd_traj = protein.get_RMSD(frame1=0)              # frame 0 vs all
        >>> rmsd_pair = protein.get_RMSD(frame1=0, frame2=42)   # frame 0 vs 42
        """

        # get the selection atoms (perform correction if required)
        selectionatoms = self.__get_selection_atoms(region=region, backbone=backbone)

        # set the reference trajectory we're working with
        ref = self.traj

        # if a second frame number was provided with which we're going to work with
        if frame2 > -1 and isinstance( frame2, int ):

            # our target is now a single (i.e. doing RMSD of two structures)
            target = self.traj.slice(frame2)
        else:

            # else we're going to carry out an RMSD comparison between *every* stride-th
            # frame and frame 1
            target = self.__get_subtrajectory(self.traj, stride)

        # return the RMSD comparison in Angstroms
        return 10*md.rmsd(target, ref, frame1, atom_indices=selectionatoms)



    # ........................................................................
    #
    def get_Q(self,
              protein_average = True,
              region = None,
              beta_const = 50.0,
              lambda_const = 1.8,
              native_contact_threshold = 4.5,
              stride = 1,
              native_state_reference_frame=0,
              weights = False):


        """Fraction-of-native-contacts order parameter :math:`Q` (Best et al.).

        Native contacts are defined from the reference frame (frame 0) as
        heavy-atom pairs separated by at least 4 residues in sequence and
        less than ``native_contact_threshold`` Angstroms in space. Each
        frame's :math:`Q` is then a sigmoid-weighted fraction of those
        contacts that remain "native-like" in that frame, per the Best,
        Hummer, Eaton formula (Best et al. 2013).

        The reference frame is hard-coded to frame 0 because preserving a
        variable reference frame became unwieldy with re-weighted ensembles.

        Parameters
        ----------
        protein_average : bool, optional
            * If True (default), return one Q value per frame (the protein
              average).
            * If False, return per-native-contact and per-residue breakdowns
              in a 5-tuple — useful for visualising which contacts are most
              persistent.
        region : list or tuple of length 2, optional
            ``[first_resid, last_resid]`` (inclusive) restricting the
            residues used to identify native contacts. ``None`` (default)
            uses the full protein.
        beta_const : float, optional
            Steepness of the sigmoid weighting (1/nm). Default 50. Do not
            change without strong justification.
        lambda_const : float, optional
            Sigmoid width factor (all-atom default 1.8). Do not change
            without strong justification.
        native_contact_threshold : float, optional
            Distance cutoff (Angstroms) that defines a contact in the
            reference frame. Default is 4.5.
        stride : int, optional
            Use every ``stride``-th frame. Default is 1. Must be 1 if
            ``weights`` is given.
        native_state_reference_frame : int, optional
            Currently ignored — the reference frame is always frame 0. Kept
            in the signature for backward compatibility.
        weights : list or np.ndarray, optional
            Per-frame weights for re-weighted analysis. Default ``False``.

        Returns
        -------
        np.ndarray OR tuple of length 5
            * If ``protein_average=True`` (default): 1D array of length
              ``ceil(n_frames/stride)`` with the per-frame Q value.
            * If ``protein_average=False``: a 5-tuple

              0. ``per_contact_fraction`` - 1D array of length
                 ``N_NATIVE_CONTACTS``; fraction of time each native
                 contact was native-like over the trajectory.
              1. ``native_contact_pairs`` - ``(N_NATIVE_CONTACTS, 2)`` int
                 array; the atom-pair definitions of native contacts.
              2. ``per_residue_contacts`` - dict keyed by ``"RESNAME-RESID"``;
                 each value is a list of fractional-native-contact scores
                 for atoms in that residue. Take the mean per key to get the
                 per-residue Q.
              3. ``ordered_residue_keys`` - the keys of (2) in residue order.
              4. ``res_res_q_matrix`` - ``(n_residues, n_residues)`` array of
                 inter-residue Q values (symmetric).

        Raises
        ------
        SSException
            If ``weights`` is given with ``stride != 1``, or any constant
            cannot be coerced to float.

        Example
        -------
        >>> q = protein.get_Q()                          # per-frame Q
        >>> q.mean()
        0.85
        >>> per_c, pairs, per_res, keys, res_res = protein.get_Q(protein_average=False)

        References
        ----------
        Best RB, Hummer G, Eaton WA. Native contacts determine protein
        folding mechanisms in atomistic simulations.
        PNAS 2013. doi:10.1073/pnas.1311599110
        """

        # SET
        native_state_frame = 0
        n_res = self.n_residues

        # less stringent weights test cos trajectory is one frame too long because we probably loaded the PDB file
        # as a frame
        #weights = self.__check_weights(weights, stride)

        if weights is not False and stride != 1:
            raise SSException('For get_Q() weights must be set for EACH frame and stride=1')

        # if we're using a subregion
        # NOTE this is WAY more elegant than the previous way of doing this but there *used* to be problems with MDTraj doing
        # things like this...
        selectionatoms = self.__get_selection_atoms(region, backbone=False, heavy=True)

        # extract out the native state frame
        native = self.traj.slice(native_state_frame)

        # get the sub-trajectory to be used
        target = self.__get_subtrajectory(self.traj, stride)

        # now align the entire trajectory to the 'native' frame
        target.superpose(target, frame=native_state_frame, atom_indices=selectionatoms)

        try:
            BETA_CONST = float(beta_const)       # in reciprocal nm (1/nm)
            LAMBDA_CONST = float(lambda_const)    # For all-atom simulations

            # Native contact threshold distance in nm (not param is passed in Angstroms but the calculation
            # happens expecting nanometers so have to update (hence divide by 10). NB: This breaks standard
            # soursop convention and we probably should rescale the various parameters so this works in
            # A instead of bm
            NATIVE_CUTOFF = float(native_contact_threshold)/10


        except ValueError as e:
            raise SSException(f'Could not convert constant into float for setting constants in get_Q().\nSee below:\n\n{e}')


        # use all pairs of atoms that are over 3 away in sequence space
        heavy_pairs = np.array(
            [(i,j) for (i,j) in combinations(selectionatoms, 2)
             if abs(native.topology.atom(i).residue.index - \
                    native.topology.atom(j).residue.index) > 3])

        # compute the distances between these pairs in the native state
        ## NB: This breaks soursop convention of converting nm->A at the site of
        ## where it's calculated.
        heavy_pairs_distances = md.compute_distances(native[0], heavy_pairs, periodic=False)[0]

        # and get the pairs s.t. the distance is less than NATIVE_CUTOFF. This returns the
        # set of interatomic residues that define the native contacts
        native_contacts = heavy_pairs[heavy_pairs_distances < NATIVE_CUTOFF]

        # now compute these distances for the whole trajectory
        r = md.compute_distances(target, native_contacts, periodic=False)

        # and recompute them for just the native state
        r0 = md.compute_distances(native[0], native_contacts, periodic=False)

        # If we're just computing the protein average then this returns the Q value for the whole protein on a per-frame basis
        if protein_average:
            if weights is not False:
                raise SSException('Reweighting for frame averaged should be done with trajectory weights OUTSIDE of SOURSOP')

            q = np.mean(expit(-BETA_CONST * (r - LAMBDA_CONST * r0)), axis=1)

            return q

        else:

            # if the analysis is to be re-weighted uses the weights here on a per-frame basis
            if weights is not False:
                q_full = expit(-BETA_CONST * (r - LAMBDA_CONST * r0)).transpose()

                q = []
                for i in q_full:

                    # note i[1:] means we ignore the first (native) frame
                    q.append(np.average(i[1:], 0, weights))

                q = np.array(q)
            else:

                # check this makes sense - aboe we do i[1:] should probably correct this to remove
                # the native state structure? Anyway, here we're averaging over every from for each
                # residue
                q = np.mean(expit(-BETA_CONST * (r - LAMBDA_CONST * r0)), axis=0)

            # get the set of unqiue atoms which are involved in native contacts
            unique_native_contact_atoms = np.unique(np.hstack((np.transpose(native_contacts)[0],np.transpose(native_contacts)[1])))

            res2at = {}
            res2res = {}

            # for each unqiue atom
            for atom in unique_native_contact_atoms:

                # determine the name of the residue it's from and update the res2at dictionary
                # if needed
                local_res = str(self.topology.atom(atom).residue)

                if local_res not in res2at:
                    res2at[local_res] = []
                    res2res[int(local_res[3:])] = local_res

                # now for every native contact pair that atom is involved in,
                # associate the fraction of the time it's native with the
                # residue in question.
                for pair_idx in range(0, len(native_contacts)):
                    if atom in native_contacts[pair_idx]:
                        res2at[local_res].append(q[pair_idx])

            # we now have a dictionary where, for each residue, we have the
            # ensemble average fraction of the simulation each atom was making
            # native contacts. Note different residues have different numbers
            # of atoms (obviously...) SO each entry in res2at is going to be
            # a variable length

            # now construct an n-res by n-rex empty matrix
            res_res_matrix = np.zeros((self.n_residues, self.n_residues))
            res_res_matrix_count = np.zeros((n_res, n_res))

            # and for each pair of atoms as identified previously, look up
            # which residue they're from, determine the q score for that pairwise
            # interaction, and increment the associated positions on the nres by
            # nres matrix, keeping count of how many such increments we make in the
            # res_res_matrix_count matrix
            for pair_idx in range(0,len(native_contacts)):

                pair = native_contacts[pair_idx]
                R1 = self.topology.atom(pair[0]).residue.index
                R2 = self.topology.atom(pair[1]).residue.index

                res_res_matrix[R1,R2] = q[pair_idx] + res_res_matrix[R1,R2]
                res_res_matrix[R2,R1] = q[pair_idx] + res_res_matrix[R2,R1]

                res_res_matrix_count[R1,R2] = 1 + res_res_matrix_count[R1,R2]
                res_res_matrix_count[R2,R1] = 1 + res_res_matrix_count[R2,R1]

            # pairwise division accounts for residues with more atoms; safe_count
            # replaces zero entries with 1 to avoid 0/0, those positions are zeroed out by np.where
            safe_count = np.where(res_res_matrix_count != 0, res_res_matrix_count, 1)
            normalized_res_matrix = np.where(res_res_matrix_count != 0, res_res_matrix / safe_count, 0.0)

            # just as a convenience, build a sorted list of the residues which makes
            # the data a bit easier to play with going forward.
            res2res_keys = list(res2res.keys())
            np.sort(res2res_keys)
            sorted_residues = []

            for lk in res2res_keys:
                sorted_residues.append(res2res[lk])

            # finally, return all the things mentioned in the function description (note we nan-to-num
            # to remove all the NaNs from the norlaized res matrix (generated by dividiving by zeor)
            return (q, native_contacts, res2at, sorted_residues, np.nan_to_num(normalized_res_matrix,0))


    # ........................................................................
    #
    #
    def get_contact_map(self, distance_thresh=5.0, mode='closest-heavy', stride=1, weights=False):
        """Inter-residue contact map and per-residue contact-order vector.

        For every pair of residues this returns the fraction of frames in
        which the chosen inter-residue distance is below
        ``distance_thresh``. The companion contact-order vector reduces the
        2D map to a per-residue summary by averaging over the appropriate
        neighbours-excluded denominator.

        ``i to i``, ``i to i+1``, and ``i to i+2`` pairs are always excluded.
        ACE / NME caps are excluded from the residue set.

        Parameters
        ----------
        distance_thresh : float, optional
            Contact threshold in Angstroms. A pair is "in contact" in a given
            frame if its distance under the chosen ``mode`` is below this
            value. Default is 5.0.
        mode : {'closest-heavy', 'ca', 'closest', 'sidechain', 'sidechain-heavy'}, optional
            How the inter-residue distance is defined (the schemes from
            ``mdtraj.compute_contacts``):

            * ``'closest-heavy'`` (default) - closest pair of heavy atoms.
            * ``'ca'`` - alpha-carbon distance.
            * ``'closest'`` - closest pair of any atoms (incl. hydrogens).
            * ``'sidechain'`` - closest pair of sidechain atoms (mdtraj >= 1.8).
            * ``'sidechain-heavy'`` - closest pair of sidechain heavy atoms
              (mdtraj >= 1.8). Raises if any residue is glycine.
        stride : int, optional
            Use every ``stride``-th frame. Default is 1. Must be 1 if
            ``weights`` is provided.
        weights : list or np.ndarray, optional
            Per-frame weights for re-weighted averaging. Default ``False``.

        Returns
        -------
        tuple of (np.ndarray, np.ndarray)
            ``(contact_map, contact_order)`` where:

            * ``contact_map`` is an ``(N, N)`` matrix of contact fractions
              in ``[0, 1]``.
            * ``contact_order`` is a 1D vector of length ``N`` giving the
              mean fractional contact count per residue, normalised by the
              number of valid partners for each residue.

        Raises
        ------
        SSException
            For ``mode='sidechain-heavy'`` when any glycine is present, or
            if ``weights`` is given with ``stride != 1``.

        Example
        -------
        >>> cmap, corder = protein.get_contact_map()
        >>> cmap.shape
        (92, 92)
        """

        ssutils.validate_keyword_option(mode, ['closest-heavy', 'ca', 'closest', 'sidechain', 'sidechain-heavy'] , 'mode')

        if weights is not False:
            if int(stride) != 1:
                raise SSException("Cannot accomodate weights if stride is not set to 1")

        # check weights are correct
        weights = self.__check_weights(weights, stride)

        # set the distance threshold to a value in nm (we use A by default)
        distance_thresh_in_nm = float(distance_thresh/10.0)

        # build a substractectory based on the stride argument
        subtraj = self.__get_subtrajectory(self.traj, stride)

        # ensure we only select main chain atoms (no termini) - NOTE, this
        # is a REALLY useful design pattern - should consider re-writing the
        # code to use this...
        mainchain_atoms = self.topology.select('(not resname "NME") and (not resname "ACE")')

        # compute the contactmap and square-form it (map per frame)
        # CMAP is a [N_FRAMES x N_RES x N_RES] array

        # sidechain-heavy fails for GLY (no heavy sidechain atoms); check explicitly
        # so the behavior is consistent regardless of mdtraj version (older versions
        # raised ValueError; 1.11+ emits a warning and uses sidechain-H instead)
        if mode == 'sidechain-heavy':
            gly_atoms = self.topology.select('resname GLY and (not resname "NME") and (not resname "ACE")')
            if len(gly_atoms) > 0:
                msg = "Failed computing contacts. This is likely because one of the residues has a glycine and\nthere are no heavy sidechain residues in glycine. Raising exception..."
                raise SSException(msg)

        try:
            CMAP_nonsquare = md.compute_contacts(subtraj.atom_slice(mainchain_atoms), scheme=mode)
        except ValueError as e:
            if str(e) == 'zero-size array to reduction operation minimum which has no identity':
                if mode == 'sidechain-heavy':
                    msg = "Failed computing contacts. This is likely because one of the residues has a glycine and\nthere are no heavy sidechain residues in glycine. Raising exception..."

                    ssio.exception_message(msg, e, with_frills=True, raise_exception=False)
                    raise SSException(msg)

            raise e

        CMAP = md.geometry.squareform(CMAP_nonsquare[0], CMAP_nonsquare[1])

        # extract the normalization factor used to compute fractional
        # contacts
        normalization_factor = np.shape(CMAP)[0]

        # build a MASK where distance is not zero (i.e. where distances were calculated
        MASK =  (CMAP[0] != 0)*1

        # if no weights...
        if weights is False:

            # for each frame set true/false if less than threshold,
            # then convert bools to ints and sum over all frames
            # and finally normalize by the normalization factor. This
            # gives us the _normalized_ contact map (i.e. each element is between 0 and 1)
            normalized_contact_map = (np.sum(1*(CMAP < distance_thresh_in_nm),0)*MASK) / float(normalization_factor)

        # else, if weights...
        else:

            # if we use weights then we multiply each frame's contact map by the weight and sum
            n_frames = CMAP.shape[0]
            normalized_contact_map = np.zeros((CMAP.shape[1],CMAP.shape[1]))
            for fid in range(0,n_frames):
                normalized_contact_map = normalized_contact_map + (np.ndarray.astype((CMAP[fid]<distance_thresh_in_nm),int)*MASK)*weights[fid]

        # we can further reduce the dimensionality to ask which residues are most involved in contacts with outher
        # residues in general (i.e. without caring about what those residues are). This gives us a normalized
        # contact order.
        n_res = int(np.shape(normalized_contact_map)[0])

        # So this is a kind of funky line, but basically because we don't calculate over distances that are
        # that are
        # i to i
        # i to i-1
        # i to i-2
        # i to i+1
        # i to i+2
        # which means for MOST residues the max possible is n_res-5, but for those at the end and start
        # there is no i-1/i-2 for the 0th residues, so the line below builds a vector that for each residues
        # calculates the TRUE max fractional contacts
        contact_order_normalization_vector = n_res - np.hstack((np.hstack(([3,4],np.repeat(5,n_res-4))),[4,3]))

        normalized_contact_order = np.sum(normalized_contact_map,0)/contact_order_normalization_vector

        return (normalized_contact_map, normalized_contact_order)



    # ........................................................................
    #
    #
    def get_clusters(self, region=None, n_clusters=10, backbone=True, stride=20):
        """Ward hierarchical clustering of conformations by RMSD.

        Builds a strided all-vs-all RMSD matrix and groups frames using
        Ward's agglomerative clustering. The number of clusters
        (``n_clusters``) must be specified up front. The chosen "centroid"
        of each cluster is the member whose summed exponential RMSD weight
        is largest (i.e. the most internally-connected frame).

        Parameters
        ----------
        region : list of length 2, optional
            ``[first_resid, last_resid]`` (inclusive) restricting the atoms
            used for the RMSD calculation. ``None`` (default) uses the full
            chain.
        n_clusters : int, optional
            Target number of clusters. Default is 10. Note that the actual
            number of clusters returned may be smaller if some labels are
            never assigned.
        backbone : bool, optional
            If True (default), use only backbone heavy atoms for the RMSD
            (much faster). If False, use every atom in the selection.
        stride : int, optional
            Use every ``stride``-th frame. Default is 20. ``stride=1`` does
            an exact all-vs-all but can be very slow.

        Returns
        -------
        tuple of length 5
            ``(cluster_members, cluster_trajs, cluster_distance_matrices,
            cluster_centroids, cluster_frames)`` where:

            * ``cluster_members`` (list of int): number of frames in each
              cluster.
            * ``cluster_trajs`` (list of ``mdtraj.Trajectory``): one
              sub-trajectory per cluster containing the assigned frames.
            * ``cluster_distance_matrices`` (list of np.ndarray): per-cluster
              all-vs-all RMSD submatrix (Angstroms).
            * ``cluster_centroids`` (list of int): for each cluster, the
              index *within that cluster* of the representative frame.
            * ``cluster_frames`` (list of np.ndarray): frame indices in the
              original trajectory that were assigned to each cluster.

        Example
        -------
        >>> members, trajs, dmats, centroids, frames = protein.get_clusters(n_clusters=5)
        >>> members
        [120, 95, 80, 60, 45]
        """

        # build an empty distance matrix
        if self.n_frames % stride == 0:
            distance_dims = int(self.n_frames/stride)
        else:
            distance_dims = int((self.n_frames/stride))+1

        distances = np.zeros((distance_dims, distance_dims))

        idx = 0

        # build an all vs. all RMSD matrix based on the parameters provided for every
        # stride-th frame
        for i in range(0, self.n_frames, stride):
            distances[idx] = self.get_RMSD(i, stride=stride, region=region, backbone=backbone)
            idx=idx+1

        # CLUSTERING
        # having computed the RMSD distance matrix we do Ward based hierachical clustering
        # and then separate out into n_clusters

        # we feed ward a redundant distance matrix
        # See: http://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.ward.html#scipy.cluster.hierarchy.ward
        linkage = scipy.cluster.hierarchy.ward(distances)

        # linkage is the hierachical clustering encoded as a linkage matrix
        labels = scipy.cluster.hierarchy.fcluster(linkage, t=n_clusters, criterion='maxclust')

        # get a subtrajectory which corresponds to the trajectory examined
        # in the all vs. all comparison (i.e. a trajectory made of every stride-th
        # frame
        subtraj = self.__get_subtrajectory(self.traj, stride)

        # if we're looking at a region further extract out ONLY the atoms
        # associated with that subregion
        if region is not None:
            selectionatoms = self.__get_selection_atoms(region, backbone)
            subtraj = subtraj.atom_slice(selectionatoms)

        # we now build n_cluster separate trajectories contaning conformations from the clustering
        cluster_trajs = []
        cluster_distance_matricies = []
        cluster_centroids = []
        cluster_members = []
        cluster_frames = []

        # regardless of how many clusters we *think* we should have, extract the number of labeles
        # we'll actually have...
        final_labels = list(set(labels))

        # for each cluster
        for i in final_labels:

            # determine the indices associated with frames which are associated
            # with the i-th cluster
            IDXs=np.where(labels == i)
            IDXs=IDXs[0]
            cluster_frames.append(IDXs)

            # record how many frames are associated with the i-th cluster
            cluster_members.append(len(IDXs))

            # create the trajectory and append to the cluster_trajectory list
            cluster_trajs.append(subtraj.slice(IDXs))

            # create the appropriate submatrix
            cluster_distances=np.zeros((len(IDXs), len(IDXs)))

            # initial distances is a subset of the all vs all RMSD distance
            # matrix which now only includes rows associated with the frames
            # from the i-th cluster
            initial_dist=distances[IDXs]

            # for each frame associated with the i-th cluster
            for k in range(0,len(IDXs)):

                # add the full all vs. all set of distances between each of the frames from
                # the i-th cluster and all the other frames from the i-th clster
                cluster_distances[k] = initial_dist[k][IDXs]

            # at this point the cluster_distances matrix is an all vs. all RMSD distance
            # matrix for all the frames in i-th cluster - this effectivly gives you a way
            # to think about how well an RMSD cluster represents those structures
            cluster_distance_matricies.append(cluster_distances)

            # find the frame closest to the cluster centroid; when std==0 all frames are
            # equidistant (single-frame or all-identical cluster) so any index is valid
            std = cluster_distances.std()
            if std == 0:
                cluster_centroids.append(0)
            else:
                cluster_centroids.append(np.exp(-1*cluster_distances / std).sum(axis=1).argmax())

        return (cluster_members, cluster_trajs, cluster_distance_matricies, cluster_centroids, cluster_frames)


    # ........................................................................
    #
    #
    def get_inter_residue_COM_distance(self, R1, R2, stride=1):
        """Per-frame distance between two residues' centres of mass.

        Distances are in Angstroms.

        Parameters
        ----------
        R1 : int
            Resid of the first residue.
        R2 : int
            Resid of the second residue.
        stride : int, optional
            Use every ``stride``-th frame. Default is 1.

        Returns
        -------
        np.ndarray
            1D array of length ``ceil(n_frames / stride)`` with the
            inter-residue COM distance for each (strided) frame, in Angstroms.

        Example
        -------
        >>> d = protein.get_inter_residue_COM_distance(5, 45)
        >>> d.mean(), d.std()
        (12.3, 1.1)
        """

        # get COM of the two residues for every stride-th frame
        COM_1 = self.get_residue_COM(R1)[::stride]
        COM_2 = self.get_residue_COM(R2)[::stride]


            
        # calculate distance
        # note 10* to get angstroms was done at the get_residue_COM level
        distances = np.linalg.norm(COM_1 - COM_2, axis=1)
        
        # old way
        #d = np.sqrt(np.square(np.transpose(COM_1)[0] - np.transpose(COM_2)[0]) + np.square(np.transpose(COM_1)[1] - np.transpose(COM_2)[1])+np.square(np.transpose(COM_1)[2] - np.transpose(COM_2)[2]))


        return distances



    # ........................................................................
    #
    #
    def get_inter_residue_COM_vector(self, R1, R2):
        """Per-frame displacement *vector* ``COM(R1) - COM(R2)`` between residue COMs.

        Unlike :meth:`get_inter_residue_COM_distance`, this returns the full
        ``(dx, dy, dz)`` separation rather than its magnitude. The convention
        is ``COM(R1) - COM(R2)``, i.e. the vector that points from residue
        ``R2``'s COM to residue ``R1``'s COM. Units are Angstroms.

        .. note:: Units changed in v0.2.0; prior to that release the vector was
            returned in nanometres.

        Parameters
        ----------
        R1 : int
            Resid of the first residue (the "head" of the vector).
        R2 : int
            Resid of the second residue (the "tail" of the vector).

        Returns
        -------
        np.ndarray
            Array of shape ``(n_frames, 3)`` giving the displacement
            ``COM(R1) - COM(R2)`` in each frame, in Angstroms.

        Example
        -------
        >>> v = protein.get_inter_residue_COM_vector(5, 45)
        >>> v.shape
        (1000, 3)
        """

        COM_1 = self.get_residue_COM(R1)
        COM_2 = self.get_residue_COM(R2)

        # note 10* to get Angstroms is done at the get_residue_COM level
        return (COM_1 - COM_2)


    # ........................................................................
    #
    #
    def get_inter_residue_atomic_distance(self, R1, R2, A1='CA', A2='CA', mode='atom', stride=1):
        """Per-frame distance between two residues, with mode-dependent semantics.

        For ``mode='atom'`` the distance is between the named atoms ``A1`` of
        residue ``R1`` and ``A2`` of residue ``R2``. For all other modes the
        ``A1``/``A2`` arguments are ignored and the distance is the one
        returned by ``mdtraj.compute_contacts`` under the chosen scheme.

        SOURSOP does no sanity checking on atom names — if a name is invalid,
        an SSException is raised with a hint about likely mistypes.

        Distance is returned in Angstroms.

        Parameters
        ----------
        R1, R2 : int
            Resids of the two residues to measure between.
        A1, A2 : str, optional
            For ``mode='atom'``, the atom names within ``R1`` and ``R2``.
            Defaults are ``'CA'``. Ignored for all other modes.
        mode : {'atom', 'ca', 'closest', 'closest-heavy', 'sidechain', 'sidechain-heavy'}, optional
            * ``'atom'`` (default) - distance between atoms ``A1`` and ``A2``.
            * ``'ca'`` - equivalent to ``mode='atom'``, ``A1=A2='CA'``.
            * ``'closest'`` - closest approach between any atom of R1 and any
              atom of R2 in each frame.
            * ``'closest-heavy'`` - same as ``'closest'`` but ignoring
              hydrogens.
            * ``'sidechain'`` - closest sidechain atom of R1 to closest
              sidechain atom of R2 (requires mdtraj >= 1.8.0).
            * ``'sidechain-heavy'`` - sidechain heavy-atom version
              (requires mdtraj >= 1.8.0). Raises for glycine residues.
        stride : int, optional
            Use every ``stride``-th frame. Default is 1.

        Returns
        -------
        np.ndarray
            1D array of length ``ceil(n_frames / stride)`` with the
            inter-residue distance for each (strided) frame, in Angstroms.

        Raises
        ------
        SSException
            If ``mode`` is invalid or if an atom name cannot be found in the
            corresponding residue.

        Example
        -------
        >>> d = protein.get_inter_residue_atomic_distance(5, 45)               # CA-CA by default
        >>> d_min = protein.get_inter_residue_atomic_distance(5, 45, mode='closest-heavy')
        """

        # check mode keyword is valid
        ssutils.validate_keyword_option(mode, ['atom', 'ca', 'closest-heavy', 'closest', 'sidechain', 'sidechain-heavy'] , 'mode')

        ## if mode is atom...
        ##

        if mode == 'atom':
            try:
                atom1 = self.__residue_atom_lookup(R1,A1)
                if len(atom1) == 0:
                    raise SSException(f'Unable to find atom [{A1}] in residue R1 ({R1})')


                TRJ_1 = self.traj.atom_slice(atom1)
                TRJ_1 = self.__get_subtrajectory(TRJ_1, stride)

                atom2 = self.__residue_atom_lookup(R2,A2)
                if len(atom2) == 0:
                    raise SSException(f'Unable to find atom [{A2}] in residue R1 ({R2})')

                TRJ_2 = self.traj.atom_slice(atom2)
                TRJ_2 = self.__get_subtrajectory(TRJ_2, stride)

                COM_1 = 10*md.compute_center_of_mass(TRJ_1)
                COM_2 = 10*md.compute_center_of_mass(TRJ_2)


                # 
                distances = np.linalg.norm(COM_1 - COM_2, axis=1)

                # old way
                #distances = np.sqrt(np.square(np.transpose(COM_1)[0] - np.transpose(COM_2)[0]) + np.square(np.transpose(COM_1)[1] - np.transpose(COM_2)[1])+np.square(np.transpose(COM_1)[2] - np.transpose(COM_2)[2]))

            except IndexError as e:
                ssio.exception_message(f"This is likely because one of [{A1}] or [{A2}] is not a valid atom type for the residue in question. Full error printed below\n{e}", e, with_frills=True)

        # if ANY of the other modes are passed
        else:
            subtraj = self.__get_subtrajectory(self.traj, stride)
            distances = 10*md.compute_contacts(subtraj, [[R1,R2]], scheme=mode, periodic=False)[0].ravel()


        return distances



    # ........................................................................
    #
    #
    def get_residue_mass(self, R1):
        """Total mass of every atom belonging to a residue.

        The mass is summed from each atom's ``element.mass`` as recorded in the
        underlying mdtraj topology, in atomic mass units (g/mol).

        Parameters
        ----------
        R1 : int
            Resid whose mass should be summed.

        Returns
        -------
        float
            Total mass of residue ``R1`` in atomic mass units.

        Example
        -------
        >>> protein.get_residue_mass(5)   # e.g. leucine
        113.16
        """

        # get the atoms associated with the resite of interest
        res_atoms = self.topology.select(f'resid {R1}')

        totalMass = 0

        for atom in res_atoms:
            totalMass = totalMass + self.topology.atom(atom).element.mass

        return totalMass


    # ........................................................................
    #
    #
    def get_asphericity(self, R1=None, R2=None, verbose=True):
        """Per-frame asphericity of the chain (or a sub-region).

        Asphericity is a dimensionless shape descriptor computed from the
        eigenvalues of the gyration tensor: it is 0 for a perfectly spherical
        distribution of mass and approaches 1 as the chain becomes
        increasingly extended/rod-like.

        See p. 65 of Andreas Vitalis' thesis (*Probing the Early Stages of
        Polyglutamine Aggregation with Computational Methods*, 2009, WashU)
        for a canonical derivation.

        Parameters
        ----------
        R1 : int, optional
            First residue of the region to consider. ``None`` (default) means
            the first residue in the sequence (including caps).
        R2 : int, optional
            Last residue of the region to consider. ``None`` (default) means
            the last residue in the sequence (including caps).
        verbose : bool, optional
            If True (default), print one status line every 500 frames during
            the gyration-tensor computation.

        Returns
        -------
        np.ndarray
            1D array of length ``n_frames`` with the per-frame asphericity.

        Example
        -------
        >>> a = protein.get_asphericity(verbose=False)
        >>> a.mean()
        0.18
        """

        # compute the gyration tensor NOTE that R1 and R2 rationalization are done in
        # the get_gyration_tensor function
        gyration_tensor_vector = self.get_gyration_tensor(R1, R2, verbose=verbose)

        # finally for each gyration tensor value compute the asphericity
        asph_vector = []
        for gyr in gyration_tensor_vector:

            # calculate the eigenvalues of the gyration tensor!
            (EIG, norm) = LA.eig(gyr)

            # finally calculate the instantanous asphericity and append to the growing vector
            asph = 1 - 3*((EIG[0]*EIG[1] + EIG[1]*EIG[2] + EIG[2]*EIG[0])/np.power(EIG[0]+EIG[1]+EIG[2],2))

            asph_vector.append(asph)

        return np.array(asph_vector)


    # ........................................................................
    #
    #
    def get_gyration_tensor(self, R1=None, R2=None, verbose=True):
        """Per-frame ``3 x 3`` gyration tensor of the chain (or a sub-region).

        The gyration tensor :math:`T_{ab} = (1/N) \\sum_i (r_{i,a} - R_a)(r_{i,b} - R_b)`
        is the second moment of mass about the centre of mass. Its
        eigenvalues are the principal moments of inertia divided by mass and
        underlie derived quantities like :meth:`get_radius_of_gyration` and
        :meth:`get_asphericity`. Units are Angstroms^2.

        Parameters
        ----------
        R1 : int, optional
            First residue of the region to consider. ``None`` (default) means
            the first residue (including caps).
        R2 : int, optional
            Last residue of the region to consider. ``None`` (default) means
            the last residue (including caps).
        verbose : bool, optional
            If True (default), print one status line every 500 frames.

        Returns
        -------
        np.ndarray
            Array of shape ``(n_frames, 3, 3)`` where ``[f]`` is the gyration
            tensor for frame ``f``.

        Example
        -------
        >>> T = protein.get_gyration_tensor(verbose=False)
        >>> T.shape
        (1000, 3, 3)
        """
        (R1,R2, _) = self.__get_first_and_last(R1,R2, withCA=False)

        all_positions_all_frames = self.traj.atom_slice(self.topology.select(f'resid {R1} to {R2}'))

        gyration_tensor_vector = []
        count = 1

        for frame in all_positions_all_frames:

            # quick status update...
            if count % 500 == 0:
                ssio.status_message(f"On frame {count} of {len(all_positions_all_frames)} [computing gyration tensor]", verbose)

            count = count + 1

            # compute the center of mass for the relevant atoms
            COM = 10*md.compute_center_of_mass(frame)

            # calculate the COM to position difference matrix. Note
            # we have to multiply the frame positions by 10 so that the COM
            # units and the frame units match up.
            ## QUICK ASSIDE: we could actually leave both in nm because the gyration tensor ends up being
            ## unitless. HOWEVER, for consistency we're following the convention of converting ANY
            ## nm-based values to angstroms at the point at which the numbers enter soursop...

            DIF = 10*frame.xyz[0] - COM

            # old slow method
            """
            T_PRE = 0.0
            # for each atom in the frame calculate the outer product (dyadic product) of
            # the difference between the position at the overall center of mass. This creates a
            # 3x3 gyration tensor
            for pos in frame.xyz[0]:
                T_PRE = T_PRE + np.outer(pos - COM, pos - COM)

            T = T_PRE/len(frame.xyz[0])
            """

            # compute the gyration tensor - this syntax is WAY faster than the above
            T_new = np.sum(np.einsum('ij...,i...->ij...',DIF,DIF),axis=0)/len(frame.xyz[0])

            gyration_tensor_vector.append(T_new)


        return np.array(gyration_tensor_vector)


    # ........................................................................
    #
    #
    def get_end_to_end_distance(self, mode='COM'):
        """Per-frame end-to-end distance between the first and last CA residues.

        Caps (ACE / NME) are excluded — the "ends" are the first and last
        residues with a CA atom, so the two modes operate on the same pair of
        residues.

        End-to-end distance is returned in Angstroms.

        Parameters
        ----------
        mode : {'CA', 'COM'}, optional
            * ``'COM'`` (default) - distance between residue centres of mass.
            * ``'CA'`` - distance between alpha-carbon atoms.

        Returns
        -------
        np.ndarray
            1D array of length ``n_frames`` with the per-frame end-to-end
            distance in Angstroms.

        Example
        -------
        >>> ee = protein.get_end_to_end_distance(mode='CA')
        >>> ee.mean(), ee.std()
        (30.1, 4.4)
        """

        ssutils.validate_keyword_option(mode, ['CA', 'COM'], 'mode')

        # extract first and last residue that contain CA
        start = self.resid_with_CA[0]
        end = self.resid_with_CA[-1]

        if mode == 'CA':
            distance = self.get_inter_residue_atomic_distance(start, end, stride=1)

        elif mode == 'COM':
            distance = self.get_inter_residue_COM_distance(start, end, stride=1)

        return distance


    # ........................................................................
    #
    #
    def get_center_of_mass(self, R1=None, R2=None):
        """Per-frame centre of mass of the protein (or a sub-region).

        Returns absolute ``(x, y, z)`` positions in the trajectory's frame of
        reference, in Angstroms.

        .. warning:: Units changed in v0.2.0; prior to that release the values
            were returned in nanometres. This is a breaking change.

        Parameters
        ----------
        R1 : int, optional
            First residue of the region. ``None`` (default) means the first
            residue including caps.
        R2 : int, optional
            Last residue of the region. ``None`` (default) means the last
            residue including caps.

        Returns
        -------
        np.ndarray
            Array of shape ``(n_frames, 3)`` giving the COM ``(x, y, z)`` for
            each frame, in Angstroms.

        Example
        -------
        >>> com = protein.get_center_of_mass()
        >>> com.shape
        (1000, 3)
        """

        (R1_new, R2_new, selection_string) = self.__get_first_and_last(R1, R2, withCA=False)

        return 10*md.compute_center_of_mass(self.traj.atom_slice(self.topology.select(selection_string)))


    # ........................................................................
    #
    #
    def get_radius_of_gyration(self, R1=None, R2=None):
        """Per-frame radius of gyration of the chain (or a sub-region).

        :math:`R_g = \\sqrt{(1/N) \\sum_i (r_i - R)^2}` is the canonical
        polymer-physics measure of overall chain size. It is computed by
        ``mdtraj.compute_rg`` on every heavy/light atom in the selection and
        returned in Angstroms.

        Parameters
        ----------
        R1 : int, optional
            First residue of the region. ``None`` (default) means the first
            residue including caps.
        R2 : int, optional
            Last residue of the region. ``None`` (default) means the last
            residue including caps.

        Returns
        -------
        np.ndarray
            1D array of length ``n_frames`` with the per-frame :math:`R_g` in
            Angstroms.

        Example
        -------
        >>> rg = protein.get_radius_of_gyration()
        >>> rg.mean()
        33.4
        """

        (_, _, selection_string) = self.__get_first_and_last(R1,R2, withCA=False)

        # in angstroms
        return 10*md.compute_rg(self.traj.atom_slice(self.topology.select(selection_string)))


    # ........................................................................
    #
    #
    def get_hydrodynamic_radius(self, R1=None, R2=None, mode='nygaard', alpha1=0.216, alpha2=4.06, alpha3=0.821, distance_mode='CA'):
        """Per-frame apparent hydrodynamic radius :math:`R_h` (in Angstroms).

        Two estimators are available:

        * ``mode='nygaard'`` - empirical Rg-to-Rh conversion of Nygaard et al.
          (2017), using the three parameters ``alpha1``, ``alpha2``, ``alpha3``.
        * ``mode='kr'`` - Kirkwood-Riseman equation (Kirkwood & Riseman 1948),
          recommended for comparison with PFG-NMR-derived :math:`R_h` values
          (Pesce et al. 2022).

        Both approximations hold for fully flexible disordered proteins and
        become less accurate for larger / more folded domains.

        Parameters
        ----------
        R1 : int, optional
            First residue of the region. ``None`` (default) means the first
            residue including caps.
        R2 : int, optional
            Last residue of the region. ``None`` (default) means the last
            residue including caps.
        mode : {'nygaard', 'kr'}, optional
            Which estimator to use. Default is ``'nygaard'``.
        alpha1, alpha2, alpha3 : float, optional
            Parameters of equation (7) of Nygaard et al. Defaults reproduce
            the published values (0.216, 4.06, 0.821). Ignored for ``mode='kr'``.
        distance_mode : {'CA', 'COM'}, optional
            For ``mode='kr'``, which inter-residue distance to feed into the
            Kirkwood-Riseman sum. Default is ``'CA'``. Ignored for
            ``mode='nygaard'``.

        Returns
        -------
        np.ndarray
            1D array of length ``n_frames`` with the per-frame :math:`R_h` in
            Angstroms.

        Example
        -------
        >>> rh_nyg = protein.get_hydrodynamic_radius()
        >>> rh_kr  = protein.get_hydrodynamic_radius(mode='kr')
        >>> rh_nyg.mean(), rh_kr.mean()
        (24.1, 23.8)

        References
        ----------
        Nygaard M, Kragelund BB, Papaleo E, Lindorff-Larsen K.
        An Efficient Method for Estimating the Hydrodynamic Radius of
        Disordered Protein Conformations. Biophys J. 2017;113:550-557.

        Kirkwood JG, Riseman J. The Intrinsic Viscosities and Diffusion
        Constants of Flexible Macromolecules in Solution. J Chem Phys.
        1948;16(6):565-573.

        Pesce F et al. Assessment of models for calculating the
        hydrodynamic radius of intrinsically disordered proteins.
        Biophysical Journal. 2022. doi:10.1016/j.bpj.2022.12.013
        """

        # check a valid mode was passed and FREAK OUT if not!
        ssutils.validate_keyword_option(mode, ['nygaard', 'kr'], 'mode')


        # if we're using the nygaard mode
        if mode == 'nygaard':

            # first compute the rg
            rg = self.get_radius_of_gyration(R1, R2)

            # precompute
            N_033 = np.power(self.n_residues, 0.33)
            N_060 = np.power(self.n_residues, 0.60)

            Rg_over_Rh = ((alpha1*(rg - alpha2*N_033)) / (N_060 - N_033)) + alpha3

            return (1/Rg_over_Rh)*rg

        # if we're using the Kirkwood-Riseman mode
        elif mode == 'kr':
            all_rij = []

            # build empty lists associated with each frame
            for _ in range(self.n_frames):
                all_rij.append([])

            # now use our efficient implementation for calculating non-redundant
            # and non-overlaping inter-residue distances for each residue over
            # every frame
            for idx in self.resid_with_CA:
                rij = self.calculate_all_CA_distances(idx, mode=distance_mode)

                # this breaks rij down into each frame
                for idx, f in enumerate(rij):

                    # extend the contribugion for each residue's set of non-redundant
                    # inverse distances
                    all_rij[idx].extend((1/f).tolist())


            # finally, take per-frame inverse of the inverse distance
            # to get Rh
            Rh = np.reciprocal(np.mean(all_rij, axis=1).astype(float))

            return Rh




    # ........................................................................
    #
    #
    def get_molecular_volume(self, **kwargs):
        """Per-frame molecular volume from the convex hull of atomic positions.

        Computes the volume of ``scipy.spatial.ConvexHull(xyz)`` for each
        frame's atomic coordinates and converts from nm^3 to Angstrom^3.

        Parameters
        ----------
        **kwargs : dict, optional
            Forwarded to ``scipy.spatial.ConvexHull`` (e.g. ``qhull_options``).

        Returns
        -------
        np.ndarray
            1D array of length ``n_frames`` with the per-frame convex-hull
            volume in Angstrom^3.

        Example
        -------
        >>> v = protein.get_molecular_volume()
        >>> v.mean()
        1.45e+05
        """
        NM_TO_ANGSTROM = 1000 # 1000 A^3 / nm^3

        # in angstroms cubued
        volumes = np.array([ConvexHull(xyz,**kwargs).volume for xyz in self.traj.xyz])*NM_TO_ANGSTROM

        return volumes


    # ........................................................................
    #
    #
    def get_t(self, R1=None, R2=None):
        r"""Per-frame dimensionless size parameter :math:`\langle t \rangle`.

        :math:`\langle t \rangle` is the Vitalis-thesis size parameter
        (Vitalis 2009) that scales the radius of gyration against the contour
        length :math:`L = 3.6 N` (Angstroms). It is dimensionless and useful
        for comparing chains of different lengths.

        Parameters
        ----------
        R1 : int, optional
            First residue of the region. ``None`` (default) means the first
            residue.
        R2 : int, optional
            Last residue of the region. ``None`` (default) means the last
            residue.

        Returns
        -------
        np.ndarray
            1D array of length ``n_frames`` with the per-frame
            :math:`\langle t \rangle`.

        Example
        -------
        >>> t = protein.get_t()
        >>> t.mean()
        2.4

        References
        ----------
        Vitalis, A. (2009). Probing the Early Stages of Polyglutamine
        Aggregation with Computational Methods (R. Pappu, ed.)
        [Ph.D. thesis from Washington University in St. Louis].
        """

        # first get the instantanoues RG
        rg = self.get_radius_of_gyration(R1, R2)

        n_res = self.n_residues
        c_length = n_res * 3.6

        # next define the exponent
        exponent = 4.0/(np.power(n_res,0.3333))

        # define a function which returns the instantaneous t value for
        # a given rg
        def inst_t(i):
            return 2.5*np.power((1.75*(i/(c_length))),exponent)

        # compile into a vectorized version
        K = np.vectorize(inst_t)

        # run over all rg
        return K(rg)


    # ........................................................................
    #
    #
    def get_angles(self, angle_name):
        """Per-frame backbone or sidechain dihedral angles for every applicable residue.

        Computes the named dihedral using the corresponding mdtraj helper
        (``compute_phi``, ``compute_psi``, ..., ``compute_chi5``) and converts
        the result to degrees. The first call for a given ``angle_name`` is
        memoised so repeated calls are essentially free.

        Parameters
        ----------
        angle_name : {'phi', 'psi', 'omega', 'chi1', 'chi2', 'chi3', 'chi4', 'chi5'}
            Dihedral to compute. ``chi3``-``chi5`` only exist for residues
            with long sidechains (Arg, Lys, Met, etc.), so the result will
            contain fewer entries than ``n_residues``.

        Returns
        -------
        list of [list_of_mdtraj_atoms, np.ndarray]
            A two-element list:

            * ``[0]`` is a list-of-lists; each sub-list holds the four
              ``mdtraj.Atom`` objects that define one dihedral.
            * ``[1]`` is an array of shape ``(n_dihedrals, n_frames)`` of the
              per-frame angles in degrees.

        Raises
        ------
        SSException
            If ``angle_name`` is not one of the eight allowed strings.

        Example
        -------
        For a protein with two arginine residues:

        >>> chi5 = protein.get_angles('chi5')
        >>> chi5[0][1]                  # atoms defining 2nd arginine's chi5
        >>> chi5[1][1]                  # per-frame chi5 angles for that residue
        """


        # check input ketword selector
        ssutils.validate_keyword_option(angle_name, ['chi1', 'chi2','chi3', 'chi4','chi5','phi', 'psi', 'omega'], 'angle_name')

        # construct the selector dictionary - enables mapping of an angle name to an internal function
        selector = {"phi":md.compute_phi,
                    "omega":md.compute_omega,
                    "psi":md.compute_psi,
                    "chi1":md.compute_chi1,
                    "chi2":md.compute_chi2,
                    "chi3":md.compute_chi3,
                    "chi4":md.compute_chi4,
                    "chi5":md.compute_chi5}

        if angle_name not in self.__all_angles:

            # compute all phi angles
            fx = selector[angle_name]
            tmp = fx(self.traj)

            # extract out atomic indices and angles
            atoms = tmp[0]
            angles = np.rad2deg(tmp[1].transpose())

            # sanity check
            if len(atoms) != len(angles):
                raise SSException('ERROR: This is a bug in how angle selection happens in get_angles - please open an Issue on GitHub')

            all_atom_names = []
            for res_atoms in atoms:
                local_atom_names = []

                for atom_idx in res_atoms:
                    local_atom_names.append(self.traj.topology.atom(atom_idx))
                all_atom_names.append(local_atom_names)

            self.__all_angles[angle_name] = [all_atom_names, angles]

        return self.__all_angles[angle_name]



    # ........................................................................
    #
    #
    def __ca_position_cache(self, stride):
        """Return a closure ``r -> 10*COM`` of residue ``r``'s CA atom over
        the strided frames, computed once per residue (per-call cache).

        This reproduces, byte-for-byte, the CA branch of
        ``get_inter_residue_atomic_distance`` (mode='atom', A1=A2='CA'):
        the same ``__residue_atom_lookup`` -> ``atom_slice`` ->
        ``__get_subtrajectory`` -> ``10*md.compute_center_of_mass``
        sequence. The only difference is the per-residue result is cached,
        so an O(n^2) pair loop performs O(n) CA-COM computations instead
        of O(n^2). Lazy caching also preserves the exact first-error
        location and message for a residue lacking a CA atom.
        """
        cache = {}

        def _ca(r):
            if r not in cache:
                atom = self.__residue_atom_lookup(r, 'CA')
                if len(atom) == 0:
                    raise SSException(f'Unable to find atom [CA] in residue R1 ({r})')
                TRJ = self.traj.atom_slice(atom)
                TRJ = self.__get_subtrajectory(TRJ, stride)
                cache[r] = 10 * md.compute_center_of_mass(TRJ)
            return cache[r]

        return _ca



    # ........................................................................
    #
    #
    def get_internal_scaling(self, R1=None, R2=None, mode='COM', mean_vals=False, stride=1, weights=False, etol=0.0000001, verbose=True):
        """Internal scaling profile: sequence separation vs. inter-residue distance.

        For every sequence separation :math:`|i-j|` from 0 to ``R2-R1``, this
        method collects the inter-residue distances of every pair at that
        separation across every frame. The result is the raw material for an
        internal-scaling plot — a foundational polymer-physics observable for
        IDPs and IDRs (Mao et al. 2013; Pappu et al. 2008).

        ACE / NME peptide caps are excluded. (Note: this differs from CAMPARI,
        which includes them.)

        Distances are in Angstroms.

        Parameters
        ----------
        R1 : int, optional
            First residue of the region. ``None`` (default) means the first
            CA-bearing residue.
        R2 : int, optional
            Last residue of the region. ``None`` (default) means the last
            CA-bearing residue.
        mode : {'COM', 'CA'}, optional
            * ``'COM'`` (default) - distances between residue centres of mass.
            * ``'CA'`` - distances between alpha-carbon atoms.
        mean_vals : bool, optional
            If True, the per-separation distances are reduced to their mean
            before being returned. If False (default), every individual
            inter-residue distance is kept (useful for plotting the full
            distribution).
        stride : int, optional
            Use every ``stride``-th frame. Default is 1. Larger values are
            recommended for long trajectories of large proteins.
        weights : list or np.ndarray, optional
            Per-frame weights for re-weighted averaging. Default ``False``
            means uniform weighting.
        etol : float, optional
            Tolerance for the weights-sum-to-one check. Default ``1e-7``.
        verbose : bool, optional
            If True (default), print one status line per sequence separation.

        Returns
        -------
        tuple of (list, list)
            ``(seq_sep, distances)`` where:

            * ``seq_sep`` is a list of integer sequence separations
              ``[0, 1, ..., R2-R1]``.
            * If ``mean_vals=True``, ``distances`` is a list of per-separation
              mean distances (one float per separation).
            * If ``mean_vals=False`` (default), ``distances`` is a list of 1D
              arrays — one array per separation containing every individual
              inter-residue distance at that separation.

        Example
        -------
        >>> sep, mean_d = protein.get_internal_scaling(mean_vals=True, verbose=False)
        >>> # plot sep vs mean_d for an internal-scaling profile

        References
        ----------
        Mao AH, Lyle N, Pappu RV. Describing sequence-ensemble
        relationships for intrinsically disordered proteins. Biochem J.
        2013;449:307-318.

        Pappu RV, Wang X, Vitalis A, Crick SL. A polymer-physics
        perspective on driving forces and mechanisms for protein
        aggregation. Arch Biochem Biophys. 2008;469:132-141.

        Notes
        -----
        Computational cost scales rapidly with chain length; consider
        increasing ``stride`` for trajectories of long proteins.
        """

        if weights is not False:
            if int(stride) != 1:
                raise SSException("For get_scaling_exponent with weights stride MUST be set to 1. If this is a HUGE deal for you please contact alex and he'll try and update the code to accomodate this, but for now we suggest creating a sub-sampled trajectory and loading that")

        # check weights are correct
        weights = self.__check_weights(weights, stride, etol)

        # check stride is ok
        self.__check_stride(stride)

        # check mode is OK
        ssutils.validate_keyword_option(mode, ['CA', 'COM'], 'mode')

        # process the R1/R2 to set the positions
        out =  self.__get_first_and_last(R1, R2, withCA = True)
        R1 = out[0]
        R2 = out[1]

        max_seq_sep = (R2 - R1) + 1

        # if chain is too short...
        if max_seq_sep < 1:
            return ([], [])

        seq_sep_distances = []
        seq_sep_vals = []

        # Hoist the per-residue CA position out of the O(n^2) pair loop.
        # get_inter_residue_atomic_distance recomputes each residue's CA
        # center-of-mass from scratch on every call; caching it per
        # residue is byte-identical but O(n) rather than O(n^2).
        if mode == 'CA':
            _ca = self.__ca_position_cache(stride)

        for seq_sep in range(0, max_seq_sep):
            ssio.status_message(f"Internal Scaling - on sequence separation {seq_sep} of {max_seq_sep-1}", verbose)

            tmp = []
            seq_sep_vals.append(seq_sep)
            for pos in range(0, max_seq_sep-seq_sep):

                # define the two positions
                A = R1 + pos
                B = R1 + pos + seq_sep

                # get the distance for every stride-th frame between those two positions using either the CA
                # mode or the COM mode
                if mode == 'CA':
                    # byte-identical to get_inter_residue_atomic_distance(A, B, stride=stride)
                    distance = np.linalg.norm(_ca(A) - _ca(B), axis=1)

                elif mode == 'COM':
                    distance = self.get_inter_residue_COM_distance(A, B, stride=stride)

                # if weights were provided subsample from the set of distances using the weights vector
                if weights is not False:
                    distance = choice(distance, len(distance), p=weights)

                tmp = np.concatenate((tmp,distance))

            seq_sep_distances.append(tmp)

        if mean_vals:
            mean_is = [np.mean(i) for i in seq_sep_distances]
            return (seq_sep_vals, mean_is)
        else:
            return (seq_sep_vals, seq_sep_distances)



    # ........................................................................
    #
    #
    def get_internal_scaling_RMS(self, R1=None, R2=None, mode='COM', stride=1, weights=False, etol=0.0000001, verbose=True):
        """Root-mean-square internal scaling profile: separation vs. RMS distance.

        Like :meth:`get_internal_scaling` but reports
        :math:`\\sqrt{\\langle r_{ij}^2 \\rangle}` for each sequence
        separation — the formally correct order parameter for polymer-physics
        analyses (Mao et al. 2013; Pappu et al. 2008). ACE / NME caps are
        excluded. Distances in Angstroms.

        Parameters
        ----------
        R1 : int, optional
            First residue of the region. ``None`` (default) means the first
            CA-bearing residue.
        R2 : int, optional
            Last residue of the region. ``None`` (default) means the last
            CA-bearing residue.
        mode : {'COM', 'CA'}, optional
            * ``'COM'`` (default) - distances between residue centres of mass.
            * ``'CA'`` - distances between alpha-carbon atoms.
        stride : int, optional
            Use every ``stride``-th frame. Default is 1.
        weights : list or np.ndarray, optional
            Per-frame weights for re-weighted averaging. Default ``False``.
        etol : float, optional
            Tolerance for the weights-sum-to-one check. Default ``1e-7``.
        verbose : bool, optional
            If True (default), print one status line per sequence separation.

        Returns
        -------
        tuple of (list, list)
            ``(seq_sep, rms_distance)`` where ``seq_sep`` is the list of
            integer separations and ``rms_distance`` is the corresponding
            list of RMS inter-residue distances in Angstroms.

        Example
        -------
        >>> sep, rms = protein.get_internal_scaling_RMS(verbose=False)
        >>> import numpy as np
        >>> # log-log fit of rms vs sep recovers the apparent scaling exponent
        >>> nu, _ = np.polyfit(np.log(sep[1:]), np.log(rms[1:]), 1)

        References
        ----------
        Mao AH, Lyle N, Pappu RV. Describing sequence-ensemble
        relationships for intrinsically disordered proteins. Biochem J.
        2013;449:307-318.

        Pappu RV, Wang X, Vitalis A, Crick SL. A polymer-physics
        perspective on driving forces and mechanisms for protein
        aggregation. Arch Biochem Biophys. 2008;469:132-141.
        """

        # compute the non RMS internal scaling behaviour
        (seq_sep_vals, seq_sep_distances) = self.get_internal_scaling(R1=R1, R2=R2, mode=mode, mean_vals=False, stride=stride, weights=weights, etol=etol, verbose=verbose)

        # calculate RMS for each distance
        mean_is = [np.sqrt(np.mean(i*i)) for i in seq_sep_distances]

        return (seq_sep_vals, mean_is)


    # ........................................................................
    #
    def get_scaling_exponent(self,
                             inter_residue_min=15,
                             end_effect=5,
                             subdivision_batch_size=20,
                             mode='COM',
                             num_fitting_points=40,
                             fraction_of_points=0.5,
                             fraction_override=False,
                             etol=0.0000001,
                             stride=1,
                             weights=False,
                             verbose=True):
        """Fit the apparent polymer scaling exponent and prefactor by loglog regression.

        Fits the homopolymer scaling relationship
        :math:`\\sqrt{\\langle r_{ij}^2 \\rangle} = A_0 |i-j|^{\\nu_{app}}`
        to the internal-scaling profile of the chain. The exponent
        :math:`\\nu_{app}` ("nu-app") reports on solvent quality and the
        prefactor :math:`A_0` captures chain stiffness / segment volume.

        Bootstrap error estimates for ``nu`` and ``A0`` come from randomly
        subdividing the trajectory into ``subdivision_batch_size`` chunks,
        fitting each, and taking the min/max across the chunks.

        .. note::
            Despite the precision of the fit, ``nu_app`` and ``A0`` are
            *qualitative* metrics subject to finite-chain-size effects, and
            are most rigorously interpreted for genuine homopolymers. For
            heteropolymers consider :meth:`get_polymer_scaled_distance_map`
            which avoids the homopolymer assumption.

        Parameters
        ----------
        inter_residue_min : int, optional
            Minimum sequence separation ``|i-j|`` used in the fit. Helps
            avoid the short-distance regime where local steric effects
            dominate. Default is 15.
        end_effect : int, optional
            Exclude pairs in which one residue is within ``end_effect``
            residues of either chain end. Default is 5.
        subdivision_batch_size : int, optional
            Target chunk size when bootstrapping fit-quality errors.
            Default is 20.
        mode : {'COM', 'CA'}, optional
            Distance type. ``'COM'`` (default) is preferred — CA-CA tends to
            inflate the apparent profile.
        num_fitting_points : int, optional
            Maximum number of evenly-spaced loglog points used in the linear
            fit. Default is 40.
        fraction_of_points : float, optional
            When ``fraction_override`` is True, or when the chain is too
            short to provide ``num_fitting_points`` points, this fraction of
            the valid sequence separations is used instead. Default is 0.5.
        fraction_override : bool, optional
            If True, always use ``fraction_of_points`` and ignore
            ``num_fitting_points``. Default is False.
        etol : float, optional
            Tolerance for the weights-sum-to-one check. Default ``1e-7``.
        stride : int, optional
            Use every ``stride``-th frame. Default is 1.
        weights : list or np.ndarray, optional
            Per-frame weights for re-weighted averaging. Default ``False``.
        verbose : bool, optional
            If True (default), print status updates during the
            internal-scaling pass.

        Returns
        -------
        tuple of length 10
            ``(best_nu, best_A0, min_nu, max_nu, min_A0, max_A0,
            redchi_fit_region, redchi_all_points, fit_region_data,
            all_points_data)`` where:

            * ``best_nu``, ``best_A0`` - point estimates from the full fit.
            * ``min_nu``..``max_A0`` - bootstrap min/max bounds.
            * ``redchi_fit_region`` - reduced chi^2 on the points used in
              the linear fit (loglog-uniform).
            * ``redchi_all_points`` - reduced chi^2 across every observed
              :math:`|i-j|`.
            * ``fit_region_data`` - shape ``(2, K)``; col 0 is sequence
              separation, col 1 is RMS distance for the fit-region points.
            * ``all_points_data`` - shape ``(3, M)``; col 0 is sequence
              separation, col 1 is observed RMS distance, col 2 is the
              best-fit prediction.

        Raises
        ------
        SSException
            If ``fraction_of_points > 1.0`` or if too few points remain to
            fit a line (``< 3``).

        Example
        -------
        >>> import numpy as np
        >>> np.random.seed(42)   # error bootstrap consumes the RNG
        >>> result = protein.get_scaling_exponent(verbose=False)
        >>> nu, A0 = result[0], result[1]
        """
                    

        # check weights are OK
        weights = self.__check_weights(weights, stride, etol)

        # check mode is OK
        ssutils.validate_keyword_option(mode, ['CA', 'COM'], 'mode')

        # compute max |i-j| distance being used...
        first = self.resid_with_CA[0]
        last = self.resid_with_CA[-1]
        max_separation = (last-first)+1

        #  if we're not using fraction override check the number of points requested makes sense given sequence length
        if not fraction_override and (max_separation - (end_effect+inter_residue_min)) < num_fitting_points:
            fraction_of_points=1.0
            ssio.warning_message(f"For scaling exponent calculation, sequence not long enough to use {num_fitting_points} points (only {max_separation - (end_effect+inter_residue_min)} valid positions once end effects and low |i-j| are accounted for), switching to using the fraction of points mode (will use {int(fraction_of_points*(max_separation - (end_effect+inter_residue_min)))} points instead)")

            fraction_override = True

        if fraction_override:
            if fraction_of_points > 1.0:
                raise SSException("Using fraction_overide to define the number of points to fit in the linear loglog analysis, but requested over 1.0 fraction (fraction_of_points must lie between >=0 and 1.0")
            # again note int to round down here
            num_fitting_points = int(fraction_of_points*(max_separation - (end_effect + inter_residue_min)))
            if num_fitting_points < 3:
                raise SSException("Less than three points - cannot fit a straight line")
            if num_fitting_points < 10:
                ssio.warning_message(f"Warning: Scaling fit has only {num_fitting_points} points - likely finite size effects!")


        # This section determines the number of subdivisions performed for error
        # bootstrapping. If we have fewer frames than we can divide the data into
        # then we just use each frame individually (although now error bootstrapping
        # is probably meaningless!
        # note integer math used here to round down - also set
        if int(self.n_frames/stride) < int(subdivision_batch_size):
            num_subdivisions_for_error = int(self.n_frames/stride)
        else:
            num_subdivisions_for_error = int(int(self.n_frames/stride) / subdivision_batch_size)

        seq_sep_vals             = []
        seq_sep_RMS_distance     = []
        seq_sep_RMS_var_distance     = []
        seq_sep_RSTDS_distance   = []
        seq_sep_subsampled_distances  = []

        # Hoist the per-residue CA position out of the O(n^2) pair loop.
        # get_inter_residue_atomic_distance recomputes each residue's CA
        # center-of-mass from scratch on every call; caching it per
        # residue is byte-identical but O(n) rather than O(n^2). The
        # hoist consumes no RNG, so the np.random error-bootstrap below
        # sees an identical stream and its output is unchanged.
        if mode == 'CA':
            _ca = self.__ca_position_cache(stride)

        # for each possible sequence separation  (|i-j| value)
        for seq_sep in range(1, max_separation):

            ssio.status_message(f"Internal Scaling - on sequence separation {seq_sep} of {max_separation}", verbose)

            tmp = []
            seq_sep_vals.append(seq_sep)

            # collect all possible average seq-sep values weighted by the weights
            for pos in range(0, max_separation-seq_sep):

                # define the two positions
                A = first + pos
                B = first + pos + seq_sep


                if mode == 'CA':
                    # byte-identical to get_inter_residue_atomic_distance(A, B, stride=stride)
                    distance = np.linalg.norm(_ca(A) - _ca(B), axis=1)

                elif mode == 'COM':
                    distance = self.get_inter_residue_COM_distance(A, B, stride=stride)

                # compute the ensemble average of the distances
                if weights is not False:
                    distance = choice(distance, len(distance), p=weights)

                tmp.extend(distance)

            tmp = np.array(tmp)

            # add mean and std vals for this sequence sep
            seq_sep_RMS_distance.append(np.sqrt(np.mean(tmp*tmp)))
            seq_sep_RSTDS_distance.append(np.sqrt(np.std(tmp*tmp)))
            seq_sep_RMS_var_distance.append(np.sqrt(np.power(np.var(tmp),2)))

            if num_subdivisions_for_error > 0:

                # note we cast this to an int to ensure subdivision_size is always
                # the value added to RMS_local_append is always the same, because
                # len(tmp) will vary with sequence separation. Basically this means
                # we take ALL the data and subidivided it into num_subdivisions_for_error
                # chunks and then use this for error calculations
                subdivision_size = int(len(tmp)/num_subdivisions_for_error)

                # get shuffled indices
                idx = np.random.permutation(list(range(0,len(tmp))))

                # split shuffled indices into $num_subdivisions_for_error sized chunks

                subdivided_idx = sstools.chunks(idx, subdivision_size)

                # finally subselect each of the randomly selected indicies
                RMS_local   = []

                for idx_set in subdivided_idx:

                    # subselect a random set of distances and compute RMS
                    RMS_local.append(np.sqrt(np.mean(tmp[idx_set]*tmp[idx_set])))

                # add distribution of values for this sequence sep
                seq_sep_subsampled_distances.append(RMS_local)

        # now sub-select the bit of the curve we actually want for the separation, distance, and distance variance data
        # note we are RE DEFINING these three variables here
        seq_sep_vals = seq_sep_vals[inter_residue_min:-end_effect]
        seq_sep_RMS_distance = seq_sep_RMS_distance[inter_residue_min:-end_effect]
        seq_sep_RMS_var_distance = seq_sep_RMS_var_distance[inter_residue_min:-end_effect]

        ## next find indices for evenly spaced points in logspace. This whole sectino
        # leads to the identification of the indices in logspaced_idx, which are the
        # list indices that will given evenly spaced points when plotted in log space
        y_data = np.log(seq_sep_vals)
        y_data_offset = y_data - y_data[0]
        interval = y_data_offset[-1]/num_fitting_points
        integer_vals = y_data_offset/interval

        logspaced_idx = []
        for i in range(0,num_fitting_points):
            [local_ix,_] = sstools.find_nearest(integer_vals, i)
            if local_ix in logspaced_idx:
                continue
            else:
                logspaced_idx.append(local_ix)

        # finally using those evenly-spaced log indices we extract out new lists
        # that have values which will be evenly spaced in logspace. Cool.
        fitting_separation = [seq_sep_vals[i] for i in logspaced_idx]
        fitting_distances  = [seq_sep_RMS_distance[i] for i in logspaced_idx]
        fitting_variance   = [seq_sep_RMS_var_distance[i] for i in logspaced_idx]

        # fit to a log/log model and extract params
        out = np.polyfit(np.log(fitting_separation), np.log(fitting_distances), 1)
        nu_best = out[0]
        R0_best = np.exp(out[1])

        ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        ### next calculated reduced chi-squared
        n_points = len(fitting_distances)
        chi2=0
        for i in range(0, n_points):
            chi2 = chi2 + (np.power(np.log(fitting_distances[i]) - nu_best*np.log(fitting_separation[i])+R0_best,2))/fitting_variance[i]

        # finally calculated reduced chi squared correcting for 2 model parameters
        reduced_chi_squared_fitting = chi2 / (n_points-2)

        full_n_points = len(seq_sep_vals)

        chi2=0

        for i in range(0, full_n_points):
            chi2 = chi2 + (np.power(np.log(seq_sep_RMS_distance[i]) - nu_best*np.log(seq_sep_vals[i])+R0_best,2))/seq_sep_RMS_var_distance[i]

        # finally calculated reduced chi squared correcting for 2 model parameters
        reduced_chi_squared_all = chi2 / (full_n_points-2)

        ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        ### Finally run the subselection protocol to subsampled

        subselected = np.array(seq_sep_subsampled_distances).transpose()

        nu_sub = []
        R0_sub = []
        for i in range(0, num_subdivisions_for_error):

            local_distances = subselected[i][inter_residue_min:-end_effect]

            OF = np.polyfit(np.log(fitting_separation), np.log([local_distances[i] for i in logspaced_idx]), 1)

            nu_sub.append(OF[0])
            R0_sub.append(np.exp(OF[1]))

        if num_subdivisions_for_error < 1:
            nu_sub.append(np.nan)
            R0_sub.append(np.nan)


        return [nu_best, R0_best, min(nu_sub), max(nu_sub), min(R0_sub), max(R0_sub),  reduced_chi_squared_fitting, reduced_chi_squared_all, np.vstack((np.array(fitting_separation),np.array(fitting_distances))), np.vstack((seq_sep_vals, seq_sep_RMS_distance, sstools.powermodel(seq_sep_vals, nu_best, R0_best)))]



    # ........................................................................
    #
    #
    def get_all_SASA(self, probe_radius=1.4, mode='residue', stride=20):
        """Solvent-accessible surface area (SASA) per residue or per atom.

        Internally uses mdtraj's Shrake-Rupley implementation (Golden-Spiral
        algorithm). Results are memoised per ``(stride, mode, probe_radius)``
        triple so repeated calls with the same parameters are free.

        SASA is returned in Angstroms^2; ``probe_radius`` is in Angstroms.
        Because the calculation is expensive the default ``stride`` is 20.

        Parameters
        ----------
        probe_radius : float, optional
            Solvent probe radius in Angstroms. Default 1.4 (water).
        mode : {'residue', 'atom', 'sidechain', 'backbone', 'all'}, optional
            Granularity:

            * ``'residue'`` (default) - per-residue SASA, shape ``(n_frames_strided, n_residues)``.
            * ``'atom'`` - per-atom SASA, shape ``(n_frames_strided, n_atoms)``.
            * ``'sidechain'`` - per-residue SASA summed over sidechain atoms
              (excluding backbone H, HA, HA2, HA3 which mdtraj's 'sidechain'
              selector erroneously includes).
            * ``'backbone'`` - per-residue SASA summed over backbone atoms
              (including the backbone hydrogens that mdtraj's 'backbone'
              selector erroneously excludes).
            * ``'all'`` - returns a 3-tuple of (residue, sidechain, backbone)
              arrays.
        stride : int, optional
            Use every ``stride``-th frame. Default is 20.

        Returns
        -------
        np.ndarray or tuple of np.ndarray
            See ``mode`` description. All arrays are in Angstroms^2.

        Example
        -------
        >>> sasa = protein.get_all_SASA(mode='residue', stride=10)
        >>> sasa.mean(axis=0)   # mean per-residue SASA across frames
        """

        ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        ## internal function
        def get_sasa_based_on_type(basis, passed_mode):
            # get all reas
            all_areas = basis[0]*100
            all_atoms = list(basis[1])

            # define initial SASA
            residue_SASA = []

            # for each residue with a CA atom...
            for i in self.resid_with_CA:

                # get the atomic indices
                if passed_mode == 'sidechain':
                    
                    # for some reason 'sidechain' selection includes the backbone hydrogen atoms??!?!
                    relevant_atom_idx = self.topology.select(f'resid {i} and {passed_mode} and (not name "H" "HA" "HA2" "HA3")')

                if passed_mode == 'backbone':
                    # for some reason 'backbone' ignores the backbone hydrogen atoms ?!?!?
                    relevant_atom_idx = self.topology.select(f'(resid {i} and {passed_mode}) or (resid {i} and name "H" "HA" "HA2" "HA3")')

                # no atoms so create an empty list
                if len(relevant_atom_idx) == 0:
                    per_res = list(np.zeros(len(all_areas.transpose()[0])))
                else:

                    # get SASA index of first atom index...
                    per_res = all_areas.transpose()[all_atoms.index(relevant_atom_idx[0])]

                    for atom in relevant_atom_idx[1:]:
                        per_res = per_res + all_areas.transpose()[all_atoms.index(atom)]

                residue_SASA.append(per_res)

            return np.array(residue_SASA).transpose()
            ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        ##
        ## START OF ACTUAL FUNCTION
        ##

        # validate input mode
        ssutils.validate_keyword_option(mode, ['residue', 'atom','backbone','sidechain','all'], 'mode')

        # build a specific memoized name for this request
        memoized_name = f'SASA_{stride}_{mode}_{probe_radius}'

        # return already computed SASA is available...
        if memoized_name in self.__SASA_saved:
            return self.__SASA_saved[memoized_name]

        # downsample based on the stride
        target = self.__get_subtrajectory(self.traj, stride)

        # 100* to convert from nm^2 to A^2, and *0.1 for probe radius to convert from A to nm
        if mode == 'residue':
            return_data = 100*md.shrake_rupley(target, mode='residue', probe_radius=probe_radius*0.1)

        if mode == 'atom':
            return_data = 100*md.shrake_rupley(target, mode='atom', probe_radius=probe_radius*0.1)

        if mode == 'sidechain' or mode == 'backbone' or mode == 'all':

            # run calculation for on the whole trajectory
            basis =  md.shrake_rupley(target, mode='atom', probe_radius=probe_radius*0.1, get_mapping=True)

            # extract sidechains
            if mode == 'sidechain' or mode == 'all':
                SC_SASA = get_sasa_based_on_type(basis, 'sidechain')

            # extract backbone
            if mode == 'backbone' or mode == 'all':
                BB_SASA = get_sasa_based_on_type(basis, 'backbone')

            if mode == 'all':
                ALL_SASA = 100*md.shrake_rupley(target, mode='residue', probe_radius=probe_radius*0.1)

            if mode == 'sidechain':
                return_data = SC_SASA

            if mode == 'backbone':
                return_data = BB_SASA

            if mode == 'all':
                return_data =  (ALL_SASA, SC_SASA, BB_SASA)

        # save for the future!!
        self.__SASA_saved[memoized_name] = return_data
        return return_data



    # ........................................................................
    #
    #
    def get_site_accessibility(self, input_list, probe_radius=1.4, mode='residue_type', stride=20):
        """Mean & std SASA for selected residue types or specific residue indices.

        Under ``mode='residue_type'``, ``input_list`` is interpreted as
        canonical three-letter residue names and every matching residue in
        the chain is summarised. Under ``mode='resid'`` the input is treated
        as a list of explicit residue indices.

        Useful for comparing accessibility of chemically-equivalent residues
        at different positions in the chain.

        Parameters
        ----------
        input_list : list
            Either:

            * three-letter residue names, e.g. ``['TRP', 'TYR', 'GLN']``
              (when ``mode='residue_type'``).
            * residue indices, e.g. ``[1, 2, 3, 4]`` (when ``mode='resid'``).
        probe_radius : float, optional
            Solvent probe radius in Angstroms. Default 1.4 (water).
        mode : {'residue_type', 'resid'}, optional
            How to interpret ``input_list``. Default ``'residue_type'``.
        stride : int, optional
            Use every ``stride``-th frame. Default is 20.

        Returns
        -------
        dict
            Keys are ``"RESNAME-RESID"`` strings; each value is
            ``[mean_SASA, std_SASA]`` in Angstroms^2.

        Raises
        ------
        SSException
            If ``input_list`` contains a residue name not in the canonical
            20 + caps set when ``mode='residue_type'``.

        Example
        -------
        >>> protein.get_site_accessibility(['LEU', 'VAL'])
        {'LEU-5': [42.1, 6.0], 'LEU-12': [55.4, 8.3], 'VAL-21': [38.0, 5.5], ...}
        """

        ## First check mode is valid and then sanity check input
        ssutils.validate_keyword_option(mode, ['residue_type', 'resid'], 'mode')


        # empty list of residues we're going to examine (resids)
        resid_list = []

        # if we're in 'residue_type' mode
        if mode == 'residue_type':

            ## First check all the passed residue types are valid...
            # set the list of valid residues to the 20 AAs + the caps
            for res in input_list:
                if res not in ALL_VALID_RESIDUE_NAMES:
                    raise SSException(f"Error: Tried to use get_site_accessibility in 'residue_type' mode but residue {res} not found in list of valid residues {ALL_VALID_RESIDUE_NAMES}")


            # get AA list and convert into an ordered list of each
            # AA residue type
            SEQS = self.get_amino_acid_sequence(numbered=False)

            # initialze the empty idx
            idx = 0

            # this cyles over each amino acid in the sequences and, if that amino acid is one of the ones
            # we're interested in adds the index of that residues to the resid_list
            for res in SEQS:
                if res in input_list:
                    resid_list.append(idx)
                idx = idx + 1


        # if we're in resid mode
        elif mode == 'resid':
            resid_list = input_list

        # next compute ALL SASA for all residues (need the full protein based SASA
        ALL_SASA = np.transpose(self.get_all_SASA(stride=stride, probe_radius=probe_radius))

        lookup = self.get_amino_acid_sequence()

        return_data = {}
        for i in resid_list:
            return_data[lookup[i]] = [np.mean(ALL_SASA[i]), np.std(ALL_SASA[i])]

        return return_data

    # ........................................................................
    #
    #
    def get_regional_SASA(self, R1, R2, probe_radius=1.4, stride=20):
        """Mean total SASA summed over a contiguous residue range.

        Calls :meth:`get_all_SASA` (whose result is memoised) and sums the
        time-averaged per-residue SASA values from ``R1`` to ``R2-1``. Note
        the SASA *must* be computed on the full protein first because atoms
        outside the region can occlude region atoms.

        SASA is returned in Angstroms^2.

        Parameters
        ----------
        R1 : int
            First residue of the region.
        R2 : int
            Last residue of the region (exclusive in the internal sum).
        probe_radius : float, optional
            Solvent probe radius in Angstroms. Default 1.4 (water).
        stride : int, optional
            Use every ``stride``-th frame. Default is 20.

        Returns
        -------
        float
            Sum over residues ``[R1, R2)`` of the time-averaged per-residue
            SASA in Angstroms^2.

        Example
        -------
        >>> protein.get_regional_SASA(5, 25)
        2150.6
        """

        # NOTE - we HAVE to compute SASA over the full ensemble to take into acount
        # atoms OUTSIDE the region getting in the way of the regional SASA
        total = self.get_all_SASA(stride=stride, probe_radius=probe_radius)

        total = np.transpose(total)

        regional_SASA = 0
        for i in range(R1, R2):
            regional_SASA = regional_SASA + np.mean(total[i])

        # return the mean sum of SASA for all atoms
        return regional_SASA


    # ........................................................................
    #
    #
    def get_sidechain_alignment_angle(self, R1, R2, sidechain_atom_1='default', sidechain_atom_2='default'):
        """Per-frame angle between two residues' sidechain orientation vectors.

        For each residue a "sidechain vector" is the CA-to-tip vector, where
        the tip atom is the canonical sidechain terminus (e.g. CB for Ala,
        CZ for Phe/Tyr/Arg). The returned angle is between those two
        vectors. Useful for asking whether two sidechains tend to point
        towards each other (small angle) or away (~180 degrees).

        Default tip atoms (from the OPLS-AA force field)::

            ALA: CB    CYS: SG    ASP: CG    GLU: CD    PHE: CZ
            HIS: NE2   ILE: CD1   LYS: NZ    LEU: CG    MET: CE
            ASN: CG    PRO: CG    GLN: CD    ARG: CZ    SER: OG
            THR: CB    VAL: CB    TRP: CG    TYR: CZ

        GLY, ACE, NME have no sidechain and raise.

        Parameters
        ----------
        R1 : int
            Resid of the first residue. Must not be GLY/ACE/NME.
        R2 : int
            Resid of the second residue. Must not be GLY/ACE/NME.
        sidechain_atom_1 : str, optional
            Override the tip atom of R1. ``'default'`` (default) uses the
            canonical OPLS-AA atom from the table above.
        sidechain_atom_2 : str, optional
            Override the tip atom of R2. Same convention as
            ``sidechain_atom_1``.

        Returns
        -------
        np.ndarray
            1D array of length ``n_frames`` with the per-frame alignment
            angle in degrees.

        Raises
        ------
        SSException
            If either residue is glycine or a cap, if a residue index is out
            of range, or if a custom sidechain atom name cannot be found.

        Example
        -------
        >>> theta = protein.get_sidechain_alignment_angle(5, 45)
        >>> theta.mean()
        82.4
        """


        # The initial section is responsible for selecting/defining the sidechain atom and catching any
        # bad inputs. It's a litte drawn out but useful for code clarity reasons.
        # note need to do this with the un-corrected residue index values, hence why it comes first
        if R1 not in self.__residue_atom_table:
            raise SSException(f'Residue index for R1: {R1:d} does not exist in the protein.')

        if R2 not in self.__residue_atom_table:
            raise SSException(f'Residue index for R2: {R2:d} does not exist in the protein.')


        if sidechain_atom_1 == 'default':
            resname_1 = self.get_amino_acid_sequence(numbered=False)[R1]
            resname_1 = sstools.fix_histadine_name(resname_1)

            try:
                sidechain_atom_1 = DEFAULT_SIDECHAIN_VECTOR_ATOMS[resname_1]

                if sidechain_atom_1 == 'ERROR':
                    raise SSException(f'Residue lacks a valid sidechain ({resname_1})')

            except KeyError:
                raise SSException(f'Cannot parse residue at position {R1} (residue name = {resname_1}) ')
        else:
            raise SSException('Unsupported sidechain atom name: "%s". Please use: "default".' % sidechain_atom_1)

        if sidechain_atom_2 == 'default':
            resname_2 = self.get_amino_acid_sequence(numbered=False)[R2]
            resname_2 = sstools.fix_histadine_name(resname_2)

            try:
                sidechain_atom_2 = DEFAULT_SIDECHAIN_VECTOR_ATOMS[resname_2]

                if sidechain_atom_2 == 'ERROR':
                    raise SSException(f'Residue lacks a valid sidechain ({resname_2})')

            except KeyError:
                raise SSException(f'Cannot parse residue at position {R2} (residue name = {resname_2}) ')
        else:
            raise SSException('Unsupported sidechain atom name: "%s". Please use: "default".' % sidechain_atom_1)

        ### At this point we have reasonable atom names defined!
        TRJ_1_SC = self.traj.atom_slice(self.topology.select(f'resid {R1} and name "{sidechain_atom_1}"'))
        TRJ_1_CA = self.traj.atom_slice(self.topology.select(f'resid {R1} and name "CA"'))

        TRJ_2_SC = self.traj.atom_slice(self.topology.select(f'resid {R2} and name "{sidechain_atom_2}"'))
        TRJ_2_CA = self.traj.atom_slice(self.topology.select(f'resid {R2} and name "CA"'))


        # compute CA-SC vector
        R1_vector = TRJ_1_SC.xyz - TRJ_1_CA.xyz
        R2_vector = TRJ_2_SC.xyz - TRJ_2_CA.xyz

        # finally compute the alignment for each frame and return
        # a vector of alignments
        nframes = R1_vector.shape[0]
        alignment=[]
        for i in range(0, nframes):

            # convert to unit vector
            V1 = R1_vector[i][0] / np.linalg.norm(R1_vector[i][0])
            V2 = R2_vector[i][0] / np.linalg.norm(R2_vector[i][0])

            # compute dot product and take arccos and do rad->deg to get angle
            # between the unit vectors
            alignment.append(np.rad2deg(np.arccos(np.dot(V1,V2))))

        return np.array(alignment)


    # ........................................................................
    #
    #
    def get_dihedral_mutual_information(self, angle_name='psi', bwidth=np.pi/5.0, stride=1, weights=False, normalize=False):
        """Mutual-information matrix between every pair of a chosen dihedral.

        For each pair of dihedrals :math:`(\\phi_i, \\phi_j)` of the same
        type, this computes the mutual information

        .. math::
            I(\\phi_i; \\phi_j) = H(\\phi_i) + H(\\phi_j) - H(\\phi_i, \\phi_j)

        where the Shannon entropies :math:`H` come from binned histograms.
        The natural baseline for interpretation is the equivalent matrix
        computed on a polymer-model reference ensemble (e.g. an
        excluded-volume simulation), since the absolute values are
        bin-width-dependent. Setting ``normalize=True`` divides each entry
        by the joint entropy, which can help with comparisons across
        ensembles of different sizes.

        Parameters
        ----------
        angle_name : {'phi', 'psi', 'omega', 'chi1'}, optional
            Dihedral to use. Default ``'psi'``. (``chi2``-``chi5`` are not
            supported by this method.)
        bwidth : float, optional
            Histogram bin width on the ``[-pi, pi]`` interval. Smaller
            widths give a finer-grained estimate but need more data.
            Default ``pi/5``.
        stride : int, optional
            Use every ``stride``-th frame. Default 1.
        weights : array_like, optional
            Per-frame weights for re-weighted MI estimation. Must satisfy
            ``len(weights) == n_frames // stride``. Default ``False``.
        normalize : bool, optional
            If True, divide each MI value by the joint entropy. Default False.

        Returns
        -------
        np.ndarray
            Symmetric ``(n_dihedrals, n_dihedrals)`` mutual-information
            matrix. The diagonal is the self-information of each dihedral.

        Raises
        ------
        SSException
            If ``bwidth`` is not in ``(0, 2*pi]`` or ``angle_name`` is not
            one of the four allowed strings.

        Example
        -------
        >>> MI = protein.get_dihedral_mutual_information(angle_name='psi')
        >>> MI.shape
        (92, 92)
        """


        ## ..................................................
        ## SAFETY FIRST!
        ##
        # verify binwidth input values
        if bwidth > 2*np.pi or not (bwidth > 0):
           raise SSException(f'The bwidth parameter must be between 2*pi and greater than 0. Received {bwidth}')

        # if stride was passed make sure it's ok
        self.__check_stride(stride)

        # if weights were passed make sure they're LEGIT!
        weights = self.__check_weights(weights, stride)

        # check input keyword selector
        ssutils.validate_keyword_option(angle_name, ['chi1', 'phi', 'psi', 'omega'], 'angle_name')

        ## ..................................................

        # define histogram bins based on passed bin width
        bins = np.arange(-np.pi, np.pi+bwidth, bwidth)

        # construct the selector dictionary
        selector = {"phi":md.compute_phi, "omega":md.compute_omega, "psi":md.compute_psi, "chi1":md.compute_chi1}

        # check the angle_name is an allowed name
        if angle_name not in list(selector.keys()):
            raise SSException(f'The variable angle_name was set to {angle_name}, which is not one of phi, omega, psi, chi1')

        # select and compute the relevant angles of the subtrajectory
        fx = selector[angle_name]
        angles = fx(self.traj[0::stride])

        # construct empty matrices
        SIZE = len(angles[0])
        MI_mat = np.zeros((SIZE,SIZE))

        # populate matrices note we only compute the upper right triangle
        # but can populate the lower half because it's a symmetrical matrix
        for i in range(0,SIZE):
            for j in range(i,SIZE):

                X = np.transpose(angles[1])[j]
                Y = np.transpose(angles[1])[i]
                MI = ssmutualinformation.calc_MI(X,Y, bins, weights,normalize)

                MI_mat[i,j] = MI
                MI_mat[j,i] = MI

        return MI_mat


    # ........................................................................
    #
    #
    def get_local_to_global_correlation(self, mode='COM', n_cycles=100, max_num_pairs=10, stride=20, weights=False, verbose=True):
        """How well finite subsets of inter-residue distances predict global Rg.

        For each ``k`` in ``[1, max_num_pairs)``, this randomly picks ``k``
        inter-residue distance pairs ``n_cycles`` times and computes the
        correlation between the chain's :math:`R_g^2` and the mean-squared
        of those ``k`` distances. Averaging the correlations across cycles
        gives a per-``k`` measure of how informative a finite local
        measurement is about the global chain dimensions.

        This is intended to inform experimental design (e.g. how many
        smFRET pairs are "enough" to pin down :math:`R_g`).

        Parameters
        ----------
        mode : {'COM', 'CA'}, optional
            Distance type. Default ``'COM'``.
        n_cycles : int, optional
            Number of random pair-resamples per ``k``. Default 100. Larger
            ``n_cycles`` reduces sampling noise; values below ~50 are not
            recommended.
        max_num_pairs : int, optional
            Upper bound on the number of pairs sampled. Default 10. Beyond
            10-20 the average correlation typically saturates close to 1.
        stride : int, optional
            Use every ``stride``-th frame. Default 20.
        weights : array_like or False, optional
            Per-frame weights for the covariance computation. Default
            ``False`` (uniform weighting).
        verbose : bool, optional
            If True (default), print one status line per ``k``.

        Returns
        -------
        tuple of length 4
            ``(raw, n_pairs, mean_corr, std_corr)`` where:

            * ``raw`` (np.ndarray of shape ``(n_cycles * (max_num_pairs - 1), 2)``):
              one row per resample; col 0 is ``k`` and col 1 is the
              :math:`R_g^2` vs. local-mean-square correlation for that
              resample.
            * ``n_pairs`` (np.ndarray): the ``k`` values used,
              ``[1, 2, ..., max_num_pairs - 1]``.
            * ``mean_corr`` (np.ndarray): mean correlation per ``k``.
            * ``std_corr`` (np.ndarray): standard deviation of correlation
              per ``k`` (assumes Gaussian distribution, may be approximate).

        Raises
        ------
        SSException
            If ``mode`` is invalid, the strided Rg length disagrees with the
            inter-residue distance length, or the installed numpy version
            cannot accept ``aweights`` in ``numpy.cov``.

        Example
        -------
        >>> import numpy as np
        >>> np.random.seed(42)   # the pair-resampling consumes the RNG
        >>> raw, ks, mean_c, std_c = protein.get_local_to_global_correlation()
        """

        weights = self.__check_weights(weights, stride)


        ssutils.validate_keyword_option(mode, ['CA', 'COM'], 'mode')

        # define first and last residue that has a CA
        start = self.resid_with_CA[0]
        end = self.resid_with_CA[-1]

        FIRST_CHECK = True

        # Hoist the per-residue CA position out of the O(n^2) pair loop.
        # get_inter_residue_atomic_distance recomputes each residue's CA
        # center-of-mass from scratch on every call; caching it per
        # residue is byte-identical but O(n) rather than O(n^2). The
        # hoist consumes no RNG, so the downstream pair-resampling sees
        # an identical stream and its output is unchanged.
        if mode == 'CA':
            _ca = self.__ca_position_cache(stride)

        # start with a sequence separation of 1
        for seq_sep in range(1, self.n_residues):
            ssio.status_message(f"On sequence separation {seq_sep}", verbose)

            for pos in range(start, end - seq_sep):

                # define the two positions
                A = pos
                B = pos + seq_sep

                # get the distance for every stride-th frame between those two positions using either the CA
                # mode or the COM mode
                if mode == 'CA':
                    # byte-identical to get_inter_residue_atomic_distance(A, B, stride=stride)
                    distance = np.linalg.norm(_ca(A) - _ca(B), axis=1)

                elif mode == 'COM':
                    distance = self.get_inter_residue_COM_distance(A, B, stride=stride)

                if FIRST_CHECK:
                    all_distances = distance
                    FIRST_CHECK = False
                else:
                    all_distances = np.vstack((all_distances, distance))


        full_rg = self.get_radius_of_gyration()
        stride_rg = full_rg[0::stride]

        if len(stride_rg) != len(all_distances[0]):
            raise SSException('Something when wrong when comparing stride-derived Rg and internal distances, this is a bug in the code...)')

        # total number of distance pairs
        n_pairs = len(all_distances)

        pair_selection_vector = np.arange(1,max_num_pairs,1)

        return_data = np.zeros((len(pair_selection_vector)*n_cycles,2))

        weights = False
        weights = np.repeat(1.0/len(stride_rg),len(stride_rg))

        idx=0
        for n_selected in pair_selection_vector:
            ssio.status_message(f"On {n_selected} pairs selected", verbose)
            for i in range(0,n_cycles):

                # select n_select different values between 0 and n_pairs
                idx_selection = np.random.randint(0, n_pairs, n_selected)

                # this next line means that for each idx_selection value (i.e. for each pair of distances)
                # we calculate the squared value of each distance, then for each conformation SUM all the pairs
                # and then for each of those n_frame summs divide by 2 * n_selected squared
                # this gives a series of values for EACH conformation (so local mean is a 1xn_configurations vector)
                local_mean_square = np.sum(np.power(all_distances[idx_selection],2),0)/(2*(n_selected*n_selected))

                # correlation between local mean square and rg_square - returns a SINGLE value

                # THIS is the covariance matrix, but we can pass a weighting factor to this making it better suited for
                # weighted ensembles

                if weights is not False:
                    try:
                        cov_matrix = np.cov(np.vstack((local_mean_square, np.power(stride_rg,2))), ddof=0, aweights=weights)
                    except TypeError:
                        # this probaby happend because the version of numpy doesn't support aweights. If this is the case try again
                        # without weights assuming weights are equal
                        if len(set(weights)) == 1:
                            cov_matrix = np.cov(np.vstack((local_mean_square, np.power(stride_rg,2))))
                        else:
                            raise SSException(f'Weights being passed to get_global_from_local but the current version of numpy ({np.version.full_version}) may not support weights via the "aweights" keyword...')

                    # this computes the upper right square from the correlation matrix from the covariance matrix using (Rij = (Cij/(sqrt(Cii - Cjj))) where i and j are
                    # 0 and 1 respectively
                    c = cov_matrix[0,1]/(np.sqrt(cov_matrix[0,0]*cov_matrix[1,1]))
                else:
                    # this is the correlation matrix and was the old way of doing this
                    c = np.corrcoef(local_mean_square, np.power(stride_rg,2))[0][1]

                return_data[idx] = [n_selected, c]
                idx=idx+1

        # leaving this here incase we want to re-introduce the 2D histogram information in later versions...
        # np.histogram2d(np.transpose(return_data)[0],np.transpose(return_data)[1], bins=[np.arange(1,n_pairs), np.arange(0,1,0.01)]))

        return  (return_data,
                 pair_selection_vector,
                 np.mean(np.reshape(np.transpose(return_data)[1], (len(pair_selection_vector),n_cycles)),1),
                 np.std(np.reshape(np.transpose(return_data)[1], (len(pair_selection_vector),n_cycles)),1))



    # ........................................................................
    #
    #
    def get_end_to_end_vs_rg_correlation(self, mode='COM'):
        """Pearson correlation between per-frame :math:`R_e^2` and :math:`R_g^2`.

        For a Gaussian chain Debye's result gives :math:`R_e^2 = 6 R_g^2`,
        but in practice the correlation is a more useful diagnostic of how
        fractal-like (or otherwise) the chain is. This method returns the
        Pearson correlation across frames between the end-to-end distance
        squared and the radius of gyration squared.

        Parameters
        ----------
        mode : {'COM', 'CA'}, optional
            Type of distance to use for the end-to-end distance.
            Default ``'COM'``.

        Returns
        -------
        float
            Pearson correlation coefficient in ``[-1, 1]``.

        Example
        -------
        >>> protein.get_end_to_end_vs_rg_correlation()
        0.92
        """

        # validate the keyword
        ssutils.validate_keyword_option(mode, ['CA', 'COM'], 'mode')

        # get end-to-end distance
        distance = self.get_end_to_end_distance(mode)

        # get radius of gyration
        full_rg = self.get_radius_of_gyration()

        # NOTE this is the same as Re^2 =  RG^2/ 6 (standard Debye result for a Gaussian chain) - basically asking
        # about fractality more than a specific chain model, because examining correlation the scalar factor doesn't
        # matter , ie really this is Rg^2 vs Re^2 - scalar is irrelevant. But, this approach is consistent with
        # the approach in the get_local_to_global_correlation() function
        local_mean_square = np.power(distance,2)/(2)
        c = np.corrcoef(local_mean_square, np.power(full_rg,2))[0][1]

        return c

    # ........................................................................
    #
    #
    def get_secondary_structure_DSSP(self, R1=None, R2=None, return_per_frame=False):
        """Per-residue secondary structure (helix / extended / coil) via DSSP.

        Calls ``mdtraj.compute_dssp`` on the chosen residue range and
        summarises the output into per-residue fractional occupancies (the
        default) or per-frame binary masks (with ``return_per_frame=True``).

        The three buckets are H (helix), E (extended/beta), C (coil). At
        every residue the three values sum to 1 (default mode).

        Parameters
        ----------
        R1 : int, optional
            First residue of the region. ``None`` (default) means the first
            CA-bearing residue.
        R2 : int, optional
            Last residue of the region. ``None`` (default) means the last
            CA-bearing residue.
        return_per_frame : bool, optional
            * False (default): return per-residue fractional occupancies.
            * True: return per-frame binary 1/0 masks.

        Returns
        -------
        tuple of length 4
            ``(resid_list, helix, extended, coil)`` where:

            * ``resid_list`` (list of int): residue indices covered.
            * If ``return_per_frame=False``: ``helix``, ``extended``, ``coil``
              are 1D arrays of length ``len(resid_list)`` with per-residue
              fractional time in each bucket.
            * If ``return_per_frame=True``: ``helix``, ``extended``, ``coil``
              are 2D arrays of shape ``(n_frames, len(resid_list))`` with
              1 / 0 entries indicating presence of that bucket in that frame.

        Example
        -------
        >>> resids, h, e, c = protein.get_secondary_structure_DSSP()
        >>> h.mean()    # mean helicity across the chain
        0.32
        """

        # build R1/R2 values
        out = self.__get_first_and_last(R1, R2, withCA=True)
        R1_real = out[0]
        R2_real = out[1]

        # select the relevant subtrajectory (out[2] is the 'resid %i to %i' where %i and %i are R1 and R2)
        ats = self.traj.atom_slice(self.topology.select(f'{out[2]}'))

        # compute DSSP over the selected subtrajectory
        dssp_data = md.compute_dssp(ats)
        reslist    = list(range(R1_real, R2_real+1))

        C_vector = []
        E_vector = []
        H_vector = []

        # if we want per-frame
        if return_per_frame is True:
            for f in dssp_data:
                C_vector.append(np.array(f=='C',dtype=int).tolist())
                E_vector.append(np.array(f=='E',dtype=int).tolist())
                H_vector.append(np.array(f=='H',dtype=int).tolist())

            return (reslist, np.array(H_vector), np.array(E_vector), np.array(C_vector))

        # note the + 1 because the R1 and R2 positions are INCLUSIVE whereas
        else:
            n_frames   = self.n_frames

            for i in range(len(reslist)):
                C_vector.append(float(sum(dssp_data.transpose()[i] == 'C'))/n_frames)
                E_vector.append(float(sum(dssp_data.transpose()[i] == 'E'))/n_frames)
                H_vector.append(float(sum(dssp_data.transpose()[i] == 'H'))/n_frames)


            return (reslist, np.array(H_vector), np.array(E_vector), np.array(C_vector))


    # ........................................................................
    #
    #
    def get_secondary_structure_BBSEG(self, R1=None, R2=None, return_per_frame=False):
        """Per-residue secondary structure via the BBSEG2 phi/psi classification.

        Classifies each (phi, psi) pair into one of 9 BBSEG2 bins using a
        lookup table distributed with CAMPARI. Implementation is in pure
        Python, so the result is version-stable across mdtraj releases.

        BBSEG2 classes:

        * 0 - unclassified
        * 1 - beta (turn/sheet)
        * 2 - PII (polyproline type II helix)
        * 3 - unusual region
        * 4 - right-handed alpha helix
        * 5 - inverse C7-equatorial (gamma-prime turn)
        * 6 - classic C7-equatorial (gamma turn)
        * 7 - 7-residue-per-turn helix
        * 8 - left-handed alpha helix

        Parameters
        ----------
        R1 : int, optional
            First residue of the region. ``None`` (default) means the first
            residue with phi/psi defined.
        R2 : int, optional
            Last residue of the region. ``None`` (default) means the last
            residue with phi/psi defined.
        return_per_frame : bool, optional
            * False (default): return per-residue fractional occupancies.
            * True: return per-frame binary masks.

        Returns
        -------
        tuple of (list, dict)
            ``(resid_list, per_class)`` where:

            * ``resid_list`` (list of int): residue indices covered.
            * ``per_class`` (dict): keys 0..8 (int), values are either
              1D fractional-occupancy arrays (length ``len(resid_list)``)
              or 2D ``(n_frames, len(resid_list))`` binary masks depending
              on ``return_per_frame``.

        Example
        -------
        >>> resids, classes = protein.get_secondary_structure_BBSEG()
        >>> classes[4].mean()    # mean right-handed-alpha occupancy across chain

        """

        # build R1/R2 values - NOTE that for BBSEG because we compute from PHI/PSI angles
        # we don't need to specify withCA becayse phi/psi are oNLY between valid peptide
        # units, so this means we automatically select the right sets of residues
        out = self.__get_first_and_last(R1, R2, withCA=False)

        # note we select RESID using withCA = True
        out_selector = self.__get_first_and_last(R1, R2, withCA=True)
        R1_real = out_selector[0]
        R2_real = out_selector[1]
        reslist    = list(range(R1_real, R2_real+1))


        # extract the phi/psi angles in degrees
        phi_data = np.degrees(md.compute_phi(self.traj.atom_slice(self.topology.select(f'{out[2]}')))[1])
        psi_data = np.degrees(md.compute_psi(self.traj.atom_slice(self.topology.select(f'{out[2]}')))[1])

        # extract the relevant information (note shape of phi_data and psi_data will be identical)
        # shape info here is (number_of_frames, number_of_residues) sized
        shape_info = np.shape(phi_data)
        all_classes = []

        # for each frame iterate through and classify each residue, building a shape_info
        # sized matrix where each elements reflects the BBSEG classification of that residue
        # in a given frame
        for f in range(0, shape_info[0]):

            # so each step through the loop we're passing two vectors, each of which
            # is nres residues long
            all_classes.append(self.__phi_psi_bbseg(phi_data[f], psi_data[f]))

        # convert to a numpy array
        all_classes = np.array(all_classes)
        return_bbseg = {}

        # if we want per-frame BBSEG secondary structure info
        if return_per_frame:

            # initialize the empty lists to be appended to
            for c in range(0, 9):
                return_bbseg[c] = []

            # cycle over each frame
            for f in all_classes:

                # cycle over each class in each frame to binarize the per-frame/per residue vector
                for c in range(0,9):
                    return_bbseg[c].append(np.array(f==c, dtype=int).tolist())

            for c in return_bbseg:
                return_bbseg[c] = np.array(return_bbseg[c])
            return (reslist, return_bbseg)
        else:

            # finally cycle through each BBSEG classification type and average
            # over each frame (shape_info[0] is number of frames)

            for c in range(0,9):
                return_bbseg[c] = np.sum((all_classes == c)*1,0)/shape_info[0]

            return (reslist, return_bbseg)


    # ........................................................................
    #
    def __phi_psi_bbseg(self, phi_vector, psi_vector):
        """Internal function that takes two equally-matched phi and psi angle
        vectors and based on the pairwise combination classified each pair of
        elements using the BBSEG2 definition. Definition was generated from the
        BBSEG2 file distributed with CAMPARI, and is encoded and stored in the
        _internal_data module.

        NOTE that because this is an internal function we do not double check
        that the phi_vector and psi_vectors are of the same length, but this
        is critical, so if this function is being called make sure this is
        true!

        Parameters
        ----------
        phi_vector :   iterable (list or numpy vector)
             ordered list of phi angles for a specific residue

        psi_vector :   iterable (list or numpy vector)
            ordered list of psu angles for a specific residue

        Returns
        -------

        classes : list
             A list of length equal to phi_vector and psi_vector that
             classifies each pair of phi/psi angles using the BBSEG2
             definition.
        """

        classes = []

        for i in range(len(phi_vector)):
            phi = phi_vector[i]
            psi = psi_vector[i]

            fixed_phi = phi-(phi%10)
            fixed_psi = psi-(psi%10)

            # following corrections for edge cases if we hit
            if fixed_phi == 180.0:
                fixed_phi = 170.0

            if fixed_psi == 180.0:
                fixed_psi = 170.0

            # classify the phi/psi values in terms of BBSEG values
            classes.append(BBSEG2[fixed_phi][fixed_psi])

        return classes


    # ........................................................................
    #
    def get_overlap_concentration(self):
        """Overlap concentration :math:`c^*` of the chain, in Molar.

        Defined as the polymer concentration at which independent chains
        begin to occupy the same volume in trans, computed by treating each
        chain as a sphere of radius equal to its mean :math:`R_g`. Useful as
        a coarse limit on dilute-regime behaviour.

        Returns
        -------
        float
            Overlap concentration :math:`c^*` in Molar.

        Example
        -------
        >>> protein.get_overlap_concentration()
        0.0042
        """

        return sspolymer.get_overlap_concentration(np.mean(self.get_radius_of_gyration()))



    # ........................................................................
    #
    #
    def get_angle_decay(self, atom1='C', atom2='N', return_all_pairs=False):

        """Bond-vector autocorrelation along the chain (persistence-length proxy).

        For each residue we form a single intra-residue vector
        :math:`\\vec{v}_i` (default ``C -> N``, the peptide bond direction).
        For every pair :math:`(i, j)` we compute
        :math:`\\langle \\hat{v}_i \\cdot \\hat{v}_j \\rangle` across
        frames; averaging over all pairs with the same sequence separation
        :math:`|i-j|` yields the decay profile. A faster initial decay
        implies a less rigid chain.

        Both ``atom1`` and ``atom2`` must exist in every residue, so the
        peptide C->N vector is essentially the only reasonable choice.

        Parameters
        ----------
        atom1 : str, optional
            Origin atom name. Default ``'C'``.
        atom2 : str, optional
            Tip atom name. Default ``'N'``.
        return_all_pairs : bool, optional
            If True, additionally return a dict of every individual
            ``"i-j"`` residue-pair correlation. Default is False.

        Returns
        -------
        np.ndarray OR tuple of (np.ndarray, dict)
            * If ``return_all_pairs=False`` (default): array of shape
              ``(n_separations, 3)`` with columns ``[separation, mean_corr,
              std_corr]``. The first row is ``[0, 1.0, 0.0]`` (self
              correlation).
            * If ``return_all_pairs=True``: ``(array, pair_dict)`` where
              ``pair_dict["i-j"]`` is the mean correlation for that
              specific residue pair.

        Example
        -------
        >>> decay = protein.get_angle_decay()
        >>> # decay[:, 0] is separation, decay[:, 1] is mean correlation
        """

        # first compute all the C-N vector for each residue

        CN_vectors = []
        CN_lengths = []

        for i in self.resid_with_CA:

            # this extracts the C->N vector for each frame for each residue
            value = np.squeeze(self.traj.atom_slice(self.__residue_atom_lookup(i, atom1)).xyz) - np.squeeze(self.traj.atom_slice(self.__residue_atom_lookup(i, atom2)).xyz)

            # CN_vectors becomes a list where each element is [3 x nframes] array where 3 is the x/y/z
            # vector coordinates. 
            CN_vectors.append(value)

            # CN_lengths extracts the ||v|| length of each vector (should be basically the same). We need
            # this for the final linalg operation later
            CN_lengths.append(np.linalg.norm(value,axis=1))

        # calculate the number of residues for which we have C->N vectors
        npos = len(CN_vectors)

        # initialize an empty dictionary - the keys for this are i-j sequence separation,
        # and values are the cos(theta) angle between pairs of vectors
        all_vals={}
        for i in range(1, npos):
            all_vals[i] = []

        # precompute || u || * || v || wich 
        length_multiplier = {}
        for i1 in range(0, npos-1):
            length_multiplier[i1] = {}
            for j1 in range(i1+1, npos):
                length_multiplier[i1][j1] = CN_lengths[i1]*CN_lengths[j1]

        # we then cycle over the non-redudant set of pairwise residues in the protein
        for i1 in range(0, npos-1):
            for j1 in range(i1+1, npos):

                # for each frame calculate (u . v) / (||u|| * ||v||)
                # where u and v are vectors and "." is the dot product between each pair. We're only
                # calculating PAIR-WISE dot product of each [x,y,z] with [x,y,z] vector, so doing
                # np.sum(Matrix*matrix) is SO SO SO much faster than anything else. We also take the
                # average to avoid storing a ton of numbers and generating these giant vectors
                # 

                all_vals[j1-i1].append(np.mean(np.sum(CN_vectors[i1]*CN_vectors[j1],axis=1)/length_multiplier[i1][j1]))

        return_matrix = []
        return_matrix.append([0, 1.0, 0.0])
        for k in all_vals:
            return_matrix.append([k, np.mean(all_vals[k]), np.std(all_vals[k])])

        # convert to a matrix at the end
        return_matrix = np.array(return_matrix)

        
        # if we each possible inter-residue autocorrelation, this generates a dictionary
        # where each unique res1-res2 pair value is represented
        if return_all_pairs:

            # minus 1 because return matrix includes the obligate self correlation which
            # is always 1
            all_pairs_dict = {}

            # for 0 to the number of residues (i.e. each row in the [nres x nres] matrix
            # note that i will be 1, 2, ... n-1 where n is number of residues in the chain

            # build the self correlation pairs first
            for i in all_vals:
                all_pairs_dict[f"{i}-{i}"] = 1.0
            all_pairs_dict[f"{i+1}-{i+1}"] = 1.0

            # i here is the |i-j| distance
            for i in all_vals:
                idx = 1
                for x in all_vals[i]:
                    all_pairs_dict[f"{idx}-{idx+i}"] = x
                    idx = idx + 1


            return (return_matrix, all_pairs_dict)
        else:
            return np.array(return_matrix)



    #oxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxo
    #
    #
    def get_local_collapse(self, window_size=10, bins=None, verbose=True):
        """Sliding-window local radius of gyration profile.

        At every starting residue ``i`` the local Rg is computed over the
        window ``[i, i + window_size]``, producing a per-position profile of
        compaction along the chain. Useful for identifying locally collapsed
        regions in otherwise extended IDPs (and vice versa).

        Parameters
        ----------
        window_size : int, optional
            Size of the sliding window in residues. Must be ``<= n_residues``.
            Default 10.
        bins : np.ndarray or list, optional
            Histogram bin edges (evenly spaced). If None (default),
            ``np.arange(0, 10, 0.01)`` is used.
        verbose : bool, optional
            If True (default), print one status line per starting residue.

        Returns
        -------
        tuple of (list, list, list, np.ndarray)
            ``(mean, std, histo, bins)`` where:

            * ``mean`` (list of float, length ``n_residues - window_size + 1``):
              mean local Rg per starting residue.
            * ``std`` (same length): standard deviation of local Rg per
              starting residue.
            * ``histo`` (list of np.ndarray, same length): histogram counts.
            * ``bins`` (np.ndarray): the bin edges actually used.

        Raises
        ------
        SSException
            If ``window_size > n_residues`` or ``bins`` is malformed.

        Example
        -------
        >>> mean, std, hist, bins = protein.get_local_collapse(window_size=8)
        >>> # plot mean vs residue index to see where the chain collapses

        Returns (legacy enumeration)
        ----------------------------
            * [0] - list of floats (length = n)
                  Each float reports on the mean Rg at a specific position
                  along the sequence.

            * [1] - list of floats (length = n)
                  Each float reports on the standard deviation of the
                  Rg distribution at a specific position along the sequence.

            * [2] - list of np.ndarrays (length = n)
                  Histogram counts  associated with the local Rg at a given
                  position along the sequence. Basically, this is the emprical
                  distribution that the standard devaitions and mean report on

            * [3] - np.ndarray (length = n)
                  Bin values for histogram counts returned in [2]
        """
        # validate bins
        if bins is None:
            bins = np.arange(0,10,0.01)
        else:
            try:
                if len(bins)  < 2:
                    raise SSException('Bins should be a numpy defined vector of values - e.g. np.arange(0,1,0.01)')
            except TypeError:
                raise SSException('Bins should be a list, vector, or numpy array of evenly spaced values')

            try:
                bins = np.array(bins, dtype=float)
            except ValueError:
                raise SSException('Passed bins could not be converted to a numpy array of floats')

            # Check whether the bins are evenly spaced. If the bins are evenly spaced, subtracting the leading
            # value of the discrete difference from itself should yield a sum of 0 (within the floating point epsilon).
            diff = np.diff(bins)
            bins_delta = diff - diff[0]
            fpe = np.finfo(diff[0].dtype).eps
            evenly_spaced = np.isclose(np.sum(bins_delta), 0, rtol=fpe)
            if not evenly_spaced:
                raise SSException(f'The spacing between bins is uneven, or you may using bins widths less than: {fpe:f}')

        n_residues = self.n_residues
        n_frames   = self.n_frames

        # check the window is an appropriate size
        if window_size > n_residues:
            raise SSException('window_size is larger than the number of residues')

        meanData = []
        stdData  = []
        histo    = []

        for i in range(window_size - 1, n_residues):

            ssio.status_message(f"On range {i}", verbose)

            # get radius of gyration (now by default is in Angstroms
            # - in previous versions we performed a conversion here)
            tmp = self.get_radius_of_gyration(i - (window_size-1), i)


            (b, c) = np.histogram(tmp, bins)
            histo.append(b)

            meanData.append(np.mean(tmp))
            stdData.append(np.std(tmp))


        return (meanData, stdData, histo, bins)











