##     _____  ____  _    _ _____   _____  ____  _____  
##   / ____|/ __ \| |  | |  __ \ / ____|/ __ \|  __ \ 
##  | (___ | |  | | |  | | |__) | (___ | |  | | |__) |
##   \___ \| |  | | |  | |  _  / \___ \| |  | |  ___/ 
##   ____) | |__| | |__| | | \ \ ____) | |__| | |     
##  |_____/ \____/ \____/|_|  \_\_____/ \____/|_|     

## Alex Holehouse (Pappu Lab and Holehouse Lab) and Jared Lalmansingh (Pappu lab)
## Simulation analysis package
## Copyright 2014 - 2024
##

import time
import mdtraj as md
import numpy as np
from .configs import *
from .ssdata  import ALL_VALID_RESIDUE_NAMES

from .ssprotein import SSProtein
from .ssexceptions import SSException
from . import ssutils
from . import ssio
from . import sstools

## Order of standard args:
## 1. stride
## 2. weights
## 3. verbose

from copy import copy
from functools import partial
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
    
    def wrapper(self, *args, **kwargs):
        if self._SSTrajectory__single_protein_traj is None:
            self._SSTrajectory__single_protein_traj = self._SSTrajectory__get_all_proteins(self.traj, self._SSTrajectory__explicit_residue_checking)
        return func(self, *args, **kwargs)
    return wrapper
    


class SSTrajectory:

    #oxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxo
    #
    #
    def __init__(self, trajectory_filename=None,
                 pdb_filename=None,
                 TRJ=None,
                 protein_grouping=None,
                 pdblead=False,
                 debug=False,
                 extra_valid_residue_names=None,
                 explicit_residue_checking=False,
                 print_warnings=False):
                 
        
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
                print('Unable to use the extra_valid_residue_names - this must be a list of strings')
        
        # first we decide if we're reading from file or from an existing trajectory
        if (trajectory_filename is None) and (pdb_filename is None):
            if TRJ is None:
                raise SSException('No input provided! Please provide ether a pdb and trajectory file OR a pre-formed traj object')
                
            # note the [:] means this is a COPY!
            self.traj = TRJ[:]
        else:
            if (trajectory_filename is None):
                raise SSException('No trajectory file provided!')
            if (pdb_filename is None):
                raise SSException('No PDB file provided!')

            # read in the raw trajectory
            self.traj = self.__readTrajectory(trajectory_filename, pdb_filename, pdblead, print_warnings)
        
        # Next, having read in the trajectory we parse out into proteins
        # extract a list of protein trajectories where each protein is assumed
        # to be in its own chain
        if protein_grouping == None:
            (self.proteinTrajectoryList, single_chain_sim) = self.__get_proteins(self.traj, debug, explicit_residue_checking=explicit_residue_checking)        
        else:
            self.proteinTrajectoryList = self.__get_proteins_by_residue(self.traj, protein_grouping, debug)

        # this is initialized to None and then gets defined using the lazy_loading_single_protein_trajectory() 
        # decorator when it's first needed
        self.__single_protein_traj = None
        

    def  __repr__(self):
        return "SSTrajectory (%s): %i proteins and %i frames" % (hex(id(self)), self.n_proteins, self.n_frames)


    def __len__(self):
        # Edited to mimic the behavior of `mdtraj` trajectory objects.
        # Originally: `return (self.n_proteins, self.n_frames)`
        return self.n_frames

    @property
    def n_frames(self):
        """
        :property: Returns the number of frames in the trajectory.

        """
        return len(self.traj)

    @property
    def n_proteins(self):
        """
        :property: Returns the number of individual proteins found.
        """
        return len(self.proteinTrajectoryList)

    @property
    def length(self):
        """
        :property: Returns a tuple with number of proteins and number of frames.
        """
        # Implemented a method that encapsulates the data output by the original `__len__` method.
        return (self.n_proteins, self.n_frames)

    @property
    def unitcell(self):
        """
        :property: Returns a list with unit cell information in Angstroms (assuming base units
        of trajectory are in nanometers)
        """
        return self.traj.unitcell_lengths[0]*10
        
        

    #oxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxo
    #
    #
    def __readTrajectory(self, trajectory_filename, pdb_filename, pdblead, print_warnings):
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
        traj =  md.load(trajectory_filename, top=pdb_filename)
                    
        # check unit cell lengths
        try:
            uc_lengths = traj.unitcell_lengths[0]

            # this is s custom warning for a specific edge-case we encounter a lot
            if (uc_lengths[0] == 0 or uc_lengths[1] == 0 or uc_lengths[2] == 0):
                if print_warnings:
                    ssio.warning_message("Trajectory file unit cell lengths are zero for at least one dimension. This is a probably a bug with an FRC generated __START.pdb file, because in the old version of CAMPARI used to do grid based FRC calculations the unit cell dimensions are not written correctly. This may cause issues but we're going to assume everything is OK for now. Check the validity of any analysis output. If you're worried, you can use the following workaround.\n\n:::: WORK AROUNDS ::::\nSimply run\n\ntrjconv -f __traj.xtc -s __START.pdb -box a b c -o frc.xtc \n\nAn then \n\ntrjconv -f frc.xtc -s __START.pdb -box a b c -o start.pdb -dump 0\n\n\nHere\n-f n__traj.xtc   : defines the trajectory file\n-s __start.pdb   : defines the pdb file used to parse the topology\n-box a b c       : defines the box lengths **in nanometers**\n-o frc.xtc       : is the name of the new trajectory file with updated box lengths\nSelect 0 (system) when asked to 'Select group for output'.The second step creates the equivalent PDB file with the header-line correctly defining the box unit cell lengths and angles. These two new files should then be used for analysis.\n\nAs an example, if my FRC simulation had a sphere radius of 100 angstroms then my correction command would look something like \n\ntrjconv -f __traj.xtc -s __START.pdb -box 20 20 20 -o frc.xtc\ntrjconv -f frc.xtc -s __START.pdb -box 20 20 20 -o start.pdb -dump 0")

        except TypeError:
            if print_warnings:
                ssio.warning_message("Warning: UnitCell lengths were not provided... This may cause issues but we're going to assume everything is OK for now...")
        
        # if pdbLead is true then load the pdb_filename as a trajectory
        # and then add it to the front (the PDB file is its own topology
        # file so no need to specificy the top= file here!
        if pdblead:
            pdbtraj = md.load(pdb_filename)
            traj = pdbtraj+traj


            # having added the PDB file now check all the unit-cells match up!
            try:
                uc_lengths_1 = traj.unitcell_lengths[0]
                uc_lengths_2 = traj.unitcell_lengths[1]
            
                if (uc_lengths_1[0] != uc_lengths_2[0]) or (uc_lengths_1[2] != uc_lengths_2[2]) or (uc_lengths_1[2] != uc_lengths_2[2]):
                    ssio.warning_message('........................\nWARNING:\nThe unit cell dimensions of the PDB file and trajectory file did not match, specifically\nPDB file = [ %s ]\nXTC file = [ %s ]\nThis may cause issues if native MDTraj untilities are used (and potentially for SOURSOP utilities that are based on these. It is not necessarily an issue, but PLEASE sanity check your outcome. To be save we reeommend editing the PDB-file untilcell dimenions to match.')

            except TypeError:
                ssio.warning_message("Warning: UnitCell lengths were not provided... This may cause issues but we're going to assume everything is OK for now...")
                                                                           
        return traj


    #oxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxo
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

                    # intialize an empty list of atoms
                    local_atoms = []

                    # get every atom in this chain
                    for atom in chain.atoms:
                        protein_atoms.append(atom.index)

                    protein_atoms.extend(local_atoms)
            else:
                
                # initialize an empty list of atoms
                local_atoms = []
                
                for res in chain.residues:

                    # for each residue in that chain ask if that residue is valid
                    if res.name in self.valid_residue_names:

                        # if yes 
                        local_atoms.extend([a.index for a in res.atoms])

                protein_atoms.extend(local_atoms)

            

        return SSProtein(trajectory.atom_slice(protein_atoms))
        

    #oxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxo
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

        # flag that gets set to true if this is a single chain simulation
        single_chain_sim = False
        
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
                        ssio.debug_message('Skipping residue %s from %s' %(chain.residue(0).name, chain))
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

        # after parsing define number of chains we're gonna have... We do this because
        # the ONLY scenario (for now) where we want to use the CG check is if we have
        # a single chain simulation. 
        for local_chain_atoms in chainAtoms:

            if single_chain_sim:
                raise SSException('Error parsing trajectory; this is a bug and this should never happen.')

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
                raise SSException('After extracting a protein subtrajectory, the first resid is not 0. This may reflect a bug, or you may not be using MDTraj 1.9.5')
                
            # if the number of atoms in the subtrajectory is the same as the 
            # the number of atoms in the FULL trajectory we know there is only
            # a single chain here - no need to pay the penalty of initializing                        
            if PT.n_atoms == trajectory.n_atoms:                                
                single_chain_sim = True
            
            # add that SSProtein to the ever-growing proteinTrajectory list. NOTE that 
            # THIS is the line that is the source of most of the slowness for trajectory
            # loading.... 
            proteinTrajectoryList.append(SSProtein(PT))
            

        if len(proteinTrajectoryList) == 0:
            ssio.warning_message('No protein chains found in the trajectory')
        
        return (proteinTrajectoryList, single_chain_sim)




    #oxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxo
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
            res_string = ''
            for r in group:
                res_string = res_string + " %i" % (int(r)) 
            
            # select atoms based on the resid string
            local_atoms = topology.select('resid %s' %(res_string))
            
            if len(local_atoms) == 0:
                ssio.warning_message('In residue group [%s ...] no residues in the trajectory were found...' %(str(group)[0:8]))

            else:
                group_atoms.append(local_atoms)

        # for each protein group that we have atomic indices
        # for cycle through and create sub-trajectories
        proteinTrajectoryList=[]
        
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
                raise SSException('After extracting a protein subtrajectory, the first resid is not 0. This may reflect a bug, or you may not be using MDTraj 1.9.5')
                
            proteinTrajectoryList.append(SSProtein(PT))

        if len(proteinTrajectoryList) == 0:
            ssio.warning_message('No protein chains found in the trajectory')


        return proteinTrajectoryList


    #oxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxo
    #
    #
    @lazy_loading_single_protein_trajectory
    def get_overall_radius_of_gyration(self):
        """
        Function which returns the per-frame OVERALL radius of gyration for 
        every  protein residue in the trajectory. For systems with multiple 
        protein chains,  all chains are combined together. For systems with 
        a single protein chain, this function offers no advantage over 
        interacting directly with the SSProtein object in the 
        ``.proteinTrajectoryList``.

        WARNING: This DOES NOT perform any PBC correction.

        Parameters
        -------------
        None

        Returns
        ----------
        np.ndarray 
            Returns a numpy array with per-frame instantaneous radius of 
            gyration
        
        """

        return self.__single_protein_traj.get_radius_of_gyration()

    #oxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxo
    #
    #
    @lazy_loading_single_protein_trajectory
    def get_overall_asphericity(self):
        """
        Function which returns the per-frame OVERALL asphericity for every
        protein residue in the trajectory. For systems with multiple protein 
        chains,  all chains are combined together. For systems with a single 
        protein chain, this function offers no advantage over interacting 
        directly with the SSProtein object in the underlying 
        ``.proteinTrajectoryList`` object.

        WARNING: This DOES NOT perform any PBC correction.

        Parameters
        -------------
        None

        Returns
        ----------
        np.ndarray 
            Returns a numpy array with per-frame instantaneous asphericity
        """

        return self.__single_protein_traj.get_asphericity()

    #oxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxo
    #
    #
    @lazy_loading_single_protein_trajectory
    def get_overall_hydrodynamic_radius(self):
        """
        Function which returns the per-frame OVERALL hydrodynamic radius for 
        everyprotein residue in the trajectory. For systems with multiple 
        protein chains,  all chains are combined together. For systems with 
        a single protein chain, this function offers no advantage over 
        interacting directly with the SSProtein object in the underlying 
        ``.proteinTrajectoryList`` object.

        Parameters
        -------------
        None

        Returns
        ----------
        np.ndarray 
            Returns a numpy array with per-frame instantaneous asphericity

        """

        return self.__single_protein_traj.get_hydrodynamic_radius()


    #oxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxo
    #
    #
    def get_interchain_distance_map(self, proteinID1, proteinID2, mode='CA', periodic=False):
        """        
        Function which returns two matrices with the mean and standard 
        deviation distances between the residues in resID1 from 
        proteinID1 and resID2 from proteinID2.
        
        This computes the (full) intramolecular distance map, where the 
        "distancemap" function computes the intermolecular distance map.

        Specifically, this allows the user to define two distinct chains
        (i.e.  an "interchain" distance map).        

        Obviously, this only makes sense if your system has two separate 
        protein objects defined, but in principle the output from::        

            TrajObj # TrajOb is an SSTrajectory object

            TrajObj.get_interchain_distance_map(0,0)

        would be the same as::

            TrajObj # TrajOb is an SSTrajectory object
        
            TrajObj.proteinTrajectoryList[0].get_distance_map()

        This is actually a useful sanity check!
        
        Parameters
        ------------

        proteinID1 : int
            The ID of the first protein of the two being considered, where
            the  ID is the proteins position in the 
            ``self.proteinTrajectoryList``  list.             

        proteinID2 : int
            The ID of the second protein of the two being considered, where 
            the ID is the proteins position in the 
            ``self.proteinTrajectoryList`` list
            
        mode : str (default = 'CA')
            String, must be one of either 'CA' or 'COM', where CA means alpha
            carbon and COM means center of mass. 

        periodic : bool (default = False)
            Flag which if distances mode is passed as anything other than 'atom' 
            then this determines if the minimum image convention should be used. 
            Note that this is only available if pdb crystal dimensions are provided, 
            and in general it's better to set this to false and center the molecule 
            first. Default = False.

        Returns
        ---------

            Tuple

                - **distanceMap** is an [n x m] numpy matrix where n and m are the number of proteinID1 residues and proteinID2 residues. Each position  in the matrix corresponds to the mean distance between those two  residues over the course of the simulation.
                - **stdMap** is an [n x m] numpy matrix where n and m are the number of proteinID1 residues and proteinID2 residues. Each position in the matrix corresponds to the standard deviation of the distances between those two residues.
        
        """

        ssutils.validate_keyword_option(mode, ['CA', 'COM'], 'mode')
        
        # get SSProtein objects for the two IDs passed (could be the same)
        P1 = self.proteinTrajectoryList[proteinID1]        
        P2 = self.proteinTrajectoryList[proteinID2]

        # create the empty distance maps
        p1_residues = P1.resid_with_CA
        p2_residues = P2.resid_with_CA
        map_shape   = (len(p1_residues), len(p2_residues))
        distanceMap = np.zeros(map_shape)
        stdMap      = np.zeros(map_shape)

        for r1 in p1_residues:
            if mode == 'COM':
                COM_1 = P1.get_residue_COM(r1)
            else:
                COM_1 = P1.get_residue_COM(r1, atom_name='CA')

            p1_index = copy(r1)
            if P1.ncap:
                p1_index = r1 - 1

            for r2 in p2_residues:

                if mode == 'COM':
                    COM_2 = P2.get_residue_COM(r2)
                else:
                    COM_2 = P2.get_residue_COM(r2, atom_name='CA')

                p2_index = copy(r2)
                if P2.ncap:
                    p2_index = r2 - 1
                
                # compute distance... 
                if periodic:
                    d = sstools.get_distance_periodic(COM_1, COM_2, self.unitcell[0], 'cube')
                else:
                    d = np.linalg.norm(COM_1 - COM_2, axis=1)                    

                    # old way
                    #d = np.sqrt(np.square(np.transpose(COM_1)[0] - np.transpose(COM_2)[0]) + np.square(np.transpose(COM_1)[1] - np.transpose(COM_2)[1])+np.square(np.transpose(COM_1)[2] - np.transpose(COM_2)[2]))
    

                
                distanceMap[p1_index, p2_index] =  np.mean(d, 0)
                stdMap[p1_index, p2_index]    =  np.std(d, 0)
                
        return (distanceMap, stdMap)


    #oxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxo
    #
    #
    def get_interchain_contact_map(self, proteinID1, proteinID2, threshold=5.0, mode='atom', A1='CA', A2='CA', periodic=False, verbose=False):
        
        """        
        Function which returns a matrix with inter-residue contact fractions 
        i.e., this returns the fraction of simulation that each residue from
        one chain is in contact with each residue from the other chain, where
        "contact" is defined as the inter-residue distance being below the
        passed threshold.

        By *default* the mode here is CA-CA distance, which, depending on the
        question, may not be what you want. See the various options under
        'mode' for alternatives.

        Note that this analysis can take some time for large trajectories. If
        this is an issue consider setting the verbose flag to True, and the
        progress will be printed. By default this is off, but it can be 
        reassuring to confirm things are running.
                
        Parameters
        ------------

        proteinID1 : int
            The ID of the first protein of the two being considered, where the 
            ID is the protein's index in the `self.proteinTrajectoryList` 
            list.             

        proteinID2 : int
            The ID of the second protein of the two being considered, where 
            the ID is the protein's inde in the 
            `self.proteinTrajectoryList` list.

        threshold : float 
            Distance that is used as the distance cutoff to define something
            as a contact or not

        mode : str (default = 'atom')
            Mode allows the user to define different modes for computing 
            atomic distance.

            The default is ``atom`` whereby a pair of atoms (A1 and A2) are 
            provided. Other options are detailed below and are identical to 
            those offered by mdtraj in compute_contacts.

            Note that if modes other than ``atom`` are used the A1 and A2 
            options are ignored.

            + ``ca`` - same as setting ``atom`` and then defining atoms \
                       1 and 2 (A1 and A2) as CA. 
                      
            + ``closest`` - closest atom associated with each of the \
                            residues, i.e. the point of closest approach \
                            between the two residues.
                                                        
            + ``closest-heavy`` - same as `'closest'`, except only non-\
                                  hydrogen atoms are considered.
                                  
            + ``sidechain`` - closest atom where that atom is in the\
                              sidechain.

            + ``sidechain-heavy`` - closest atom where that atom is\
                                    in the sidechain and is heavy.

        A1 : str (default = 'CA')
            Atom name of the atom in R1 we're looking at. 

        A2 : str (default = 'CA')
            Atom name of the atom in R2 we're looking at. 

        periodic : bool (default = False)
            Flag which if distances mode is passed as anything other than 'atom' 
            then this determines if the minimum image convention should be used. 
            Note that this is only available if pdb crystal dimensions are provided, 
            and in general it's better to set this to false and  center the molecule 
            first. Default = False. 

        verbose : bool
            Flag which, if set to true, will print each iteraction through the outer
            loop of each residue in the protein as the distances are calculated.
            Default = False.

        Returns
        ---------
        np.ndarray [n x m] 
            Returns an matrix with inter-residue contact fractions reported at each
            intersection.
        
        """

        # get number of residues/bases for the two proteins
        n_res_P1 = self.proteinTrajectoryList[proteinID1].n_residues
        n_res_P2 = self.proteinTrajectoryList[proteinID2].n_residues

        # 
        all_contact_fractions = []

        # cycle over each residue in protein 1
        for p1_res_idx in range(0, n_res_P1):
            if verbose:
                print(f'On {p1_res_idx} of {n_res_P1}')
    
            tmp = []

            # cycle over each residue in protein 2
            for p2_res_idx in range(0, n_res_P2):

                # calculate distances between unique residues in P1 and P2
                distances = self.get_interchain_distance(proteinID1, proteinID2, p1_res_idx, p2_res_idx, A1=A1, A2=A2, mode=mode, periodic=periodic)

                
                # find number of distances below the threshold and calculat as a 
                # a fraction of all frames
                contact_fraction = np.sum(distances<threshold)/self.n_frames

                # add that fraction to the temporary 
                tmp.append(contact_fraction)
    
            all_contact_fractions.append(tmp)   

        return np.array(all_contact_fractions)



    #oxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxo
    #
    #
    def get_interchain_distance(self, proteinID1, proteinID2, R1, R2, A1='CA', A2='CA', mode='atom', periodic=False):
        """
        Function which returns the distance between two specific atoms on 
        two residues, or between two residues based on mdtraj's atom 
        selection mode rules (discussed below). Required input are protein        
        ID selectors and the resid being used. Resids should be used as 
        would be normally used for the SSProtein objects associated with 
        proteinID1 and proteinID2.
        
        For inter-atomic distances, atoms are selected from the passed 
        residue and their 'name' field from the topology selection language 
        (e.g. "CA", "CB" "NH" etc). By default CA atoms are used, but one can 
        define any residue of interest. We do not perform any sanity checking 
        on the atom name - this gets really hard - so have an explicit 
        try/except block which will warn you that you've probably 
        selected an illegal atom name from the residues.

        For inter-residue distances the associated rules are defined by the 
        'mode' selector. By default mode is set to 'atom', which means the 
        variables A1 and A2 are used (with CA as default) to define inter-
        residue distance. However, if one of the other modes are used the 
        A1/A2 parameters are ignored and alternative rules for computing 
        inter-residue distance are used. These modes are detailed below.

        Distance is returned in Angstroms.

        Parameters
        ----------

        proteinID1 : int
            Index of the first protein of interest

        proteinID2 : int
            Index of the second protein of interest

        R1 : int
            Residue index of first residue

        R2 : int
            Residue index of second residue

        A1 : str (default = 'CA')
            Atom name of the atom in R1 we're looking at. 

        A2 : str (default = 'CA')
            Atom name of the atom in R2 we're looking at. 

        mode : str (default = 'atom')
            Mode allows the user to define different modes for computing atomic 
            distance.

            The default is ``atom`` whereby a pair of atoms (A1 and A2) are 
            provided. Other options are detailed below and are identical to 
            those offered by mdtraj in compute_contacts.

            Note that if modes other than ``atom`` are used the A1 and A2 
            options are ignored.

            + ``ca`` - same as setting ``atom`` and then defining atoms \
                       1 and 2 (A1 and A2) as CA. 
                      
            + ``closest`` - closest atom associated with each of the \
                            residues, i.e. the point of closest approach \
                            between the two residues.
                                                        
            + ``closest-heavy`` - same as `'closest'`, except only non-\
                                  hydrogen atoms are considered.
                                  
            + ``sidechain`` - closest atom where that atom is in the\
                              sidechain.

            + ``sidechain-heavy`` - closest atom where that atom is\
                                    in the sidechain and is heavy.
                                    

        periodic : bool (default = False)
            Flag which if distances mode is passed as anything other than 'atom'
            then this determines if the minimum image convention should be used.
            Note that this is only available if pdb crystal dimensions are
            provided, and in general it's better to set this to false and
           center the molecule first. Default = False.

        Returns
        -----------
        np.array
            Returns a 1D numpy array with the distance-per-frame betwee the specified residues

        """

        # check mode keyword is valid
        allowed_modes = [ 'atom', 'ca', 'closest', 'closest-heavy', 'sidechain', 'sidechain-heavy' ]
        if mode not in allowed_modes:
            raise SSException("Provided mode keyword must be one of 'atom', 'ca', 'closest', 'closest-heavy', 'sidechain', or 'sidechain-heavy'. Provided keyword was [%s]" % (mode))

        # get SSProtein objects for the two IDs passed (could be the same)        
        try:
            P1 = self.proteinTrajectoryList[proteinID1]
            P2 = self.proteinTrajectoryList[proteinID2]

        except IndexError as e:
            raise SSException('In get_interchain_distance(): When selecting protein indices %i and %i at least one of these was out of range (indices are from 0...%i)' % (proteinID1, proteinID2, len(self.proteinTrajectoryList)-1))
            

        # next build a new trajectory that contains ONLY the two residues selected
        local_atoms1 = P1.topology.select('resid %i' % (R1))
        local_atoms2 = P2.topology.select('resid %i' % (R2))

        if len(local_atoms1) == 0:
            raise SSException("In get_interchain_distance(): When selecting resid %i from proteinID1 found no atoms" %(R1))

        if len(local_atoms2) == 0:
            raise SSException("In get_interchain_distance(): When selecting resid %i from proteinID2 found no atoms" %(R2))
            
        subtraj_p1 = P1.traj.atom_slice(local_atoms1)
        subtraj_p2 = P2.traj.atom_slice(local_atoms2)
        
        # this is now a subtrajectory which in principle contains just two residues. We can check
        # this to ensure that the trajectory has exactly 2 residues
        full_subtraj = subtraj_p1.stack(subtraj_p2)
        full_subtraj_residues = [i for i in full_subtraj.topology.residues]
        if len(full_subtraj_residues) != 2:
            raise SSException("In get_interchain_distance(): When passed in two residues (R1=%i, R2=%i) in proteins %i and %i found multiple residues (%i)...these resids could not be found " %(R1, R2, proteinID1, proteinID2, len(full_subtraj_residues)))

        
        # if we're looking at a specific pair of atoms (note we use resid 0 and 1 because we KNOW this trajectory only has 2 residues and we know R1 is 0 and R2 is 1
        if mode == 'atom':

            atom1 = full_subtraj.topology.select('resid 0 and name "%s"' % A1)
            if len(atom1) != 1:
                raise SSException("In get_interchain_distance() when selecting atom %s from residue %i in protein %i no atoms were found " % (A1, R1,  proteinID1))                

            COM_1 = 10*md.compute_center_of_mass(full_subtraj.atom_slice(atom1))

            atom2 = full_subtraj.topology.select('resid 1 and name "%s"' % A2)
            if len(atom2) != 1:
                raise SSException("In get_interchain_distance() when selecting atom %s from residue %i in protein %i no atoms were found " % (A2, R2,  proteinID2))                

            COM_2 = 10*md.compute_center_of_mass(full_subtraj.atom_slice(atom2))
            
            # finally compute distances. Use minimum image convention if the periodic keyword is passed
            if periodic:
                distances = sstools.get_distance_periodic(COM_1, COM_2, self.unitcell[0], 'cube')
                
            else:

                # revised way
                distances = np.linalg.norm(COM_1 - COM_2, axis=1)

                # old way
                # distances = np.sqrt(np.square(np.transpose(COM_1)[0] - np.transpose(COM_2)[0]) + np.square(np.transpose(COM_1)[1] - np.transpose(COM_2)[1])+np.square(np.transpose(COM_1)[2] - np.transpose(COM_2)[2]))

        else:

            # TODO: Documentation missing!
            # use the compute_contacts() function from mdtraj, multiplying by 10 because this will
            # by default give you numbers that... 
            distances = 10*md.compute_contacts(full_subtraj, [[0, 1]], scheme=mode, periodic=periodic)[0].ravel()

        return distances
      

def __load_trajectory(trj_filename, top_filename, **kwargs):
    """ Private helper function for loading trajectories in parallel

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
        raise RuntimeError(f"Failed to load SSTrajectory for trajectory '{trj_filename}' with topology '{top_filename}'. "
                           f"Perhaps you supplied the wrong topology?: {e}")

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

    # Assume the same topology file for all trajectories if only one is provided
    if len(top_filenames) < len(trj_filenames) and len(top_filenames) == 1:
        top_filenames = [top_filenames] * len(trj_filenames)

    # Partially apply load_trajectory with topology file path and kwargs
    partial_load = partial(__load_trajectory, **kwargs)

    # Parallelize load with **kwargs
    with Pool(processes=n_procs) as pool:
        trjs = pool.starmap(partial_load, zip(trj_filenames, top_filenames))

    return trjs