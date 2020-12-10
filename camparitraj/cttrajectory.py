"""
cttrajectory.py

This is where some stuff will be described

"""

##
##                                       _ _              _ 
##   ___ __ _ _ __ ___  _ __   __ _ _ __(_) |_ _ __ __ _ (_)
##  / __/ _` | '_ ` _ \| '_ \ / _` | '__| | __| '__/ _` || |
## | (_| (_| | | | | | | |_) | (_| | |  | | |_| | | (_| || |
##  \___\__,_|_| |_| |_| .__/ \__,_|_|  |_|\__|_|  \__,_|/ |
##                     |_|                             |__/ 
##
##
## Alex Holehouse (Pappu Lab)
## Simulation analysis package
## Copyright 2014 - 2019
##

import mdtraj as md
import numpy as np
from .configs import *
from .ctdata  import ALL_VALID_RESIDUE_NAMES

from .ctprotein import CTProtein
from .ctexceptions import CTException
from . import ctutils
from . import ctio

def testfunct():
    """
    This is a test function for documentation - do we build from source?
    """
    print('YEAH OK')




class CTTrajectory:
    """
    CTrajectory class that holds a single simulation trajectory object. 
    

    """

    #oxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxo
    #
    #
    def __init__(self, trajectory_filename=None, pdb_filename=None, TRJ=None, protein_grouping=None, pdblead=False, debug=False):
        """
        CAMPARITraj trajectory object initializer. 

        CAMPARITraj will, by default, extract out the protein component from 
        your trajectory automatically, which lets you ask questions about the
        protein only (i.e. without salt ions getting in the way).
        
        Note that by default the mechanism by which individual proteins are
        identified is by cycling through the unique chains and determining 
        if they are protein or not. You can also provide manual grouping
        via the protein_grouping option, which lets you define which
        residues should make up an individual protein. This can be useful
        if you have multiple proteins associated with the same chain, which
        happens in CAMPARI if you have more than 26 separate chains (i.e.
        every protein after the 26th is the 'Z' chain).
               
        Parameters
        ----------        
        trajectory_filename : str
            Filename which contains the trajectory file of interest. Normally \
            this is `__traj.xtc` or `__traj.dcd`.

        pdb_filename : str
            Filename which contains the pdb file associated with the trajectory \
            of interest. Normally this is `__START.pdb`.

        TRJ : mdtraj.Trajectory
            It is sometimes useful to re-defined a trajectory and create a new CTTraj \
            object from that trajectory. This could be done by writing that new trajectory \
            to file, but this is extremely slow due to the I/O impact of reading/writing \
            from disk. If an mdtraj trajectory objected is passed, this is used as the \
            new trajectory from which the CTTrajectory object is constructed. 

            Default = None

        protein_grouping : list of lists of ints
            Lets you manually define protein groups to be considered independently.
        
            Default = None

        pdblead : bool
            Lets you set the PDB file (which is normally ONLY used as a topology \
            file) to be the first frame of the trajectory. This is useful when \
            the first PDB file holds some specific reference information which \
            you want to use (e.g. RMSD or Q). 

            Default = False

        debug : book
            Prints warning/help information to help debug weird stuff during initial trajectory read-in. 
            Default = False.
        """
        
        # first we decide if we're reading from file or from an existing trajectory
        if (trajectory_filename is None) and (pdb_filename is None):
            if TRJ is None:
                raise CTException('No input provided! Please provide ether a pdb and trajectory file OR a pre-formed traj object')
                
            # note the [:] means this is a COPY!
            self.traj = TRJ[:]
        else:
            if (trajectory_filename is None):
                raise CTException('No trajectory file provided!')
            if (pdb_filename is None):
                raise CTException('No PDB file provided!')

            # read in the raw trajectory
            self.traj = self.__readTrajectory(trajectory_filename, pdb_filename, pdblead)


        # Next, having read in the trajectory we parse out into proteins
        # extract a list of protein trajectories where each protein is assumed
        # to be in its own chain
        if protein_grouping == None:
            (self.proteinTrajectoryList, self.resid_offset_list, self.atom_offset_list) = self.__get_proteins(self.traj, debug)        
        else:
            (self.proteinTrajectoryList, self.resid_offset_list, self.atom_offset_list)  = self.__get_proteins_by_residue(self.traj, protein_grouping, debug)

        
        self.num_proteins = len(self.proteinTrajectoryList)
        self.n_frames = len(self.traj)


    def  __repr__(self):
        return "CTTrajectory (%s): %i proteins and %i frames" % (hex(id(self)), self.num_proteins, self.n_frames)

    def __len__(self):
        return (self.num_proteins, self.n_frames)



    #oxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxo
    #
    #
    def __readTrajectory(self, trajectory_filename, pdb_filename, pdblead):
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
        ------------------
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
   
        """

        # straight up read the trajectory first using mdtraj's awesome
        # trajectory reading facility
        traj =  md.load(trajectory_filename, top=pdb_filename)
                    
        # check unit cell lengths
        try:
            uc_lengths = traj.unitcell_lengths[0]

            # this is s custom warning for a specific edge-case we encounter a lot
            if (uc_lengths[0] == 0 or uc_lengths[1] == 0 or uc_lengths[2] == 0):
                ctio.warning_message("Trajectory file unit cell lengths are zero for at least one dimension. This is a probably a bug with an FRC generated __START.pdb file, because in the old version of CAMPARI used to do grid based FRC calculations the unit cell dimensions are not written correctly. This may cause issues but we're going to assume everything is OK for now. Check the validity of any analysis output. If you're worried, you can use the following workaround.\n\n:::: WORK AROUNDS ::::\nSimply run\n\ntrjconv -f __traj.xtc -s __START.pdb -box a b c -o frc.xtc \n\nAn then \n\ntrjconv -f frc.xtc -s __START.pdb -box a b c -o start.pdb -dump 0\n\n\nHere\n-f n__traj.xtc   : defines the trajectory file\n-s __start.pdb   : defines the pdb file used to parse the topology\n-box a b c       : defines the box lengths **in nanometers**\n-o frc.xtc       : is the name of the new trajectory file with updated box lengths\nSelect 0 (system) when asked to 'Select group for output'.The second step creates the equivalent PDB file with the header-line correctly defining the box unit cell lengths and angles. These two new files should then be used for analysis.\n\nAs an example, if my FRC simulation had a sphere radius of 100 angstroms then my correction command would look something like \n\ntrjconv -f __traj.xtc -s __START.pdb -box 20 20 20 -o frc.xtc\ntrjconv -f frc.xtc -s __START.pdb -box 20 20 20 -o start.pdb -dump 0")


        except TypeError:
            ctio.warning_message("Warning: UnitCell lengths were not provided... This may cause issues but we're going to assume everything is OK for now...")
        

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
                    ctio.warning_message('........................\nWARNING:\nThe unit cell dimensions of the PDB file and trajectory file did not match, specifically\nPDB file = [ %s ]\nXTC file = [ %s ]\nThis may cause issues if native MDTraj untilities are used (and potentially for CTraj utilities that are based on these. It is not necessarily an issue, but PLEASE sanity check your outcome. To be save we reeommend editing the PDB-file untilcell dimenions to match.')

            except TypeError:
                ctio.warning_message("Warning: UnitCell lengths were not provided... This may cause issues but we're going to assume everything is OK for now...")
                                   
                                        
        return traj



    #oxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxo
    #
    #
    def __get_proteins(self, trajectory, debug):
        """
        Internal function that takes an MDTraj trajectory and returns a list of mdtraj 
        trajectory objects corresponding to each protein in the system, ASSUMING that 
        each protein is in its own chain.
        
        The way this works is to cycle through each chain, identify if that chain
        contains protein or not, and if it does grab all the atoms in that chain
        and perform an atomslice using those atoms on the main trajectory.
        

        Parameters
        -----------
        trajectory : mdtraj.Trajectory
            An already parsed trajectory object (i.e. checked for CAMPARI-
            relevant defects such as unitcell issues etc)

        Returns
        ---------
        tuple :
            Returns a tuple with three lists:
        
            proteinTrajectoryList - contains a list of 0 or more CTProtein objcts        
            resid_offset_list     - contains a list of 0 or more integers which are 
                                    resid offset values
            atom_offset_list      - contains a list of 0 or more integers which are
                                    atom offset values
 
            Note all three lists must be the same length (by definition)
        
        """

        
        atom_offset_list  = []
        resid_offset_list = []

        # extract full system topology
        topology = trajectory.topology

        chainAtoms = []
        
        # for each chain in this toplogy determine if the 
        # first residue is protein or not
        for chain in topology.chains:

            # if the first residue in the chain is protein
            # note that a formic acid cap ('FOR') is not recognized as protein
            # so we include an edgecase here for that
            if chain.residue(0).name in ALL_VALID_RESIDUE_NAMES:

                # intialize an empty list of atoms
                local_atoms = []

                # save that index value associated with
                # the first residue in the chain. 
                resid_offset_list.append(chain.residue(0).index)

                for atom in chain.atoms:
                    local_atoms.append(atom.index)

                chainAtoms.append(local_atoms)

                # save first atom associated with this chain for offset
                atom_offset_list.append(local_atoms[0])

                
                ## Code below was the old way of builing the atom sets which
                ## as of mdtraj 1.9.3 the above seems more efficient
                """
                # for each residue in the chain
                for residue in chain.residues:
                                        
                    # for each atom in the residue
                    for atom in residue.atoms:

                        # append all those atoms to our list of atom indicies
                        atoms.append(atom.index)

                # now add all the atom indices for this protein
                # chain as a single list to the super-list of 
                # chain atoms
                chainAtoms.append(atoms)
                """
            else:
                if debug:
                    ctio.debug_message('Skipping residue %s from %s' %(chain.residue(0).name, chain))


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

            # gets the resid offset in a way that is ensures internal
            # consistency for the CTProtein object
            resid_offset = PT.topology.chain(0).residue(0).index
            
            # add that trajectory, along the index value associated with
            # the resid offset 
            proteinTrajectoryList.append(CTProtein(PT, resid_offset))

        if len(proteinTrajectoryList) == 0:
            ctio.warning_message('No protein chains found in the trajectory')

        return (proteinTrajectoryList, resid_offset_list, atom_offset_list)




    #oxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxo
    #
    #
    def __get_proteins_by_residue(self, trajectory, residue_grouping, debug):
        """
        Internal function which returns a list of mdtraj trajectory objects corresponding 
        to each protein where we *explicitly* define the residues in each protein.

        Unlike the `__get_proteins()` function, which doesn't require any manual
        input in identifying the proteins, here we provide a list of groups, where 
        each group is the set of residues associated with a protein.

        The way this works is to cycle through each group, and for each residue in 
        each group grabs all the atoms and uses these to carry out an atomslice 
        on the full trajectory. 

        Parameters
        -----------
        trajectory : mdtraj.Trajectory
            An already parsed trajectory object (i.e. checked for CAMPARI-
            relevant defects such as unitcell issues etc)

        residue_grouping : list of list of integers
            Must be a list containing one or more lists, where each internal list
            contains a set of monotonically increasing residues (which correspond
            to the full protein trajectory). In other words, each sub-list defines
            a single protein. The integer indexing here - importantly - uses the 
            CAMPARITraj internal residue indexing, meaning that indexing begins at 0 from
            the first residue in the PDB file.

        """
        
        atom_offset_list  = []
        resid_offset_list = []

        group_atoms = []

        # extract full system topology
        topology = trajectory.topology
        
        # for each chain in this toplogy determine if the 
        # first residue is protein or not        
        for group in residue_grouping:
            local_atoms = [] 
            
            for res_index in group:

                # get residue object
                residue = topology.residue(res_index)
                resid_offset_list.append(residue.index)
                                
                # for each atom in the residue
                for atom in residue.atoms:                    
                    # append all those atoms to our list of atom indicies
                    local_atoms.append(atom.index)

            # save first atom associated with this group for offset
            atom_offset_list.append(local_atoms[0])

            # now add all the atom indices for this protein
            # chain as a single list to the superlist of 
            # chain atoms
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
            # consistency for the CTProtein object
            resid_offset = PT.topology.chain(0).residue(0).index

            proteinTrajectoryList.append(CTProtein(PT, resid_offset))

        if len(proteinTrajectoryList) == 0:
            ctio.warning_message('No protein chains found in the trajectory')


        return (proteinTrajectoryList, resid_offset_list, atom_offset_list)


    #oxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxo
    #
    #
    def get_interchain_distance_map(self, proteinID1, proteinID2, resID1=None, resID2=None):
        """        
        Function which returns two matrices with the mean and standard deviation distances
        between the residues in resID1 from proteinID1 and resID2 from proteinID2

        This computes the (full) intramolecular distance map, where the "distancemap"
        function computes the intermolecular distance map.

        Specifically, this allows the user to define two distinct chains (i.e. an "interchain"
        distance map).

        Obviously this only makes sense if your system has two separate protein objects
        defined, but in principle the output from:

            `intra_chain_distance_Map(0,0)`

        would be the same as:
        
            `proteinTrajectoryList[0].get_distance_map()`

        This is actually a useful sanity check!
        
        Parameters
        ----------

        proteinID1 : int
            The ID of the first protein of the two being considered, where the ID is the proteins position in the
            `self.proteinTrajectoryList` list. 

        proteinID2 : int
            The ID of the second protein of the two being considered, where the ID is the proteins position in the
            `self.proteinTrajectoryList` list

        resID1 : list of integers, default=None
            Is the list of residues from protein 1 we're considering. If this is left as None (default), then it 
            is assumed that all residues in proteinID1 should be used

        resID2 : list of integers, default=None
            Is the list of residues from protein 2 we're considering.If this is left as None (default), then it 
            is assumed that all residues in proteinID2 should be used

        Returns
        -------
        tuple : tuple containing `distanceMap` and `STDMap`

        distanceMap : numpy matrix
            Is an [n x m] matrix where n and m are the number of proteinID1 residues
            and proteinID2 residues. Each position in the matrix corresponds to the
            mean distance between those two residues over the course of the simulation.

        stdMap : numpy matrix
            Is an [n x m] matrix where n and m are the number of proteinID1 residues
            and proteinID2 residues. Each position in the matrix corresponds to the
            standard devaiation associated with the distances between those two
            residues.

        
        """
        
        # get CTProtein objects for the two IDs passed (could be the same)
        P1 = self.proteinTrajectoryList[proteinID1]        
        P2 = self.proteinTrajectoryList[proteinID2]        

        P1_atom_offset = self.atom_offset_list[proteinID1]
        P2_atom_offset = self.atom_offset_list[proteinID2]

        # get all the raw CA atom indices we're gonna be working with
        CA_p1_raw = P1.get_multiple_CA_index(resID1)
        CA_p2_raw = P2.get_multiple_CA_index(resID2)

        # atom indices in subtrajectories (i.e. in the trajectories found in proteinTrajectory)
        # always start from 0 from their own trajectory! To get around this and relate atom
        # numbers in subtrajectories back to the FULL trajectory we have to add the atomic
        # offset associated with each subtrajectory
        CA_p1 = []
        CA_p2 = []
        for atom in CA_p1_raw:
            CA_p1.append(atom + P1_atom_offset)

        for atom in CA_p2_raw:
            CA_p2.append(atom + P2_atom_offset)
                
        # create the empty distance maps
        distanceMap = np.zeros([len(CA_p1),len(CA_p2),])
        stdMap      = np.zeros([len(CA_p1),len(CA_p2),])
        
        # calculate the FULL distance map (so 2N rather than N time).
        # note this is actually the non-redundant map, because only in the limit
        # of perfect sampling is it necesessarily true that
        # 
        # P1-R5 ::: P2-R10 == P1-R10 :::: P2-R5
        #
        #
        
        count = 0
        for CA1 in CA_p1:
            pairs = []
            for CA2 in CA_p2:

                pairs.append([CA1, CA2])
                
            data = 10*md.compute_distances(self.traj, np.array(pairs))
            print("On residue %i (protein 1) out of %i" % (count+1, len(CA_p1)))

            # calculate mean and standard deviation over the trajectory data
            mean_data = np.mean(data,0)
            std_data = np.std(data,0)

            # assign rows in the matrix
            distanceMap[count] = mean_data
            stdMap[count]      = std_data

            # increment the outercounter
            count = count + 1

        return (distanceMap, stdMap)


    #oxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxoxoxoxoxoxoxoxoxoxooxoxo
    #
    #
    def get_interchain_distance(self, proteinID1, proteinID2, R1, R2, A1='CA', A2='CA', stride=1, mode='atom'):
        """
        Function which returns the distance between two specific atoms on two residues, or between 
        two residues based on mdtraj' atomselection mode rules (discussed below). Required input are protein
        ID selectors and the resid being used. Resids should be used as would be normally used for the CTProtein
        objects associated with proteinID1 and proteinID2.

        For inter-atomic distances, atoms are selected from the passed residue and their 'name' field from the topology 
        selection language (e.g. "CA", "CB" "NH" etc). By default CA atoms are used, but one can define
        any residue of interest. We do not perform any sanity checking on the atom name - this gets 
        really hard - so have an explicit try/except block which will warn you that you've probably 
        selected an illegal atom name from the residues.

        For inter-residue distances the associated rules are defined by the 'mode' selector. By default
        mode is set to 'atom', which means the variables A1 and A2 are used (with CA as default) to define
        inter-residue distance. However, if one of the other modes are used the A1/A2 parameters are ignored
        and alternative rules for computing inter-residue distance are used. These modes are detailed below.

        Distance is returned in Angstroms.

        Parameters
        ----------
        R1 : int
            Residue index of first residue

        R2 : int
            Residue index of second residue

        A1 : str
            Atom name of the atom in R1 we're looking at. Default = 'CA'

        A2 : str
            Atom name of the atom in R2 we're looking at. Default='CA'

        stride : int
            Defines the spacing between frames to compare with - i.e. take every $stride-th frame.
            Setting `stride=1` would mean every frame is used, which would mean you're doing an
            all vs. all comparisons, which would be ideal BUT may be slow. Default = 1

        mode : str
            Mode allows the user to define different modes for computing atomic distance.

            The default is 'atom' whereby a pair of atoms (A1 and A2) are provided. Other options
            are detailed below and are identical to those offered by mdtraj in compute_contacts.

            Note that if modes other than 'atom' are used the A1 and A2 options are ignored.

            + `'ca'` - same as setting `'atom'` and `A1='CA'` and `A2='CA'`, this uses the C-alpha atoms.
            + `'closest'` - closest atom associated with each of the residues, i.e. the is the point \
                            of closest approach between the two residues.
            + `'closest-heavy'` - same as `'closest'`, except only non-hydrogen atoms are considered.
            + `'sidechain'` - closest atom where that atom is in the sidechain. Note this requires \
                              mdtraj version 1.8.0 or higher.
            + `'sidechain-heavy'` - closest atom where that atom is in the sidechain and is heavy. \
                                    Note this requires mdtraj version 1.8.0 or higher.

        Returns
        -----------
        np.array
            Returns a 1D numpy array with the distance-per-frame betwee the specified residues

        """

        # check mode keyword is valid
        if mode in [ 'closest', 'ca',  'closest-heavy',  'sidechain', 'sidechain-heavy',  'atom']:
            pass
        else:
            raise CTException("Provided mode keyword must be one of 'closest', 'ca', 'closest-heavy', 'sidechain', sidechain-heavy', or 'atom'. Provided keyword was [%s]" % (mode))

        # get CTProtein objects for the two IDs passed (could be the same)
        P1 = self.proteinTrajectoryList[proteinID1]        
        P2 = self.proteinTrajectoryList[proteinID2]        

        # get resids that correspond to the FULL trajectory
        R1_re_referenced = R1 + self.resid_offset_list[proteinID1]
        R2_re_referenced = R2 + self.resid_offset_list[proteinID2]

        try:
            # if atom mode was used
            if mode == 'atom':
                if A1 == 'CA' and A2 == 'CA':
                    if float(stride) != float(1): 
                        subtraj = self.traj.slice(list(range(0, self.n_frames, stride)))
                    else:
                        subtraj = self.traj

                    distances = 10*md.compute_contacts(subtraj, [[R1_re_referenced, R2_re_referenced]], scheme='ca')[0].ravel()
                    
                else:
                    
                    atom1 = P1._CTProtein__residue_atom_lookup(R1,A1)
                    if len(atom1) == 0:
                        raise CTException('Unable to find atom [%s] in residue R1 (%i)' % (A1, R1))

                    TRJ_1 = self.traj.atom_slice(atom1)
                    TRJ_1 = TRJ_1.slice(list(range(0, self.traj.n_frames, stride)))

                    
                    atom2 = P2._CTProtein__residue_atom_lookup(R2,A2)
                    if len(atom2) == 0:
                        raise CTException('Unable to find atom [%s] in residue R1 (%i)' % (A2, R2))

                    TRJ_2 = self.traj.atom_slice(atom1)
                    TRJ_2 = TRJ_2.slice(list(range(0, self.traj.n_frames, stride)))

                    COM_1 = md.compute_center_of_mass(TRJ_1)
                    COM_2 = md.compute_center_of_mass(TRJ_2)


                    distances = 10*np.sqrt(np.square(np.transpose(COM_1)[0] - np.transpose(COM_2)[0]) + np.square(np.transpose(COM_1)[1] - np.transpose(COM_2)[1])+np.square(np.transpose(COM_1)[2] - np.transpose(COM_2)[2]))

        except IndexError as e:
            
            print("<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@")
            print("")
            print("This is likely because one of [%s] or [%s] is not a valid atom type for the residue in question. Full error printed below" %( A1,A2))
            print("")
            print("<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@")
            print(e)
            raise e

        # parse any of the allowed modes in compute_contacts (see http://mdtraj.org/1.8.0/api/generated/mdtraj.compute_contacts.html
        # for more details!)
        if mode == 'closest' or mode == 'ca' or mode == 'closest-heavy' or mode == 'sidechain' or mode == 'sidechain-heavy':
            if float(stride) != float(1): 
                subtraj = self.traj.slice(list(range(0, self.traj.n_frames, stride)))
            else:
                subtraj = self.traj

            
            try:
                distances = 10*md.compute_contacts(subtraj, [[R1,R2]], scheme=mode)[0].ravel()
            except ValueError as e:
                print(e)
                raise CTException('Your current version of mdtraj does not support [%s] - please update mdtraj to 1.8.0 or later to facilitate support. Alternatively this may be because residue %i or %i is not parsed correctly by mdtraj' % (mode, R1, R2)) 

                                
        # note 10* to get Angstroms
        return distances
      

