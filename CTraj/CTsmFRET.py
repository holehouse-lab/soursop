##
################################################
##  ,-----.,--------.                 ,--.    ##
## '  .--./'--.  .--',--.--. ,--,--.  `--'    ##
## |  |       |  |   |  .--'' ,-.  |  ,--.    ##
## '  '--'\   |  |   |  |   \ '-'  |  |  |    ##
##  `-----'   `--'   `--'    `--`--'.-'  /    ##
##                                  '---'     ##
################################################
##
## Alex Holehouse (Pappu Lab)
## Original code and idea for smFRET based analysis by 
## Kiersten Ruff (Pappu Lab)
## Simulation analysis package
## Copyright 2014 - 2018
##

import mdtraj as md
import numpy as np
import scipy.optimize
import os, sys
import random
from itertools import chain

from .CTExceptions import CTsmFRET_Exception
from . import CTUtils

###<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
###
### DYE CONFIGURATION INFO
###
###<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
##
## The following consection holds global variables specific to dyes used by CTsmFRET.
##
##



SUPPORTED_DYES=['A488','A594']
DYE_RESNAMES     = {'A488':'AX4', 'A594':'AM5'} # making of dye name to dye resname as found in PDB file
DYE_LIGHT_SOURCE = {'A488':'C19', 'A594':'C19'} # defines the atom name associated with the photon point source in the dye #PhotoPhysicsLOLs
                                                # AGAIN we don't actually use this, but worth having for reference 
COCOFRET_LOG='COCOFRET_STATUS.log'
MAXCORES=1

CS_BOND = 0.181
CCS_ANG = 108.6
CSC_ANG = 109.5
MAXDISTANG = 15.0
MAXDIST = 0.2

# Refractive index of water (units = unitless)
REFRACTIVE_INDEX = 1.334

# Quantum yeild of dye pairs (units = unitless)
QUANTUM_YEILD = {'A488-A594' : 0.92}

# Overlap integral of dyes M(-1) cm^(-1) nm^4
OVERLAP_INTEGRAL = {'A488-A594' : 2.4e15}

# VALUE from Sebastian for work with ben
DEFAULT_R0_VALS = {'A488-A594':5.900} 


# use MDtraj selection language to define the first ring in the dye - NOTE that for A488 and A594 these rings are actually
# the same BUT we don't want to need to guarentee that so we explicitly define them here. For other dyes you'd need to
# define the first bulky structure to avoid overly-extreme rejection when placing dyes... 
# NB: Actually we don't use this but I've kept it anyway...
FIRST_RING_RES = {'A488': 'name C62 or name H62 or name C63 or name C64 or name N65 or name C66 or name 067 or name O69',
                  'A594': 'name C62 or name H62 or name C63 or name C64 or name N65 or name C66 or name O67 or name O69 '}

# Define the atom IDX in the dyes that corresponds to the first ring
# NOTE these index assuming the FIRST residue after the SG is 0 (!!) - i.e. SG is 'NOT' part of the
# dye
FIRST_RING_IDX  = {'A488': [0, 1, 2, 67, 68, 69, 70, 71], 'A594' : [0, 1, 2, 93, 94, 95, 96, 97]}
DYE_LIGHT_SOURCE_IDX = {'A488': 19, 'A594': 17} # defines the IDX on the dye where the photon source is. # NOTE these index assuming the 
                                                # FIRST residue after the SG is 0 (!!) - i.e. SG is 'NOT' part of the dye

# These two positions define the positions of the two atoms that make up the points which define the 
# dipole in the dyes. For clarity I've included the atom names. The index values here assume that the 0th index
# is the FIRST non sulphur atom in the dye (i.e. using the same indexing scheme as used in the FIRST_RING_IDX
# dictionary above
# A488 = C16 and C21
# A594 = C88 and C73
# The dipole moment defined by these moments goes straight through the LIGHT_SOURCE_IDX atom along the bottom row of the dye
DYE_TRANSITION_DIPOLE_IDX = {'A488':[31,14], 'A594':[26,11]}


# General notes on this code all dye index values used are internal to each dye. i.e. 0 is either the start of the 
# linker (the SG) or the residue after SG
#
#
#
#
#
#
#
#
#
# 
# 
#
#
#
#
#         
#


###<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>


class CTsmFRET:
    """
    
    CTsmFRET is a class that takes a CTProtein objecs, adds FRET dyes to posistions, and calculates
    the FRET efficienciy. 

    This is, in many ways, a purely functional class, but given it's very specific goal
    it is separated out into its own class in the interest of more robust modularity.

    There are a number of internal functions, but it is the public facing functions which are 
    those relevant for users. These are outlined below:


    > run_dye_sampling
      First generation algorithm for placing dyes
    
    
    > run_dye_sampling_independent_placement
      Second generation where dyes are placed independently first    

    > get_efficiencies_from_distances
      Takes output from run_dye_sampling or run_dye_sampling_independent_placement and returns 
      the corresponding FRET efficiencies.

    > get_FRET_effiencies
      Combines run_dye_sampling_independent_placement and get_efficiencies_from_distances into a
      single function that generates FRET efficiencies. 

    > run_COCOFRET
      Third generation approach, where the number of dye distances calculated
      is set in terms of relative convergence of the efficiency value. This is the recommended 
      approach right now. A fourth generation approach is being devloped which takes dye orientation
      into account, and thereafter a fifth generation approach which includes dye dynamics and lifetime
      information.

    """

    # ........................................................................
    #                 
    def __init__(self, CTProteinObject, dye1_name, dye1_pos_w_offset, dye2_name, dye2_pos_w_offset, R0_VALUES=DEFAULT_R0_VALS, verbose=True):
        """
        Initialization function for creating a CTsmFRET object. The resulting object is used to 
        run smFRET dye placement on.

        ........................................
        OPTIONS 
        ........................................
        keyword [type] {default} 
        Description
        ........................................
        
        CTProteinObject [CTProtein-derived object]  
        CTProtein object extracted from a CTTraj objected. The CTProtein object is the main object
        that most protein-based analysis is performed over in CTraj

        dye1_name [string]
        Name of the first dye we're adding (note dye names are from a restricted vocabulary of 
        'A488' or 'A594', as these are the only two dyes supported right now.

        dye1_pos_w_offset [int]
        Residue position where the dye should be added including offset.

        dye2_name [string]
        Name of the first dye we're adding (note dye names are from a restricted vocabulary of 
        'A488' or 'A594', as these are the only two dyes supported right now.

        R0_values [key value dictionary] { {'A488-A594':5.900,'A594-A488':5.900} }
        Default R0 values to be used for dyes. Note that by default the values used are those
        provided by the Schuler group, which are SLIGHTLY different to those provided by the
        Haran group, which are {'A488-A594':5.637,'A594-A488':5.637}. 

        verbose [bool] {True}
        Print status messages during analysis
        

        """

        # NOTE: NO CHANGES SHOULD BE MADE TO THE CTPO by the object - this should be treated 
        # as a read only object by CTsmFRET (no such explicit control in Python)
        self.CTPO = CTProteinObject        

        if dye1_name in SUPPORTED_DYES:
            self.D1 = self.__load_dye_library(dye1_name)
        else:
            raise CTsmFRET_Exception('Unknown dye name used [%s]. Please ensure dyename is one of [%s]' % (dye1_name, SUPPORTED_DYES))

        if dye2_name in SUPPORTED_DYES:
            self.D2 = self.__load_dye_library(dye2_name)
        else:
            raise CTsmFRET_Exception('Unknown dye name used [%s]. Please ensure dyename is one of [%s]' % (dye2_name, SUPPORTED_DYES))


        ## Check R0 values correspond to the dye names
        try:
            self.__get_dye_dye_value(R0_VALUES, dye1_name, dye2_name)
        except KeyError:
            # if this happens the dyenames passed were fine but the R0 values couldn't be accessed
            raise CTsmFRET_Exception('The R0 value to be used (%s) could not be parsed. Recal this should be of the format {"dye1name-dye2name" : R0}. Note that the order of dye1 and dye2 does not matter' % (R0_VALUES))
            
            

        self.D1_idx = dye1_pos_w_offset
        self.D2_idx = dye2_pos_w_offset
        self.dye1_name = dye1_name
        self.dye2_name = dye2_name
        self.verbose = bool(verbose)
        
        # check R0_values make sense

        self.R0_VALS = R0_VALUES

        ## Verify the idx values provided make sense (specifically need a CA, CB and N)
        for atomname in ['N','CA','CB']:
            for resid in [self.D1_idx, self.D2_idx]:
                try:
                    self.CTPO.topology.select('resid %i and name %s' % (resid, atomname))[0]
                except IndexError:
                    raise CTsmFRET_Exception("Residue [%i (%s)] lacks atom [%s] - please try a different residue ID" % (resid, str(self.CTPO.topology.select('resid %i'%resid )), atomname))


        # numpy libraries will try and use as many cores as available - set to use 1 so we don't overclock the nodes!
        try:
            CTUtils.mkl_set_num_threads(MAXCORES)
            print("MKL will use %i threads" % (CTUtils.mkl_get_max_threads()))
        except Exception as e:
            print("MKL libraries not available [%s]" % str(e))
            pass

        # REMEMBER THIS IS IN A^6 NOT NM^6 
        self.R0_pow6_prefactor = self.__compute_R0_pow6_prefactor(self.dye1_name, self.dye2_name)
            

    def __message(self, msg):
        """
        Statau 

        """
        if self.verbose:
            print(msg)

                    
    # ........................................................................
    #                 
    def __get_dye_dye_value(self, relevant_dict, D1_name, D2_name):
        """
        Internal function that you can pass a pair of dynames to an a dye-dye dictionary (i.e.
        of format {'A488-A584': <VALUE OF INTEREST>} and it returns value of interest irrepsective
        of which dye is D1 or D2. Basically removes order dependence on the dyes for the user!

        Note this casts the return value to a float - so not appropriate for strings!!

        """

        try:
            return float(relevant_dict['%s-%s' % (D1_name, D2_name)])
        except KeyError:
            return float(relevant_dict['%s-%s' % (D2_name, D1_name)])
                

    # ........................................................................
    #                 
    def __load_dye_library(self, dyename):
        """
        Internal function that loads in the rotomer libraries from CTraj. If more rotomer libraries are
        needed this function should be updated
        
        """

        # Alexa488
        if dyename == 'A488':
            pdb  = "%s/data/dye_libraries/AF488.pdb" %(os.path.dirname(sys.modules['CTraj'].__file__))
            traj = "%s/data/dye_libraries/AF488_rotamers.dcd" %(os.path.dirname(sys.modules['CTraj'].__file__))
        elif dyename == 'A594':
            pdb  = "%s/data/dye_libraries/AF594m.pdb" %(os.path.dirname(sys.modules['CTraj'].__file__))
            traj = "%s/data/dye_libraries/AF594m_rotamers.dcd" %(os.path.dirname(sys.modules['CTraj'].__file__))
        else:
            raise CTsmFRET_Exception('Unexpected dye name passed... [%s]'%dyename)

        # read in the rotomer trajectory
        rotomer_traj = md.load(traj, top=pdb)

        # SANITY CHECKS BELOW...        
        if len(rotomer_traj.topology.select('name SG')) != 1:
            raise CTsmFRET_Exception('Did not find a single unique SG atom in the structure - required for dye alignment')

        return rotomer_traj


    # ........................................................................
    #                     
    def __compute_R0_pow6_prefactor(self, D1name, D2name):
        """
        Internal function that returns the R0^6 prefactor term, which when multiplied by kappa^2 and the 6th root
        is taken returns the R0 distance. This is based on equation 13.6 (page 369) from
        Principles of Fluorescent Spectroscopy.

        NOTE - this gives you a value in Angstroms^6, so to get R0 in nm need to divide by 10

        """


        Q = self.__get_dye_dye_value(QUANTUM_YEILD, D1name, D2name)
        J = self.__get_dye_dye_value(OVERLAP_INTEGRAL, D1name, D2name)
        
        return (8.79e-5)*np.power(REFRACTIVE_INDEX,-4)*Q*J
            


    # ........................................................................
    #                     
    def run_dye_sampling(self, ntrials=100, dye_solvation_shell=5.0, write_DyeTraj=False, framerange=[0,-1]):
        """

        Dye sampling approach: version 1.
        Returns a nested list of length equal to the number of protein conformations examined, where each 
        element is a list of dye-dye distances calculated for that conformation. 

        This method adds dyes and calculates the mean FRET efficiency in the following way:

        for each frame in the simulation 
           for ntrials per frame
             Dye 1 is placed at position 1 on the protein
             Dye 2 is placed at position 2 on the protein
             Dye-protein and dye-dye clashes are checked
             If all good then the dye-dye distance is calculated 

        This is equivalent to Kiersten's original method.

        ........................................
        OPTIONS 
        ........................................
        keyword [type] {default} 
        Description
        ........................................
        
        ntrials [int] {100}
        Number of trials to try and fit dyes associated with each protein conformation

        dye_solvation_shell [float] {5.0}
        Solvation shell around majority of the dye in Angstroms. NOTE that the initial 
        ring that connects the dye to protein has a more permissive solvation shell of 
        2.0  Angstroms.

        write_DyeTraj [bool] {False}
        Boolean flag which defines if we should write a PDB/XTC pair for EACH protein
        conformation with the dyes explicitly placed. Useful for visualization and 
        sanity checking but slows things down a lot so general not worth switching
        to on. Output files are generated in the directory the analysis is being 
        run from.

        framerange [list] [0,-1]
        Default is to use all frames. Otherwise, a  start and end frame must be 
        provided as a list/tuple of length 2. Allows for analysis of subregions
        or for testing things out.              
        
        """

        # input validation
        self.__validate_framerange(framerange)

        # set the dye solvation shell in nm - i.e. default is 0.5 nm (5 angstroms)
        dyeshell_in_nm = dye_solvation_shell * 0.1
        
        # define the frames to use
        frame_start = framerange[0]
        if framerange[1] == -1:
            frame_end = self.CTPO.n_frames
        else:            
            frame_end   = framerange[1]
        
        # lists for the the PROTEIN atom positions where dyes are going to be
        # placed
        D1_idx = []        
        D2_idx = []

        CP = self.CTPO
        
        # extract the atomic positions associated with the protein where we're going                 
        # to place the dye
        D1_idx.append(CP.topology.select('resid %i and name N'  % self.D1_idx)[0])
        D1_idx.append(CP.topology.select('resid %i and name CA' % self.D1_idx)[0])
        D1_idx.append(CP.topology.select('resid %i and name CB' % self.D1_idx)[0])

        D2_idx.append(CP.topology.select('resid %i and name N'  % self.D2_idx)[0])
        D2_idx.append(CP.topology.select('resid %i and name CA' % self.D2_idx)[0])
        D2_idx.append(CP.topology.select('resid %i and name CB' % self.D2_idx)[0])

        # get the position of the SG atom in all rotomers of the dye
        IDX_of_SG_from_dye1 = self.D1.topology.select('name SG')[0]
        D1_SG_positions     = np.transpose(self.D1.xyz,(1,0,2))[IDX_of_SG_from_dye1]

        IDX_of_SG_from_dye2 = self.D2.topology.select('name SG')[0]
        D2_SG_positions     = np.transpose(self.D2.xyz,(1,0,2))[IDX_of_SG_from_dye2]

        # get the index of the start and end of the 'dye' parts of the dyes
        D1_start  = self.D1.topology.select('resname %s' % (DYE_RESNAMES[self.dye1_name]))[0]
        D1_end    = self.D1.topology.select('resname %s' % (DYE_RESNAMES[self.dye1_name]))[-1]

        D2_start  = self.D2.topology.select('resname %s' % (DYE_RESNAMES[self.dye2_name]))[0]
        D2_end    = self.D2.topology.select('resname %s' % (DYE_RESNAMES[self.dye2_name]))[-1]

        # define the dye-based index values of the first ring(s) in each of dyes we're a bit more
        # permisisve with clashes associated with these rings
        D1_ring = np.array(FIRST_RING_IDX[self.dye1_name])
        D2_ring = np.array(FIRST_RING_IDX[self.dye2_name])

        # protein atoms to be compared against dyes for steric clashes (excludes the
        # host residue sidechain). This basically defines the list of all protein atoms
        # OTHER than the atoms from the host residue
        D1_protein_atom_IDs = self.__non_host_atomic_ids( self.D1_idx)
        D2_protein_atom_IDs = self.__non_host_atomic_ids( self.D2_idx)
                
        success_count_per_frame = np.zeros((self.CTPO.n_frames))
        dye_dye_distance_per_frame = []
        for frame in range(frame_start, frame_end):
            
            # D1 base residue positions
            D1_N  = CP.traj.xyz[frame][D1_idx[0]]
            D1_CA = CP.traj.xyz[frame][D1_idx[1]]
            D1_CB = CP.traj.xyz[frame][D1_idx[2]]

            # D2 base residue positions
            D2_N  = CP.traj.xyz[frame][D2_idx[0]]
            D2_CA = CP.traj.xyz[frame][D2_idx[1]]
            D2_CB = CP.traj.xyz[frame][D2_idx[2]]
            
            # set the empty dye-dye distance list for this frame
            dye_dye_distance = []
            topology_built = False
                        
            print("On frame %i [total success = %i]" % (frame, np.sum(success_count_per_frame)))

            for trial in range(0, ntrials):

                # set
                clash=False
            
                ## ----------------------------------------------------------------------------------------
                ## STAGE 1 - determine the position of the SG atom based on the position
                ##           of the N, CB and CA atoms in the host residues
                pSG_D1 = self.__get_posSG(D1_N, D1_CA, D1_CB)
                pSG_D2 = self.__get_posSG(D2_N, D2_CA, D2_CB)


                ## ----------------------------------------------------------------------------------------
                ## STAGE 2 - randomly add the dyes to the molecule, selecting a random dihedral and a random
                ##           rotomer
                transDyePos_D1 = self.__add_dye(D1_SG_positions, pSG_D1, D1_CB, self.D1, D1_start, D1_end)
                transDyePos_D2 = self.__add_dye(D2_SG_positions, pSG_D2, D2_CB, self.D2, D2_start, D2_end)


                # collect all the dye atom positions together (this is the 3D position of all the atoms
                # associated with the dye) assuming the first position in the dye comes after the SG atom
                # this means allDyePos is a 3D array of length num_of_atoms_in_D1 + number_of_atoms_in_D2
                allDyePos = np.vstack((transDyePos_D1, transDyePos_D2))

                
                # construct a list of atom indices associated with the first bulky ring in the dyes. We're going
                # to be less restrictive in terms of distance clashes with this ring. The offset is applied so the
                # ring_idx indicies apply directly to the positions defined in allDyePos
                D2_ring_offset = D2_ring + len(transDyePos_D1) 
                ring_idx = np.vstack((D1_ring, D2_ring_offset))

                # for now for EACH atom in the dye... (idx = general counter to keep on track where we are)
                idx = 0
                num_of_atoms_in_D1 = len(transDyePos_D1)
                for atom in allDyePos:
                                    
                    # if we're on an atom from dye 1 don't want to worry about clashes with the host residue
                    if idx < num_of_atoms_in_D1:
                        trans = np.transpose(CP.traj.xyz[frame][D1_protein_atom_IDs])                       
                    # if we're on an atom from dye 2 don't want to worry about clashes with the host residue
                    else:
                        trans = np.transpose(CP.traj.xyz[frame][D2_protein_atom_IDs])
                    
                    # distance between the dye atom and all (relevant) atoms in the protein
                    dist = np.sqrt(np.square(atom[0] - trans[0]) + np.square(atom[1] - trans[1]) + np.square(atom[2] - trans[2]))

                    # if the current dye atom is in first ring (of either dye) the use a 2 angstrom distance
                    # cutoff
                    if idx in ring_idx:
                        cutoff=0.2
                    # else use a 5 angstrom distance cutoff (i.e. dyes are well solvated)
                    else:
                        cutoff=dyeshell_in_nm

                    # if, based on these cutoffs, we clash then reject
                    if np.sum(dist < cutoff) > 0:
                        #print "DYE-PROTEIN clash on atom %i of trial %i [cutoff = %3.3f" %(idx, trial,cutoff)
                        clash=True
                        break

                    # if we're on a D1 atom then check for clashes with any D2 positions (the D2 check is implictly
                    # done because each comparison is reciprocal - i.e. if a D1 atom doesn't hit ANY D2 atom then
                    # no D2 atom hits that same D1 atom. 
                    if idx < num_of_atoms_in_D1:
                        trans = np.transpose(transDyePos_D2)
                        dist = np.sqrt(np.square(atom[0] - trans[0]) + np.square(atom[1] - trans[1]) + np.square(atom[2] - trans[2]))
                        if np.sum(dist < cutoff) > 0:
                            #print "DYE-DYE clash on atom %i of trial %i" %(idx, trial)
                            clash=True
                            break

                    # incremement atom idx
                    idx=idx+1

                # if we hit a clash then carry on, nothing to see, try again
                if clash:
                    continue
                
                # increment the success counter
                success_count_per_frame[frame] = success_count_per_frame[frame]+1

                # calculate distance between dye light sources - this is SUPER explicit here
                # just for clarity!
                x_D1 = transDyePos_D1[DYE_LIGHT_SOURCE_IDX[self.dye1_name]][0]
                y_D1 = transDyePos_D1[DYE_LIGHT_SOURCE_IDX[self.dye1_name]][1]
                z_D1 = transDyePos_D1[DYE_LIGHT_SOURCE_IDX[self.dye1_name]][2]

                x_D2 = transDyePos_D2[DYE_LIGHT_SOURCE_IDX[self.dye2_name]][0]
                y_D2 = transDyePos_D2[DYE_LIGHT_SOURCE_IDX[self.dye2_name]][1]
                z_D2 = transDyePos_D2[DYE_LIGHT_SOURCE_IDX[self.dye2_name]][2]
            
                # calculate dye-dye distance based on distance between light-source atom IDX
                dye_dye_distance.append(np.sqrt(np.square(x_D1 - x_D2) + np.square(y_D1 - y_D2) + np.square(z_D1 - z_D2)))
                
                # if we're generating DyeTraj files (trajectory files with the dyes placed)
                if write_DyeTraj:

                    # always need a topology (pdb) file, so if one has already been built just add the current dye information 
                    # to the growing DyeTraj vector
                    if topology_built:
                        prot_pos = CP.traj.xyz[frame]
                        all_pos  = np.vstack((prot_pos, allDyePos))
                    
                        DyeTraj.xyz = np.concatenate((DyeTraj.xyz,np.array([all_pos])),axis=0)
                        DyeTraj.time = np.zeros(DyeTraj.xyz.shape[0])

                    # if topology not yet build then build this
                    else:
                        topology_built=True
                        protein_top = CP.topology
                        D1_top      = self.D1.topology.subset(self.D1.topology.select('resname %s' % (DYE_RESNAMES[self.dye1_name])))
                        D2_top      = self.D2.topology.subset(self.D2.topology.select('resname %s' % (DYE_RESNAMES[self.dye2_name])))
                        
                        full_top = protein_top.join(D1_top)
                        full_top = full_top.join(D2_top)

                        prot_pos = CP.traj.xyz[frame]
                        all_pos  = np.vstack((prot_pos, allDyePos))
                        DyeTraj  = md.Trajectory(all_pos, full_top)
                        DyeTraj.time = np.zeros(DyeTraj.xyz.shape[0])
                                                            
            # end of frame
            dye_dye_distance_per_frame.append(dye_dye_distance)

            if write_DyeTraj:
                if len(dye_dye_distance) > 0:                    
                    DyeTraj.save_xtc('DyeTraj_frame_%i.xtc'%frame)
                    DyeTraj.save_pdb('DyeTraj_frame_%i.pdb'%frame)


        return dye_dye_distance_per_frame



    # ........................................................................
    #                 
    def run_dye_sampling_independent_placement(self, ntrials=100, dye_solvation_shell=5.0, write_DyeTraj=False, framerange=[0,-1]):
        """

        Dye sampling approach: version 2.
        Returns a nested list of length equal to the number of protein conformations examined, where each 
        element is a list of dye-dye distances calculated for that conformation. 

        This method adds dyes and calculates the mean FRET efficiency in the following way:

        for each frame in the simulation 
           for ntrials per frame
             Dye 1 is placed at position 1 on the protein
             Dye1 - protein clashes are checked and clashes are discarded
   
           for ntrials per frame              
             Dye 2 is placed at position 2 on the protein
             Dye 2 - protein clashes are checked

           for all dye pairs that did not clash with protein we check for
           dye-dye clash and if none then append distance

        Because the two dyes are now decoupled this leads to an ENORMOUS increase in the 
        acceptance probability of having two dyes without an increase to computational
        cost.

        ........................................
        OPTIONS 
        ........................................
        keyword [type] {default} 
        Description
        ........................................
        
        ntrials [int] {100}
        Number of trials to try and fit dyes associated with each protein conformation

        dye_solvation_shell [float] {5.0}
        Solvation shell around majority of the dye in Angstroms. NOTE that the initial 
        ring that connects the dye to protein has a more permissive solvation shell of 
        2.0  Angstroms.

        write_DyeTraj [bool] {False}
        Boolean flag which defines if we should write a PDB/XTC pair for EACH protein
        conformation with the dyes explicitly placed. Useful for visualization and 
        sanity checking but slows things down a lot so general not worth switching
        to on. Output files are generated in the directory the analysis is being 
        run from.

        framerange [list] [0,-1]
        Default is to use all frames. Otherwise, a  start and end frame must be 
        provided as a list/tuple of length 2. Allows for analysis of subregions
        or for testing things out.              
        
        """       

        # input validation
        self.__validate_framerange(framerange)
        
        # define the frames to use
        frame_start = framerange[0]
        if framerange[1] == -1:
            frame_end = self.CTPO.n_frames
        else:            
            frame_end   = framerange[1]
        

        # set the dye solvation shell in nm
        dyeshell_in_nm = dye_solvation_shell * 0.1
        
        D1_idx = []        
        D2_idx = []

        CP = self.CTPO
        
        # extract the positions                 
        D1_idx.append(CP.topology.select('resid %i and name N'  % self.D1_idx)[0])
        D1_idx.append(CP.topology.select('resid %i and name CA' % self.D1_idx)[0])
        D1_idx.append(CP.topology.select('resid %i and name CB' % self.D1_idx)[0])

        D2_idx.append(CP.topology.select('resid %i and name N'  % self.D2_idx)[0])
        D2_idx.append(CP.topology.select('resid %i and name CA' % self.D2_idx)[0])
        D2_idx.append(CP.topology.select('resid %i and name CB' % self.D2_idx)[0])

        # get the position of the SG atom in all rotomers of the dye
        IDX_of_SG_from_dye1 = self.D1.topology.select('name SG')[0]
        D1_SG_positions     = np.transpose(self.D1.xyz,(1,0,2))[IDX_of_SG_from_dye1]

        IDX_of_SG_from_dye2 = self.D2.topology.select('name SG')[0]
        D2_SG_positions     = np.transpose(self.D2.xyz,(1,0,2))[IDX_of_SG_from_dye2]

        # get the index of the start and end of the 'dye' parts of the dyes
        D1_start  = self.D1.topology.select('resname %s' % (DYE_RESNAMES[self.dye1_name]))[0]
        D1_end    = self.D1.topology.select('resname %s' % (DYE_RESNAMES[self.dye1_name]))[-1]

        D2_start  = self.D2.topology.select('resname %s' % (DYE_RESNAMES[self.dye2_name]))[0]
        D2_end    = self.D2.topology.select('resname %s' % (DYE_RESNAMES[self.dye2_name]))[-1]

        # define the dye-based index values of the first ring(s) in each of dyes
        D1_ring = np.array(FIRST_RING_IDX[self.dye1_name])
        D2_ring = np.array(FIRST_RING_IDX[self.dye2_name])

        # protein atoms to be compared against dyes for steric clashes (excludes the
        # host residue sidechain)
        D1_protein_atom_IDs = self.__non_host_atomic_ids( self.D1_idx)
        D2_protein_atom_IDs = self.__non_host_atomic_ids( self.D2_idx)
                
        success_count_per_frame = np.zeros((self.CTPO.n_frames))
        dye_dye_distance_per_frame = []
        for frame in range(frame_start, frame_end):
            
            # D1 base residue positions
            D1_N  = CP.traj.xyz[frame][D1_idx[0]]
            D1_CA = CP.traj.xyz[frame][D1_idx[1]]
            D1_CB = CP.traj.xyz[frame][D1_idx[2]]

            # D2 base residue positions
            D2_N  = CP.traj.xyz[frame][D2_idx[0]]
            D2_CA = CP.traj.xyz[frame][D2_idx[1]]
            D2_CB = CP.traj.xyz[frame][D2_idx[2]]
            
            # set the empty dye-dye distance list for this frame
            dye_dye_distance = []
            topology_built = False

            #
            good_D1 = []
            good_D2 = []
                        
            print("On frame %i [total success = %i]" % (frame, len(np.array(list(chain.from_iterable(dye_dye_distance_per_frame))))))

            ## ----------------------------------------------------------------------------------------
            ## STAGE 1 - place D1 and figure out any protein-dye clashes
            print("Stage 1")
            for trial in range(0, ntrials):                
                clash=False
            
                pSG_D1 = self.__get_posSG(D1_N, D1_CA, D1_CB)
                transDyePos_D1 = self.__add_dye(D1_SG_positions, pSG_D1, D1_CB, self.D1, D1_start, D1_end)

                idx = 0
                for atom in transDyePos_D1:
                    trans = np.transpose(CP.traj.xyz[frame][D1_protein_atom_IDs])  
                    dist = np.sqrt(np.square(atom[0] - trans[0]) + np.square(atom[1] - trans[1]) + np.square(atom[2] - trans[2]))

                    # define solvation shell ( 2 Angstrom for first ring in dye)
                    if idx in D1_ring:
                        cutoff=0.2
                    else:
                        cutoff=dyeshell_in_nm

                    if np.sum(dist < cutoff) > 0:
                        #print "DYE-PROTEIN clash on atom %i of trial %i [cutoff = %3.3f" %(idx, trial,cutoff)
                        clash=True
                        break
                    idx=idx+1

                # if we get here managed to place the dye
                if not clash:
                    good_D1.append(transDyePos_D1)


            ## ----------------------------------------------------------------------------------------
            ## STAGE 2 - place D2 and figure out any protein-dye clashes
            print("Stage 2")
            for trial in range(0, ntrials):                
                clash=False
            
                pSG_D2 = self.__get_posSG(D2_N, D2_CA, D2_CB)
                transDyePos_D2 = self.__add_dye(D2_SG_positions, pSG_D2, D2_CB, self.D2, D2_start, D2_end)

                idx = 0
                for atom in transDyePos_D2:
                    trans = np.transpose(CP.traj.xyz[frame][D2_protein_atom_IDs])  
                    dist = np.sqrt(np.square(atom[0] - trans[0]) + np.square(atom[1] - trans[1]) + np.square(atom[2] - trans[2]))

                    # define solvation shell ( 2 Angstrom for first ring in dye)
                    if idx in D2_ring:
                        cutoff=0.2
                    else:
                        cutoff=dyeshell_in_nm

                    if np.sum(dist < cutoff) > 0:
                        #print "DYE-PROTEIN clash on atom %i of trial %i [cutoff = %3.3f" %(idx, trial,cutoff)
                        clash=True
                        break
                    idx=idx+1

                # if we get here managed to place the dye
                if not clash:
                    good_D2.append(transDyePos_D2)


            ## ----------------------------------------------------------------------------------------
            ## STAGE 3 - check for dye-dye clashes and get distances where not observed
            print("Stage 3")
            topology_built = False
            for transDyePos_D1 in good_D1:
                for transDyePos_D2 in good_D2:
                        
                    clash=False
                    for atom in transDyePos_D1:
                        trans = np.transpose(transDyePos_D2)
                        dist = np.sqrt(np.square(atom[0] - trans[0]) + np.square(atom[1] - trans[1]) + np.square(atom[2] - trans[2]))
                            
                        if np.sum(dist < dyeshell_in_nm) > 0:
                            clash=True
                            break
                        
                    if not clash:
                        x_D1 = transDyePos_D1[DYE_LIGHT_SOURCE_IDX[self.dye1_name]][0]
                        y_D1 = transDyePos_D1[DYE_LIGHT_SOURCE_IDX[self.dye1_name]][1]
                        z_D1 = transDyePos_D1[DYE_LIGHT_SOURCE_IDX[self.dye1_name]][2]

                        x_D2 = transDyePos_D2[DYE_LIGHT_SOURCE_IDX[self.dye2_name]][0]
                        y_D2 = transDyePos_D2[DYE_LIGHT_SOURCE_IDX[self.dye2_name]][1]
                        z_D2 = transDyePos_D2[DYE_LIGHT_SOURCE_IDX[self.dye2_name]][2]
            
                        dye_dye_distance.append(np.sqrt(np.square(x_D1 - x_D2) + np.square(y_D1 - y_D2) + np.square(z_D1 - z_D2)))
                            
                        # if we're generating DyeTraj files...
                        if write_DyeTraj:
                            
                            # combine all dye positions
                            allDyePos = np.vstack((transDyePos_D1, transDyePos_D2))

                            if topology_built:                                
                                prot_pos = CP.traj.xyz[frame]
                                all_pos  = np.vstack((prot_pos, allDyePos))
                    
                                DyeTraj.xyz = np.concatenate((DyeTraj.xyz,np.array([all_pos])),axis=0)
                                DyeTraj.time = np.zeros(DyeTraj.xyz.shape[0])
                            
                            else:
                                topology_built=True
                                protein_top = CP.topology
                                D1_top      = self.D1.topology.subset(self.D1.topology.select('resname %s' % (DYE_RESNAMES[self.dye1_name])))
                                D2_top      = self.D2.topology.subset(self.D2.topology.select('resname %s' % (DYE_RESNAMES[self.dye2_name])))
                            
                                full_top = protein_top.join(D1_top)
                                full_top = full_top.join(D2_top)
                                
                                prot_pos = CP.traj.xyz[frame]
                                all_pos  = np.vstack((prot_pos, allDyePos))
                                DyeTraj  = md.Trajectory(all_pos, full_top)
                                DyeTraj.time = np.zeros(DyeTraj.xyz.shape[0])
                            
            # end of frame
            dye_dye_distance_per_frame.append(dye_dye_distance)
            if write_DyeTraj:
                if len(dye_dye_distance) > 0:                    
                    DyeTraj.save_xtc('DyeTraj_frame_%i.xtc'%frame)
                    DyeTraj[0].save_pdb('DyeTraj_frame_%i.pdb'%frame)


        return dye_dye_distance_per_frame



    # ........................................................................
    #                 
    def get_efficiencies_from_distances(self, dye_dye_distance_per_frame):
        """
        Function that takes the output from either

        run_dye_sampling_independent_placement
    
        OR
        
        run_dye_sampling

        And returns a single list where each element is 2 place tuple of frame index
        and FRET efficiency for EACH distance passed in. 
        

        ........................................
        OPTIONS 
        ........................................
        keyword [type] {default} 
        Description
        ........................................

        dye_dye_distance_per_frame [nested list]
        
        Formatted as the output from run_dye_sampling_independent_placement or
        run_dye_sampling - i.e. a list equal to the number of frames, where each
        element contains a number of dye-dye distances. Elements can (and often
        will) be empty.
        
        """

        R0 = self.__get_dye_dye_value(self.R0_VALS, self.dye1_name, self.dye2_name)
        idx = 0
        EFF=[]
        for i in dye_dye_distance_per_frame:
            if len(i) > 0:
                for dis in i:
                    EFF.append([idx, self.__get_efficiency_from_distance(dis, R0)])
            idx=idx+1
        return EFF
              

      
    # ........................................................................
    #                                                         
    def get_efficiencies(self, ntrials=100, dye_solvation_shell=5.0, write_DyeTraj=False, framerange=[0,-1]):
        """
        Function that simply combines run_dye_sampling_independent_placement and get_FRET_effiencies. Directly
        generates FRET efficiencies, returning a 2D vector of length Nx2 where N = number of dye distances, and
        for each dye distance the simulation frame is element 0 and the FRET efficiency is element 1.

        ........................................
        OPTIONS 
        ........................................
        keyword [type] {default} 
        Description
        ........................................
        
        ntrials [int] {100}
        Number of trials to try and fit dyes associated with each protein conformation

        dye_solvation_shell [float] {5.0}
        Solvation shell around majority of the dye in Angstroms. NOTE that the initial 
        ring that connects the dye to protein has a more permissive solvation shell of 
        2.0  Angstroms.

        write_DyeTraj [bool] {False}
        Boolean flag which defines if we should write a PDB/XTC pair for EACH protein
        conformation with the dyes explicitly placed. Useful for visualization and 
        sanity checking but slows things down a lot so general not worth switching
        to on. Output files are generated in the directory the analysis is being 
        run from.

        framerange [list] [0,-1]
        Default is to use all frames. Otherwise, a  start and end frame must be 
        provided as a list/tuple of length 2. Allows for analysis of subregions
        or for testing things out.  
               
        """
        DDDPF = self.run_dye_sampling_independent_placement(ntrials, dye_solvation_shell, write_DyeTraj, framerange)
        return self.get_efficiencies_from_distances(DDDPF)
   

    # ........................................................................
    #                
    def run_COCOFRET(self, ntrials=100, 
                     bailthresh=10, 
                     FRETerror_thresh=0.05, 
                     dye_solvation_shell=5.0, 
                     base_connection_shell=2.0,
                     min_dye_set = 200,
                     min_individual_dye=20, 
                     max_individual_dye=40, 
                     framerange=[0,-1], 
                     orientational=False, 
                     write_DyeTraj=False,                      
                     logfreq=-1, 
                     logfile=COCOFRET_LOG):
        """

        Dye sampling approach: version 5

        COCOFRET (COnformational COnvergence of FRET) is a more general algorithm for the reconstruction of FRET efficiencies.

        This method adds dyes and calculates the mean FRET efficiency in the following way:
        _________________________________________________________________________________________

        for each frame in the simulation 

           while not converged 
             
             for ntrials per frame              
               Dye 1 is placed at position 1 on the protein
               Dye1 - protein clashes are checked and clashes are discarded
   
             for ntrials per frame              
               Dye 2 is placed at position 2 on the protein
               Dye 2 - protein clashes are checked

             for all dye pairs that did not clash with protein we check for
             dye-dye clash and if none then calculate FRET effiency and dye-dye
             distance associated with that conformation. This means for each protein
             conformation we generate an ensemble of compatible dye configurations.

             compute mean FRET efficiency and standard error on FRET efficiency
             if standard error is below FRETerror_thresh then converged, else 
             not converged and repeat. Also increment a histogram of all observed 
             FRET efficiencies and dye-dye distances. 

           save ALL dye_dye_distances per frame, mean_FRET per frame, error per frame, and the 
           overall dye-dye and FRET effiency histograms.
        _________________________________________________________________________________________


        Different dye distances converge at different rates, so by ensuring a uniform convergence
        criterion, COCOFRET means we can be equally confident in different FRET effiencies from
        different protein conformations. 

        There are a few other criterion in this algorithm that deserve some comment:

        1) ntrials = number of trails for each individual dye placement. As this increases the speed
                     will dramatically drop, and often a convergent number of pairs can be generated on
                     the first trial so we recommend keeping this < 100.                    

        2) bailthresh = how many times we should try to add dyes (ntrials per dye) until giving up.
                        the interplay between ntrials and bailthreshold may depend on a situation 
                        by situation basis. In general we recommend ntrails=100 and bailthresh is
                        5 to 10. 

        3) FRETerror_thresh = error tollerance (std(values)/sqrt(number of values)). Suggest this remains at
                              0.005

        4) min_dye_set = minimum number of unique dye pairs. Generall 200 is a good number. Note COCOFRET will
                         check that this number is possible given the max_individual_dye and will adjust 
                         the max_individuall_dye in the case that an impossible combination is provided
        
        5) min_individual_dye = To ensure we get a true ensemble this is the minimum number of individual
                                conformations for each dye we require

        6) max_individual_dye = This value sets the maximum number of unique dye conformations for each of
                                the dyes needed. Recall the number of unique pairs is the product of the 
                                the number of conformations for dye1 and dye2. This max provides a method
                                to substantially reduce computational inefficiency by avoiding a scenario where
                                we slowly find conformations of dye1 but rapidly find conformations of dye2,
                                leading to a scenario where we have (say) 10 * 2000. The larger this value the
                                more accurate the derived values, but the more computationally intensive the
                                procedure. The COCOFRET algorithm will dynamically increase this value IF the
                                error is above the FRETerror_thresh and we reach this max.
             
        
        The return value from run_COCOFRET gives a 4-element tuple

        Element 0 is equivalent to the return from the other two run_dye_sampling functions - i.e. a list of N elements
        where N is the number of frames analyzed. The element at each frame contains the number of dye-dye distances 
        from the dye ensemble associated with that frame.
        
        Element 1 is a list where each element contains six sub-elements: [frame idx, mean FRET efficiency, error, 
        number of distances used, mean distance, standard deviation on distance distribution, acceptance fraction]
        e.g.

        [[0, 6.77e-02, 2.76-03, 1170,  9.29, 6.86-01, 0.0558],
         [1, 2.12-01,  4.87-03, 5190,  7.54, 8.33-01, 0.1201],
         [2, 7.21-01,  4.59-03, 16200, 4.86, 9.51-01, 0.2002],
         [3, 8.05-01,  4.52-03, 13300, 4.39, 9.81-01, 0.1854],
         [6, 4.23-01,  4.83-03, 13300, 6.28, 8.11-01, 0.2140],
         ...
         ]

        In this example, frames, 0, 1,2,3 and 6 generated reproducible converged results, while frames 
        4 and 5 did not form converged values. 

        Element 2 list of N elements where N is the number of frames analyzed. The element at each frame 
        contains the number of dye-dye transfer effiencies from the dye ensemble associated with that frame.
        This could be computed directly for a fixed-orientation analysis but when R0 depends on the dye-dye
        orientation these efficiency values capture that.

        Element 3 contains the a pair of bins/counts for the global effiency histogram and the global distance histogram between the dyes

        ........................................
        OPTIONS 
        ........................................
        keyword [type] {default} 
        Description
        ........................................
        
        ntrials [int] {100}
        Number of trials to try and fit dyes associated with each protein conformation (NB: number
        of traisl PER dye)

        bailthresh [int] {10}
        Number of attemps to make in case dyes cannot be fit - avoids infinite loop. If there are
        conformations which cannot be fit will make ntrais * bailthresh attemps, so consider that
        as a potential bottleneck for compact conformations.

        FRETerror_thresh [float] {0.05}:
        Standard error tollerance allowed from distribution of FRET efficiencies (not distances
        but efficiencies). 

        dye_solvation_shell [float] {5.0}
        Solvation shell around majority of the dye in Angstroms. NOTE that the initial 
        ring that connects the dye to protein has by default a more permissive solvation 
        shell of 2.0 Angstroms.
        
        base_connection_shell [float] {5.0}
        Solvation shell around the initial ring in the dye in Angstroms. This is set to the same
        value as the dye_solvation_shell to maintain backwards compatibility

        min_individual_dye [int] {20}
        Minimum number of dyes found for each dye. 5 is an absolute lowest minimum and the
        value should probably be more like 10 or 20.

        max_individual_dye [int] {40}
        Minimum number of dyes found for each dye. The default of 300 is means that
        in the most extreme cases - where adding one dye is easy and the other is challenging -
        the number of unique dye pairs will be ~1500 (5x300) based on default values. However,
        if the min_individual_dye is increased this value can also probably be decreased. Note
        that if convergence is challenging this value will be dynamically increased on a case-
        by-case automatically to ensure error convergence is equivalent for all analyses.

        write_DyeTraj [bool] {False}
        Boolean flag which defines if we should write a PDB/XTC pair for EACH protein
        conformation with the dyes explicitly placed. Useful for visualization and 
        sanity checking but slows things down a lot so general not worth switching
        to on. Output files are generated in the directory the analysis is being 
        run from.

        logfreq [int] {-1}
        Default is off. Frequency with which logfile information is written to disk.
        Useful when performing largescale analysis.        

        framerange [list] [0,-1]
        Default is to use all frames. Otherwise, a  start and end frame must be 
        provided as a list/tuple of length 2. Allows for analysis of subregions
        or for testing things out.  

        orientational [bool] {False}
        Flag which determines if a dye-dye orientation should be used. 

        logfile [string] {COCOFRET_LOG}
        Name of the output logfile (default is "COCOFRET_STATUS.log")
                
        """

        ## -----------------------------------------------------------------------------------------
        ## PARSE VARIOUS OPTIONS here
        if ntrials < 1:
            raise CTsmFRET_Exception('Set ntrials to a positive integer >= 1')

        if bailthresh < 1:
            raise CTsmFRET_Exception('Set bailthresh to a positive integer >= 1')

        if FRETerror_thresh <= 0:
            raise CTsmFRET_Exception('Set FRETError_thresh to a positive value > 0')

        if dye_solvation_shell < 0:
            raise CTsmFRET_Exception('Set dye_solvation_shell to a positive value > 0. Note 0 imples no excluded volume!')

        if base_connection_shell < 0:
            raise CTsmFRET_Exception('Set base_connection_shell to a positive value > 0. Note 0 imples no excluded volume!')                        

        if min_dye_set < 1:
            raise CTsmFRET_Exception('Set min_dye_set to a positive integer >= 1')

        if not len(framerange) == 2:
            raise CTsmFRET_Exception('framerange must be a tuple or list containing two elements definining the first and last frame over which the analysis will be performed') 
            
        if framerange[1] <= framerange[0] and not (framerange[1] == -1):
            raise CTsmFRET_Exception('The first element in framerange must be less than the second element in framerange')

        if min_individual_dye >= max_individual_dye:
            raise CTsmFRET_Exception('max_individual_dye must be greater than min_individual_dye (min = %i and max = %i)' % (min_individual_dye, max_individual_dye))

        # WARNING and auto-correct if necessary...
        if max_individual_dye*max_individual_dye < min_dye_set:
            print('WARNING: With the provided max_individual_dye value (%i) even if both dyes hit this target the dye_dye distance set would be less that min_dye_set (%i). Increasing max_indvidual_dye to %i' % (max_individual_dye, min_dye_set, int(np.sqrt(min_dye_set)+10)))
            max_individual_dye=int(np.sqrt(min_dye_set)+10)
            
        if logfreq > 0:
            self.__zero_COCOFRET_log(logfile)
            logout = True
        else:
            logout = False
            
        # define the simulation frames over which we're going to analyze
        frame_start = framerange[0]
        if framerange[1] == -1:
            frame_end = self.CTPO.n_frames
        else:            
            frame_end   = framerange[1]
                                        
        # set the dye solvation shell and base connection shell in nm
        dyeshell_in_nm = dye_solvation_shell * 0.1
        base_connection_in_nm = base_connection_shell * 0.1

        D1_idx = []        
        D2_idx = []

        CP = self.CTPO
        R0 = self.__get_dye_dye_value(self.R0_VALS, self.dye1_name, self.dye2_name)
        
        # extract the positions where we're going to add the dyes to                
        D1_idx.append(CP.topology.select('resid %i and name N'  % self.D1_idx)[0])
        D1_idx.append(CP.topology.select('resid %i and name CA' % self.D1_idx)[0])
        D1_idx.append(CP.topology.select('resid %i and name CB' % self.D1_idx)[0])

        D2_idx.append(CP.topology.select('resid %i and name N'  % self.D2_idx)[0])
        D2_idx.append(CP.topology.select('resid %i and name CA' % self.D2_idx)[0])
        D2_idx.append(CP.topology.select('resid %i and name CB' % self.D2_idx)[0])

        # get the position of the SG atom in all rotomers of the dye
        IDX_of_SG_from_dye1 = self.D1.topology.select('name SG')[0]
        D1_SG_positions     = np.transpose(self.D1.xyz,(1,0,2))[IDX_of_SG_from_dye1]

        IDX_of_SG_from_dye2 = self.D2.topology.select('name SG')[0]
        D2_SG_positions     = np.transpose(self.D2.xyz,(1,0,2))[IDX_of_SG_from_dye2]

        # get the index of the start and end of the 'dye' parts of the dyes - the dye rotomer libraries
        # include a connecting CYS residue, so this extracts out the first/last residue excluding anything
        # from that cysteine residue
        D1_start  = self.D1.topology.select('resname %s' % (DYE_RESNAMES[self.dye1_name]))[0]
        D1_end    = self.D1.topology.select('resname %s' % (DYE_RESNAMES[self.dye1_name]))[-1]

        D2_start  = self.D2.topology.select('resname %s' % (DYE_RESNAMES[self.dye2_name]))[0]
        D2_end    = self.D2.topology.select('resname %s' % (DYE_RESNAMES[self.dye2_name]))[-1]

        # define the dye-based index values of the first ring(s) in each of dyes
        D1_ring = np.array(FIRST_RING_IDX[self.dye1_name])
        D2_ring = np.array(FIRST_RING_IDX[self.dye2_name])

        # protein atoms to be compared against dyes for steric clashes (excludes the
        # host residue sidechain)
        D1_protein_atom_IDs = self.__non_host_atomic_ids( self.D1_idx)
        D2_protein_atom_IDs = self.__non_host_atomic_ids( self.D2_idx)
        
        success_count_per_frame = np.zeros((self.CTPO.n_frames))
        dye_dye_distance_per_frame = []
        dye_dye_transfer_efficiency_per_frame = []
        mean_and_error_on_FRET_per_frame = []

        # buld empty vectors that will store efficiency and distance histogram data
        efficiency_bins = np.arange(0,1.0025,0.0025)
        efficiency_histogram = np.zeros(len(efficiency_bins)-1)
        
        distance_bins = np.arange(0,((self.CTPO.get_numberOfResidues()*0.36)+1),0.01)
        distance_histogram = np.zeros(len(distance_bins)-1)
        
        for frame in range(frame_start, frame_end):
            print(frame)

            if logout and frame % logfreq == 0:
                self.__write_COCOFRET_status(frame_start, frame, logfreq, dye_dye_distance_per_frame, dye_dye_transfer_efficiency_per_frame, logfile)

            # D1 base residue positions
            D1_N  = CP.traj.xyz[frame][D1_idx[0]]
            D1_CA = CP.traj.xyz[frame][D1_idx[1]]
            D1_CB = CP.traj.xyz[frame][D1_idx[2]]

            # D2 base residue positions
            D2_N  = CP.traj.xyz[frame][D2_idx[0]]
            D2_CA = CP.traj.xyz[frame][D2_idx[1]]
            D2_CB = CP.traj.xyz[frame][D2_idx[2]]
            
            # set the empty dye-dye distance list for this frame
            dye_dye_distance = []
            dye_dye_transfer_efficiency = []
            topology_built = False
                        
            # set the exit flag (this is set to true if a solution cannot be found
            # or the FRET efficiency reaches convergence)
            convergence_reached = False

            # set the total number of dyepairs identified so far to 0. If we continously
            # fail to find new dye pairs then this triggers a failure 
            total_num_dyepairs_identified = 0

            # set the bailcount incrementer, which keeps track of how many iterations
            # we go through without finding a new pair of dyes
            bailcount = 0
            itercount = 0

            # variable the counts the number of times we try and place a dye
            total_add_count = 0
            total_good_add_count = 0
            
            # lists that hold the overall good sets of dye positions for dye 1 and dye 2
            overall_good_D1 = []
            overall_good_D2 = []
            
            while not convergence_reached:

                # list that hold this iterations set of good dye positions
                good_D1 = []
                good_D2 = []

                # incremenet the bailcounter
                bailcount=bailcount+1
                itercount=itercount+1


                ## ----------------------------------------------------------------------------------------
                ## STAGE 1 - place D1 and figure out any protein-dye clashes
                D1clashcount=0
                for trial in range(0, ntrials):                

                    # if we've already hit the max number required for this dye placement
                    # then just escape and move on to the next dye...
                    if len(overall_good_D1) >= max_individual_dye:
                        break

                    total_add_count=total_add_count+1
                    
                    # initialize the clash variable to false
                    clash=False
            
                    # first we set a randomly positioned of a gamma sulphur and make sure angles and bond
                    # lengths are ideal
                    pSG_D1 = self.__get_posSG(D1_N, D1_CA, D1_CB)

                    # then we add a dye using this gamma sulphur as a reference, by selecting a random rotomer
                    # and place at a rangdom angle off from that gamma sulphur. The returned values (transDyePos_D1)
                    # contain ONLY the positions of the dye, so the gamma suphur atom is never actually written
                    # written out to the trajectory file
                    transDyePos_D1 = self.__add_dye(D1_SG_positions, pSG_D1, D1_CB, self.D1, D1_start, D1_end)

                    idx = 0                    
                    for atom in transDyePos_D1:

                        trans = np.transpose(CP.traj.xyz[frame][D1_protein_atom_IDs])  
                        dist = np.sqrt(np.square(atom[0] - trans[0]) + np.square(atom[1] - trans[1]) + np.square(atom[2] - trans[2]))

                        # define solvation shell ( 2 Angstrom for first ring in dye)
                        if idx in D1_ring:
                            cutoff=base_connection_in_nm
                        else:
                            cutoff=dyeshell_in_nm

                        if np.sum(dist < cutoff) > 0:
                            #print "D1 DYE-PROTEIN clash on atom %i of trial %i [cutoff = %3.3f" %(idx, trial,cutoff)
                            D1clashcount=D1clashcount+1
                            clash=True                            
                            break
                        idx=idx+1

                    # if we get here managed to place the dye
                    if not clash:
                        good_D1.append(transDyePos_D1)
                
                ## ----------------------------------------------------------------------------------------
                ## STAGE 2 - place D2 and figure out any protein-dye clashes
                D2clashcount=0
                for trial in range(0, ntrials):   

                    if len(overall_good_D2) >= max_individual_dye:
                        break
             
                    total_add_count=total_add_count+1
                    clash=False                
                    pSG_D2 = self.__get_posSG(D2_N, D2_CA, D2_CB)
                    transDyePos_D2 = self.__add_dye(D2_SG_positions, pSG_D2, D2_CB, self.D2, D2_start, D2_end)

                    idx = 0
                    for atom in transDyePos_D2:

                        trans = np.transpose(CP.traj.xyz[frame][D2_protein_atom_IDs])  
                        dist = np.sqrt(np.square(atom[0] - trans[0]) + np.square(atom[1] - trans[1]) + np.square(atom[2] - trans[2]))

                        # define solvation shell ( 2 Angstrom for first ring in dye)
                        if idx in D2_ring:
                            cutoff=base_connection_in_nm
                        else:
                            cutoff=dyeshell_in_nm

                        if np.sum(dist < cutoff) > 0:
                            #print "D2 DYE-PROTEIN clash on atom %i of trial %i [cutoff = %3.3f" %(idx, trial,cutoff)
                            clash=True
                            D2clashcount=D2clashcount+1
                            break
                            idx=idx+1

                    # if we get here managed to place the dye
                    if not clash:
                        good_D2.append(transDyePos_D2)

                ## ----------------------------------------------------------------------------------------
                ## STAGE 3 - check for dye-dye clashes and get distances where not observed
                
                # set the total number of D1 and D2 positions currently found
                D1_total = len(good_D1)+len(overall_good_D1)
                D2_total = len(good_D2)+len(overall_good_D2)

                #print "D1 total = %i D2 total = %i" %(D1_total, D2_total)

                if D1_total > min_individual_dye and D2_total > min_individual_dye:
                    
                    # if this is the first time we gotta back-calculate everything we'd been saving up until now..
                    if len(dye_dye_distance) == 0:
                        (dye_dye_distance, dye_dye_transfer_efficiency) = self.__check_clash_and_compute_distances_and_effiencies(overall_good_D1, overall_good_D2, orientational, dyeshell_in_nm, R0)

                    # first calculate the distances/efficiencies for all the new dye pairs
                    [set_1_distances, set_1_efficiences] = self.__check_clash_and_compute_distances_and_effiencies(good_D1, good_D2, orientational, dyeshell_in_nm, R0)

                    # then compute all the old D1 positions with all the new D2 positions
                    [set_2_distances, set_2_efficiences] = self.__check_clash_and_compute_distances_and_effiencies(overall_good_D1, good_D2, orientational, dyeshell_in_nm, R0)

                    # finally compute all the old D2 positions with all the new D1 positions
                    [set_3_distances, set_3_efficiences] = self.__check_clash_and_compute_distances_and_effiencies(good_D1, overall_good_D2, orientational, dyeshell_in_nm, R0)

                    # having computed these extend the old values
                    dye_dye_distance.extend(set_1_distances)
                    dye_dye_distance.extend(set_2_distances)
                    dye_dye_distance.extend(set_3_distances)

                    dye_dye_transfer_efficiency.extend(set_1_efficiences)
                    dye_dye_transfer_efficiency.extend(set_2_efficiences)
                    dye_dye_transfer_efficiency.extend(set_3_efficiences)

                # regardless of if we computed distances or not, add those newly found dye positions to the overall good lists (note
                # 'good' refers to the fact we don't have dye-protein clashes, but says nothing about dye-dye clashes as that'll vary
                # on a pair to pair basis!)                
                overall_good_D1.extend(good_D1)
                overall_good_D2.extend(good_D2)

                # compute the maximum possible number of dye-dye pairs assuming NONE of the dyes 
                # clash (used to see if we're finding new dye conformations)
                max_possible_pairs = len(overall_good_D1)*len(overall_good_D2)          

                # if we haven't encountered a new dye-dye pair in $bailthres iterations
                # then its time to abort!
                if bailcount == bailthresh and max_possible_pairs == total_num_dyepairs_identified:
                    convergence_reached=True
                    self.__message("Frame %i: unable to reached convergence after %i iterations " % (frame, bailcount))
                    dye_dye_distance = [] # ensures we don't save non-converged distances

                # if we haven't found an addition pair carry on and don't reset the bailcount variable
                # (i.e. jump to the start of the loop again)
                if max_possible_pairs == total_num_dyepairs_identified:
                    continue
                             
                # if we have found one or more additional pairs then *do* reset bailcount 
                # and update the total_num_dyepairs_identified
                bailcount = 0
                total_num_dyepairs_identified = max_possible_pairs

                # finally, if we haven't actually computed any distances (yet) then continue
                # and try again (i.e. no mean FRET to calculate)
                if len(dye_dye_distance) == 0:
                    continue
                    
                # calculate mean fret efficienct based on all dye dye distances
                new_MFRET =  np.mean(dye_dye_transfer_efficiency)
                FRETerror = np.std(dye_dye_transfer_efficiency)/np.sqrt(len(dye_dye_transfer_efficiency))
                
                # evaluate if the FRETerror is low enough and we've found enough pairs and enough of each dye
                if FRETerror < FRETerror_thresh and len(dye_dye_distance) > min_dye_set and (len(overall_good_D1) > min_individual_dye) and (len(overall_good_D2) > min_individual_dye):
                    self.__message("Frame %i: reached convergence after %i iterations (mean = %3.4f (+/- %3.4f) on %i pairs [%i and %i]" % (frame, itercount, new_MFRET, FRETerror, len(dye_dye_distance), len(overall_good_D1), len(overall_good_D2)))
                    convergence_reached=True

                # The following double-if statement gives us a way to dynamically increase the upper limit on the number of each dye 
                # conformation we allow. This 
                # 
                # if we're not converged due to error on the FRET efficiency
                if FRETerror < FRETerror_thresh:

                    # and if one of our sets has maxed out...
                    if D1_total == max_individual_dye or D2_total == max_individual_dye:
                        # then increase
                        self.__message('Dynamically increasing the max_individual_dye from %i to %i' % (max_individual_dye, max_individual_dye+20))
                        max_individual_dye=max_individual_dye+20

                # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
                #  this is the end of the loop associated with a single conformation
                # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
 

            # WE GO HERE once convergence is reached
            # end of frame !!!!!
            dye_dye_distance_per_frame.append(dye_dye_distance)
            dye_dye_transfer_efficiency_per_frame.append(dye_dye_transfer_efficiency)
            
            # ------------------------------------
            # if we bailed out then no need to do anything else - we're not actually
            # saving the value, so skip to the next conformation
            if bailcount == bailthresh:
                continue

            # ELSE we now update all the info!
            # else save distance
            mean_and_error_on_FRET_per_frame.append([frame, new_MFRET, FRETerror, len(dye_dye_distance), np.mean(dye_dye_distance), np.std(dye_dye_distance), (len(overall_good_D1)+len(overall_good_D2))/float(total_add_count)])

            # update histograms
            efficiency_histogram = efficiency_histogram + np.histogram(dye_dye_transfer_efficiency, efficiency_bins)[0]
            distance_histogram = distance_histogram + np.histogram(dye_dye_distance, distance_bins)[0]

            if write_DyeTraj:
                self.__dyetraj_out(frame, overall_good_D1, overall_good_D2, DYE_RESNAMES, CP)

        return (dye_dye_distance_per_frame, mean_and_error_on_FRET_per_frame, dye_dye_transfer_efficiency_per_frame, [[efficiency_histogram, efficiency_bins[:-1]], [distance_histogram, distance_bins[:-1]]])


    def __check_clash_and_compute_distances_and_effiencies(self, position_set_1, position_set_2, orientational, dyeshell_in_nm, R0):
        """
        position_set_1 and position_set_2 are equal-length lists containing the positions of dyes. For the pairwise set of possible dyes in sets 1 and 2
        this function computes if (A) those dyes are overlapping and (B) if not computes the dye-dye distances AND the FRET efficiency in either the 
        standard of orientation-dependent manner.

        """

        local_dye_dye_distance          = []
        local_dye_dye_transfer_efficiency = []

        for transDyePos_D1 in position_set_1:
            for transDyePos_D2 in position_set_2:
                        
                clash=False

                # for each each atom in D1 see if we hit any atoms in D2
                for atom in transDyePos_D1:
                    trans = np.transpose(transDyePos_D2)
                    dist = np.sqrt(np.square(atom[0] - trans[0]) + np.square(atom[1] - trans[1]) + np.square(atom[2] - trans[2]))
                            
                    # as soon as we clash break and set Clash to true
                    if np.sum(dist < dyeshell_in_nm) > 0:
                        clash=True
                        break
                        
                # if no dye-dye clashes were observed...
                if not clash:

                    # then calculate interdye distance
                    x_D1 = transDyePos_D1[DYE_LIGHT_SOURCE_IDX[self.dye1_name]][0]
                    y_D1 = transDyePos_D1[DYE_LIGHT_SOURCE_IDX[self.dye1_name]][1]
                    z_D1 = transDyePos_D1[DYE_LIGHT_SOURCE_IDX[self.dye1_name]][2]
                        
                    x_D2 = transDyePos_D2[DYE_LIGHT_SOURCE_IDX[self.dye2_name]][0]
                    y_D2 = transDyePos_D2[DYE_LIGHT_SOURCE_IDX[self.dye2_name]][1]
                    z_D2 = transDyePos_D2[DYE_LIGHT_SOURCE_IDX[self.dye2_name]][2]
                            
                    # actually calculate the distance
                    local_distance = np.sqrt(np.square(x_D1 - x_D2) + np.square(y_D1 - y_D2) + np.square(z_D1 - z_D2))
                            
                    # add that distance to the growing dye_dye_distance vector
                    local_dye_dye_distance.append(local_distance)

                    # if we're calculating a dye-orientation dependent kappa factor...
                    if orientational:
                        (lk2, local_R0) = self.__compute_orientational_dependent_R0(transDyePos_D1, transDyePos_D2)
                        #print "%1.6f, %1.6f, %1.6f, %1.6f" % (self.__get_efficiency_from_distance(local_distance, local_R0), local_distance, local_R0, lk2)
                    else:
                        local_R0 = R0
                            
                    # and calculate the transfer efficiency using the R0 as determined by the orientational
                    # switch
                    local_dye_dye_transfer_efficiency.append(self.__get_efficiency_from_distance(local_distance, local_R0))

        return (local_dye_dye_distance, local_dye_dye_transfer_efficiency)
        
        


    def __dyetraj_out(self, frame, overall_good_D1, overall_good_D2, DYE_RESNAMES, CP):
        """
        Function that writes a set of dye positions to disk...

        """

        # could equally be if overall_good_D2 was 0...
        if len(overall_good_D1) == 0:                    
            return None

        # first build the topology using the 0th position from both lists
        allDyePos = np.vstack((overall_good_D1[0], overall_good_D2[0]))
        protein_top = CP.topology.__deepcopy__()
        D1_top      = self.D1.topology.subset(self.D1.topology.select('resname %s' % (DYE_RESNAMES[self.dye1_name])))
        D2_top      = self.D2.topology.subset(self.D2.topology.select('resname %s' % (DYE_RESNAMES[self.dye2_name])))

        ## SANITY CHECK: Useful to ensure index values used internally match up as expected
        ## Means the values printed to screens should match the values in the PDB file (by a value of 10x)
        ## - tested and validated during development at several different points!

        """
        print "D1 lightsource position: %s" %(overall_good_D1[0][DYE_LIGHT_SOURCE_IDX[self.dye1_name]])
        print "D1 transition dipole atoms: %s -- %s" % (overall_good_D1[0][DYE_TRANSITION_DIPOLE_IDX[self.dye1_name][0]],overall_good_D1[0][DYE_TRANSITION_DIPOLE_IDX[self.dye1_name][1]])
        print "D1 first ring positions:"
        for i in FIRST_RING_IDX[self.dye1_name]:
            print overall_good_D1[0][i]


        print "D2 lightsource position: %s" %(overall_good_D2[0][DYE_LIGHT_SOURCE_IDX[self.dye2_name]])
        print "D2 transition dipole atoms: %s -- %s" % (overall_good_D2[0][DYE_TRANSITION_DIPOLE_IDX[self.dye2_name][0]],overall_good_D2[0][DYE_TRANSITION_DIPOLE_IDX[self.dye2_name][1]])
        print ""
        print "D2 first ring positions:"
        for i in FIRST_RING_IDX[self.dye2_name]:
            print overall_good_D2[0][i]
        """


        #########################################3
        # OLD MECHANISM whereby the dyes were added as seperate chains. Left in case this becomes a desired behaviour under some circumstances...
        #full_top = protein_top.join(D1_top)
        #full_top = full_top.join(D2_top)
        #########################################3


        #########################################3
        ## NEW MECHANISM whereby dyes are added to the main chain
        # add dye1 to the protein topolog
        newres = protein_top.add_residue(D1_top.residue(0).name, protein_top.chain(0))
        for atom in D1_top.residue(0).atoms:
            protein_top.add_atom(atom.name, atom.element, newres)

        # add dye2
        newres = protein_top.add_residue(D2_top.residue(0).name, protein_top.chain(0))
        for atom in D2_top.residue(0).atoms:
            protein_top.add_atom(atom.name, atom.element, newres)

        full_top = protein_top
        #########################################3

        
        prot_pos = CP.traj.xyz[frame]
        all_pos  = np.vstack((prot_pos, allDyePos))
        DyeTraj  = md.Trajectory(all_pos, full_top)
        DyeTraj.time = np.zeros(DyeTraj.xyz.shape[0])
        
        # now using the position at index 1 from D1 write all the D2 positions 
        # out        
        transDyePos_D1 = overall_good_D1[1]
        for transDyePos_D2 in overall_good_D2[2:]:
            
            allDyePos = np.vstack((transDyePos_D1, transDyePos_D2))
            prot_pos = CP.traj.xyz[frame]
            all_pos  = np.vstack((prot_pos, allDyePos))
                    
            DyeTraj.xyz = np.concatenate((DyeTraj.xyz,np.array([all_pos])),axis=0)
            DyeTraj.time = np.zeros(DyeTraj.xyz.shape[0])

        # finally do the same using the dye position at index 1 for dye 2
        # and cyle over 2-end indx positions for dye 1
        transDyePos_D2 = overall_good_D2[1]
        for transDyePos_D1 in overall_good_D1[2:]:
            
            allDyePos = np.vstack((transDyePos_D1, transDyePos_D2))
            prot_pos = CP.traj.xyz[frame]
            all_pos  = np.vstack((prot_pos, allDyePos))
            
            # this is a fairly expensive operation...
            DyeTraj.xyz = np.concatenate((DyeTraj.xyz,np.array([all_pos])),axis=0)
            DyeTraj.time = np.zeros(DyeTraj.xyz.shape[0])
     
        # and finally write to disk..
        DyeTraj.save_xtc('DyeTraj_frame_%i.xtc'%frame)
        DyeTraj[0].save_pdb('DyeTraj_frame_%i.pdb'%frame)


    # ........................................................................
    #
    def __zero_COCOFRET_log(self, logfile):
        """
        Wipes a previous log file so all subsequent writes append to an existing file

        """
        with open(logfile, 'w') as fh:
            fh.write('')


    # ........................................................................
    #
    def __write_COCOFRET_status(self, frame_start, frame, logfreq, dye_dye_distance_per_frame, dye_dye_transfer_efficiency_per_frame, logfile):
        """
        Logging function for writing current COCOFRET status. Useful for non-interactive jobs
        so status is written to a local log file


        """
        if frame == frame_start:
            return 

        # set end (e) and start (s) indices
        e=frame - frame_start
        s=(frame-logfreq) - frame_start
        
        with open(logfile, 'a') as fh:            
            for i in range(s,e):
                if len(dye_dye_distance_per_frame[i]) ==0:
                    fh.write('%i, 0, -1, -1\n' % (i+frame_start))
                else:                            
                    fh.write('%i, %i, %4.4f, %4.4f\n' %( i+frame_start, len(dye_dye_distance_per_frame[i]), np.mean(dye_dye_distance_per_frame[i]), np.mean(dye_dye_transfer_efficiency_per_frame[i])))


        
        

    # ........................................................................
    #
    def __compute_kappa_squared(self, tDP_1, tDP_2):
        """
        Compute the dye-dependent configurational parameter kappa.

        Sanity checks done: the direction of D1_dm and D2_dm doesn't matter (as well
        it shouldn't, but nice to actually know this was double checked!)

        """
        # 3 unit vectors to compute...

        # r vector is the vector that goes from donor to acceptor (2 is acceptor, 1 is donor)
        r_vec    = np.array(tDP_2[DYE_LIGHT_SOURCE_IDX[self.dye2_name]] - tDP_1[DYE_LIGHT_SOURCE_IDX[self.dye1_name]])

        # D1_m is the vector along the transition dipole of the donor (note direction gets normalized out)
        D1_dm    = np.array(tDP_1[DYE_TRANSITION_DIPOLE_IDX[self.dye1_name][0]] - tDP_1[DYE_TRANSITION_DIPOLE_IDX[self.dye1_name][1]])

        # D2_m is the vector along the transition dipole of the acceptor (note direction gets normalized out)
        D2_dm    = np.array(tDP_2[DYE_TRANSITION_DIPOLE_IDX[self.dye2_name][0]] - tDP_2[DYE_TRANSITION_DIPOLE_IDX[self.dye2_name][1]])

        # convert all vectors into unit vectors
        r_1  = r_vec / np.linalg.norm(r_vec)
        D1_1   = D1_dm / np.linalg.norm(D1_dm)
        D2_1   = D2_dm / np.linalg.norm(D2_dm)

        # compute kappa 
        k = np.dot(D1_1, D2_1) - (3 * np.dot(r_1, D1_1) * np.dot(r_1, D2_1) )

        # return kappa squared
        return k*k
        


    def __compute_orientational_dependent_R0(self, tDP_1, tDP_2):
        """
        Calculates an orientational dependent R0 value using the dye transition dipoles to compute the
        orientational dependence factor. 
        
        Return value is in units of nm.


        """
        k2 = self.__compute_kappa_squared(tDP_1, tDP_2)
        
        """
        print k2
        if k2 > 4:
            raise Exception
        """


        # ******************************************************
        # note the 0.1 converts angstroms to nm - this is super
        # important!!!!!! Don't mess up units, kids!
        # ******************************************************
        LR0 = 0.1 * np.power(k2*self.R0_pow6_prefactor, 1/6.0)


        # note we return both of these here PURELY for simpler debugging!
        return (k2, LR0)
        
        

    # ........................................................................
    #
    def __get_posSG(self, N, CA, CB):
        """
        Set position of gamma sulphur so both distance and angles are minimized to give ideal bond 
        length and angles

        """

        """

        # old implementation that doesnt precompute anything, but kept around to more explicitly show the math required to get
        # the geometry right... Kiersten figured this out, because she's the MVP.

        def posfunct(x):            
            #
            #Function to minimize bond length and angle to match ideal values as closely as possible
            #after making some initial guess for the value. Function is defined here such that the various
            #values are initialized each time
            
            
            A = np.abs(np.cos(np.deg2rad(dih)) - np.dot(np.cross(CA-N, CB-CA) / np.linalg.norm(np.cross(CA - N, CB - CA)), np.cross(CB - CA, x - CB)/np.linalg.norm(np.cross(CB - CA, x- CB))))
            # cos(ccs_ang) - ([CA - CB] . [x - CB]) / |CA - CB| * |x - DB|
            B = np.abs(np.cos(np.deg2rad(CCS_ANG)) - np.dot(CA - CB, x - CB) / ( np.linalg.norm(CA - CB) * np.linalg.norm(x - CB)))
            C = np.abs(np.linalg.norm(x - CB) - CS_BOND)
            return A + B + C
            ## ----------------------------------------------------------------------------------------            
        """

        def posfunct_OPT(x):            
            """
            Function to minimize bond length and angle to match ideal values as closely as possible
            after making some initial guess for the value. Function is defined here such that the various
            values are initialized WITHIN the __get_posSG() function as a closure.
            
            """
            A = np.abs(A1 - np.dot(A2, np.cross(CB_MIN_CA, x - CB)/np.linalg.norm(np.cross(CB_MIN_CA, x- CB))))
            B = np.abs(B1 - np.dot(CA_MIN_CB, x - CB) / ( B2 * np.linalg.norm(x - CB)))
            C = np.abs(np.linalg.norm(x - CB) - CS_BOND)

            return A + B + C
            ## ----------------------------------------------------------------------------------------


        distTest = 0.5
        angTest  = 5


        A2 = np.cross(CA-N, CB-CA) / np.linalg.norm(np.cross(CA - N, CB - CA))
        B1 = np.cos(np.deg2rad(CCS_ANG))
        B2 = np.linalg.norm(CA - CB)
        CB_MIN_CA = CB - CA
        CA_MIN_CB = CA - CB
        

        while (distTest > 0.001) or (angTest > 0.01):                
            # initial random dihedral angle
            dih = -180.0 + 360*random.random()                                

            # randomly select an x/y/z value within MAXDIST from the CB atom
            x0 = CB + -MAXDIST+(2*MAXDIST)*np.random.random(3)                

            # optimize             
            A1 = np.cos(np.deg2rad(dih))

            # THIS LINE is where we spend ~90% of all COCOFRET's CPU time... If anyone can
            # find a way to optimize this even a little that would be SWELL
            """
            print "x = %s" % str(x0)
            print "A1 = %s" % str(A1)
            print "A2 = %s" % str(A2)
            print "B1 = %s" % str(B1)
            print "B2 = %s" % str(B2)
            print "CB_MIN_CA = %s" % str(CB_MIN_CA)
            print "CA_MIN_CB = %s" % str(CA_MIN_CB)
            """

            posSG = scipy.optimize.minimize(posfunct_OPT, [x0], method='Nelder-Mead').x


            distTest = abs(np.linalg.norm(posSG - CB) - CS_BOND);
            angTest  = abs(np.cos(np.deg2rad(CCS_ANG)) - np.dot(CA - CB, posSG - CB) / ( np.linalg.norm(CA - CB)*np.linalg.norm(posSG - CB)))

        return posSG


    # ........................................................................
    #
    def __get_mean_efficiency_from_distances(self, distance_list, R0):
        """
        Internal function that, given and R0 and a distance returns the instantaneous FRET
        efficiency. Note both dis and R0 should be in nm.

        """
        
        AVs = []
        for i in distance_list:            
            AVs.append(self.__get_efficiency_from_distance(i,R0))
        return (np.mean(AVs), np.std(AVs)/np.sqrt(len(AVs)))
                       

    # ........................................................................
    #
    def __get_efficiency_from_distance(self, dis, R0):
        """
        Internal function that, given and R0 and a distance returns the instantaneous FRET
        efficiency. Note both dis and R0 should be in nm

        """
        return 1/(1+np.power((float(dis)/R0),6.0))



    # ........................................................................
    #
    def __add_dye(self, dye_SG_positions, posSG, CB, dye_traj, dye_start, dye_end):
        """
        Add a dye from the rotomer library to a specific position. Returns the
        atomic positions of the dye starting with the position next excluding the 
        after the gamma sulphur.        
        
        """

        diffAng=40
        while diffAng>MAXDISTANG:

            # random rotomer
            rotomer_id = np.random.randint(0, len(dye_traj))
                
            # get the position of the gamma sulphur atom (which we pre-computed)
            SGPos  = dye_SG_positions[rotomer_id]

            # get the positions of all the atoms in the actual dye - NOTE this
            # 1) Includes the gamma sulphur atom (SGPos)
            # 2) We must offset the final residue so it's inclusive
            DyePos = np.vstack((SGPos, dye_traj.xyz[rotomer_id][dye_start:dye_end+1]))

            # randomly select a 3D rotation through which we're gonna rotate the dye
            # and build the associated rotation matrices
            angX = np.random.randint(0,360)
            angY = np.random.randint(0,360)
            angZ = np.random.randint(0,360)

            sinX = np.sin(np.deg2rad(angX))
            cosX = np.cos(np.deg2rad(angX))
            Rx   = np.array([[1,0,0], [0, cosX,-sinX],[0, sinX, cosX]])
                
            sinY = np.sin(np.deg2rad(angY))
            cosY = np.cos(np.deg2rad(angY))
            Ry   = np.array([[cosY, 0, sinY], [0,1,0], [-sinY, 0, cosY]])

            sinZ = np.sin(np.deg2rad(angZ))
            cosZ = np.cos(np.deg2rad(angZ))
            Rz   = np.array([[cosZ, -sinZ, 0], [sinZ, cosZ, 0], [0, 0, 1]])
                
            # this is Rx * Ry * Rz in terms of dot product 
            R = Rx.dot(Ry).dot(Rz)
            
            # Perform rotation of the dye
            rotDyePos = np.transpose(np.dot(R,DyePos.transpose()))

            # now determine and perform translation offset - note the 
            # transDyePos now doesn't have the SG atom associated with
            # it
            transDyeAmount = posSG - rotDyePos[0]
            transDyePos = rotDyePos[1:] + transDyeAmount

            # finally we're going to see how tight the angle is
            v1 = CB - posSG
            v2 = transDyePos[0] - posSG

            angv = 180.0/np.pi*(np.arccos(np.dot(v1,v2)/(np.linalg.norm(v1)*np.linalg.norm(v2))))

            diffAng = abs(angv - CSC_ANG)

        return transDyePos

            

    # ........................................................................
    #
    def __non_host_atomic_ids(self, resid):
        """
        Function which returns the atomic IDX of all atoms EXCEPT those in
        residue $resid
        
        """

        # get atoms in residue sidechain (backbone clashes still bad!)
        ignorelist = self.CTPO.topology.select('resid %i and sidechain' % resid)

        # return list of all IDXs between 0 and num of atoms EXCLUDING those
        return [i for i in range(0, self.CTPO.topology.n_atoms) if i not in ignorelist]



    #################
    ##
    ## Series of functions used to validate user input and catch edge-cases
    ##
    #################

    # ........................................................................
    #
    def __validate_framerange(self, framerange):
        """
        Validate the framerange passed is reasonable

        """
        
        if len(framerange) != 2:
            raise CTsmFRET_Exception('Framerange passed must be a list of 2 elements defining the first and last frame. Passed [%s]' %(str(framerange)))
            
        if framerange[0] < 0:
            raise CTsmFRET_Exception('Framerange first element less than 0 (first frame) [%s]- invalid' % (framerange[0]))
            
        if framerange[1] > self.CTPO.n_frames:
            raise CTsmFRET_Exception('Framerange second element is beyond the number of frames in the trajectory [%s] - invalid' %(framerange[1]))


        

        
                

                
                                                                                   
                            
                    
                    

        

    
