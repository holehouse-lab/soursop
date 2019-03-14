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
## Simulation analysis package
## Copyright 2014 - 2018
##

import numpy
import ctypes


##
## This is included the force numpy to use defined number of cores. For
## some of the linear algebra routines numpy will default to using as many
## cores as it can get its greedy little hands on - this function allows that 
## thirst to be quenched...
##

def mkl_set_num_threads(cores):
    mkl_rt = ctypes.CDLL('libmkl_rt.so')
    mkl_get_max_threads = mkl_rt.mkl_get_max_threads
    mkl_rt.mkl_set_num_threads(ctypes.byref(ctypes.c_int(cores)))


# usage
#mkl_set_num_threads(1)
#print mkl_get_max_threads() 
