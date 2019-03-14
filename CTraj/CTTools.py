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

## CTTools contains random python functions which provide
## some useful functionality 
##
import numpy as np

def chunks(l, n):
    """Yield successive n-sized chunks from l."""

    maxval = len(l) - len(l)%n

    for i in range(0, maxval, n):
        yield l[i:i + n]

def fix_histadine_name(name):
    """
    Corrects the histadine residue name which can be HIE, HID 
    or HIP and unifies all of these to HIS
    """
    
    if name == 'HIE' or name == 'HID' or name == 'HIP':
        return 'HIS'
    else:
        return name


def find_nearest(array, target):
    """
    Find the value nearest to a target in an array, returns a tuple
    with the value and the tuple

    """
    idx = (np.abs(array-target)).argmin()
    return (idx,array[idx])


def powermodel(X, nu, R0):            
    """
    Function that computes the resulting power-model values

    """
    return R0*np.power(X,nu)            
