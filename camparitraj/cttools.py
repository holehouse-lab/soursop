"""
cttools contains misc. numerical python functions which provide some useful functionality 

"""
##
##                                       _ _              _ 
##   ___ __ _ _ __ ___  _ __   __ _ _ __(_) |_ _ __ __ _ (_)
##  / __/ _` | '_ ` _ \| '_ \ / _` | '__| | __| '__/ _` || |
## | (_| (_| | | | | | | |_) | (_| | |  | | |_| | | (_| || |
##  \___\__,_|_| |_| |_| .__/ \__,_|_|  |_|\__|_|  \__,_|/ |
##                     |_|                             |__/ 
##
## Alex Holehouse (Pappu Lab and Holehouse Lab)
## Simulation analysis package
## Copyright 2014 - 2019
##

import numpy as np

# ........................................................................
#
def chunks(l, n):
    """Yield successive n-sized chunks from l."""


    maxval = int(round(len(l) - len(l)%n))

    for i in range(0, maxval, n):
        yield l[i:i + n]


# ........................................................................
#
def fix_histadine_name(name):
    """
    Corrects the histadine residue name which can be HIE, HID 
    or HIP and unifies all of these to HIS
    """
    
    if name == 'HIE' or name == 'HID' or name == 'HIP':
        return 'HIS'
    else:
        return name


# ........................................................................
#
def find_nearest(array, target):
    """
    Find the value nearest to a target in an array, returns a tuple
    with the value and the tuple

    """
    idx = (np.abs(array-target)).argmin()
    return (idx,array[idx])


# ........................................................................
#
def powermodel(X, nu, R0):            
    """
    Function that computes the resulting power-model values

    """
    return R0*np.power(X,nu)            
