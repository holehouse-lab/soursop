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
## Copyright 2014 - 2010
##

import numpy
import ctypes
from . ctexceptions import CTException


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

def validate_keyword_option(keyword, allowed_vals, keyword_name, error_message=None):
    """
    Helper function that checks a passed keyword is only one of a set of possible
    valid keywords

    Parameters
    -----------
    keyword : str
        The actual passed keyword value

    allowed_vals : list of str
        A list of possible keywords

    keyword_name : str
        the name of the keyword as the user would select it in the function call

    error_message : str
        Allows the user to pass a custom error message


    Returns
    --------
    None

        No return value, but raises ctexceptions.CTException if ``keyword `` is not
        found in the allowed_vals list
           
    """


    if keyword not in allowed_vals:
        if error_message is None:
            raise CTException('Keyword %s passed value [%s], but this is not valid.\nMust be one of :%s'  % (keyword_name, keyword, ", ".join(allowed_vals)))

