##     _____  ____  _    _ _____   _____  ____  _____  
##   / ____|/ __ \| |  | |  __ \ / ____|/ __ \|  __ \ 
##  | (___ | |  | | |  | | |__) | (___ | |  | | |__) |
##   \___ \| |  | | |  | |  _  / \___ \| |  | |  ___/ 
##   ____) | |__| | |__| | | \ \ ____) | |__| | |     
##  |_____/ \____/ \____/|_|  \_\_____/ \____/|_|     

## Alex Holehouse (Pappu Lab and Holehouse Lab) and Jared Lalmansingh (Pappu lab)
## Simulation analysis package
## Copyright 2014 - 2022
##

"""
cttools contains misc. numerical python functions which provide some useful functionality 

"""


import numpy as np

# ........................................................................
#
def chunks(l, n):
    """Yield successive n-sized chunks from l.

    Parameters
    ----------

    l : list
        The list from which chunks of size `n` will be selected.

    n : int
        The size of the chunks to select from an input list, `l`.
    """


    maxval = int(round(len(l) - len(l)%n))

    for i in range(0, maxval, int(n)):
        yield l[i:i + n]


# ........................................................................
#
def fix_histadine_name(name):
    """
    Corrects the histadine residue name which can be HIE, HID 
    or HIP and unifies all of these to HIS.

    Parameters
    ----------

    name : str
        The input residue name (3-letter code) that may be transformed if it is Histidine-related.

    Returns
    -------
    The transformed name of the Histidine residue, otherwise the input residue is returned.
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
    with the value and the tuple.

    Parameters
    ----------

    array : np.array
        The array that will be searched for a target value.

    target: a value with dtype of `array`
        The search value to locate within `array`.

    Returns
    -------

    The first index and value of the target located within the input array.

    """
    idx = (np.abs(array-target)).argmin()
    return (idx,array[idx])


# ........................................................................
#
def powermodel(X, nu, R0):            
    """
    Function that computes the resulting power-model values.

    Parameters
    ----------

    X : array or numeric value
        An array or value whose power-model will be calculated.

    nu : int or float
        The exponent which the input array is raised to.

    R0: int or float
        A scaling constant.

    Returns
    -------
    array or numeric value

    """
    return R0*np.power(X,nu)            
