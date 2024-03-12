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
from typing import Union, List
import os
from natsort import natsorted
import pathlib

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

def find_trajectory_files(
    root_dir: Union[str, pathlib.Path],
    num_replicates: int,
    exclude_dirs: Union[None, List] = ["eq", "FULL"],
    traj_name: str = "__traj.xtc",
    top_name: str = "__START.pdb",
):
    """
    This function assembles the list of trajectory and topology paths
    for a set of simulations in a given directory tree.

    Parameters
    ----------
    root_dir : Union[str, pathlib.path]
        Filepath or list of file paths
    num_replicates : int, optional
        Number of replicas to gather trajectories from child directories [1,num_replicates+1], by default None
    exclude_dirs : Union[None,list], optional
        List of directory names you want to exclude, by default ["eq", "FULL"]
        Set to None to include "eq" and "FULL" or specify list.
    traj_name : str, optional
        trajectory filename, by default "__traj.xtc"
    top_name : str, optional
        topology filename, by default "__START.pdb"
    """
    if exclude_dirs is None:
        exclude_dirs = []

    traj_files = []
    start_files = []
    dir_dict = {}
    for dirpath, dirnames, filenames in os.walk(os.path.abspath(root_dir)):
        if any(d in exclude_dirs for d in dirnames):
            # exclude directories in exclude_dirs
            dirnames[:] = [d for d in dirnames if d not in exclude_dirs]
        if dirpath.endswith(tuple(str(i) for i in range(0, num_replicates + 1))):
            # extract the parent directory name
            parent_dirname = os.path.basename(os.path.dirname(dirpath))
            traj_file = os.path.join(dirpath, f"{traj_name}")
            start_file = os.path.join(dirpath, f"{top_name}")
            if os.path.isfile(traj_file) and os.path.isfile(start_file):
                if parent_dirname not in dir_dict:
                    # create a new list for the current parent directory name
                    dir_dict[parent_dirname] = []
                dir_dict[parent_dirname].append((traj_file, start_file))

    for parent_dirname in natsorted(dir_dict.keys()):
        # sort the list of files for the current parent directory name
        sorted_files = natsorted(dir_dict[parent_dirname])
        for traj_file, start_file in sorted_files:
            traj_files.append(traj_file)
            start_files.append(start_file)
    return traj_files, start_files

