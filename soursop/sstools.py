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

"""
sstools - miscellaneous numerical helper functions used across SOURSOP.

This module collects small, dependency-light utility routines (list
chunking, residue-name normalisation, nearest-value search, the polymer
power-law model, minimum-image distances and trajectory-file discovery)
that are shared by the higher-level analysis modules. None of these
functions hold state; they operate purely on their arguments.
"""

import numpy as np
from soursop.ssexceptions import SSException
from soursop import ssutils
from typing import Union, List
import os
from natsort import natsorted
import pathlib


# ........................................................................
#
def chunks(l, n):
    """Yield successive ``n``-sized chunks from a list.

    Splits ``l`` into consecutive, non-overlapping sublists of length
    ``n``. If ``len(l)`` is not an exact multiple of ``n`` the trailing
    remainder (fewer than ``n`` elements) is discarded, so every yielded
    chunk is guaranteed to contain exactly ``n`` elements. The list is
    consumed lazily as the generator is iterated.

    Parameters
    ----------
    l : list
        The list to be divided into chunks.

    n : int
        The size of each chunk. Must be a positive integer.

    Yields
    ------
    list
        Successive length-``n`` sublists of ``l``, in order.

    Example
    -------
    >>> list(chunks([1, 2, 3, 4, 5], 2))
    [[1, 2], [3, 4]]
    """

    maxval = int(round(len(l) - len(l) % n))

    for i in range(0, maxval, int(n)):
        yield l[i : i + n]


# ........................................................................
#
def fix_histadine_name(name):
    """Normalise protonation-specific histidine names to ``'HIS'``.

    Molecular dynamics force fields commonly encode the protonation state
    of histidine in the residue name (``HIE``, ``HID`` or ``HIP``). For
    sequence- and identity-based analysis these should all be treated as a
    single residue type. This function maps any of those three variants
    onto the canonical ``'HIS'`` and returns every other residue name
    unchanged.

    Parameters
    ----------
    name : str
        A three-letter residue name. May be a histidine protonation
        variant (``'HIE'``, ``'HID'``, ``'HIP'``) or any other code.

    Returns
    -------
    str
        ``'HIS'`` if ``name`` is a histidine protonation variant,
        otherwise ``name`` returned unmodified.

    Example
    -------
    >>> fix_histadine_name('HIE')
    'HIS'
    >>> fix_histadine_name('ALA')
    'ALA'
    """

    if name == "HIE" or name == "HID" or name == "HIP":
        return "HIS"
    else:
        return name


# ........................................................................
#
def find_nearest(array, target):
    """Find the array element closest in value to a target.

    Performs a brute-force search for the element of ``array`` whose
    absolute difference from ``target`` is smallest, returning both the
    index of that element and the element itself. If several elements are
    equidistant from ``target`` the first (lowest-index) match is
    returned.

    Parameters
    ----------
    array : np.ndarray
        Numeric array to be searched. Must support element-wise
        subtraction with ``target``.

    target : scalar
        The value to locate. Should share a comparable dtype with
        ``array``.

    Returns
    -------
    tuple of (int, scalar)
        A 2-tuple ``(index, value)`` where ``index`` is the position of
        the nearest element within ``array`` and ``value`` is that
        element.

    Example
    -------
    >>> import numpy as np
    >>> find_nearest(np.array([0.0, 1.0, 2.0, 3.0]), 2.2)
    (2, 2.0)
    """
    idx = (np.abs(array - target)).argmin()
    return (idx, array[idx])


# ........................................................................
#
def powermodel(X, nu, R0):
    """Evaluate the polymer-scaling power law ``R0 * X**nu``.

    Implements the standard polymer-physics scaling relationship used
    throughout SOURSOP - for example relating an inter-residue distance to
    sequence separation via ``R = R0 * |i - j|**nu``. This is a thin
    wrapper around ``numpy.power`` and broadcasts over array input.

    Parameters
    ----------
    X : np.ndarray or float
        Independent variable(s), typically a sequence separation. May be a
        scalar or any array-like accepted by ``numpy.power``.

    nu : int or float
        The scaling exponent (e.g. the apparent Flory scaling exponent).

    R0 : int or float
        The scaling prefactor.

    Returns
    -------
    np.ndarray or float
        ``R0 * X**nu``, with the same shape as ``X``.

    Example
    -------
    >>> powermodel(4, 0.5, 5.5)
    11.0
    """
    return R0 * np.power(X, nu)


# ------------------------------------------------------------------
#
def get_distance_periodic(distance1, distance2, box_size, box_shape="cube"):
    """Pairwise point distances under the minimum-image convention.

    Computes, for two equal-length arrays of 3D coordinates, the distance
    between each corresponding pair of points using the minimum-image
    convention for a cubic periodic box. This corrects distances for
    particles that are nearest across a periodic boundary. The
    implementation is a straightforward per-pair Python loop and is not
    performance-optimised; for a non-periodic system the result is
    identical to the naive Euclidean distance.

    Parameters
    ----------
    distance1 : np.ndarray
        Array of length ``n`` in which each element is the x/y/z position
        of a particle / atom / bead.

    distance2 : np.ndarray
        Array of length ``n`` in which each element is the x/y/z position
        of a particle / atom / bead. Must be the same length as
        ``distance1``.

    box_size : int or float
        Edge length of the cubic periodic box (in the same units as the
        coordinates).

    box_shape : {'cube'}, optional
        Geometry of the periodic cell used for the minimum-image
        convention. Currently only ``'cube'`` is supported. Default
        ``'cube'``.

    Returns
    -------
    list of float
        Inter-position distances (length ``n``) under the minimum-image
        convention.

    Raises
    ------
    soursop.ssexceptions.SSException
        If ``distance1`` and ``distance2`` differ in length, or if an
        unsupported ``box_shape`` is passed.

    Example
    -------
    >>> import numpy as np
    >>> a = np.array([[0.1, 0.1, 0.1]])
    >>> b = np.array([[9.9, 9.9, 9.9]])
    >>> get_distance_periodic(a, b, 10.0)        # nearest across the boundary
    [0.34641016151377546]
    """

    if len(distance1) != len(distance2):
        raise SSException(
            "The two distance vectors in get_distance_periodic() must be the same length"
        )

    ssutils.validate_keyword_option(box_shape, ["cube"], "box_shape")

    if len(distance1) != len(distance2):
        raise SSException(
            "The two distance vectors in get_distance_periodic() must be the same length"
        )

    # should rewrite in cython at some point...
    distances = []
    for idx in range(len(distance1)):
        c1 = distance1[idx]
        c2 = distance2[idx]

        delta = np.abs(c1 - c2)

        # apply minimum image convention
        delta = np.where(delta > 0.5 * box_size, box_size - delta, delta)

        distances.append(np.linalg.norm(delta))

    return distances


# ------------------------------------------------------------------
#
def find_trajectory_files(
    root_dir: Union[str, pathlib.Path],
    num_replicates: int,
    exclude_dirs: Union[None, List] = ["eq", "FULL"],
    traj_name: str = "__traj.xtc",
    top_name: str = "__START.pdb",
):
    """Discover matching trajectory/topology file pairs in a directory tree.

    Walks ``root_dir`` recursively and collects paired trajectory and
    topology files for a set of replicate simulations. A directory is
    treated as a replicate when its name ends in a digit in the range
    ``[0, num_replicates]`` and it contains both a ``traj_name`` and a
    ``top_name`` file. Results are grouped by parent-directory name and
    returned in natural-sorted order so replicate ordering is stable and
    human-intuitive (e.g. ``rep2`` before ``rep10``).

    Parameters
    ----------
    root_dir : str or pathlib.Path
        Root of the directory tree to search.

    num_replicates : int
        Highest replicate index to gather; child directories whose name
        ends in a digit in ``[0, num_replicates]`` are considered.

    exclude_dirs : list or None, optional
        Directory names to prune from the walk. Pass ``None`` to search
        everything (including ``"eq"`` and ``"FULL"``). Default
        ``["eq", "FULL"]``.

    traj_name : str, optional
        Trajectory filename to look for in each replicate directory.
        Default ``"__traj.xtc"``.

    top_name : str, optional
        Topology filename to look for in each replicate directory.
        Default ``"__START.pdb"``.

    Returns
    -------
    tuple of (list of str, list of str)
        A 2-tuple ``(traj_files, start_files)`` of equal length, where
        ``traj_files[i]`` and ``start_files[i]`` are the absolute paths to
        a matched trajectory and its topology, ordered by natural-sorted
        parent-directory name then natural-sorted file path.

    Example
    -------
    >>> trajs, tops = find_trajectory_files('/data/sim_set', num_replicates=5)
    >>> len(trajs) == len(tops)
    True
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
