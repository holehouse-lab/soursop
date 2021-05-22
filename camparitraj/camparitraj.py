"""
camparitraj module

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
## Copyright 2014 - 2021
##

from ._version import get_versions


def version():
    """
    Prints the full version and the git hash of the current commit

    Returns
    --------
    None 
        But prints version and git revision 
    """


    versions = get_versions()
    __version__ = versions['version']
    __git_revision__ = versions['full-revisionid']
    
    print(__version__)
    print(__git_revision__)


def version_full():
    """
    Returns the full camparitraj version as a string

    Returns
    --------
    str 
        camparitraj version
    """
    versions = get_versions()
    return versions['version']


def version_git_revision():
    """
    Returns the full git revision id as a string

    Returns
    --------
    str 
        camparitraj git revision id
    """

    versions = get_versions()
    return versions['full-revisionid']
    

if __name__ == "__main__":
    # Do something if this file is invoked on its own - e.g. print the version.
    print(version_full())
