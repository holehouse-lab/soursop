##     _____  ____  _    _ _____   _____  ____  _____  
##   / ____|/ __ \| |  | |  __ \ / ____|/ __ \|  __ \ 
##  | (___ | |  | | |  | | |__) | (___ | |  | | |__) |
##   \___ \| |  | | |  | |  _  / \___ \| |  | |  ___/ 
##   ____) | |__| | |__| | | \ \ ____) | |__| | |     
##  |_____/ \____/ \____/|_|  \_\_____/ \____/|_|     

## Alex Holehouse (Pappu Lab and Holehouse Lab) and Jared Lalmansing (Pappu lab)
## Simulation analysis package
## Copyright 2014 - 2022
##


"""
camparitraj module

"""

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
