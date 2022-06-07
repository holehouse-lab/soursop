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

import os

# Add imports here
from .soursop import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions


# code that allows access to the data directory
_ROOT = os.path.abspath(os.path.dirname(__file__))
def get_data(path):
    return os.path.join(_ROOT, 'data', path)


def get_version():
    return "%s - %s" % (str(__version__), str(__git_revision__))

from . import _version
__version__ = _version.get_versions()['version']
