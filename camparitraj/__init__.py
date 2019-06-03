"""
camparitraj __ini__ file
A short description of the project.
"""
import os

# Add imports here
from .camparitraj import *

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
    return "%s - %s" % (str(__version__), str(_git_revision))
