"""
CTraj
Analysis package for all-atom simulations of proteins, with a specific focus on intrinsically disordered proteins.
"""

# Add imports here
from .CTraj import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
