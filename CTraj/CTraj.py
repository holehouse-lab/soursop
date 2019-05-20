"""
CTraj.py
Analysis package for all-atom simulations of proteins, with a specific focus on intrinsically disordered proteins.

Handles the primary functions
"""
from ._version import get_versions

def version(with_attribution=True):
    """
    Placeholder function to show example docstring (NumPy format)

    Replace this function and doc string for your own project

    Parameters
    ----------
    with_attribution : bool, Optional, default: True
        Set whether or not to display who the quote is from

    Returns
    -------
    quote : str
        Compiled string including quote and optional attribution
    """


    versions = get_versions()
    __version__ = versions['version']
    __git_revision__ = versions['full-revisionid']
    
    print(__version__)
    print(__git_revision__)


if __name__ == "__main__":
    # Do something if this file is invoked on its own
    print(canvas())
