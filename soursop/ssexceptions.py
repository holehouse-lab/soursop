##     _____  ____  _    _ _____   _____  ____  _____
##   / ____|/ __ \| |  | |  __ \ / ____|/ __ \|  __ \
##  | (___ | |  | | |  | | |__) | (___ | |  | | |__) |
##   \___ \| |  | | |  | |  _  / \___ \| |  | |  ___/
##   ____) | |__| | |__| | | \ \ ____) | |__| | |
##  |_____/ \____/ \____/|_|  \_\_____/ \____/|_|

## Alex Holehouse (Pappu Lab and Holehouse Lab) and Jared Lalmansing (Pappu lab)
## Simulation analysis package
## Copyright 2014 - 2026
##


"""
Exception class where all the cool exception stuff happens. Empty exceptions that inherit from
the Exceptions

"""

import warnings


# ........................................................................
#
class SSException(Exception):
    """
    Exception class for raising custom exceptions

    """

    pass


# ........................................................................
#
class notYetImplementedException(Exception):
    """
    Exception for functionality not yet implemented

    """

    pass


# ........................................................................
#
class SoursopWarning(UserWarning):
    """
    Warning category used for every non-fatal advisory SOURSOP emits.

    It subclasses :class:`UserWarning`, so existing filters on
    ``UserWarning`` keep working, but it also lets you target SOURSOP's
    advisories specifically, e.g.
    ``warnings.filterwarnings('error', category=SoursopWarning)``.

    """

    pass


# ........................................................................
#
def SSWarning(string, stacklevel=2):
    """
    Emit a non-fatal SOURSOP warning.

    This is a thin wrapper around :func:`warnings.warn` that always uses the
    :class:`SoursopWarning` category. By default the warning is attributed
    to the SOURSOP function that called ``SSWarning`` (rather than to this
    helper), so a ``module='soursop'`` filter matches it.

    Parameters
    ----------
    string : str
        The warning message.

    stacklevel : int, optional
        Passed to :func:`warnings.warn`. Default is 2 (the caller of this
        function).

    Returns
    -------
    None

    """
    warnings.warn(string, SoursopWarning, stacklevel=stacklevel)
