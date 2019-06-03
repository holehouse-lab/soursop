"""
Exception class where all the cool exception stuff happens. Empty exceptions that inherit from 
the Exceptions

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
## Copyright 2014 - 2019
##

import warnings



# ........................................................................
#
class CTException(Exception):
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
def CTWarning(string):
    """
    Custom function to display non-fatal warnings.

    """
    warnings.warn(string)
    
