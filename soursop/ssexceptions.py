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
def SSWarning(string):
    """
    Custom function to display non-fatal warnings.

    """
    warnings.warn(string)
    
