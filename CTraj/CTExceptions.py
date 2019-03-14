##
################################################
##  ,-----.,--------.                 ,--.    ##
## '  .--./'--.  .--',--.--. ,--,--.  `--'    ##
## |  |       |  |   |  .--'' ,-.  |  ,--.    ##
## '  '--'\   |  |   |  |   \ '-'  |  |  |    ##
##  `-----'   `--'   `--'    `--`--'.-'  /    ##
##                                  '---'     ##
################################################
##
## Alex Holehouse (Pappu Lab)
## Simulation analysis package
## Copyright 2014 - 2018
##
##

import warnings

class CTException(Exception):
    """
    Exception class for raising custom exceptions

    """
    pass

class notYetImplementedException(Exception):
    """
    Exception for functionality not yet implemented

    """
    pass


class CTsmFRET_Exception(Exception):
    """
    Exception class for raising exceptions associated with
    the smFRET analysis

    """
    pass


def CTWarning(string):
    """
    Custom function to display non-fatal warnings.

    """
    warnings.warn(string)
    
