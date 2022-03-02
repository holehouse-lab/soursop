"""
ssio contains heloer functions for reporting information

"""

##
##
## Alex Holehouse (Pappu Lab and Holehouse Lab)
## Simulation analysis package
## Copyright 2014 - 2022
##


def warning_message(msg, with_frills=False):
    if with_frills:
        print("<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@")

    print("")
    print("WARNING: %s"%(msg))
    print("")
    if with_frills:
        print("<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@")

def exception_message(msg, exception, with_frills=False, raise_exception=True):

    if with_frills:
        print("\n<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@")

    print("")
    print("ERROR: %s"%(msg))
    print("")
    if with_frills:
        print("<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@<>@")

    if raise_exception:
        raise exception


def debug_message(msg):
    print("")
    print("DEBUG: %s"%(msg))
    print("")


def status_message(msg, verbose):
    if verbose:
        print("STATUS: %s"%(msg))

