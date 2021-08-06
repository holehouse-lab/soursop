"""
configs contain all global configuration.

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
## Copyright 2014 - 2021
##
import multiprocessing as mp
import platform
import tempfile


MAXCORES = mp.cpu_count()  # Added to resolve a reference with `sstrajectory.SSTrajectory.__init__`
DEBUGGING = False


# See: https://stackoverflow.com/questions/847850/cross-platform-way-of-getting-temp-directory-in-python
TMP_DIR = tempfile.gettempdir()
if platform.system().lower() == 'darwin':
    TMP_DIR = '/tmp'



