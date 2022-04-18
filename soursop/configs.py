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
configs contain all global configuration.

"""

import multiprocessing as mp
import platform
import tempfile


MAXCORES = mp.cpu_count()  # Added to resolve a reference with `sstrajectory.SSTrajectory.__init__`
DEBUGGING = False


# See: https://stackoverflow.com/questions/847850/cross-platform-way-of-getting-temp-directory-in-python
TMP_DIR = tempfile.gettempdir()
if platform.system().lower() == 'darwin':
    TMP_DIR = '/tmp'



