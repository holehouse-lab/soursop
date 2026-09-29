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

import os

# code that allows access to the data directory
_ROOT = os.path.abspath(os.path.dirname(__file__))

# Generate _version.py if missing and in the Read the Docs environment. The
# path is resolved relative to this file (it used to depend on the current
# working directory), and a source checkout without _version.py falls back to
# the version recorded by soursop.soursop rather than failing to import.
if os.getenv("READTHEDOCS") == "True" and not os.path.isfile(
    os.path.join(_ROOT, "_version.py")
):
    import versioningit

    __version__ = versioningit.get_version(os.path.dirname(_ROOT))
else:
    try:
        from soursop._version import __version__
    except ImportError:  # pragma: no cover - only before the build step runs
        from soursop.soursop import __version__

# The git revision is derived from the version string written by versioningit
# (see soursop.soursop.version_git_revision). Guarded so a partially-built or
# source checkout still imports cleanly.
try:
    from soursop.soursop import version_git_revision as _version_git_revision

    __git_revision__ = _version_git_revision()
except Exception:  # pragma: no cover - defensive fallback
    __git_revision__ = "unknown"


def get_data(path):
    return os.path.join(_ROOT, "data", path)


def get_version():
    return "%s - %s" % (str(__version__), str(__git_revision__))
