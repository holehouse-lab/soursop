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
soursop module

"""

# The package version is written by versioningit into soursop/_version.py at
# build time (see [tool.versioningit.write] in pyproject.toml). Fall back to the
# installed-distribution metadata, then to "unknown", so that importing this
# module never fails in a source checkout where _version.py has not yet been
# generated (e.g. on Read the Docs).
try:
    from ._version import __version__
except ImportError:  # pragma: no cover - only hit before the build step runs
    try:
        from importlib.metadata import version as _dist_version

        __version__ = _dist_version("soursop")
    except Exception:
        __version__ = "unknown"


def version_git_revision():
    """Return the abbreviated git revision recorded in the version string.

    versioningit encodes the VCS revision in the PEP 440 local-version segment
    for any build made *after* a release tag (the configured format is
    ``{base_version}+{distance}.{vcs}{rev}``), e.g. ``2.0.0+5.g1a2b3c4`` or
    ``2.0.0+5.g1a2b3c4.dirty``. A clean, exactly-tagged release carries no local
    segment, so no revision is recorded.

    Returns
    -------
    str
        The git revision hash (without the leading ``g``), or ``"unknown"`` if
        the current build does not embed one (e.g. a tagged release).
    """
    if "+" not in __version__:
        return "unknown"

    local_segment = __version__.split("+", 1)[1]
    for token in local_segment.split("."):
        # versioningit prefixes the git revision with "g" (the VCS marker).
        if token.startswith("g") and len(token) > 1:
            return token[1:]
    return "unknown"


def version():
    """Print the full version and the git revision of the current build.

    Returns
    -------
    None
        Prints the version string and the git revision to stdout.
    """
    print(__version__)
    print(version_git_revision())


def version_full():
    """Return the full soursop version as a string.

    Returns
    -------
    str
        The soursop version string (as produced by versioningit).
    """
    return __version__


if __name__ == "__main__":
    # Do something if this file is invoked on its own - e.g. print the version.
    print(version_full())
