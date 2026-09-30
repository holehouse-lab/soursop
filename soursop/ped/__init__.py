"""
soursop.ped - access to the Protein Ensemble Database (PED).

PED (https://proteinensemble.org) is the main public repository of
conformational ensembles of intrinsically disordered proteins. This
sub-package wraps its REST API so that PED ensembles can be found and loaded
straight into SOURSOP.

Quick start::

    from soursop import ped

    # find ensembles
    ped.search("alpha-synuclein")             # -> ['PED00024', ...]
    ped.search(uniprot="P37840", ensembles=True)
                                              # -> ['PED00024e001', ...]

    # look at an entry
    entry = ped.get_entry("PED00001")
    print(entry.title, entry.identifiers)

    # load an ensemble as an SSTrajectory (one frame per conformer)
    traj = ped.load_ensemble("PED00001e001")
    rg = traj.proteinTrajectoryList[0].get_radius_of_gyration()

    # weighted ensembles
    w = ped.get_ensemble_weights("PED00001e001")   # None if unweighted

Identifiers are either an entry (``"PED00001"``) or an ensemble within an
entry (``"PED00001e001"``); every function that needs an ensemble accepts
either the full form or an entry plus ``ensemble_id``.

Every request goes over HTTPS to :data:`PED_API_URL`; failures raise
:class:`PEDError`, a subclass of :class:`~soursop.ssexceptions.SSException`.
"""

from ._http import DEFAULT_TIMEOUT, PED_API_URL, PEDError
from .ensemble import (
    ENSEMBLE_ASSETS,
    TABULAR_ASSETS,
    download_all_data,
    download_ensemble,
    get_ensemble_asset,
    get_ensemble_weights,
    load_ensemble,
)
from .entry import (
    PEDChainStats,
    PEDEnsembleSummary,
    PEDEntry,
    get_entry,
    list_ensembles,
)
from .identifiers import full_identifier, normalise_ensemble_id, parse_identifier
from .search import search, search_entries

__all__ = [
    "PED_API_URL",
    "DEFAULT_TIMEOUT",
    "PEDError",
    "search",
    "search_entries",
    "get_entry",
    "list_ensembles",
    "load_ensemble",
    "download_ensemble",
    "get_ensemble_weights",
    "get_ensemble_asset",
    "download_all_data",
    "ENSEMBLE_ASSETS",
    "TABULAR_ASSETS",
    "PEDEntry",
    "PEDEnsembleSummary",
    "PEDChainStats",
    "parse_identifier",
    "full_identifier",
    "normalise_ensemble_id",
]
