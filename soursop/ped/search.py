"""
Searching PED.

:func:`search` returns identifiers and :func:`search_entries` returns full
:class:`~soursop.ped.entry.PEDEntry` records. Both take the same criteria,
which map onto the query parameters of PED's ``/entries`` endpoint, and both
page through the results automatically.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

from ..ssexceptions import SSException
from . import _http
from .entry import PEDEntry

#: Our keyword -> PED's query parameter.
_CRITERIA = {
    "query": "free_text",
    "protein_name": "protein_name",
    "uniprot": "uniprot_acc",
    "term": "term",
    "author": "data_owner",
    "publication": "publication_identifier",
    "publication_title": "publication_html",
    "cross_ref": "cross_ref",
}

#: Entries requested per page when paging through search results.
DEFAULT_PAGE_SIZE = 100


def search_entries(
    query: Optional[str] = None,
    *,
    protein_name: Optional[str] = None,
    uniprot: Optional[str] = None,
    term: Optional[str] = None,
    author: Optional[str] = None,
    publication: Optional[str] = None,
    publication_title: Optional[str] = None,
    cross_ref: Optional[str] = None,
    max_results: Optional[int] = None,
    page_size: int = DEFAULT_PAGE_SIZE,
    timeout: float = _http.DEFAULT_TIMEOUT,
) -> List[PEDEntry]:
    """Search PED and return the matching entries.

    Every criterion given must match (they are passed to PED together). With
    no criteria at all, every entry in PED is returned - over 25,000 entries
    (mostly AlphaFlex and MD ensembles), which takes several minutes; use
    ``max_results`` or a criterion to narrow it down.

    Parameters
    ----------
    query : str, optional
        Free-text search across the entry (PED's ``free_text``), e.g.
        ``"alpha-synuclein"`` or ``"SAXS"``.
    protein_name : str, optional
        Protein name, e.g. ``"Sic1"``.
    uniprot : str, optional
        UniProt accession, e.g. ``"P37840"``.
    term : str, optional
        IDPO ontology term (a measurement or ensemble-generation method),
        e.g. ``"NMR"``.
    author : str, optional
        PED's "data owner" field (the depositor), which is not the same as
        the publication's author list: searching for a paper's author often
        finds nothing. To find entries by author, use ``query`` instead.
    publication : str, optional
        PubMed ID or DOI of the associated publication.
    publication_title : str, optional
        Publication title.
    cross_ref : str, optional
        Database or experimental cross-reference, e.g. a DisProt or BMRB
        identifier. At the time of writing PED answers this filter with a
        server error (HTTP 500); a free-text ``query`` for the identifier
        (e.g. ``query="DP00631"``) works.
    max_results : int, optional
        Stop after this many entries. Default ``None`` (all matches).
    page_size : int, optional
        Entries requested per call to PED. Default 100.
    timeout : float, optional
        Seconds to wait for each response. Default 120.

    Returns
    -------
    list of PEDEntry
        The matching entries, in the order PED returns them.

    Raises
    ------
    SSException
        If ``max_results`` or ``page_size`` is not a positive integer.
    PEDError
        If PED cannot be reached or returns an error.

    Example
    -------
    >>> from soursop import ped
    >>> for entry in ped.search_entries(uniprot="P37840"):
    ...     print(entry)
    """
    for label, value in (("max_results", max_results), ("page_size", page_size)):
        if value is None and label == "max_results":
            continue
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise SSException(f"{label} must be a positive integer; received {value!r}")

    supplied = {
        "query": query,
        "protein_name": protein_name,
        "uniprot": uniprot,
        "term": term,
        "author": author,
        "publication": publication,
        "publication_title": publication_title,
        "cross_ref": cross_ref,
    }
    params: Dict[str, Any] = {
        _CRITERIA[key]: value for key, value in supplied.items() if value is not None
    }

    records: List[Dict[str, Any]] = []
    offset = 0
    while True:
        limit = (
            page_size
            if max_results is None
            else min(page_size, max_results - len(records))
        )
        page = _http.get_json(
            "/entries", {**params, "limit": limit, "offset": offset}, timeout=timeout
        )
        if isinstance(page, list):
            # tolerate a bare list of entries (no pagination wrapper)
            records.extend(page)
            break
        items = page.get("result") or []
        records.extend(items)
        offset += len(items)
        count = page.get("count")
        if (
            len(items) == 0
            or (max_results is not None and len(records) >= max_results)
            or (isinstance(count, int) and offset >= count)
        ):
            break

    if max_results is not None:
        records = records[:max_results]
    return [PEDEntry.from_json(record) for record in records]


def search(
    query: Optional[str] = None,
    *,
    ensembles: bool = False,
    protein_name: Optional[str] = None,
    uniprot: Optional[str] = None,
    term: Optional[str] = None,
    author: Optional[str] = None,
    publication: Optional[str] = None,
    publication_title: Optional[str] = None,
    cross_ref: Optional[str] = None,
    max_results: Optional[int] = None,
    page_size: int = DEFAULT_PAGE_SIZE,
    timeout: float = _http.DEFAULT_TIMEOUT,
) -> List[str]:
    """Search PED and return the identifiers of the matches.

    Takes the same criteria as :func:`search_entries`; see there for their
    meaning.

    Parameters
    ----------
    query : str, optional
        Free-text search (PED's ``free_text``).
    ensembles : bool, optional
        If False (default) return entry identifiers (``"PED00001"``); if
        True return one identifier per ensemble (``"PED00001e001"``,
        ``"PED00001e002"``, ...), ready to pass to :func:`load_ensemble`.
    protein_name, uniprot, term, author, publication, publication_title, cross_ref : str, optional
        Further criteria (see :func:`search_entries`).
    max_results : int, optional
        Maximum number of *entries* to fetch. Default ``None`` (all).
    page_size : int, optional
        Entries requested per call to PED. Default 100.
    timeout : float, optional
        Seconds to wait for each response. Default 120.

    Returns
    -------
    list of str
        Entry identifiers, or ensemble identifiers if ``ensembles=True``.

    Raises
    ------
    SSException
        If ``max_results`` or ``page_size`` is invalid.
    PEDError
        If PED cannot be reached or returns an error.

    Example
    -------
    >>> from soursop import ped
    >>> ped.search("sic1")
    ['PED00001', ...]
    >>> ped.search("sic1", ensembles=True)
    ['PED00001e001', 'PED00001e002', 'PED00001e003', ...]
    """
    entries = search_entries(
        query,
        protein_name=protein_name,
        uniprot=uniprot,
        term=term,
        author=author,
        publication=publication,
        publication_title=publication_title,
        cross_ref=cross_ref,
        max_results=max_results,
        page_size=page_size,
        timeout=timeout,
    )
    if ensembles:
        return [identifier for entry in entries for identifier in entry.identifiers]
    return [entry.entry_id for entry in entries]
