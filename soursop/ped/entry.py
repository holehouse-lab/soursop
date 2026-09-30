"""
PED entries and their ensembles.

:class:`PEDEntry` is a light, typed view of the JSON record PED returns for an
entry (from :func:`get_entry` or a search). It exposes the fields you usually
want - title, authors, UniProt accessions, protein names and the ensembles
it contains - and keeps the full JSON as ``raw`` for anything else.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple, Union

from . import _http
from .identifiers import normalise_ensemble_id, parse_identifier

if TYPE_CHECKING:  # pragma: no cover
    from ..sstrajectory import SSTrajectory


def _unique(values: List[str]) -> Tuple[str, ...]:
    """Unique, non-empty strings in first-seen order."""
    seen: List[str] = []
    for value in values:
        if value and value not in seen:
            seen.append(value)
    return tuple(seen)


def _as_float(value: Any) -> Optional[float]:
    """``float(value)``, or ``None`` if it is missing or not a number."""
    try:
        return None if value is None else float(value)
    except (TypeError, ValueError):
        return None


@dataclass(frozen=True)
class PEDChainStats:
    """What PED reports about one chain of one ensemble.

    Attributes
    ----------
    chain_name : str
        Chain identifier in the ensemble's PDB file.
    chain_type : str or None
        Chain type, e.g. ``"protein"``.
    sequence : str or None
        One-letter sequence of the chain in the ensemble.
    rg_mean : float or None
        Mean radius of gyration, if PED reports it per chain.
    entropy_dssp_mean : float or None
        Mean per-residue DSSP entropy, if PED reports it per chain.
    relative_asa_mean : float or None
        Mean relative accessible surface area, if PED reports it per chain.
    """

    chain_name: str
    chain_type: Optional[str] = None
    sequence: Optional[str] = None
    rg_mean: Optional[float] = None
    entropy_dssp_mean: Optional[float] = None
    relative_asa_mean: Optional[float] = None


@dataclass(frozen=True)
class PEDEnsembleSummary:
    """One ensemble within a PED entry.

    Attributes
    ----------
    entry_id : str
        The entry the ensemble belongs to, e.g. ``"PED00001"``.
    ensemble_id : str
        The ensemble within the entry, e.g. ``"e001"``.
    n_models : int or None
        Number of conformers (models) in the ensemble.
    chains : tuple of PEDChainStats
        The chains in the ensemble (names, types, sequences).
    only_ca : bool or None
        True if the ensemble only has C-alpha coordinates (a coarse-grained
        ensemble, which SOURSOP treats as one bead per residue).
    ped_statistics : dict
        PED's own summary of the ensemble (``analysis.summary_dataframe``),
        e.g. ``rg_mean``, ``end_to_end_mean`` and ``flory_exponent``.
        Distances here are in **nanometres**, as PED reports them. Empty if
        PED does not provide a summary (search results never include one).
    """

    entry_id: str
    ensemble_id: str
    n_models: Optional[int]
    chains: Tuple[PEDChainStats, ...]
    only_ca: Optional[bool] = None
    ped_statistics: Dict[str, Any] = field(default_factory=dict, compare=False)

    @property
    def identifier(self) -> str:
        """Full ensemble identifier, e.g. ``"PED00001e001"``."""
        return f"{self.entry_id}{self.ensemble_id}"

    def load(self, **kwargs: Any) -> "SSTrajectory":
        """Download this ensemble and return it as an ``SSTrajectory``.

        Parameters
        ----------
        **kwargs
            Passed to :func:`soursop.ped.load_ensemble` (e.g. ``cache_dir``,
            or keyword arguments for ``SSTrajectory``).

        Returns
        -------
        soursop.sstrajectory.SSTrajectory
            The ensemble, one frame per conformer.
        """
        from .ensemble import load_ensemble

        return load_ensemble(self.entry_id, self.ensemble_id, **kwargs)


@dataclass(frozen=True)
class PEDEntry:
    """A PED entry: its metadata and the ensembles it contains.

    Attributes
    ----------
    entry_id : str
        Entry identifier, e.g. ``"PED00001"``.
    title : str
        Entry title.
    version : int or None
        Entry version (only reported by :func:`get_entry`, not by searches).
    creation_date : str or None
        Date the entry was created, as reported by PED.
    authors : tuple of str
        Author names.
    uniprot_accessions : tuple of str
        UniProt accessions of the constructs' fragments.
    protein_names : tuple of str
        Descriptions of the constructs' fragments (usually protein names).
    publication : str or None
        PubMed ID or DOI of the associated publication.
    ontology_terms : tuple of str
        Names of the IDPO terms attached to the entry (measurement and
        ensemble-generation methods).
    cross_references : tuple of (str, str)
        ``(database, identifier)`` pairs, from both the entry and
        experimental cross-references (e.g. ``("disprot", "DP00631")``).
    ensembles : tuple of PEDEnsembleSummary
        The ensembles in this entry.
    raw : dict
        The full JSON record returned by PED.
    """

    entry_id: str
    title: str
    version: Optional[int]
    creation_date: Optional[str]
    authors: Tuple[str, ...]
    uniprot_accessions: Tuple[str, ...]
    protein_names: Tuple[str, ...]
    publication: Optional[str]
    ontology_terms: Tuple[str, ...]
    cross_references: Tuple[Tuple[str, str], ...]
    ensembles: Tuple[PEDEnsembleSummary, ...]
    raw: Dict[str, Any] = field(repr=False, compare=False)

    @classmethod
    def from_json(cls, data: Dict[str, Any]) -> "PEDEntry":
        """Build a :class:`PEDEntry` from PED's JSON record for an entry.

        Missing fields are tolerated (search results omit some of them).

        Parameters
        ----------
        data : dict
            One entry, as returned by ``/entries/{id}`` or in the ``result``
            list of ``/entries``.

        Returns
        -------
        PEDEntry
            The parsed entry.

        Raises
        ------
        PEDError
            If the record has no ``entry_id``.
        """
        if not isinstance(data, dict) or not data.get("entry_id"):
            raise _http.PEDError(
                f"PED returned an entry record without an entry_id: {data!r}"
            )
        entry_id = str(data["entry_id"]).upper()
        description = data.get("description") or {}

        fragments = [
            fragment
            for chain in data.get("construct_chains") or []
            for fragment in chain.get("fragments") or []
        ]
        cross_refs = [
            (str(ref.get("db", "")), str(ref.get("id", "")))
            for key in ("entry_cross_reference", "experimental_cross_reference")
            for ref in description.get(key) or []
        ]

        ensembles = [
            _parse_ensemble(entry_id, ensemble)
            for ensemble in data.get("ensembles") or []
            if isinstance(ensemble, dict) and ensemble.get("ensemble_id")
        ]

        version = data.get("version")
        return cls(
            entry_id=entry_id,
            title=str(description.get("title") or ""),
            version=int(version) if isinstance(version, (int, float)) else None,
            creation_date=data.get("creation_date"),
            authors=_unique(
                [str(a.get("name") or "") for a in description.get("authors") or []]
            ),
            uniprot_accessions=_unique(
                [str(f.get("uniprot_acc") or "") for f in fragments]
            ),
            protein_names=_unique([str(f.get("description") or "") for f in fragments]),
            publication=description.get("publication_identifier") or None,
            ontology_terms=_unique(
                [
                    str(t.get("name") or "")
                    for t in description.get("ontology_terms") or []
                ]
            ),
            cross_references=tuple(ref for ref in cross_refs if ref[0] and ref[1]),
            ensembles=tuple(ensembles),
            raw=data,
        )

    @property
    def ensemble_ids(self) -> Tuple[str, ...]:
        """Ensemble identifiers within this entry, e.g. ``("e001", "e002")``."""
        return tuple(e.ensemble_id for e in self.ensembles)

    @property
    def identifiers(self) -> Tuple[str, ...]:
        """Full ensemble identifiers, e.g. ``("PED00001e001", "PED00001e002")``."""
        return tuple(e.identifier for e in self.ensembles)

    def ensemble(self, ensemble_id: Union[str, int]) -> PEDEnsembleSummary:
        """Return one of this entry's ensembles.

        Parameters
        ----------
        ensemble_id : str or int
            The ensemble, e.g. ``"e001"`` or ``1``.

        Returns
        -------
        PEDEnsembleSummary
            The requested ensemble.

        Raises
        ------
        SSException
            If the entry has no such ensemble.
        """
        wanted = normalise_ensemble_id(ensemble_id)
        for summary in self.ensembles:
            if summary.ensemble_id == wanted:
                return summary
        from ..ssexceptions import SSException

        raise SSException(
            f"{self.entry_id} has no ensemble {wanted}; its ensembles are {list(self.ensemble_ids)}"
        )

    def load(
        self, ensemble_id: Optional[Union[str, int]] = None, **kwargs: Any
    ) -> "SSTrajectory":
        """Download one of this entry's ensembles as an ``SSTrajectory``.

        Parameters
        ----------
        ensemble_id : str or int, optional
            Which ensemble to load. May be omitted if the entry contains
            only one ensemble.
        **kwargs
            Passed to :func:`soursop.ped.load_ensemble`.

        Returns
        -------
        soursop.sstrajectory.SSTrajectory
            The ensemble, one frame per conformer.
        """
        from .ensemble import load_ensemble

        if ensemble_id is None:
            ensemble_id = _only_ensemble(self)
        return load_ensemble(self.entry_id, ensemble_id, **kwargs)

    def __str__(self) -> str:
        proteins = ", ".join(self.protein_names) or "unknown protein"
        return (
            f"{self.entry_id}: {self.title} [{proteins}; "
            f"{len(self.ensembles)} ensemble(s): {', '.join(self.ensemble_ids)}]"
        )


def _parse_ensemble(entry_id: str, ensemble: Dict[str, Any]) -> PEDEnsembleSummary:
    """Parse one element of an entry's ``ensembles`` list.

    PED's live API puts chain names, sequences, the C-alpha-only flag and its
    own summary statistics under ``ensemble_details``; the API description
    shows per-chain statistics under ``chains``. Both are read if present.
    """
    details = ensemble.get("ensemble_details")
    if not isinstance(details, dict):
        details = {}

    chains: Dict[str, PEDChainStats] = {}
    for chain in details.get("chains") or []:
        if isinstance(chain, dict) and chain.get("chain_name") is not None:
            name = str(chain["chain_name"])
            chains[name] = PEDChainStats(
                chain_name=name,
                chain_type=chain.get("chain_type"),
                sequence=chain.get("sequence"),
            )
    for chain in ensemble.get("chains") or []:
        if not isinstance(chain, dict) or chain.get("chain_name") is None:
            continue
        name = str(chain["chain_name"])
        known = chains.get(name, PEDChainStats(chain_name=name))
        chains[name] = PEDChainStats(
            chain_name=name,
            chain_type=known.chain_type,
            sequence=known.sequence,
            rg_mean=_as_float(chain.get("rg_mean")),
            entropy_dssp_mean=_as_float(chain.get("entropy_dssp_mean")),
            relative_asa_mean=_as_float(chain.get("relative_asa_mean")),
        )

    models = ensemble.get("models", details.get("models"))
    only_ca = details.get("only_CA")
    return PEDEnsembleSummary(
        entry_id=entry_id,
        ensemble_id=normalise_ensemble_id(str(ensemble["ensemble_id"])),
        n_models=int(models) if isinstance(models, (int, float)) else None,
        chains=tuple(chains.values()),
        only_ca=bool(only_ca) if isinstance(only_ca, bool) else None,
        ped_statistics=_ped_statistics(details),
    )


def _ped_statistics(details: Dict[str, Any]) -> Dict[str, Any]:
    """PED's summary row for the ensemble, or ``{}`` if there is none."""
    analysis = details.get("analysis")
    if not isinstance(analysis, dict):
        return {}
    rows = analysis.get("summary_dataframe")
    if isinstance(rows, str):
        try:
            rows = json.loads(rows)
        except json.JSONDecodeError:
            return {}
    if isinstance(rows, dict):
        rows = [rows]
    if not isinstance(rows, list):
        return {}
    candidates = [row for row in rows if isinstance(row, dict)]
    for row in candidates:
        if row.get("ensemble_code") == "input_ensemble":
            return dict(row)
    return dict(candidates[0]) if candidates else {}


def _only_ensemble(entry: PEDEntry) -> str:
    """The ensemble of a single-ensemble entry; raise if there is not exactly one."""
    from ..ssexceptions import SSException

    if len(entry.ensembles) == 1:
        return entry.ensembles[0].ensemble_id
    if len(entry.ensembles) == 0:
        raise SSException(f"PED reports no ensembles for {entry.entry_id}")
    raise SSException(
        f"{entry.entry_id} contains {len(entry.ensembles)} ensembles "
        f"({', '.join(entry.identifiers)}); say which one to use, e.g. "
        f"'{entry.identifiers[0]}' or ensemble_id='{entry.ensemble_ids[0]}'"
    )


def get_entry(identifier: str, timeout: float = _http.DEFAULT_TIMEOUT) -> PEDEntry:
    """Fetch the metadata for a PED entry.

    Parameters
    ----------
    identifier : str
        Entry identifier (``"PED00001"``). A full ensemble identifier
        (``"PED00001e001"``) is also accepted; its ensemble part is ignored.
    timeout : float, optional
        Seconds to wait for PED. Default 120.

    Returns
    -------
    PEDEntry
        The entry, including the list of ensembles it contains.

    Raises
    ------
    PEDError
        If the entry does not exist (``status == 404``) or PED cannot be
        reached.
    SSException
        If ``identifier`` is malformed.

    Example
    -------
    >>> from soursop import ped
    >>> entry = ped.get_entry("PED00001")
    >>> entry.identifiers
    ('PED00001e001', 'PED00001e002', 'PED00001e003')
    """
    entry_id, _ = parse_identifier(identifier)
    data = _http.get_json(f"/entries/{entry_id}", timeout=timeout)
    return PEDEntry.from_json(data)


def list_ensembles(
    identifier: str, timeout: float = _http.DEFAULT_TIMEOUT
) -> Tuple[PEDEnsembleSummary, ...]:
    """List the ensembles in a PED entry.

    Parameters
    ----------
    identifier : str
        Entry identifier, e.g. ``"PED00001"``.
    timeout : float, optional
        Seconds to wait for PED. Default 120.

    Returns
    -------
    tuple of PEDEnsembleSummary
        One per ensemble, with its identifier, number of models and PED's
        per-chain summary statistics.

    Raises
    ------
    PEDError
        If the entry does not exist or PED cannot be reached.
    """
    return get_entry(identifier, timeout=timeout).ensembles
