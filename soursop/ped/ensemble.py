"""
Downloading PED ensembles and turning them into SOURSOP objects.

The main entry point is :func:`load_ensemble`, which downloads an ensemble's
multi-model PDB file and returns it as an
:class:`~soursop.sstrajectory.SSTrajectory` with one frame per conformer.
:func:`download_ensemble` just saves the PDB file, :func:`get_ensemble_weights`
fetches per-conformer weights for weighted ensembles, and
:func:`get_ensemble_asset` / :func:`download_all_data` give access to PED's
other per-ensemble assets (its own DSSP, Rg, validation data and so on).
"""

from __future__ import annotations

import gzip
import io
import json
import os
import re
import tarfile
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple, Union

import numpy as np

from ..ssexceptions import SSException
from . import _http
from .entry import _only_ensemble, get_entry
from .identifiers import parse_identifier

if TYPE_CHECKING:  # pragma: no cover
    from ..sstrajectory import SSTrajectory

#: Per-ensemble assets PED serves at ``/entries/{id}/ensembles/{ens}/{asset}``.
ENSEMBLE_ASSETS = (
    "ensemble-pdb",
    "ensemble-sample",
    "dssp-consensus",
    "dssp-data",
    "cb-dev",
    "clash",
    "geom-angles",
    "geom-bonds",
    "geom-carbons",
    "ramachandran",
    "rotamer",
    "gyration-chains",
    "gyration-global",
    "dmax-data",
    "dmax-global",
    "weights",
)

#: Assets PED can return as JSON or CSV (the ``response_format`` option).
TABULAR_ASSETS = (
    "dssp-consensus",
    "dssp-data",
    "cb-dev",
    "clash",
    "geom-angles",
    "geom-bonds",
    "geom-carbons",
    "ramachandran",
    "rotamer",
    "gyration-chains",
    "gyration-global",
    "dmax-data",
    "dmax-global",
)

PathLike = Union[str, "os.PathLike[str]"]


def _resolve(
    identifier: str, ensemble_id: Optional[Union[str, int]], timeout: float
) -> Tuple[str, str]:
    """Entry and ensemble for a request, looking the ensemble up if needed.

    If no ensemble is given, the entry must contain exactly one.
    """
    entry_id, ensemble = parse_identifier(identifier, ensemble_id)
    if ensemble is None:
        ensemble = _only_ensemble(get_entry(entry_id, timeout=timeout))
    return entry_id, ensemble


def _extract_pdb(payload: bytes, label: str) -> bytes:
    """Return the PDB text inside a downloaded ``ensemble-pdb`` payload.

    PED's API description does not pin down the container, so this accepts a
    plain PDB file, a gzip-compressed one, or a (possibly compressed) tar
    archive holding a single PDB file.

    Parameters
    ----------
    payload : bytes
        The downloaded bytes.
    label : str
        Identifier used in error messages.

    Returns
    -------
    bytes
        The PDB file contents.

    Raises
    ------
    PEDError
        If the payload is not a PDB file, is mmCIF, or is an archive that
        does not hold exactly one PDB file.
    """
    data = payload
    if data[:2] == b"\x1f\x8b":
        try:
            data = gzip.decompress(data)
        except (OSError, EOFError) as e:
            raise _http.PEDError(
                f"{label}: could not decompress the download ({e})"
            ) from None

    if len(data) > 262 and data[257:262] == b"ustar":
        with tarfile.open(fileobj=io.BytesIO(data)) as tar:
            members = [
                m
                for m in tar.getmembers()
                if m.isfile() and m.name.lower().endswith((".pdb", ".pdb.gz", ".ent"))
            ]
            if len(members) != 1:
                names = [m.name for m in tar.getmembers() if m.isfile()]
                raise _http.PEDError(
                    f"{label}: expected one PDB file in the downloaded archive, found {names}"
                )
            handle = tar.extractfile(members[0])
            if handle is None:  # pragma: no cover - isfile() above
                raise _http.PEDError(
                    f"{label}: could not read {members[0].name} from the archive"
                )
            return _extract_pdb(handle.read(), label)

    if data.lstrip()[:5] == b"data_":
        raise _http.PEDError(
            f"{label}: PED returned the ensemble in mmCIF format, which SOURSOP cannot read yet"
        )
    if re.search(rb"^(ATOM  |HETATM|MODEL )", data, re.MULTILINE) is None:
        raise _http.PEDError(
            f"{label}: the downloaded ensemble does not look like a PDB file"
        )
    return data


def download_ensemble(
    identifier: str,
    ensemble_id: Optional[Union[str, int]] = None,
    *,
    directory: PathLike = ".",
    overwrite: bool = False,
    timeout: float = _http.DEFAULT_TIMEOUT,
) -> str:
    """Download a PED ensemble as a multi-model PDB file.

    The file is saved as ``<directory>/<PEDxxxxxeNNN>.pdb``, one model per
    conformer. If that file already exists it is reused rather than
    downloaded again, unless ``overwrite=True``.

    Parameters
    ----------
    identifier : str
        Ensemble identifier (``"PED00001e001"``), or an entry identifier
        (``"PED00001"``) together with ``ensemble_id``. An entry identifier
        on its own is accepted if the entry contains only one ensemble.
    ensemble_id : str or int, optional
        Ensemble within the entry (``"e001"`` or ``1``).
    directory : str or os.PathLike, optional
        Directory to save into (created if needed). Default: the current
        directory.
    overwrite : bool, optional
        Download again even if the file already exists. Default False.
    timeout : float, optional
        Seconds to wait for PED. Default 120.

    Returns
    -------
    str
        Path to the PDB file.

    Raises
    ------
    PEDError
        If the ensemble does not exist, PED cannot be reached, or the
        download is not a PDB file.
    SSException
        If the identifier is malformed or ambiguous.

    Example
    -------
    >>> from soursop import ped
    >>> ped.download_ensemble("PED00001e001", directory="ped_data")
    'ped_data/PED00001e001.pdb'
    """
    entry_id, ensemble = _resolve(identifier, ensemble_id, timeout)
    label = f"{entry_id}{ensemble}"
    folder = Path(directory)
    path = folder / f"{label}.pdb"
    if path.exists() and not overwrite:
        return str(path)

    payload = _http.get_bytes(
        f"/entries/{entry_id}/ensembles/{ensemble}/ensemble-pdb", timeout=timeout
    )
    pdb = _extract_pdb(payload, label)

    folder.mkdir(parents=True, exist_ok=True)
    # write to a temporary file first so an interrupted download never
    # leaves a truncated PDB behind that a later call would reuse
    handle, tmp_name = tempfile.mkstemp(dir=folder, prefix=f".{label}.", suffix=".tmp")
    try:
        with os.fdopen(handle, "wb") as out:
            out.write(pdb)
        os.replace(tmp_name, path)
    except BaseException:
        if os.path.exists(tmp_name):
            os.remove(tmp_name)
        raise
    return str(path)


def load_ensemble(
    identifier: str,
    ensemble_id: Optional[Union[str, int]] = None,
    *,
    cache_dir: Optional[PathLike] = None,
    timeout: float = _http.DEFAULT_TIMEOUT,
    **kwargs: Any,
) -> "SSTrajectory":
    """Download a PED ensemble and return it as an ``SSTrajectory``.

    Each conformer in the ensemble becomes one frame. The protein chains are
    then available as usual through ``proteinTrajectoryList``.

    Parameters
    ----------
    identifier : str
        Ensemble identifier (``"PED00001e001"``), or an entry identifier
        (``"PED00001"``) together with ``ensemble_id``. An entry identifier
        on its own is accepted if the entry contains only one ensemble.
    ensemble_id : str or int, optional
        Ensemble within the entry (``"e001"`` or ``1``).
    cache_dir : str or os.PathLike, optional
        If given, the PDB file is kept in this directory and reused by later
        calls, so each ensemble is only downloaded once. If ``None``
        (default) it is downloaded to a temporary directory that is deleted
        once the ensemble has been read in.
    timeout : float, optional
        Seconds to wait for PED. Default 120.
    **kwargs
        Passed to :class:`~soursop.sstrajectory.SSTrajectory` (e.g.
        ``protein_grouping`` or ``extra_valid_residue_names``).

    Returns
    -------
    soursop.sstrajectory.SSTrajectory
        The ensemble.

    Raises
    ------
    PEDError
        If the ensemble does not exist, PED cannot be reached, or the
        download is not a PDB file.
    SSException
        If the identifier is malformed or ambiguous.

    Notes
    -----
    Some PED ensembles are *weighted*: their conformers are not equally
    likely. Fetch the weights with :func:`get_ensemble_weights` and pass them
    to SOURSOP's analysis methods via ``weights=``.

    Example
    -------
    >>> from soursop import ped
    >>> traj = ped.load_ensemble("PED00001e001")
    >>> protein = traj.proteinTrajectoryList[0]
    >>> protein.get_radius_of_gyration().mean()
    """
    reserved = {"trajectory_filename", "pdb_filename", "TRJ"} & set(kwargs)
    if reserved:
        raise SSException(
            f"load_ensemble() sets {sorted(reserved)} itself; pass only other SSTrajectory options"
        )

    if cache_dir is not None:
        path = download_ensemble(
            identifier, ensemble_id, directory=cache_dir, timeout=timeout
        )
        return _read_pdb_ensemble(path, kwargs)

    with tempfile.TemporaryDirectory(prefix="soursop_ped_") as tmp:
        path = download_ensemble(
            identifier, ensemble_id, directory=tmp, timeout=timeout
        )
        return _read_pdb_ensemble(path, kwargs)


def _read_pdb_ensemble(path: str, kwargs: Dict[str, Any]) -> "SSTrajectory":
    """Read a multi-model PDB file (trajectory and topology in one) into SOURSOP."""
    from ..sstrajectory import SSTrajectory

    # SSTrajectory itself is not type-annotated
    trajectory: SSTrajectory = SSTrajectory(path, path, **kwargs)  # type: ignore[no-untyped-call]
    return trajectory


def get_ensemble_asset(
    identifier: str,
    ensemble_id: Optional[Union[str, int]] = None,
    *,
    asset: str,
    response_format: str = "json",
    only_features: Optional[bool] = None,
    timeout: float = _http.DEFAULT_TIMEOUT,
) -> Any:
    """Fetch one of PED's per-ensemble assets.

    PED's API description lists per-ensemble analysis files (DSSP, radius of
    gyration and Dmax, Ramachandran and rotamer statistics, clash and
    geometry checks). At the time of writing PED no longer serves most of
    them separately (e.g. the gyration and Dmax files return 404); its
    summary statistics are available as ``PEDEnsembleSummary.ped_statistics``
    from :func:`get_entry` instead.

    Parameters
    ----------
    identifier : str
        Ensemble identifier (``"PED00001e001"``), or an entry identifier
        together with ``ensemble_id``.
    ensemble_id : str or int, optional
        Ensemble within the entry.
    asset : str
        One of :data:`ENSEMBLE_ASSETS`, e.g. ``"gyration-chains"``,
        ``"dssp-consensus"`` or ``"weights"``.
    response_format : {'json', 'csv'}, optional
        Format for the tabular assets (:data:`TABULAR_ASSETS`). Default
        ``'json'``. Ignored for the others.
    only_features : bool, optional
        PED's ``only_features`` option (``dssp-consensus`` only).
    timeout : float, optional
        Seconds to wait for PED. Default 120.

    Returns
    -------
    Any
        Decoded JSON for a tabular asset with ``response_format='json'``;
        the text for ``response_format='csv'``; otherwise the raw bytes.

    Raises
    ------
    SSException
        If ``asset`` or ``response_format`` is not recognised.
    PEDError
        If the asset does not exist for this ensemble or PED cannot be
        reached.

    Example
    -------
    >>> from soursop import ped
    >>> rg = ped.get_ensemble_asset("PED00001e001", asset="gyration-chains")
    """
    if asset not in ENSEMBLE_ASSETS:
        raise SSException(
            f"Unknown PED asset {asset!r}; choose one of {list(ENSEMBLE_ASSETS)}"
        )
    if response_format not in ("json", "csv"):
        raise SSException(
            f"response_format must be 'json' or 'csv'; received {response_format!r}"
        )

    entry_id, ensemble = _resolve(identifier, ensemble_id, timeout)
    path = f"/entries/{entry_id}/ensembles/{ensemble}/{asset}"
    params: Dict[str, _http.QueryValue] = {"only_features": only_features}

    if asset in TABULAR_ASSETS:
        params["response_format"] = response_format
        if response_format == "json":
            return _http.get_json(path, params, timeout=timeout)
        return _http.get_bytes(path, params, timeout=timeout).decode(
            "utf-8", errors="replace"
        )
    return _http.get_bytes(path, params, timeout=timeout)


def _parse_weights(payload: bytes, label: str) -> np.ndarray:
    """Per-conformer weights from PED's ``weights`` asset.

    Accepts JSON (a list of numbers, a list of records with a ``weight``
    field, or an object with a ``weights`` list) or delimited text with the
    weight in the last column of each line (comment and header lines are
    skipped).
    """
    text = payload.decode("utf-8", errors="replace").strip()
    values: List[float] = []
    try:
        data = json.loads(text)
    except json.JSONDecodeError:
        data = None

    if data is not None:
        if isinstance(data, dict):
            data = data.get("weights", data.get("result"))
        if isinstance(data, list):
            for item in data:
                value = item.get("weight") if isinstance(item, dict) else item
                if isinstance(value, bool) or not isinstance(value, (int, float, str)):
                    raise _http.PEDError(
                        f"{label}: could not read the weights PED returned"
                    )
                try:
                    values.append(float(value))
                except ValueError:
                    raise _http.PEDError(
                        f"{label}: could not read the weights PED returned"
                    ) from None
    else:
        for line in text.splitlines():
            tokens = [t for t in re.split(r"[,;\s]+", line.strip()) if t]
            if not tokens or tokens[0].startswith("#"):
                continue
            try:
                values.append(float(tokens[-1]))
            except ValueError:
                continue  # header line

    if len(values) == 0:
        raise _http.PEDError(f"{label}: PED returned no weights")
    weights = np.asarray(values, dtype=np.float64)
    if not np.all(np.isfinite(weights)) or np.any(weights < 0):
        raise _http.PEDError(
            f"{label}: PED returned invalid (negative or non-finite) weights"
        )
    return weights


def get_ensemble_weights(
    identifier: str,
    ensemble_id: Optional[Union[str, int]] = None,
    *,
    normalise: bool = True,
    timeout: float = _http.DEFAULT_TIMEOUT,
) -> Optional[np.ndarray]:
    """Fetch the per-conformer weights of a weighted PED ensemble.

    Parameters
    ----------
    identifier : str
        Ensemble identifier (``"PED00001e001"``), or an entry identifier
        together with ``ensemble_id``.
    ensemble_id : str or int, optional
        Ensemble within the entry.
    normalise : bool, optional
        If True (default) the weights are scaled to sum to 1, ready to pass
        straight to SOURSOP's ``weights=`` arguments.
    timeout : float, optional
        Seconds to wait for PED. Default 120.

    Returns
    -------
    numpy.ndarray or None
        One weight per conformer, in model order; ``None`` if PED has no
        weights for this ensemble (i.e. its conformers are equally
        weighted).

    Raises
    ------
    PEDError
        If PED cannot be reached, or returns weights that cannot be read.

    Example
    -------
    >>> from soursop import ped
    >>> w = ped.get_ensemble_weights("PED00001e001")
    >>> traj = ped.load_ensemble("PED00001e001")
    >>> rg = traj.proteinTrajectoryList[0].get_radius_of_gyration(weights=w)
    """
    entry_id, ensemble = _resolve(identifier, ensemble_id, timeout)
    try:
        payload = _http.get_bytes(
            f"/entries/{entry_id}/ensembles/{ensemble}/weights", timeout=timeout
        )
    except _http.PEDError as e:
        if e.status == 404:
            return None
        raise
    if len(payload.strip()) == 0:
        return None
    weights = _parse_weights(payload, f"{entry_id}{ensemble}")
    if normalise:
        total = float(np.sum(weights))
        if not total > 0:
            raise _http.PEDError(f"{entry_id}{ensemble}: PED's weights sum to zero")
        weights = weights / total
    return weights


def download_all_data(
    identifier: str,
    ensemble_id: Optional[Union[str, int]] = None,
    *,
    directory: PathLike = ".",
    response_format: str = "json",
    overwrite: bool = False,
    timeout: float = _http.DEFAULT_TIMEOUT,
) -> str:
    """Download every PED asset for an ensemble as a single ``.tar.gz`` file.

    Parameters
    ----------
    identifier : str
        Ensemble identifier (``"PED00001e001"``), or an entry identifier
        together with ``ensemble_id``.
    ensemble_id : str or int, optional
        Ensemble within the entry.
    directory : str or os.PathLike, optional
        Directory to save into (created if needed). Default: the current
        directory.
    response_format : {'json', 'csv'}, optional
        Format of the tabular files inside the archive. Default ``'json'``.
    overwrite : bool, optional
        Download again even if the file already exists. Default False.
    timeout : float, optional
        Seconds to wait for PED. Default 120.

    Returns
    -------
    str
        Path to ``<directory>/<PEDxxxxxeNNN>_all_data.tar.gz``.

    Raises
    ------
    SSException
        If ``response_format`` is not recognised.
    PEDError
        If the ensemble does not exist or PED cannot be reached.
    """
    if response_format not in ("json", "csv"):
        raise SSException(
            f"response_format must be 'json' or 'csv'; received {response_format!r}"
        )
    entry_id, ensemble = _resolve(identifier, ensemble_id, timeout)
    folder = Path(directory)
    path = folder / f"{entry_id}{ensemble}_all_data.tar.gz"
    if path.exists() and not overwrite:
        return str(path)
    payload = _http.get_bytes(
        f"/entries/{entry_id}/ensembles/{ensemble}/download-all-data/",
        {"response_format": response_format},
        timeout=timeout,
    )
    folder.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)
    return str(path)
