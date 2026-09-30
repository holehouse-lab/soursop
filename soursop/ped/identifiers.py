"""
Parsing and normalising PED identifiers.

PED uses two kinds of identifier:

* an **entry** identifier, ``PED`` followed by five digits (``PED00001``),
  which groups one or more ensembles deposited together; and
* an **ensemble** identifier, ``e`` followed by three digits (``e001``),
  which is only unique within its entry. The two are often written together
  as a single string, e.g. ``PED00001e001``.
"""

from __future__ import annotations

import re
from typing import Optional, Tuple, Union

from ..ssexceptions import SSException

_ENTRY_RE = re.compile(r"^PED\d{5,}$", re.IGNORECASE)
_FULL_RE = re.compile(r"^(PED\d{5,})(E\d{3,})$", re.IGNORECASE)
_ENSEMBLE_RE = re.compile(r"^E?(\d+)$", re.IGNORECASE)


def normalise_ensemble_id(ensemble_id: Union[str, int]) -> str:
    """Normalise an ensemble identifier to PED's ``eNNN`` form.

    Parameters
    ----------
    ensemble_id : str or int
        The ensemble, as ``"e001"``, ``"E1"``, ``"1"`` or ``1``.

    Returns
    -------
    str
        The ensemble identifier in PED's form, e.g. ``"e001"``.

    Raises
    ------
    SSException
        If ``ensemble_id`` is not a positive integer or ``e``-prefixed
        integer.

    Example
    -------
    >>> normalise_ensemble_id(3)
    'e003'
    """
    if isinstance(ensemble_id, bool):
        raise SSException(f"Invalid PED ensemble identifier: {ensemble_id!r}")
    if isinstance(ensemble_id, int):
        number = ensemble_id
    else:
        match = _ENSEMBLE_RE.match(str(ensemble_id).strip())
        if match is None:
            raise SSException(
                f"Invalid PED ensemble identifier {ensemble_id!r}; expected e.g. 'e001' or 1"
            )
        number = int(match.group(1))
    if number < 1:
        raise SSException(f"PED ensemble numbers start at 1; received {ensemble_id!r}")
    return f"e{number:03d}"


def parse_identifier(
    identifier: str, ensemble_id: Optional[Union[str, int]] = None
) -> Tuple[str, Optional[str]]:
    """Split a PED identifier into its entry and ensemble parts.

    Accepts either an entry identifier (``"PED00001"``), optionally with a
    separate ``ensemble_id``, or a full ensemble identifier
    (``"PED00001e001"``). Case and surrounding whitespace are ignored.

    Parameters
    ----------
    identifier : str
        ``"PED00001"`` or ``"PED00001e001"``.
    ensemble_id : str or int, optional
        Ensemble within the entry (``"e001"``, ``1``, ...). If
        ``identifier`` already names an ensemble, this must agree with it.

    Returns
    -------
    tuple of (str, str or None)
        ``(entry_id, ensemble_id)``, e.g. ``("PED00001", "e001")``. The
        ensemble part is ``None`` if no ensemble was given.

    Raises
    ------
    SSException
        If the identifier is malformed, or it names a different ensemble
        from ``ensemble_id``.

    Example
    -------
    >>> parse_identifier("ped00001E001")
    ('PED00001', 'e001')
    >>> parse_identifier("PED00001", 2)
    ('PED00001', 'e002')
    """
    if not isinstance(identifier, str):
        raise SSException(f"PED identifiers are strings; received {identifier!r}")
    text = identifier.strip()

    full = _FULL_RE.match(text)
    if full is not None:
        entry = full.group(1).upper()
        embedded = normalise_ensemble_id(full.group(2))
        if ensemble_id is not None and normalise_ensemble_id(ensemble_id) != embedded:
            raise SSException(
                f"{identifier!r} names ensemble {embedded}, but ensemble_id={ensemble_id!r} was also given"
            )
        return entry, embedded

    if _ENTRY_RE.match(text) is None:
        raise SSException(
            f"Invalid PED identifier {identifier!r}; expected an entry such as "
            "'PED00001' or an ensemble such as 'PED00001e001'"
        )
    entry = text.upper()
    if ensemble_id is None:
        return entry, None
    return entry, normalise_ensemble_id(ensemble_id)


def full_identifier(entry_id: str, ensemble_id: Union[str, int]) -> str:
    """Join an entry and ensemble into a single identifier (``PED00001e001``).

    Parameters
    ----------
    entry_id : str
        Entry identifier, e.g. ``"PED00001"``.
    ensemble_id : str or int
        Ensemble identifier, e.g. ``"e001"`` or ``1``.

    Returns
    -------
    str
        The combined identifier.

    Raises
    ------
    SSException
        If either part is malformed.
    """
    entry, ensemble = parse_identifier(entry_id, ensemble_id)
    return f"{entry}{ensemble}"
