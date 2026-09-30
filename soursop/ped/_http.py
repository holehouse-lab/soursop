"""
Low-level HTTP access to the PED REST API.

Every request SOURSOP makes to PED goes through :func:`_request`, so this is
the one place that deals with URLs, timeouts and network errors. It uses only
the Python standard library (no extra dependency). The public functions in
:mod:`soursop.ped` build on :func:`get_json` and :func:`get_bytes`.

The API is described at https://proteinensemble.org/api (OpenAPI spec,
API version 5.0).
"""

from __future__ import annotations

import json
import urllib.error
import urllib.parse
import urllib.request
from typing import Any, Mapping, Optional, Union

from ..ssexceptions import SSException

#: Base URL of the PED REST API.
PED_API_URL = "https://deposition.proteinensemble.org/v1"

#: Default time (seconds) to wait for PED to respond. Ensemble downloads can
#: be tens of MB, so this is deliberately generous.
DEFAULT_TIMEOUT = 120.0

QueryValue = Union[str, int, bool, None]


class PEDError(SSException):
    """Raised when a request to the Protein Ensemble Database fails.

    Subclasses :class:`~soursop.ssexceptions.SSException`, so catching
    ``SSException`` catches PED errors too.

    Parameters
    ----------
    message : str
        Human-readable description of the problem.
    status : int or None, optional
        HTTP status code returned by PED, if the server responded at all
        (e.g. ``404`` for an entry or asset that does not exist). ``None``
        when PED could not be reached.
    """

    def __init__(self, message: str, status: Optional[int] = None) -> None:
        super().__init__(message)
        self.status = status


def _user_agent() -> str:
    """User-Agent header identifying SOURSOP (and its version) to PED."""
    try:
        from .. import __version__
    except ImportError:  # pragma: no cover - only in a broken install
        __version__ = "unknown"
    return f"soursop/{__version__} (+https://github.com/holehouse-lab/soursop)"


def build_url(path: str, params: Optional[Mapping[str, QueryValue]] = None) -> str:
    """Build a full PED API URL from a path and optional query parameters.

    Parameters
    ----------
    path : str
        Endpoint path relative to :data:`PED_API_URL`, e.g. ``"/entries"``.
    params : mapping, optional
        Query parameters. Entries whose value is ``None`` are dropped, and
        booleans are sent as ``true`` / ``false``.

    Returns
    -------
    str
        The full URL.
    """
    url = PED_API_URL.rstrip("/") + "/" + path.lstrip("/")
    if params:
        query = {
            key: (str(value).lower() if isinstance(value, bool) else str(value))
            for key, value in params.items()
            if value is not None
        }
        if query:
            url = url + "?" + urllib.parse.urlencode(query)
    return url


def _request(url: str, timeout: float) -> bytes:
    """GET ``url`` and return the raw response body.

    Parameters
    ----------
    url : str
        Full URL to fetch.
    timeout : float
        Seconds to wait for the server.

    Returns
    -------
    bytes
        The response body.

    Raises
    ------
    PEDError
        If PED responds with an HTTP error (``status`` is set), or cannot be
        reached at all (``status`` is ``None``).
    """
    request = urllib.request.Request(
        url, headers={"User-Agent": _user_agent(), "Accept": "*/*"}
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            body: bytes = response.read()
            return body
    except urllib.error.HTTPError as e:
        if e.code == 404:
            raise PEDError(
                f"PED has no such resource (HTTP 404): {url}", status=404
            ) from None
        raise PEDError(
            f"PED returned HTTP {e.code} ({e.reason}) for {url}", status=e.code
        ) from None
    except (urllib.error.URLError, TimeoutError, OSError) as e:
        reason = getattr(e, "reason", e)
        raise PEDError(
            f"Could not reach PED at {url} ({reason}). Check your internet "
            "connection; the PED server may also be temporarily unavailable."
        ) from None


def get_bytes(
    path: str,
    params: Optional[Mapping[str, QueryValue]] = None,
    timeout: float = DEFAULT_TIMEOUT,
) -> bytes:
    """GET a PED endpoint and return the raw response body.

    Parameters
    ----------
    path : str
        Endpoint path relative to :data:`PED_API_URL`.
    params : mapping, optional
        Query parameters (``None`` values are dropped).
    timeout : float, optional
        Seconds to wait for PED. Default :data:`DEFAULT_TIMEOUT`.

    Returns
    -------
    bytes
        The response body.

    Raises
    ------
    PEDError
        If the request fails.
    """
    return _request(build_url(path, params), timeout)


def get_json(
    path: str,
    params: Optional[Mapping[str, QueryValue]] = None,
    timeout: float = DEFAULT_TIMEOUT,
) -> Any:
    """GET a PED endpoint and decode the JSON response.

    Parameters
    ----------
    path : str
        Endpoint path relative to :data:`PED_API_URL`.
    params : mapping, optional
        Query parameters (``None`` values are dropped).
    timeout : float, optional
        Seconds to wait for PED. Default :data:`DEFAULT_TIMEOUT`.

    Returns
    -------
    Any
        The decoded JSON (usually a dict).

    Raises
    ------
    PEDError
        If the request fails or the response is not valid JSON.
    """
    url = build_url(path, params)
    body = _request(url, timeout)
    try:
        return json.loads(body)
    except (json.JSONDecodeError, UnicodeDecodeError) as e:
        raise PEDError(
            f"PED returned a response that is not valid JSON for {url}: {e}"
        ) from None
