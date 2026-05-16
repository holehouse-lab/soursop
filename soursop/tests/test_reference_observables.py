"""Regression tests against committed reference observable pickles.

For every ``(name, resolution)`` entry in
:data:`soursop.tests.build_reference.build_references.TRAJECTORIES`, this
module loads the stored reference pickle from
``soursop/data/test_data/<name>_<resolution>_reference.pkl`` and re-runs the
same builder used to produce it. Every leaf in the dictionary tree must match;
any drift means an observable's numerical behaviour has changed and either the
code change is incorrect, or the references need to be rebuilt by re-running
``build_references.py`` (only do this when a behavioural change is deliberate
and well-understood).

The tests intentionally reuse :func:`build_reference` itself rather than
duplicating the call signatures, so the test stays in sync with the build
script automatically.
"""

from __future__ import annotations

import math
import pickle
from typing import Any, Iterable

import numpy as np
import pytest

import soursop
from soursop.sstrajectory import SSTrajectory

from soursop.tests.build_reference.build_references import (
    TRAJECTORIES,
    build_reference,
)


# Loose enough to absorb FP-level drift across mdtraj/numpy/BLAS versions,
# tight enough to flag any real behavioural regression. Matches the
# tolerance used by the existing distance-map regression tests in this repo.
RTOL = 1e-5
ATOL = 1e-4

# Top-level meta keys that are environment-dependent and therefore exempt
# from the leaf comparison. Other meta keys (R1, R2, save_stride, ...) are
# compared because they parameterise the observables themselves.
META_VERSION_KEYS = {
    'soursop_version',
    'mdtraj_version',
    'numpy_version',
    'python_version',
}


# -- helpers ------------------------------------------------------------------


def _spec_id(spec: dict) -> str:
    return f"{spec['name']}_{spec['resolution']}"


def _load_protein(name: str, resolution: str):
    pdb = soursop.get_data(f'test_data/{name}_{resolution}.pdb')
    xtc = soursop.get_data(f'test_data/{name}_{resolution}.xtc')
    return SSTrajectory(xtc, pdb).proteinTrajectoryList[0]


def _load_stored_reference(name: str, resolution: str) -> dict:
    path = soursop.get_data(f'test_data/{name}_{resolution}_reference.pkl')
    with open(path, 'rb') as f:
        return pickle.load(f)


def _walk_leaves(d: dict, prefix: str = '') -> Iterable[tuple[str, Any]]:
    for key, value in d.items():
        path = f"{prefix}.{key}" if prefix else str(key)
        if isinstance(value, dict):
            yield from _walk_leaves(value, path)
        else:
            yield path, value


def _get_by_path(d: dict, path: str) -> Any:
    """Walk dotted ``path`` through nested dicts. Raises ``KeyError`` on miss.

    Dict keys are stringified in the path, but the actual key type might be
    ``int`` (BBSEG ``per_class`` uses 0..8 as integer keys). Try the string
    lookup first and fall back to ``int`` when the literal string isn't a
    member.
    """
    obj: Any = d
    for part in path.split('.'):
        if not isinstance(obj, dict):
            raise KeyError(f"cannot descend into {type(obj).__name__} at '{part}'")
        if part in obj:
            obj = obj[part]
            continue
        try:
            int_key = int(part)
        except ValueError:
            raise KeyError(part)
        if int_key in obj:
            obj = obj[int_key]
            continue
        raise KeyError(part)
    return obj


def _is_skipped_meta_path(path: str) -> bool:
    if not path.startswith('meta.'):
        return False
    return path.split('.', 1)[1] in META_VERSION_KEYS


def _compare_leaf(path: str, expected: Any, actual: Any) -> str | None:
    """Return None on match or a one-line error message on mismatch."""
    if expected is None:
        if actual is None:
            return None
        return f"[{path}] expected None, got {type(actual).__name__}={actual!r}"

    if isinstance(expected, np.ndarray):
        if not isinstance(actual, np.ndarray):
            return f"[{path}] expected ndarray, got {type(actual).__name__}"
        if expected.shape != actual.shape:
            return f"[{path}] shape {expected.shape} != {actual.shape}"
        if expected.dtype.kind in ('f', 'c'):
            try:
                np.testing.assert_allclose(actual, expected, rtol=RTOL, atol=ATOL)
            except AssertionError as e:
                return f"[{path}] allclose failed: {str(e).splitlines()[0]}"
            return None
        if not np.array_equal(actual, expected):
            return f"[{path}] integer/bool/object array values differ"
        return None

    if isinstance(expected, (bool, np.bool_)):
        if bool(actual) != bool(expected):
            return f"[{path}] bool {expected!r} != {actual!r}"
        return None

    if isinstance(expected, (int, np.integer)):
        if int(actual) != int(expected):
            return f"[{path}] int {expected!r} != {actual!r}"
        return None

    if isinstance(expected, float):
        if not math.isfinite(expected) and not math.isfinite(actual):
            if math.isnan(expected) == math.isnan(actual):
                return None
            return f"[{path}] non-finite float mismatch {expected!r} vs {actual!r}"
        if not math.isclose(actual, expected, rel_tol=RTOL, abs_tol=ATOL):
            return f"[{path}] float {expected!r} != {actual!r}"
        return None

    if isinstance(expected, str):
        if actual != expected:
            return f"[{path}] string {expected!r} != {actual!r}"
        return None

    if isinstance(expected, (list, tuple)):
        if not isinstance(actual, type(expected)):
            return (
                f"[{path}] expected {type(expected).__name__}, "
                f"got {type(actual).__name__}"
            )
        if len(actual) != len(expected):
            return f"[{path}] len {len(expected)} != {len(actual)}"
        for i, (e, a) in enumerate(zip(expected, actual)):
            err = _compare_leaf(f"{path}[{i}]", e, a)
            if err is not None:
                return err
        return None

    # Fallback for anything else (e.g., numpy scalars not caught above).
    if isinstance(expected, np.floating):
        return _compare_leaf(path, float(expected), float(actual))
    if isinstance(expected, np.integer):
        return _compare_leaf(path, int(expected), int(actual))

    if actual != expected:
        return f"[{path}] {type(expected).__name__} {expected!r} != {actual!r}"
    return None


# -- fixtures ------------------------------------------------------------------


@pytest.fixture(scope='module')
def stored_references() -> dict:
    """Map of ``(name, resolution) -> stored reference dict`` (from pickle)."""
    return {
        (spec['name'], spec['resolution']): _load_stored_reference(
            spec['name'], spec['resolution'],
        )
        for spec in TRAJECTORIES
    }


@pytest.fixture(scope='module')
def fresh_references() -> dict:
    """Map of ``(name, resolution) -> freshly-built reference dict``.

    Each entry is built by re-running :func:`build_reference` against the
    on-disk trajectory. Built once per pytest module to amortise the SASA
    and stochastic-fit cost.
    """
    out: dict = {}
    for spec in TRAJECTORIES:
        key = (spec['name'], spec['resolution'])
        protein = _load_protein(*key)
        out[key] = build_reference(
            protein, spec['name'], spec['resolution'], spec['R1'], spec['R2'],
        )
    return out


# -- tests ---------------------------------------------------------------------


@pytest.mark.parametrize('spec', TRAJECTORIES, ids=_spec_id)
def test_top_level_categories_match(spec, stored_references, fresh_references):
    """The set of top-level categories must be identical."""
    key = (spec['name'], spec['resolution'])
    stored_cats = set(stored_references[key].keys())
    fresh_cats = set(fresh_references[key].keys())
    assert fresh_cats == stored_cats, (
        f"top-level category set diverged for {_spec_id(spec)}: "
        f"only in stored = {stored_cats - fresh_cats}, "
        f"only in fresh = {fresh_cats - stored_cats}"
    )


@pytest.mark.parametrize('spec', TRAJECTORIES, ids=_spec_id)
def test_leaf_keys_match(spec, stored_references, fresh_references):
    """The set of leaf paths must be identical (no missing or extra keys)."""
    key = (spec['name'], spec['resolution'])
    stored_paths = {p for p, _ in _walk_leaves(stored_references[key])}
    fresh_paths = {p for p, _ in _walk_leaves(fresh_references[key])}
    only_stored = stored_paths - fresh_paths
    only_fresh = fresh_paths - stored_paths
    assert not only_stored and not only_fresh, (
        f"leaf-path set diverged for {_spec_id(spec)}:\n"
        f"  only in stored ({len(only_stored)}): {sorted(only_stored)[:10]}\n"
        f"  only in fresh  ({len(only_fresh)}): {sorted(only_fresh)[:10]}"
    )


def _collect_leaf_params() -> tuple[list, list[str]]:
    """Enumerate every (spec, leaf_path) pair to parametrize per-observable tests.

    Pytest parametrize collection happens at module-import time, so we walk the
    stored references right here. If a reference pickle is missing we skip its
    trajectory silently — the structural tests above will surface that case
    with a more useful error.
    """
    params: list = []
    ids: list[str] = []
    for spec in TRAJECTORIES:
        try:
            stored = _load_stored_reference(spec['name'], spec['resolution'])
        except FileNotFoundError:
            continue
        for path, _ in _walk_leaves(stored):
            if _is_skipped_meta_path(path):
                continue
            params.append((spec, path))
            ids.append(f"{_spec_id(spec)}-{path}")
    return params, ids


_LEAF_PARAMS, _LEAF_IDS = _collect_leaf_params()


@pytest.mark.parametrize('spec,path', _LEAF_PARAMS, ids=_LEAF_IDS)
def test_observable_value(spec, path, stored_references, fresh_references):
    """One test per observable per trajectory: stored value == freshly recomputed value."""
    key = (spec['name'], spec['resolution'])
    expected = _get_by_path(stored_references[key], path)
    try:
        actual = _get_by_path(fresh_references[key], path)
    except KeyError:
        # Structural divergence is reported by test_leaf_keys_match; skip
        # here so we don't flood the report with redundant failures.
        pytest.skip(
            f"path '{path}' is in stored reference but missing in fresh build; "
            "see test_leaf_keys_match for the structural diff"
        )
    err = _compare_leaf(path, expected, actual)
    if err is not None:
        raise AssertionError(err)
