"""Unit tests for the inlined pure-math helpers in bin/transcript_qc_processing.py.

The script imports scanpy (and other heavy deps) at module top level, which is
not available in the test environment, so we cannot ``import`` it directly.
Instead we parse the source with ``ast``, extract the two pure functions by
name, and ``exec`` only those function definitions into a controlled namespace
that provides ``np`` (numpy) and ``pd`` (pandas). This exercises the REAL
inlined code without triggering the module-level heavy imports.
"""

import ast
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

SCRIPT_PATH = (
    Path(__file__).resolve().parent.parent / "bin" / "transcript_qc_processing.py"
)

FUNCTIONS_TO_EXTRACT = ("calculate_noise_bound", "estimate_min_mols_per_cell")


def _load_inlined_functions():
    """Extract the two target functions from the script source and exec them.

    Returns a namespace dict containing the compiled functions plus ``np``/``pd``.
    """
    source = SCRIPT_PATH.read_text()
    tree = ast.parse(source)

    namespace = {"np": np, "pd": pd}
    found = {}
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name in FUNCTIONS_TO_EXTRACT:
            segment = ast.get_source_segment(source, node)
            assert segment is not None, f"could not extract source for {node.name}"
            exec(segment, namespace)
            found[node.name] = namespace[node.name]

    missing = set(FUNCTIONS_TO_EXTRACT) - set(found)
    assert not missing, f"functions not found in script: {missing}"
    return namespace


@pytest.fixture(scope="module")
def funcs():
    ns = _load_inlined_functions()
    return {name: ns[name] for name in FUNCTIONS_TO_EXTRACT}


# ---------------------------------------------------------------------------
# estimate_min_mols_per_cell
# ---------------------------------------------------------------------------


def test_estimate_min_mols_returns_int(funcs):
    estimate_min_mols_per_cell = funcs["estimate_min_mols_per_cell"]
    rng = np.random.default_rng(0)
    # An obviously-high distribution centred around ~500 molecules/cell.
    data = rng.normal(loc=500, scale=50, size=5000).clip(min=1)
    result = estimate_min_mols_per_cell(data)
    assert isinstance(result, int)
    assert result >= 10


def test_estimate_min_mols_respects_min_value_floor(funcs):
    estimate_min_mols_per_cell = funcs["estimate_min_mols_per_cell"]
    # A tiny, low-count distribution: the computed threshold would be well
    # below the floor, so the floor must win.
    data = np.array([0, 0, 1, 1, 2])
    result = estimate_min_mols_per_cell(data, min_value=10)
    assert result == 10


def test_estimate_min_mols_custom_min_value(funcs):
    estimate_min_mols_per_cell = funcs["estimate_min_mols_per_cell"]
    data = np.array([0, 1, 2, 3])
    # A higher floor is respected on a distribution whose estimate is tiny.
    assert estimate_min_mols_per_cell(data, min_value=42) == 42


def test_estimate_min_mols_high_distribution_exceeds_floor(funcs):
    estimate_min_mols_per_cell = funcs["estimate_min_mols_per_cell"]
    rng = np.random.default_rng(1)
    # Bimodal: a big population of high-count cells makes the mode high enough
    # that the returned value comfortably exceeds the min_value floor.
    high = rng.normal(loc=1000, scale=30, size=10000).clip(min=1)
    result = estimate_min_mols_per_cell(high, min_value=10)
    assert result > 10


# ---------------------------------------------------------------------------
# calculate_noise_bound
# ---------------------------------------------------------------------------


def test_calculate_noise_bound_empty_series(funcs):
    calculate_noise_bound = funcs["calculate_noise_bound"]
    result = calculate_noise_bound(pd.Series([], dtype=float))
    assert result == (0, 0)


def test_calculate_noise_bound_nonempty_series(funcs):
    calculate_noise_bound = funcs["calculate_noise_bound"]
    rng = np.random.default_rng(2)
    # Per-feature molecule counts for negative-control probes.
    counts = pd.Series(rng.integers(low=5, high=500, size=200))
    lb, ub = calculate_noise_bound(counts)
    assert lb > 0
    assert ub > 0
    assert lb < ub


def test_calculate_noise_bound_quantile_widens_bounds(funcs):
    calculate_noise_bound = funcs["calculate_noise_bound"]
    rng = np.random.default_rng(3)
    counts = pd.Series(rng.integers(low=10, high=1000, size=300))
    lb_narrow, ub_narrow = calculate_noise_bound(counts, quant=0.90)
    lb_wide, ub_wide = calculate_noise_bound(counts, quant=0.999)
    # A higher quantile pushes the bounds further apart.
    assert ub_wide > ub_narrow
    assert lb_wide < lb_narrow
