"""Public package-root import surface (no mandatory alias)."""

from __future__ import annotations

import pytest

pytest.importorskip("pyeonclient")


def test_direct_imports_from_package_root():
    from pyeonclient import (  # noqa: F401
        Matter,
        NEB,
        NebSpec,
        Parameters,
        PathInit,
        append_timing,
        make_backend,
        pot_registry_total_force_calls,
        write_neb_results,
        write_potcall_summary,
        write_minimization_results,
        steady_clock_now,
        io_ok,
        neb_write_results,
        list_backends,
        NEBStatus,
    )
    assert callable(make_backend)
    assert callable(write_neb_results)
    assert callable(pot_registry_total_force_calls)
    assert "lj" in list_backends()
