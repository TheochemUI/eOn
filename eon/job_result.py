"""In-process JobResult consumption.

Cluster and HPC workers still ship ``results.dat``. The in-process path
attaches a JobResult mapping (scalars plus ConFrame geometries) and does
not require a ``StringIO`` of that file. ``parse_results`` stays the
reader for legacy file workers.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any, Mapping

_pkg_src = Path(__file__).resolve().parents[1] / "packages" / "eon-schema" / "src"
if _pkg_src.is_dir() and str(_pkg_src) not in sys.path:
    sys.path.insert(0, str(_pkg_src))

from eon_schema.jobs import job_result_legacy_dict  # noqa: E402


def results_mapping(result: Mapping[str, Any]) -> dict:
    """Legacy results dict from a job record.

    Prefers ``job_result``. Falls back to ``parse_results`` when the
    record only has ``results.dat`` (file communicators).
    """
    jr = result.get("job_result")
    if isinstance(jr, Mapping):
        return job_result_legacy_dict(jr)
    blob = result.get("results.dat")
    if blob is None:
        existing = result.get("results")
        if isinstance(existing, Mapping) and "termination_reason" in existing:
            return dict(existing)
        raise KeyError("results.dat")
    from eon import fileio as io

    return io.parse_results(blob)
