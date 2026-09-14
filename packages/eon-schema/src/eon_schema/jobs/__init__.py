"""Job request/result L0 schema and results.dat adapters.

Authoring home: monorepo ``schema/eon_job_result.capnp``.
This package vendors a copy for PyPI installs.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Mapping, Optional

_PKG = Path(__file__).resolve().parent
_VEND_CAPNP = _PKG / "eon_job_result.capnp"
_MONO_CAPNP = Path(__file__).resolve().parents[5] / "schema" / "eon_job_result.capnp"


def job_result_capnp_path() -> Path:
    """Path to JobResult Cap'n Proto schema (vendored, else monorepo)."""
    if _VEND_CAPNP.is_file():
        return _VEND_CAPNP
    if _MONO_CAPNP.is_file():
        return _MONO_CAPNP
    raise FileNotFoundError("eon_job_result.capnp not found (package or monorepo)")


def results_dat_to_dict(text: str) -> Dict[str, Any]:
    """Parse classic ``results.dat`` lines into a dict (legacy adapter)."""
    results: Dict[str, Any] = {}
    for line in text.splitlines():
        parts = line.split()
        if len(parts) < 2:
            continue
        key = parts[1]
        raw = parts[0]
        if "." in raw or "e" in raw.lower() or "E" in raw:
            try:
                results[key] = float(raw)
                continue
            except ValueError:
                pass
        try:
            results[key] = int(raw)
        except ValueError:
            results[key] = raw
    return results


def dict_to_results_dat(data: Mapping[str, Any]) -> str:
    """Serialize a dict to classic ``results.dat`` text (legacy adapter)."""
    lines = []
    for key, val in data.items():
        if isinstance(val, float):
            lines.append(f"{val:.12e} {key}")
        else:
            lines.append(f"{val} {key}")
    return "\n".join(lines) + ("\n" if lines else "")


def job_result_scalars_from_results_dat(text: str) -> Dict[str, Any]:
    """Map results.dat keys into JobResult-oriented scalar names.

    Geometries are not present in results.dat; load .con separately into
    Geometry / ConFrame.
    """
    d = results_dat_to_dict(text)
    out: Dict[str, Any] = {
        "status_code": d.get("termination_reason", 0),
        "status_text": d.get("termination_reason_text", ""),
        "job_type": d.get("job_type", ""),
        "potential_type": d.get("potential_type", ""),
        "random_seed": d.get("random_seed", -1),
        "potential_energy": d.get("potential_energy", d.get("Energy", 0.0)),
        "potential_energy_saddle": d.get("potential_energy_saddle", 0.0),
        "potential_energy_reactant": d.get("potential_energy_reactant", 0.0),
        "potential_energy_product": d.get("potential_energy_product", 0.0),
        "barrier_reactant_to_product": d.get("barrier_reactant_to_product", 0.0),
        "barrier_product_to_reactant": d.get("barrier_product_to_reactant", 0.0),
        "prefactor_reactant_to_product": d.get("prefactor_reactant_to_product", 0.0),
        "prefactor_product_to_reactant": d.get("prefactor_product_to_reactant", 0.0),
        "displacement_saddle_distance": d.get("displacement_saddle_distance", 0.0),
        "force_calls": {
            "total": d.get("total_force_calls", d.get("force_calls", 0)),
            "minimization": d.get("force_calls_minimization", 0),
            "saddle": d.get("force_calls_saddle", 0),
            "prefactors": d.get("force_calls_prefactors", 0),
            "neb": d.get("force_calls_neb", 0),
        },
        "wall_time_seconds": d.get("time_seconds", 0.0),
        "user_time_seconds": d.get("user_time", 0.0),
        "system_time_seconds": d.get("system_time", 0.0),
    }
    if "simulation_time" in d:
        out["simulation_time"] = d["simulation_time"]
        out["md_temperature"] = d.get("md_temperature", 0.0)
        out["has_dynamics"] = True
    return out


__all__ = [
    "job_result_capnp_path",
    "results_dat_to_dict",
    "dict_to_results_dat",
    "job_result_scalars_from_results_dat",
]
