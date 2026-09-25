"""Job request/result L0 schema and results.dat adapters.

Authoring home: monorepo ``schema/eon_job_result.capnp``.
This package vendors a copy for PyPI installs.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Mapping, Optional

from .trajectory import (
    TrajectoryManifestError,
    geometry_digest,
    landfold_consume,
    landfold_embedding_inputs,
    landfold_fes_inputs,
    potential_digest,
    rgpot_identity,
    trajectory_manifest,
    trajectory_manifest_dumps,
    trajectory_manifest_loads,
)

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


_SNAKE_TO_WIRE = {
    "job_id": "jobId",
    "job_type": "jobType",
    "status_code": "statusCode",
    "status_text": "statusText",
    "potential_type": "potentialType",
    "random_seed": "randomSeed",
    "potential_energy": "potentialEnergy",
    "potential_energy_saddle": "potentialEnergySaddle",
    "potential_energy_reactant": "potentialEnergyReactant",
    "potential_energy_product": "potentialEnergyProduct",
    "barrier_reactant_to_product": "barrierReactantToProduct",
    "barrier_product_to_reactant": "barrierProductToReactant",
    "prefactor_reactant_to_product": "prefactorReactantToProduct",
    "prefactor_product_to_reactant": "prefactorProductToReactant",
    "displacement_saddle_distance": "displacementSaddleDistance",
    "simulation_time": "simulationTime",
    "md_temperature": "mdTemperature",
    "has_dynamics": "hasDynamics",
    "wall_time_seconds": "wallTimeSeconds",
    "user_time_seconds": "userTimeSeconds",
    "system_time_seconds": "systemTimeSeconds",
    "client_version": "clientVersion",
    "landfold_artifacts": "landfoldArtifacts",
}
_ENGINE_SNAKE_TO_WIRE = {
    "engine_id": "engineId",
    "protocol_family": "protocolFamily",
    "protocol_major": "protocolMajor",
    "protocol_minor": "protocolMinor",
    "abi_major": "abiMajor",
    "abi_minor": "abiMinor",
    "layout_revision": "layoutRevision",
    "build_identity": "buildIdentity",
}
_ENGINE_WIRE_TO_SNAKE = {v: k for k, v in _ENGINE_SNAKE_TO_WIRE.items()}
_ARTIFACT_SNAKE_TO_WIRE = {
    "source_run_id": "sourceRunId",
    "input_digest": "inputDigest",
    "engine_compatibility": "engineCompatibility",
}
_ARTIFACT_WIRE_TO_SNAKE = {v: k for k, v in _ARTIFACT_SNAKE_TO_WIRE.items()}
_WIRE_TO_SNAKE = {v: k for k, v in _SNAKE_TO_WIRE.items()}


def _engine_to_wire(data: Mapping[str, Any]) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    for key, val in data.items():
        out[_ENGINE_SNAKE_TO_WIRE.get(key, key)] = val
    return out


def _engine_from_wire(data: Mapping[str, Any]) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    for key, val in data.items():
        out[_ENGINE_WIRE_TO_SNAKE.get(key, key)] = val
    return out


def _artifact_to_wire(data: Mapping[str, Any]) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    for key, val in data.items():
        dest = _ARTIFACT_SNAKE_TO_WIRE.get(key, key)
        if dest == "engineCompatibility" and isinstance(val, Mapping):
            out[dest] = _engine_to_wire(val)
        else:
            out[dest] = val
    return out


def _artifact_from_wire(data: Mapping[str, Any]) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    for key, val in data.items():
        dest = _ARTIFACT_WIRE_TO_SNAKE.get(key, key)
        if dest == "engine_compatibility" and isinstance(val, Mapping):
            out[dest] = _engine_from_wire(val)
        else:
            out[dest] = val
    return out


def _landfold_artifacts_to_wire(val: Any) -> list:
    if not isinstance(val, list):
        return []
    return [_artifact_to_wire(item) if isinstance(item, Mapping) else item for item in val]


def _landfold_artifacts_from_wire(val: Any) -> list:
    if not isinstance(val, list):
        return []
    return [_artifact_from_wire(item) if isinstance(item, Mapping) else item for item in val]


def job_result_to_wire(data: Mapping[str, Any]) -> Dict[str, Any]:
    """CamelCase dict matching JobResult field names."""
    wire: Dict[str, Any] = {}
    for key, val in data.items():
        if key == "force_calls" and isinstance(val, Mapping):
            wire["forceCalls"] = {
                "total": int(val.get("total", 0)),
                "minimization": int(val.get("minimization", 0)),
                "saddle": int(val.get("saddle", 0)),
                "prefactors": int(val.get("prefactors", 0)),
                "neb": int(val.get("neb", 0)),
            }
            continue
        if key in ("landfold_artifacts", "landfoldArtifacts"):
            wire["landfoldArtifacts"] = _landfold_artifacts_to_wire(val)
            continue
        dest = _SNAKE_TO_WIRE.get(key, key)
        wire[dest] = val
    return wire


def job_result_from_wire(data: Mapping[str, Any]) -> Dict[str, Any]:
    """Snake_case dict matching job_result_scalars_from_results_dat."""
    out: Dict[str, Any] = {}
    for key, val in data.items():
        if key == "forceCalls" and isinstance(val, Mapping):
            out["force_calls"] = {
                "total": int(val.get("total", 0)),
                "minimization": int(val.get("minimization", 0)),
                "saddle": int(val.get("saddle", 0)),
                "prefactors": int(val.get("prefactors", 0)),
                "neb": int(val.get("neb", 0)),
            }
            continue
        if key == "landfoldArtifacts":
            out["landfold_artifacts"] = _landfold_artifacts_from_wire(val)
            continue
        dest = _WIRE_TO_SNAKE.get(key, key)
        out[dest] = val
    return out


_ENGINE_TEXT = ("schema", "engineId", "protocolFamily", "buildIdentity")
_ENGINE_INT = (
    "protocolMajor",
    "protocolMinor",
    "abiMajor",
    "abiMinor",
    "layoutRevision",
)


def _fill_engine_compatibility(node: Any, data: Mapping[str, Any]) -> None:
    for name in _ENGINE_TEXT:
        val = data.get(name)
        if isinstance(val, str):
            setattr(node, name, val)
    for name in _ENGINE_INT:
        if name in data and data[name] is not None:
            setattr(node, name, int(data[name]))


def _fill_landfold_artifacts(msg: Any, artifacts: list) -> None:
    nodes = msg.init("landfoldArtifacts", len(artifacts))
    for index, art in enumerate(artifacts):
        if not isinstance(art, Mapping):
            continue
        node = nodes[index]
        for name in ("schema", "sourceRunId", "inputDigest"):
            val = art.get(name)
            if isinstance(val, str):
                setattr(node, name, val)
        compat = art.get("engineCompatibility")
        if isinstance(compat, Mapping):
            _fill_engine_compatibility(node.engineCompatibility, compat)


def job_result_dumps(data: Mapping[str, Any]) -> bytes:
    """Encode a JobResult dict. Uses pycapnp when installed, else JSON."""
    wire = job_result_to_wire(data)
    try:
        import capnp  # type: ignore
    except ImportError:
        import json

        return json.dumps(wire, separators=(",", ":")).encode("utf-8")
    schema = capnp.load(str(job_result_capnp_path()))
    msg = schema.JobResult.new_message()
    for key, val in wire.items():
        if key == "forceCalls" and isinstance(val, Mapping):
            fc = msg.forceCalls
            for sub, sval in val.items():
                if hasattr(fc, sub):
                    setattr(fc, sub, sval)
            continue
        if key == "landfoldArtifacts" and isinstance(val, list):
            _fill_landfold_artifacts(msg, val)
            continue
        if hasattr(msg, key):
            try:
                setattr(msg, key, val)
            except Exception:
                pass
    return msg.to_bytes_packed()


def job_result_loads(blob: bytes) -> Dict[str, Any]:
    """Decode bytes from job_result_dumps."""
    if blob[:1] == b"{":
        import json

        return job_result_from_wire(json.loads(blob.decode("utf-8")))
    try:
        import capnp  # type: ignore
    except ImportError as exc:
        raise ValueError("packed JobResult needs pycapnp") from exc
    schema = capnp.load(str(job_result_capnp_path()))
    msg = schema.JobResult.from_bytes_packed(blob)
    wire = msg.to_dict()
    return job_result_from_wire(wire)


__all__ = [
    "job_result_capnp_path",
    "results_dat_to_dict",
    "dict_to_results_dat",
    "job_result_scalars_from_results_dat",
    "job_result_to_wire",
    "job_result_from_wire",
    "job_result_dumps",
    "job_result_loads",
    "TrajectoryManifestError",
    "geometry_digest",
    "rgpot_identity",
    "potential_digest",
    "trajectory_manifest",
    "trajectory_manifest_dumps",
    "trajectory_manifest_loads",
    "landfold_embedding_inputs",
    "landfold_fes_inputs",
    "landfold_consume",
]
