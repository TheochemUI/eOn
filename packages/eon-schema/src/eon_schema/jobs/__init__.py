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


def job_result_legacy_dict(data: Mapping[str, Any]) -> Dict[str, Any]:
    """Historical results.dat keys for a JobResult mapping.

    Geometries (ConFrame / Geometry) are not serialized. Cluster and HPC
    adapters call :func:`job_result_to_results_dat`; in-process callers use
    this dict and never build a ``results.dat`` buffer.
    """
    fc = data.get("force_calls") or {}
    if not isinstance(fc, Mapping):
        fc = {}
    out: Dict[str, Any] = {}

    def put(src: str, dest: str) -> None:
        if src in data and data[src] is not None:
            out[dest] = data[src]

    put("status_code", "termination_reason")
    put("status_text", "termination_reason_text")
    put("job_type", "job_type")
    put("potential_type", "potential_type")
    if data.get("random_seed", -1) not in (None, -1):
        out["random_seed"] = data["random_seed"]
    if "total" in fc:
        out["total_force_calls"] = fc["total"]
    elif "force_calls" in data and not isinstance(data["force_calls"], Mapping):
        out["total_force_calls"] = data["force_calls"]
    for src, dest in (
        ("minimization", "force_calls_minimization"),
        ("saddle", "force_calls_saddle"),
        ("prefactors", "force_calls_prefactors"),
        ("neb", "force_calls_neb"),
        ("dephase", "force_calls_dephase"),
        ("dynamics", "force_calls_dynamics"),
        ("refine", "force_calls_refine"),
        ("sampling", "force_calls_sampling"),
    ):
        if src in fc:
            out[dest] = fc[src]
    for src, dest in (
        ("potential_energy", "potential_energy"),
        ("potential_energy_saddle", "potential_energy_saddle"),
        ("potential_energy_reactant", "potential_energy_reactant"),
        ("potential_energy_product", "potential_energy_product"),
        ("barrier_reactant_to_product", "barrier_reactant_to_product"),
        ("barrier_product_to_reactant", "barrier_product_to_reactant"),
        ("prefactor_reactant_to_product", "prefactor_reactant_to_product"),
        ("prefactor_product_to_reactant", "prefactor_product_to_reactant"),
        ("displacement_saddle_distance", "displacement_saddle_distance"),
        ("wall_time_seconds", "time_seconds"),
        ("user_time_seconds", "user_time"),
        ("system_time_seconds", "system_time"),
    ):
        if src in data and data[src] is not None:
            out[dest] = data[src]
    if data.get("has_dynamics"):
        put("simulation_time", "simulation_time")
        put("md_temperature", "md_temperature")
    extras = data.get("extras") or []
    if isinstance(extras, Mapping):
        for key, val in extras.items():
            out[str(key)] = val
    else:
        for item in extras:
            if isinstance(item, Mapping) and "key" in item:
                out[str(item["key"])] = item.get("value")
    return out


def job_result_to_results_dat(data: Mapping[str, Any]) -> str:
    """Serialize a JobResult to classic ``results.dat`` text.

    For cluster and HPC adapters only. The in-process path keeps the
    JobResult object and does not call this.
    """
    return dict_to_results_dat(job_result_legacy_dict(data))





_ENGINE_STAMP_FIELDS = (
    "compatibility_schema",
    "engine_id",
    "compatibility_engine_protocol_family",
    "compatibility_engine_protocol_major",
    "compatibility_engine_protocol_minor",
    "compatibility_engine_abi_major",
    "compatibility_engine_abi_minor",
    "compatibility_engine_layout_revision",
    "engine_build_identity",
)


def _engine_compatibility(d: Mapping[str, Any]) -> Optional[Dict[str, Any]]:
    """Strict eon.compatibility.v1 object. Unknown schemas stay un-upgraded."""
    if d.get("compatibility_schema") != "eon.compatibility.v1":
        return None
    if any(key not in d for key in _ENGINE_STAMP_FIELDS):
        return None
    stamp: Dict[str, Any] = {
        "schema": d["compatibility_schema"],
        "engineId": d["engine_id"],
        "protocolFamily": d["compatibility_engine_protocol_family"],
        "protocolMajor": d["compatibility_engine_protocol_major"],
        "protocolMinor": d["compatibility_engine_protocol_minor"],
        "abiMajor": d["compatibility_engine_abi_major"],
        "abiMinor": d["compatibility_engine_abi_minor"],
        "layoutRevision": d["compatibility_engine_layout_revision"],
        "buildIdentity": d["engine_build_identity"],
    }
    for source, target in (
        ("compatibility_readcon_spec_version", "readconSpecVersion"),
        ("compatibility_readcon_min_version", "readconMinVersion"),
        ("compatibility_eon_schema_min_version", "eonSchemaMinVersion"),
        ("compatibility_rgpycrumbs_min_version", "rgpycrumbsMinVersion"),
        ("compatibility_chemparseplot_min_version", "chemparseplotMinVersion"),
        ("rgpot_name", "rgpotName"),
        ("rgpot_version", "rgpotVersion"),
    ):
        if source in d:
            stamp[target] = d[source]
    return stamp


def _eindir_from_results(d: Mapping[str, Any]) -> Optional[Dict[str, Any]]:
    if "optimizer_eindir_abi_major" not in d:
        return None
    return {
        "abi_major": d.get("optimizer_eindir_abi_major", 0),
        "abi_minor": d.get("optimizer_eindir_abi_minor", 0),
        "objective_layout": d.get("optimizer_eindir_objective_layout", 0),
        "objective_size": d.get("optimizer_eindir_objective_size", 0),
        "objective_align": d.get("optimizer_eindir_objective_align", 0),
        "dlpack_major": d.get("optimizer_eindir_dlpack_major", 0),
        "dlpack_minor": d.get("optimizer_eindir_dlpack_minor", 0),
        "features": d.get("optimizer_eindir_features", 0),
    }


def _provenance_from_results(d: Mapping[str, Any]) -> Dict[str, Any]:
    eindir = _eindir_from_results(d)
    optimizer: Dict[str, Any] = {
        "schema": d.get("optimizer_provenance_schema", ""),
        "backend": d.get("optimizer_backend", ""),
        "xts_abi": {
            "major": d.get("optimizer_xts_abi_major", 0),
            "minor": d.get("optimizer_xts_abi_minor", 0),
            "layout": d.get("optimizer_xts_abi_layout", 0),
        },
        "has_eindir": eindir is not None,
    }
    if eindir is not None:
        optimizer["eindir"] = eindir
    rgpot_version = d.get("rgpot_version", "")
    rgpot = {
        "schema": d.get("rgpot_schema", ""),
        "name": d.get("rgpot_name", d.get("potential_type", "")),
        "version": rgpot_version if isinstance(rgpot_version, str) else str(rgpot_version),
    }
    compatibility: Dict[str, Any] = {
        "schema": d.get("compatibility_schema", ""),
        "readcon": {
            "spec_version": d.get("compatibility_readcon_spec_version", 0),
            "min_version": d.get("compatibility_readcon_min_version", ""),
        },
        "eon_schema": {"min_version": d.get("compatibility_eon_schema_min_version", "")},
        "rgpycrumbs": {"min_version": d.get("compatibility_rgpycrumbs_min_version", "")},
        "chemparseplot": {
            "min_version": d.get("compatibility_chemparseplot_min_version", "")
        },
        "engine": {
            "id": d.get("engine_id", d.get("potential_type", "")),
            "version": d.get("engine_version"),
            "build_identity": d.get("engine_build_identity"),
        },
        "rgpot": {"name": rgpot["name"], "version": rgpot["version"]},
    }
    stamp = _engine_compatibility(d)
    if stamp is not None:
        compatibility["engine_compatibility"] = stamp
    return {"optimizer": optimizer, "compatibility": compatibility, "rgpot": rgpot}


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
    out.update(_provenance_from_results(d))
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
        if key == "optimizer" and isinstance(val, Mapping):
            xts = val.get("xts_abi") or {}
            eindir = val.get("eindir") or {}
            wire["optimizer"] = {
                "schema": val.get("schema", ""),
                "backend": val.get("backend", ""),
                "xtsAbiMajor": int(xts.get("major", 0)),
                "xtsAbiMinor": int(xts.get("minor", 0)),
                "xtsAbiLayout": int(xts.get("layout", 0)),
                "hasEindir": bool(val.get("has_eindir", False)),
                "eindir": {
                    "abiMajor": int(eindir.get("abi_major", 0)),
                    "abiMinor": int(eindir.get("abi_minor", 0)),
                    "objectiveLayout": int(eindir.get("objective_layout", 0)),
                    "objectiveSize": int(eindir.get("objective_size", 0)),
                    "objectiveAlign": int(eindir.get("objective_align", 0)),
                    "dlpackMajor": int(eindir.get("dlpack_major", 0)),
                    "dlpackMinor": int(eindir.get("dlpack_minor", 0)),
                    "features": int(eindir.get("features", 0)),
                },
            }
            continue
        if key == "rgpot" and isinstance(val, Mapping):
            wire["rgpot"] = {
                "schema": val.get("schema", ""),
                "name": val.get("name", ""),
                "version": val.get("version", ""),
            }
            continue
        if key == "compatibility" and isinstance(val, Mapping):
            stamp = val.get("engine_compatibility")
            if isinstance(stamp, Mapping):
                wire["compatibility"] = dict(stamp)
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
        if key == "optimizer" and isinstance(val, Mapping):
            eindir = val.get("eindir") or {}
            opt = {
                "schema": val.get("schema", ""),
                "backend": val.get("backend", ""),
                "xts_abi": {
                    "major": int(val.get("xtsAbiMajor", 0)),
                    "minor": int(val.get("xtsAbiMinor", 0)),
                    "layout": int(val.get("xtsAbiLayout", 0)),
                },
                "has_eindir": bool(val.get("hasEindir", False)),
            }
            if opt["has_eindir"]:
                opt["eindir"] = {
                    "abi_major": int(eindir.get("abiMajor", 0)),
                    "abi_minor": int(eindir.get("abiMinor", 0)),
                    "objective_layout": int(eindir.get("objectiveLayout", 0)),
                    "objective_size": int(eindir.get("objectiveSize", 0)),
                    "objective_align": int(eindir.get("objectiveAlign", 0)),
                    "dlpack_major": int(eindir.get("dlpackMajor", 0)),
                    "dlpack_minor": int(eindir.get("dlpackMinor", 0)),
                    "features": int(eindir.get("features", 0)),
                }
            out["optimizer"] = opt
            continue
        if key == "rgpot" and isinstance(val, Mapping):
            out["rgpot"] = {
                "schema": val.get("schema", ""),
                "name": val.get("name", ""),
                "version": val.get("version", ""),
            }
            continue
        if key == "compatibility" and isinstance(val, Mapping):
            out["compatibility"] = {"engine_compatibility": dict(val)}
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
        if isinstance(val, Mapping):
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
    "job_result_legacy_dict",
    "job_result_to_results_dat",
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
