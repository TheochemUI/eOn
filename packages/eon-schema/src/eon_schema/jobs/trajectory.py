"""Run-level trajectory manifest for Landfold.

Geometries are hashed as little-endian IEEE-754 bytes. Coordinates are stored
as ``float.hex`` text so a JSON round trip does not change the digest.
``results.dat`` lines are not an input.
"""
from __future__ import annotations

import hashlib
import json
import math
import struct
from typing import Any, Dict, List, Mapping, Optional, Sequence

TRAJECTORY_SCHEMA = "eon.trajectory.v1"
RGPOT_SCHEMA = "eon.rgpot.v1"
LANDFOLD_SCHEMA = "eon.landfold.v1"
EMBEDDING_SCHEMA = "landfold.embedding.input.v1"
FES_SCHEMA = "landfold.fes.input.v1"

# eOn PotType tokens whose force kernel is an rgpot type. Names match
# Potential.cpp (RgpotAdapter / makeRgpot), not results.dat free text.
_KERNELS = {
    "lj": "LJPot",
    "ljcluster": "LJClusterPot",
    "morse_pt": "MorsePot",
    "cuh2": "CuH2Pot",
    "tip4p_h": "WaterHPot",
    "eam_al": "EAMAlPot",
    "edip": "EDIPPot",
    "fehe": "FeHePot",
    "lenosky_si": "LenoskyPot",
    "sw_si": "SWPot",
    "tersoff_si": "TersoffPot",
    "dftd3": "D3Pot",
    "dftd4": "D4Pot",
    "zbl": "ZBLPot",
    "mopac": "MOPACPot",
    "expr": "ExprPot",
    "rgpot": "RgpotPot",
}
_BACKEND_ALIASES = {
    "nwchem": "nwchemc",
    "nwchemc": "nwchemc",
    "cpmd": "cpmdc",
    "cpmdc": "cpmdc",
    "mta": "metatomic",
    "metatomic": "metatomic",
    "gfn": "xtb",
    "xtb": "xtb",
}
class TrajectoryManifestError(ValueError):
    """The manifest is not a Landfold trajectory input."""


def _reject_results_dat(value: Any) -> None:
    if isinstance(value, str):
        raise TrajectoryManifestError(
            "trajectory manifests are not built from results.dat lines"
        )


def _f64(value: Any, what: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TrajectoryManifestError(f"{what} must be a finite number")
    number = float(value)
    if not math.isfinite(number):
        raise TrajectoryManifestError(f"{what} must be finite")
    return number


def _hex_f64(value: Any, what: str) -> str:
    if isinstance(value, str):
        try:
            number = float.fromhex(value)
        except ValueError as exc:
            raise TrajectoryManifestError(f"{what} is not float.hex") from exc
        if not math.isfinite(number):
            raise TrajectoryManifestError(f"{what} must be finite")
        return float(number).hex()
    return _f64(value, what).hex()


def _from_hex(text: str, what: str) -> float:
    try:
        number = float.fromhex(text)
    except ValueError as exc:
        raise TrajectoryManifestError(f"{what} is not float.hex") from exc
    if not math.isfinite(number):
        raise TrajectoryManifestError(f"{what} must be finite")
    return number


def _u32(n: int) -> bytes:
    return struct.pack("<I", n)


def _texts(values: Sequence[str]) -> bytes:
    blob = _u32(len(values))
    for item in values:
        raw = item.encode("utf-8")
        blob += _u32(len(raw)) + raw
    return blob


def geometry_bytes(geometry: Mapping[str, Any]) -> bytes:
    """Canonical little-endian encoding of one Geometry."""
    _reject_results_dat(geometry)
    if not isinstance(geometry, Mapping):
        raise TrajectoryManifestError("geometry must be a mapping")
    raw_positions = geometry.get("positions")
    if not isinstance(raw_positions, Sequence) or isinstance(raw_positions, (str, bytes)):
        raise TrajectoryManifestError("positions must be a sequence")
    positions = [_f64(v, "position") for v in raw_positions]
    if not positions or len(positions) % 3 != 0:
        raise TrajectoryManifestError("positions length must be a positive multiple of 3")
    n_atoms = len(positions) // 3
    blob = b"EONGEO1" + _u32(len(positions)) + struct.pack("<" + "d" * len(positions), *positions)

    box = geometry.get("box")
    if box is None:
        blob += b"\x00"
    else:
        if not isinstance(box, Sequence) or isinstance(box, (str, bytes)) or len(box) != 9:
            raise TrajectoryManifestError("box must have length 9")
        values = [_f64(v, "box") for v in box]
        blob += b"\x01" + struct.pack("<9d", *values)

    atomic = geometry.get("atomic_numbers")
    if atomic is None:
        blob += _u32(0)
    else:
        if not isinstance(atomic, Sequence) or isinstance(atomic, (str, bytes)):
            raise TrajectoryManifestError("atomic_numbers must be a sequence")
        if len(atomic) != n_atoms:
            raise TrajectoryManifestError("atomic_numbers length must match positions")
        numbers = []
        for z in atomic:
            if isinstance(z, bool) or not isinstance(z, int) or not 0 <= z <= 65535:
                raise TrajectoryManifestError("atomic number out of uint16 range")
            numbers.append(z)
        blob += _u32(len(numbers)) + struct.pack("<" + "H" * len(numbers), *numbers)

    frozen = geometry.get("frozen")
    if frozen is None:
        blob += _u32(0)
    else:
        if not isinstance(frozen, Sequence) or isinstance(frozen, (str, bytes)):
            raise TrajectoryManifestError("frozen must be a sequence")
        if len(frozen) != n_atoms:
            raise TrajectoryManifestError("frozen length must match positions")
        flags = []
        for flag in frozen:
            if not isinstance(flag, bool):
                raise TrajectoryManifestError("frozen flags must be bool")
            flags.append(1 if flag else 0)
        blob += _u32(len(flags)) + bytes(flags)

    masses = geometry.get("masses")
    if masses is None:
        blob += _u32(0)
    else:
        if not isinstance(masses, Sequence) or isinstance(masses, (str, bytes)):
            raise TrajectoryManifestError("masses must be a sequence")
        if len(masses) != n_atoms:
            raise TrajectoryManifestError("masses length must match positions")
        values = [_f64(v, "mass") for v in masses]
        blob += _u32(len(values)) + struct.pack("<" + "d" * len(values), *values)

    symbols = geometry.get("symbols")
    if symbols is None:
        blob += _u32(0)
    else:
        if not isinstance(symbols, Sequence) or isinstance(symbols, (str, bytes)):
            raise TrajectoryManifestError("symbols must be a sequence")
        if len(symbols) != n_atoms:
            raise TrajectoryManifestError("symbols length must match positions")
        text = []
        for sym in symbols:
            if not isinstance(sym, str):
                raise TrajectoryManifestError("symbols must be strings")
            text.append(sym)
        blob += _texts(text)

    has_energy = bool(geometry.get("has_energy", "energy" in geometry and geometry.get("energy") is not None))
    if not has_energy:
        blob += b"\x00"
    else:
        blob += b"\x01" + struct.pack("<d", _f64(geometry.get("energy"), "energy"))
    return blob


def geometry_digest(geometry: Mapping[str, Any]) -> str:
    """sha256 of :func:`geometry_bytes`."""
    digest = hashlib.sha256(geometry_bytes(geometry)).hexdigest()
    return "sha256:" + digest


def _canon_config(value: Any) -> Any:
    if isinstance(value, float):
        if not math.isfinite(value):
            raise TrajectoryManifestError("potential config must be finite")
        return {"__f64__": value.hex()}
    if isinstance(value, bool) or value is None or isinstance(value, (str, int)):
        if isinstance(value, int) and not isinstance(value, bool) and abs(value) > 2**53:
            raise TrajectoryManifestError("potential config integer is not exact")
        return value
    if isinstance(value, Mapping):
        return {str(k): _canon_config(value[k]) for k in sorted(value, key=str)}
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        return [_canon_config(item) for item in value]
    raise TrajectoryManifestError("potential config has an unsupported value")


def _normalize_potential_type(potential_type: str) -> str:
    token = potential_type.strip().lower()
    if token not in _KERNELS:
        raise TrajectoryManifestError(
            f"{potential_type!r} is not an rgpot potential identity"
        )
    return token


def rgpot_identity(
    potential_type: str,
    *,
    backend: Optional[str] = None,
    revision: Optional[str] = None,
    config: Optional[Mapping[str, Any]] = None,
) -> Dict[str, Any]:
    """Structured rgpot kernel identity. Non-rgpot PotType names are rejected."""
    if not isinstance(potential_type, str) or not potential_type.strip():
        raise TrajectoryManifestError("potential_type is required")
    token = _normalize_potential_type(potential_type)
    if token == "rgpot":
        if not isinstance(backend, str) or not backend.strip():
            raise TrajectoryManifestError("rgpot identity requires a backend")
        resolved = _BACKEND_ALIASES.get(backend.strip().lower())
        if resolved is None:
            raise TrajectoryManifestError(f"unknown rgpot backend {backend!r}")
    else:
        if backend not in (None, ""):
            raise TrajectoryManifestError("backend is only valid for potential_type rgpot")
        resolved = ""
    rev = "" if revision is None else str(revision).strip().lower()
    if rev and (len(rev) != 40 or any(c not in "0123456789abcdef" for c in rev)):
        raise TrajectoryManifestError("rgpot revision must be 40 hex digits")
    identity: Dict[str, Any] = {
        "schema": RGPOT_SCHEMA,
        "provider": "rgpot",
        "potential_type": token,
        "kernel": _KERNELS[token],
        "backend": resolved,
        "revision": rev,
    }
    if config:
        identity["config"] = _canon_config(config)
    return identity


def _potential_bytes(identity: Mapping[str, Any]) -> bytes:
    payload = json.dumps(identity, sort_keys=True, separators=(",", ":"))
    return b"EONRGP1" + payload.encode("utf-8")


def potential_digest(identity: Mapping[str, Any]) -> str:
    """sha256 of the canonical rgpot identity."""
    return "sha256:" + hashlib.sha256(_potential_bytes(identity)).hexdigest()


def _frame_geometry(frame: Mapping[str, Any]) -> Dict[str, Any]:
    geo: Dict[str, Any] = {}
    positions = frame.get("positions", frame.get("positions_hex"))
    if positions is None:
        raise TrajectoryManifestError("frame is missing positions")
    geo["positions"] = [_from_hex(_hex_f64(v, "position"), "position") for v in positions]
    if "box" in frame or "box_hex" in frame:
        box = frame.get("box", frame.get("box_hex"))
        if box is not None:
            geo["box"] = [_from_hex(_hex_f64(v, "box"), "box") for v in box]
    if frame.get("atomic_numbers") is not None:
        geo["atomic_numbers"] = list(frame["atomic_numbers"])
    if frame.get("frozen") is not None:
        geo["frozen"] = [bool(v) if isinstance(v, bool) else v for v in frame["frozen"]]
    if "masses" in frame or "masses_hex" in frame:
        masses = frame.get("masses", frame.get("masses_hex"))
        if masses is not None:
            geo["masses"] = [_from_hex(_hex_f64(v, "mass"), "mass") for v in masses]
    if frame.get("symbols") is not None:
        geo["symbols"] = list(frame["symbols"])
    if frame.get("has_energy") or "energy" in frame or "energy_hex" in frame:
        energy = frame.get("energy", frame.get("energy_hex"))
        if energy is not None and frame.get("has_energy", True):
            geo["has_energy"] = True
            geo["energy"] = _from_hex(_hex_f64(energy, "energy"), "energy")
    return geo


def _stored_frame(index: int, frame: Mapping[str, Any]) -> Dict[str, Any]:
    if not isinstance(frame, Mapping):
        raise TrajectoryManifestError("each frame must be a mapping")
    geo = _frame_geometry(frame)
    n_atoms = len(geo["positions"]) // 3
    stored: Dict[str, Any] = {
        "index": index,
        "frame_id": str(frame.get("frame_id", index)),
        "role": str(frame.get("role", "image")),
        "n_atoms": n_atoms,
        "digest": geometry_digest(geo),
        "positions_hex": [v.hex() for v in geo["positions"]],
    }
    if "box" in geo:
        stored["box_hex"] = [v.hex() for v in geo["box"]]
    if "atomic_numbers" in geo:
        stored["atomic_numbers"] = [int(z) for z in geo["atomic_numbers"]]
    if "frozen" in geo:
        stored["frozen"] = [bool(v) for v in geo["frozen"]]
    if "masses" in geo:
        stored["masses_hex"] = [v.hex() for v in geo["masses"]]
    if "symbols" in geo:
        stored["symbols"] = [str(s) for s in geo["symbols"]]
    if geo.get("has_energy"):
        stored["has_energy"] = True
        stored["energy_hex"] = geo["energy"].hex()
    else:
        stored["has_energy"] = False
    return stored


def _ordered_geometry_digest(frame_digests: Sequence[str]) -> str:
    blob = b"EONTRAJ1\n" + ("\n".join(frame_digests) + "\n").encode("utf-8")
    return "sha256:" + hashlib.sha256(blob).hexdigest()


def trajectory_manifest(
    frames: Sequence[Mapping[str, Any]],
    potential: Mapping[str, Any] | str,
    *,
    source_run_id: str,
    length_unit: str = "angstrom",
    energy_unit: str = "eV",
) -> Dict[str, Any]:
    """Build an ``eon.trajectory.v1`` manifest from geometries and an rgpot identity."""
    _reject_results_dat(frames)
    if isinstance(potential, str):
        identity = rgpot_identity(potential)
    elif isinstance(potential, Mapping):
        if potential.get("schema") == RGPOT_SCHEMA and potential.get("kernel"):
            identity = rgpot_identity(
                str(potential.get("potential_type", "")),
                backend=potential.get("backend") or None,
                revision=potential.get("revision") or None,
                config=potential.get("config"),
            )
        else:
            identity = rgpot_identity(
                str(potential.get("potential_type", "")),
                backend=potential.get("backend"),
                revision=potential.get("revision"),
                config=potential.get("config"),
            )
    else:
        raise TrajectoryManifestError("potential must be an rgpot identity")
    if not isinstance(source_run_id, str) or not source_run_id.strip():
        raise TrajectoryManifestError("source_run_id is required")
    if not isinstance(length_unit, str) or not length_unit.strip():
        raise TrajectoryManifestError("length_unit is required")
    if not isinstance(energy_unit, str) or not energy_unit.strip():
        raise TrajectoryManifestError("energy_unit is required")
    if not isinstance(frames, Sequence) or isinstance(frames, (str, bytes)) or not frames:
        raise TrajectoryManifestError("at least one frame is required")
    stored = [_stored_frame(i, frame) for i, frame in enumerate(frames)]
    manifest = {
        "schema": TRAJECTORY_SCHEMA,
        "source_run_id": source_run_id,
        "length_unit": length_unit,
        "energy_unit": energy_unit,
        "potential": identity,
        "potential_digest": potential_digest(identity),
        "geometry_digest": _ordered_geometry_digest([f["digest"] for f in stored]),
        "frames": stored,
    }
    return _verify(manifest)


def trajectory_manifest_dumps(manifest: Mapping[str, Any]) -> str:
    """JSON text of a verified manifest."""
    checked = _verify(manifest)
    return json.dumps(checked, sort_keys=True, separators=(",", ":")) + "\n"


def trajectory_manifest_loads(text: str) -> Dict[str, Any]:
    """Parse manifest JSON. ``results.dat`` text is rejected."""
    if not isinstance(text, str):
        raise TrajectoryManifestError("manifest text must be a string")
    stripped = text.lstrip()
    if not stripped.startswith("{"):
        raise TrajectoryManifestError(
            "trajectory manifests are JSON, not results.dat lines"
        )
    try:
        data = json.loads(stripped)
    except json.JSONDecodeError as exc:
        raise TrajectoryManifestError("manifest is not JSON") from exc
    if not isinstance(data, dict):
        raise TrajectoryManifestError("manifest must be a JSON object")
    return _verify(data)


def _verify(manifest: Mapping[str, Any]) -> Dict[str, Any]:
    if not isinstance(manifest, Mapping):
        raise TrajectoryManifestError("manifest must be a mapping")
    if manifest.get("schema") != TRAJECTORY_SCHEMA:
        raise TrajectoryManifestError("manifest schema is not eon.trajectory.v1")
    frames = manifest.get("frames")
    if not isinstance(frames, list) or not frames:
        raise TrajectoryManifestError("manifest has no frames")
    stored = []
    for index, frame in enumerate(frames):
        rebuilt = _stored_frame(index, frame)
        claimed = frame.get("digest")
        if claimed != rebuilt["digest"]:
            raise TrajectoryManifestError(
                f"frame {index} digest does not match its geometry"
            )
        if frame.get("frame_id") is not None:
            rebuilt["frame_id"] = str(frame["frame_id"])
        if frame.get("role") is not None:
            rebuilt["role"] = str(frame["role"])
        if int(frame.get("index", index)) != index:
            raise TrajectoryManifestError("frame index is not the trajectory order")
        if int(frame.get("n_atoms", rebuilt["n_atoms"])) != rebuilt["n_atoms"]:
            raise TrajectoryManifestError("frame n_atoms does not match positions")
        stored.append(rebuilt)
    identity = manifest.get("potential")
    if not isinstance(identity, Mapping):
        raise TrajectoryManifestError("manifest is missing an rgpot identity")
    checked_id = rgpot_identity(
        str(identity.get("potential_type", "")),
        backend=identity.get("backend") or None,
        revision=identity.get("revision") or None,
        config=identity.get("config"),
    )
    if checked_id.get("kernel") != identity.get("kernel"):
        raise TrajectoryManifestError("rgpot kernel does not match potential_type")
    if identity.get("provider") not in (None, "rgpot"):
        raise TrajectoryManifestError("potential provider is not rgpot")
    if identity.get("schema") not in (None, RGPOT_SCHEMA):
        raise TrajectoryManifestError("potential schema is not eon.rgpot.v1")
    digest = potential_digest(checked_id)
    claimed_pot = manifest.get("potential_digest")
    if claimed_pot is not None and claimed_pot != digest:
        raise TrajectoryManifestError("potential digest does not match the rgpot identity")
    geo = _ordered_geometry_digest([f["digest"] for f in stored])
    claimed_geo = manifest.get("geometry_digest")
    if claimed_geo is not None and claimed_geo != geo:
        raise TrajectoryManifestError("geometry digest does not match the frames")
    source = manifest.get("source_run_id")
    if not isinstance(source, str) or not source.strip():
        raise TrajectoryManifestError("source_run_id is required")
    length_unit = manifest.get("length_unit", "angstrom")
    energy_unit = manifest.get("energy_unit", "eV")
    if not isinstance(length_unit, str) or not length_unit.strip():
        raise TrajectoryManifestError("length_unit is required")
    if not isinstance(energy_unit, str) or not energy_unit.strip():
        raise TrajectoryManifestError("energy_unit is required")
    return {
        "schema": TRAJECTORY_SCHEMA,
        "source_run_id": source,
        "length_unit": length_unit,
        "energy_unit": energy_unit,
        "potential": checked_id,
        "potential_digest": digest,
        "geometry_digest": geo,
        "frames": stored,
    }


def _coordinates(frame: Mapping[str, Any]) -> List[float]:
    return [_from_hex(v, "position") for v in frame["positions_hex"]]


def landfold_embedding_inputs(manifest: Mapping[str, Any]) -> Dict[str, Any]:
    """Frame-ordered coordinates for a Landfold embedding."""
    checked = _verify(manifest)
    counts = {frame["n_atoms"] for frame in checked["frames"]}
    if len(counts) != 1:
        raise TrajectoryManifestError("embedding requires one atom count for every frame")
    return {
        "schema": EMBEDDING_SCHEMA,
        "source_run_id": checked["source_run_id"],
        "length_unit": checked["length_unit"],
        "n_atoms": checked["frames"][0]["n_atoms"],
        "frame_ids": [frame["frame_id"] for frame in checked["frames"]],
        "digests": [frame["digest"] for frame in checked["frames"]],
        "coordinates": [_coordinates(frame) for frame in checked["frames"]],
        "geometry_digest": checked["geometry_digest"],
        "potential": checked["potential"],
        "potential_digest": checked["potential_digest"],
    }


def landfold_fes_inputs(manifest: Mapping[str, Any]) -> Dict[str, Any]:
    """Per-frame energies for a Landfold free-energy sample."""
    checked = _verify(manifest)
    if any(not frame["has_energy"] for frame in checked["frames"]):
        raise TrajectoryManifestError("FES inputs need an energy on every frame")
    return {
        "schema": FES_SCHEMA,
        "source_run_id": checked["source_run_id"],
        "energy_unit": checked["energy_unit"],
        "frame_ids": [frame["frame_id"] for frame in checked["frames"]],
        "digests": [frame["digest"] for frame in checked["frames"]],
        "energies": [_from_hex(frame["energy_hex"], "energy") for frame in checked["frames"]],
        "geometry_digest": checked["geometry_digest"],
        "potential": checked["potential"],
        "potential_digest": checked["potential_digest"],
    }


def landfold_consume(manifest: Mapping[str, Any]) -> Dict[str, Any]:
    """Verify the manifest and return embedding plus FES inputs.

    FES is omitted when any frame has no energy. Embedding still requires a
    single atom count.
    """
    embedding = landfold_embedding_inputs(manifest)
    checked = _verify(manifest)
    fes: Optional[Dict[str, Any]]
    if all(frame["has_energy"] for frame in checked["frames"]):
        fes = landfold_fes_inputs(checked)
    else:
        fes = None
    return {
        "schema": LANDFOLD_SCHEMA,
        "source_run_id": checked["source_run_id"],
        "geometry_digest": checked["geometry_digest"],
        "potential": checked["potential"],
        "potential_digest": checked["potential_digest"],
        "embedding": embedding,
        "fes": fes,
    }
