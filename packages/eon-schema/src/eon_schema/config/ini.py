"""INI helpers for eOn ``config.ini`` (shared by eon-akmc tooling + rgpycrumbs).

Priority surface for packages that **write** or **hydrate** job configs without
importing the full eon-akmc runtime:

* :func:`write_ini` — case-preserving dict → ``config.ini``
* :func:`defaults_from_catalog` / :func:`hydrate_ini` — L0 defaults + overrides
* :func:`unknown_ini_keys` — catch drift on covered L0 sections
* :func:`model_to_ini_section` — one L1 pydantic model → INI option map

L1 field names sometimes carry a section prefix (e.g. ``dimer_improved`` vs
INI ``improved``). Known renames live in :data:`INI_FIELD_ALIASES`.
"""

from __future__ import annotations

import configparser
from pathlib import Path
from typing import Any, Mapping, MutableMapping, Optional, Sequence, Union

from eon_schema.ssot import load_catalog

PathLike = Union[str, Path]
IniSections = dict[str, dict[str, Any]]

# (INI section, L1 pydantic field name) → INI option name
INI_FIELD_ALIASES: dict[tuple[str, str], str] = {
    ("Dimer", "dimer_improved"): "improved",
    ("Dimer", "dimer_max_iterations"): "max_iterations",
    ("Dimer", "dimer_rotations_min"): "rotations_min",
    ("Dimer", "dimer_torque_min"): "torque_min",
    ("Dimer", "dimer_torque_max"): "torque_max",
    ("Dimer", "dimer_remove_rotation"): "remove_rotation",
}

# L1 class name → INI section title (ConfigParser section)
MODEL_INI_SECTION: dict[str, str] = {
    "MainConfig": "Main",
    "PotentialConfig": "Potential",
    "OptimizerConfig": "Optimizer",
    "StructureComparisonConfig": "Structure Comparison",
    "ProcessSearchConfig": "Process Search",
    "SaddleSearchConfig": "Saddle Search",
    "DimerConfig": "Dimer",
    "NudgedElasticBandConfig": "Nudged Elastic Band",
    "LanczosConfig": "Lanczos",
    "DavidsonConfig": "Davidson",
    "ARTnConfig": "ARTn",
    "IRAConfig": "IRA",
    "PrefactorConfig": "Prefactor",
    "DynamicsConfig": "Dynamics",
    "ParallelReplicaConfig": "Parallel Replica",
    "HyperdynamicsConfig": "Hyperdynamics",
    "AKMCConfig": "AKMC",
    "BasinHoppingConfig": "Basin Hopping",
    "CommunicatorConfig": "Communicator",
    "PathsConfig": "Paths",
    "DebugConfig": "Debug",
    "KDBConfig": "KDB",
    "RecyclingConfig": "Recycling",
    "CoarseGrainingConfig": "Coarse Graining",
    "HessianConfig": "Hessian",
    "ServeConfig": "Serve",
    "RefineConfig": "Refine",
    "Metatomic": "Metatomic",
    "XTBPot": "XTBPot",
    "ZBLPot": "ZBLPot",
    "SocketNWChemPot": "SocketNWChemPot",
    "RgpotPot": "RgpotPot",
    "ASE_NWCHEM": "ASE_NWCHEM",
    "ASE_ORCA": "ASE_ORCA",
    "LBFGSConfig": "LBFGS",
    "CGConfig": "CG",
    "FIREConfig": "FIRE",
    "QuickMinConfig": "QuickMin",
    "SDConfig": "SD",
    "XtsciConfig": "Xtsci",
    "GPRDimerConfig": "GPR Dimer",
    "BGSDConfig": "BGSD",
    "DistributedReplicaConfig": "Distributed Replica",
    "AMSConfig": "AMS",
    "AMSIOConfig": "AMS_IO",
    "AMSEnvConfig": "AMSEnv",
}


def format_ini_value(value: Any) -> str:
    """Format a Python value for eOn ``config.ini`` (lowercase bools)."""
    if isinstance(value, bool):
        return "true" if value else "false"
    if value is None:
        return ""
    return str(value)


def write_ini(
    path: PathLike,
    sections: Mapping[str, Mapping[str, Any]],
    *,
    dir_ok: bool = True,
) -> Path:
    """Write ``sections`` to ``config.ini``, preserving option case.

    If *path* is an existing directory (or *dir_ok* and path has no suffix),
    write ``path/config.ini``. Values are stringified via :func:`format_ini_value`.
    """
    out = Path(path)
    if out.is_dir() or (dir_ok and out.suffix == "" and not out.exists()):
        out.mkdir(parents=True, exist_ok=True)
        out = out / "config.ini"
    elif out.is_dir():
        out = out / "config.ini"
    parser = configparser.ConfigParser()
    parser.optionxform = str  # type: ignore[method-assign, assignment]
    for section, options in sections.items():
        if not parser.has_section(section):
            parser.add_section(section)
        for key, val in options.items():
            parser.set(section, str(key), format_ini_value(val))
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as fh:
        parser.write(fh)
    return out


def read_ini(path: PathLike) -> IniSections:
    """Read ``config.ini`` into a nested dict (case-preserving keys)."""
    parser = configparser.ConfigParser()
    parser.optionxform = str  # type: ignore[method-assign, assignment]
    p = Path(path)
    if p.is_dir():
        p = p / "config.ini"
    read = parser.read(p)
    if not read:
        raise FileNotFoundError(f"config.ini not found: {path}")
    out: IniSections = {}
    for section in parser.sections():
        out[section] = dict(parser.items(section))
    return out


def _scalar_fields(section: str) -> list[dict[str, Any]]:
    cat = load_catalog()
    data = cat.get("sections", {}).get(section, {})
    return [f for f in data.get("fields", []) if not f.get("nested")]


def defaults_from_catalog(
    section: Optional[str] = None,
    *,
    include_flat_aliases: bool = True,
) -> IniSections | dict[str, Any]:
    """L0 defaults for one section (flat dict) or all covered sections.

    Keys use **snake** names from the Cap'n Proto catalog (same as INI for
    covered groups). Optimizer flat aliases (``lbfgs_memory``, …) are included
    when *include_flat_aliases* is true.
    """
    cat = load_catalog()
    if section is not None:
        out: dict[str, Any] = {
            f["snake"]: f["default"]
            for f in _scalar_fields(section)
            if f.get("default") is not None
        }
        if include_flat_aliases and section == "Optimizer":
            for a in cat.get("flat_aliases", []):
                if a.get("yaml_section") == "Optimizer" and a.get("default") is not None:
                    out[a["flat_key"]] = a["default"]
        for path, val in cat.get("server_yaml_default_overrides", {}).items():
            sec, key = path.split(".", 1)
            if sec == section:
                out[key] = val
        return out

    all_secs: IniSections = {}
    for sec in cat.get("sections", {}):
        # skip nested Optimizer.* sub-structs as top-level INI sections
        if sec.startswith("Optimizer."):
            continue
        d = defaults_from_catalog(sec, include_flat_aliases=include_flat_aliases)
        assert isinstance(d, dict)
        if d:
            all_secs[sec] = d
    return all_secs


def hydrate_ini(
    user: Mapping[str, Mapping[str, Any]],
    *,
    base_sections: Optional[Sequence[str]] = None,
) -> IniSections:
    """Merge L0 catalog defaults with *user* overrides (user wins).

    If *base_sections* is None, hydrate every top-level catalog section that
    has defaults, then layer user sections (including unknown sections).
    """
    if base_sections is None:
        base = defaults_from_catalog()
        assert isinstance(base, dict)
        # type: defaults_from_catalog() without section returns IniSections
        merged: IniSections = {
            sec: dict(opts) for sec, opts in base.items()  # type: ignore[arg-type]
        }
    else:
        merged = {}
        for sec in base_sections:
            d = defaults_from_catalog(sec)
            assert isinstance(d, dict)
            merged[sec] = dict(d)

    for sec, opts in user.items():
        bucket = merged.setdefault(sec, {})
        for k, v in opts.items():
            bucket[k] = v
    return merged


def allowed_keys(section: str) -> set[str]:
    """Known option names for a covered L0 section (snake + flat aliases)."""
    keys = {f["snake"] for f in _scalar_fields(section)}
    cat = load_catalog()
    if section == "Optimizer":
        for a in cat.get("flat_aliases", []):
            if a.get("yaml_section") == "Optimizer":
                keys.add(a["flat_key"])
    return keys


def unknown_ini_keys(
    settings: Mapping[str, Mapping[str, Any]],
    *,
    covered_only: bool = True,
) -> list[str]:
    """Return ``section.key`` paths not in the L0 catalog.

    Sections absent from the catalog are skipped when *covered_only* is True
    (so hand-written pot blocks like ``SocketNWChemPot`` are not flagged).
    """
    cat = load_catalog()
    known_sections = {
        s for s in cat.get("sections", {}) if not s.startswith("Optimizer.")
    }
    bad: list[str] = []
    for sec, opts in settings.items():
        if sec not in known_sections:
            if not covered_only:
                bad.append(f"{sec}.*")
            continue
        allowed = allowed_keys(sec)
        # Also accept nested Optimizer.* only via flat aliases under Optimizer
        for key in opts:
            if key not in allowed:
                bad.append(f"{sec}.{key}")
    return bad


def model_to_ini_section(
    model: Any,
    section: Optional[str] = None,
    *,
    exclude_none: bool = True,
) -> dict[str, str]:
    """Dump one L1 pydantic model to an INI option map (string values).

    Applies :data:`INI_FIELD_ALIASES` for known L1→INI renames.
    """
    name = type(model).__name__
    sec = section or MODEL_INI_SECTION.get(name)
    if not sec:
        raise ValueError(
            f"no INI section mapping for {name!r}; pass section= explicitly"
        )
    if hasattr(model, "model_dump"):
        data = model.model_dump()
    elif isinstance(model, Mapping):
        data = dict(model)
    else:
        raise TypeError(f"expected pydantic model, got {type(model)!r}")

    out: dict[str, str] = {}
    for key, val in data.items():
        if exclude_none and val is None:
            continue
        ini_key = INI_FIELD_ALIASES.get((sec, key), key)
        out[ini_key] = format_ini_value(val)
    return out


def models_to_ini(
    *models: Any,
    extra: Optional[Mapping[str, Mapping[str, Any]]] = None,
) -> IniSections:
    """Combine several L1 models into an INI section dict.

    Later models / *extra* override earlier keys on the same section.
    """
    sections: IniSections = {}
    for model in models:
        name = type(model).__name__
        sec = MODEL_INI_SECTION.get(name)
        if not sec:
            raise ValueError(f"no INI section mapping for {name!r}")
        opts = model_to_ini_section(model, sec)
        bucket = sections.setdefault(sec, {})
        bucket.update(opts)
    if extra:
        for sec, opts in extra.items():
            bucket = sections.setdefault(sec, {})
            for k, v in opts.items():
                bucket[k] = format_ini_value(v) if not isinstance(v, str) else v
    return sections


def write_models_ini(
    path: PathLike,
    *models: Any,
    extra: Optional[Mapping[str, Mapping[str, Any]]] = None,
    validate: bool = False,
) -> Path:
    """:func:`models_to_ini` then :func:`write_ini`.

    If *validate*, raise :class:`ValueError` when covered L0 sections contain
    unknown keys.
    """
    sections = models_to_ini(*models, extra=extra)
    if validate:
        bad = unknown_ini_keys(sections, covered_only=True)
        if bad:
            raise ValueError(f"unknown INI keys for covered sections: {bad}")
    return write_ini(path, sections)


__all__ = [
    "INI_FIELD_ALIASES",
    "MODEL_INI_SECTION",
    "IniSections",
    "allowed_keys",
    "defaults_from_catalog",
    "format_ini_value",
    "hydrate_ini",
    "model_to_ini_section",
    "models_to_ini",
    "read_ini",
    "unknown_ini_keys",
    "write_ini",
    "write_models_ini",
]
