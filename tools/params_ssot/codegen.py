#!/usr/bin/env python3
"""Generate consumer projections from Cap'n Proto L0 SSoT.

Authoring SSoT: schema/eon_params.capnp (monorepo root).

The PyPI split ``eon-schema`` vendors a copy under
``packages/eon-schema/src/eon_schema/ssot/`` for standalone installs;
codegen refreshes that catalog JSON in lockstep. Fat release tarballs
(``git archive`` of the full monorepo) still ship everything for
conda-forge ``eon-feedstock`` — layout cleanup is fine; the release
contract is one fat archive + optional split PyPI packages.

A struct tagged ``# project: ini,json`` is also projected into the INI
reader, the JSON reader and writer, and ``apply_ssot_defaults``. The tag
and the alias, overlay, and fallback lines above that struct are the
adapter declaration. Structs without the tag stay in the handwritten
adapters.

Outputs (relative to repo root):
  schema/eon_params_catalog.json
  packages/eon-schema/.../ssot/eon_params_catalog.json  (vendored for split)
  eon/_params_ssot_catalog.py
  include/eon/generated/ParametersSSOTDefaults.h
  include/eon/generated/ParametersSSOTFieldIndex.inc
  include/eon/generated/ParametersSSOTIni.inc
  include/eon/generated/ParametersSSOTJson.inc
  include/eon/generated/ParametersSSOTApply.inc
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
CAPNP = REPO / "schema" / "eon_params.capnp"

SECTION_MAP = {
    "MainOptions": "Main",
    "PotentialOptions": "Potential",
    "StructureComparisonOptions": "Structure Comparison",
    "ProcessSearchOptions": "Process Search",
    "OptimizerOptions": "Optimizer",
    "OptimizerLbfgsOptions": "Optimizer.LBFGS",
    "OptimizerCgOptions": "Optimizer.CG",
    "OptimizerQuickminOptions": "Optimizer.Quickmin",
    "OptimizerSdOptions": "Optimizer.SD",
    "OptimizerXtsciOptions": "Optimizer.Xtsci",
}

# Nested SSoT path → flat config.ini / config.yaml option name (historical)
# Used only for Optimizer nested knobs that clients write under [Optimizer] flat
# or [LBFGS]/[CG] sections with prefixed keys.
FLAT_ALIASES: dict[str, str] = {
    "Optimizer.LBFGS.memory": "lbfgs_memory",
    "Optimizer.LBFGS.inverse_curvature": "lbfgs_inverse_curvature",
    "Optimizer.LBFGS.max_inverse_curvature": "lbfgs_max_inverse_curvature",
    "Optimizer.LBFGS.auto_scale": "lbfgs_auto_scale",
    "Optimizer.LBFGS.angle_reset": "lbfgs_angle_reset",
    "Optimizer.LBFGS.distance_reset": "lbfgs_distance_reset",
    "Optimizer.LBFGS.curvature": "lbfgs_curvature",
    "Optimizer.LBFGS.project_rigid": "lbfgs_project_rigid",
    "Optimizer.LBFGS.secant": "lbfgs_secant",
    "Optimizer.LBFGS.precon": "lbfgs_precon",
    "Optimizer.LBFGS.h0": "lbfgs_h0",
    "Optimizer.LBFGS.accept": "lbfgs_accept",
    "Optimizer.LBFGS.extra_updates": "lbfgs_extra_updates",
    "Optimizer.LBFGS.cautious_eps": "lbfgs_cautious_eps",
    "Optimizer.LBFGS.cautious_alpha": "lbfgs_cautious_alpha",
    "Optimizer.LBFGS.precon_A": "lbfgs_precon_A",
    "Optimizer.LBFGS.precon_mu": "lbfgs_precon_mu",
    "Optimizer.LBFGS.precon_rcut": "lbfgs_precon_rcut",
    "Optimizer.LBFGS.step": "lbfgs_step",
    "Optimizer.CG.no_overshooting": "cg_no_overshooting",
    "Optimizer.CG.knock_out_max_move": "cg_knock_out_max_move",
    "Optimizer.CG.line_search": "cg_line_search",
    "Optimizer.CG.line_converged": "cg_line_converged",
    "Optimizer.CG.line_search_max_iter": "cg_max_iter_line_search",
    "Optimizer.CG.max_iter_before_reset": "cg_max_iter_before_reset",
    "Optimizer.Quickmin.steepest_descent": "qm_steepest_descent",
    "Optimizer.SD.alpha": "sd_alpha",
    "Optimizer.SD.two_point": "sd_two_point",
    "Optimizer.Xtsci.method": "xtsci_method",
    "Optimizer.Xtsci.qn_step": "xtsci_qn_step",
    "Optimizer.Xtsci.precon": "xtsci_precon",
    "Optimizer.Xtsci.accept": "xtsci_accept",
    "Optimizer.Xtsci.highs": "xtsci_highs",
    "Optimizer.Xtsci.manifold": "xtsci_manifold",
}

FIELD_SNAKE = {
    "randomSeed": "random_seed",
    "writeLog": "write_log",
    "iniFilename": "ini_filename",
    "conFilename": "con_filename",
    "finiteDifference": "finite_difference",
    "maxForceCalls": "max_force_calls",
    "removeNetForce": "remove_net_force",
    "writeConForces": "write_con_forces",
    "mpiPollPeriod": "mpi_poll_period",
    "lammpsLogging": "lammps_logging",
    "lammpsThreads": "lammps_threads",
    "emtRasmussen": "emt_rasmussen",
    "logPotential": "log_potential",
    "extPotPath": "ext_pot_path",
    "potentialsPath": "potentials_path",
    "distanceDifference": "distance_difference",
    "neighborCutoff": "neighbor_cutoff",
    "checkRotation": "check_rotation",
    "indistinguishableAtoms": "indistinguishable_atoms",
    "energyDifference": "energy_difference",
    "removeTranslation": "remove_translation",
    "useCovalent": "use_covalent",
    "covalentScale": "covalent_scale",
    "bruteNeighbors": "brute_neighbors",
    "minimizeFirst": "minimize_first",
    "minimizationOffset": "minimization_offset",
    "optMethod": "opt_method",
    "convergenceMetric": "convergence_metric",
    "maxIterations": "max_iterations",
    "maxMove": "max_move",
    "convergedForce": "converged_force",
    "timeStep": "time_step",
    "maxTimeStep": "max_time_step",
    "inverseCurvature": "inverse_curvature",
    "maxInverseCurvature": "max_inverse_curvature",
    "autoScale": "auto_scale",
    "angleReset": "angle_reset",
    "distanceReset": "distance_reset",
    "projectRigid": "project_rigid",
    "extraUpdates": "extra_updates",
    "cautiousEps": "cautious_eps",
    "cautiousAlpha": "cautious_alpha",
    "preconA": "precon_A",
    "preconMu": "precon_mu",
    "preconRcut": "precon_rcut",
    "noOvershooting": "no_overshooting",
    "knockOutMaxMove": "knock_out_max_move",
    "lineSearch": "line_search",
    "lineConverged": "line_converged",
    "lineSearchMaxIter": "line_search_max_iter",
    "maxIterBeforeReset": "max_iter_before_reset",
    "steepestDescent": "steepest_descent",
    "twoPoint": "two_point",
    "schemaVersion": "schema_version",
}


def to_snake(name: str) -> str:
    if name in FIELD_SNAKE:
        return FIELD_SNAKE[name]
    s1 = re.sub(r"(.)([A-Z][a-z]+)", r"\1_\2", name)
    return re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", s1).lower()


# Comment directives are lowercase keys. Prose comments in the schema are
# sentences and do not match.
_DIRECTIVE_RE = re.compile(r"^#\s*([a-z][a-z0-9-]*)\s*:\s*(.*?)\s*$")
_STRUCT_RE = re.compile(r"^struct\s+(\w+)\s*\{$")
_FIELD_RE = re.compile(
    r"^(\w+)\s+@(\d+)\s*:\s*([A-Za-z0-9_.]+)(?:\s*=\s*([^;]+))?;\s*$"
)
_IDENT_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
_CAST_RE = re.compile(r"^static_cast<[A-Za-z_][A-Za-z0-9_:<>\s]*>$")

_GETTER = {
    "Bool": "GetBoolean",
    "Text": "Get",
    "Float64": "GetReal",
    "Float32": "GetReal",
    "Int32": "GetInteger",
    "Int64": "GetInteger",
    "UInt32": "GetInteger",
    "UInt64": "GetInteger",
    "Int16": "GetInteger",
    "UInt16": "GetInteger",
}
_CPP_TYPE = {
    "Bool": "bool",
    "Text": "std::string",
    "Float64": "double",
    "Float32": "float",
    "Int32": "int",
    "Int64": "long",
    "UInt32": "unsigned",
    "UInt64": "unsigned long",
    "Int16": "short",
    "UInt16": "unsigned short",
}


def _ident(value: str, what: str) -> str:
    if not _IDENT_RE.match(value):
        raise SystemExit(f"refusing non-identifier {what}: {value!r}")
    return value


def _parse_overlay(text: str) -> dict:
    parts: dict[str, str] = {}
    for tok in text.split():
        if "=" not in tok:
            raise SystemExit(f"overlay token needs key=value: {tok!r}")
        key, val = tok.split("=", 1)
        parts[key] = val
    if "section" not in parts:
        raise SystemExit(f"overlay needs section=: {text!r}")
    overlay: dict = {"section": _ident(parts["section"], "overlay section")}
    if "when" in parts:
        kind, _, rest = parts["when"].partition(":")
        if kind != "backend" or not rest:
            raise SystemExit(f"overlay when must be backend:names, got {parts['when']!r}")
        names = [_ident(n.strip(), "backend") for n in rest.split(",") if n.strip()]
        if not names:
            raise SystemExit(f"overlay when has no backends: {parts['when']!r}")
        overlay["when"] = "backend:" + ",".join(names)
    if "fields" in parts and "map" in parts:
        raise SystemExit("overlay uses fields= or map=, not both")
    if "fields" in parts:
        overlay["fields"] = [
            _ident(n.strip(), "overlay field")
            for n in parts["fields"].split(",")
            if n.strip()
        ]
    elif "map" in parts:
        mapped = []
        for item in parts["map"].split(","):
            member, sep, key = item.partition(":")
            if not sep:
                raise SystemExit(f"overlay map item needs member:key: {item!r}")
            mapped.append(
                {
                    "member": _ident(member.strip(), "overlay member"),
                    "key": _ident(key.strip(), "overlay key"),
                }
            )
        overlay["map"] = mapped
    else:
        raise SystemExit(f"overlay needs fields= or map=: {text!r}")
    return overlay


def _consume_directives(pending: list[tuple[str, str]]) -> dict:
    out: dict = {}
    overlays = []
    for key, val in pending:
        if key == "overlay":
            overlays.append(_parse_overlay(val))
        elif key in {"ini-order", "overlay-order"}:
            out[key] = [
                _ident(n.strip(), key) for n in val.split(",") if n.strip()
            ]
        elif key == "project":
            kinds = [n.strip() for n in val.split(",") if n.strip()]
            unknown = [n for n in kinds if n not in {"ini", "json"}]
            if unknown or not kinds:
                raise SystemExit(f"project must be ini,json, got {val!r}")
            out["project"] = ",".join(kinds)
        elif key == "cxx-cast":
            if not _CAST_RE.match(val):
                raise SystemExit(f"refusing cxx-cast {val!r}")
            out["cxx-cast"] = val
        elif key == "ini-when":
            if not re.fullmatch(r"potential=[A-Za-z0-9_]+", val):
                raise SystemExit(f"ini-when must be potential=ENUM, got {val!r}")
            out["ini-when"] = val
        elif key in {"section", "accessor", "json", "ini-fallback-member"}:
            out[key] = _ident(val, key)
        else:
            raise SystemExit(f"unknown schema directive {key}")
    if overlays:
        out["overlays"] = overlays
    return out


def parse_default(raw: str | None, type_name: str):
    if raw is None:
        return None
    raw = raw.strip()
    if type_name == "Bool":
        return raw.lower() == "true"
    if type_name in ("Float64", "Float32"):
        return float(raw)
    if type_name in ("Int64", "Int32", "UInt64", "UInt32", "Int16", "UInt16"):
        return int(raw)
    if type_name == "Text":
        if raw.startswith('"') and raw.endswith('"'):
            return raw[1:-1]
        return raw
    return raw


def _strip_trailing_comment(line: str) -> str:
    if "#" not in line or line.lstrip().startswith("#"):
        return line.strip()
    head = line.split("#", 1)[0]
    if '"' in head:
        return line.strip()
    return head.strip()


def parse_capnp(text: str) -> dict:
    """Parse structs. Comment directives bind to the next struct or field.

    A blank line or a prose comment drops a pending directive, so a tag
    has to sit on the declaration it describes.
    """
    structs: dict[str, dict] = {}
    pending: list[tuple[str, str]] = []
    current: str | None = None

    def take() -> dict:
        nonlocal pending
        meta = _consume_directives(pending)
        pending = []
        return meta

    for raw in text.splitlines():
        stripped = raw.strip()
        if not stripped:
            pending.clear()
            continue
        if stripped.startswith("#"):
            matched = _DIRECTIVE_RE.match(stripped)
            if matched:
                pending.append((matched.group(1), matched.group(2)))
            else:
                pending.clear()
            continue
        code = _strip_trailing_comment(raw)
        struct_m = _STRUCT_RE.match(code)
        if struct_m:
            if current is not None:
                raise SystemExit(f"nested struct is not supported: {struct_m.group(1)}")
            current = struct_m.group(1)
            structs[current] = {"fields": [], "directives": take()}
            continue
        if code == "}":
            if current is None:
                raise SystemExit("closing brace outside a struct")
            current = None
            pending.clear()
            continue
        field_m = _FIELD_RE.match(code)
        if field_m and current is not None:
            fname, ord_, typ, default = field_m.groups()
            nested = typ[0].isupper() and typ.endswith("Options") and default is None
            meta = take()
            field = {
                "name": fname,
                "ordinal": int(ord_),
                "type": typ,
                "default": None if nested else parse_default(default, typ),
                "snake": to_snake(fname),
            }
            if nested:
                field["nested"] = True
            if "ini-order" in meta:
                field["ini_order"] = meta["ini-order"]
            if "overlay-order" in meta:
                field["overlay_order"] = meta["overlay-order"]
            if "ini-fallback-member" in meta:
                field["ini_fallback_member"] = meta["ini-fallback-member"]
            if "cxx-cast" in meta:
                field["cxx_cast"] = meta["cxx-cast"]
            structs[current]["fields"].append(field)
            continue
        pending.clear()
    if current is not None:
        raise SystemExit(f"unclosed struct {current}")
    return structs


def _catalog_field(field: dict) -> dict:
    entry = {
        "snake": field["snake"],
        "capnp": field["name"],
        "ordinal": field["ordinal"],
        "type": field["type"],
        "default": field["default"],
    }
    if field.get("nested"):
        entry["nested"] = True
    if field.get("ini_order"):
        entry["ini_order"] = list(field["ini_order"])
    if field.get("overlay_order"):
        entry["overlay_order"] = list(field["overlay_order"])
    if field.get("ini_fallback_member"):
        entry["ini_fallback_member"] = field["ini_fallback_member"]
    if field.get("cxx_cast"):
        entry["cxx_cast"] = field["cxx_cast"]
    return entry


def build_catalog(structs: dict) -> dict:
    sections = {}
    for sname, struct in structs.items():
        if sname == "EonParameters":
            continue
        directives = struct["directives"]
        sec = directives.get("section") or SECTION_MAP.get(sname)
        if not sec:
            continue
        entry = {
            "struct": sname,
            "fields": [_catalog_field(f) for f in struct["fields"]],
        }
        if "project" in directives:
            if "accessor" not in directives:
                raise SystemExit(f"{sname} project needs accessor:")
            entry["project"] = [
                part.strip() for part in directives["project"].split(",") if part.strip()
            ]
            entry["accessor"] = directives["accessor"]
            if "json" in directives:
                entry["json_name"] = directives["json"]
            if "ini-when" in directives:
                entry["ini_when"] = directives["ini-when"]
            if directives.get("overlays"):
                entry["overlays"] = directives["overlays"]
        sections[sec] = entry

    # Flat aliases for Optimizer nested → historical yaml/ini option names
    flat_aliases = []
    for path, flat in FLAT_ALIASES.items():
        # path like Optimizer.LBFGS.memory
        parts = path.split(".")
        sec = ".".join(parts[:-1])  # Optimizer.LBFGS
        snake = parts[-1]
        default = None
        if sec in sections:
            for f in sections[sec]["fields"]:
                if f["snake"] == snake and not f.get("nested"):
                    default = f["default"]
                    break
        flat_aliases.append(
            {
                "flat_key": flat,
                "ssot_path": path,
                "yaml_section": "Optimizer",
                "default": default,
            }
        )

    # Server config.yaml dispatcher defaults that intentionally differ from
    # client runtime Parameters defaults for the same SSoT field.
    server_yaml_default_overrides = {
        "Main.job": "akmc",  # eon-server default job; client default process_search
        "Potential.ext_pot_path": "ext_pot",  # yaml relative path without ./
    }

    return {
        "source": "schema/eon_params.capnp",
        "schema_version": 1,
        "sections": sections,
        "flat_aliases": flat_aliases,
        "server_yaml_default_overrides": server_yaml_default_overrides,
    }


def emit_python(catalog: dict) -> str:
    body = json.dumps(catalog, indent=2, sort_keys=True)
    py_body = (
        body.replace(": true", ": True")
        .replace(": false", ": False")
        .replace(": null", ": None")
    )
    return (
        '"""AUTO-GENERATED from schema/eon_params.capnp — do not edit.\n\n'
        "Regenerate: python tools/params_ssot/codegen.py\n"
        '"""\n'
        "from __future__ import annotations\n\n"
        f"CATALOG = {py_body}\n"
    )


def emit_cpp_header(catalog: dict) -> str:
    lines = [
        "// AUTO-GENERATED from schema/eon_params.capnp — do not edit.",
        "// Regenerate: python tools/params_ssot/codegen.py",
        "#pragma once",
        "#include <cstdint>",
        "#include <string_view>",
        "namespace eonc::params_ssot {",
        "struct GeneratedDefaults {",
    ]

    def cpp_val(v, typ: str) -> str:
        if v is None:
            return "/* nested */"
        if typ == "Bool":
            return "true" if v else "false"
        if typ == "Text":
            return f'std::string_view{{"{v}"}}'
        if typ in ("Float64", "Float32"):
            return f"{float(v)}"
        return str(int(v))

    for sec, data in catalog["sections"].items():
        safe = re.sub(r"[^A-Za-z0-9]+", "_", sec)
        for f in data["fields"]:
            if f.get("nested") or f["default"] is None:
                continue
            cname = f"{safe}_{f['snake']}".upper()
            lines.append(
                f"  static constexpr auto {cname} = {cpp_val(f['default'], f['type'])};"
            )
    lines.append("};")
    lines.append("} // namespace eonc::params_ssot")
    lines.append("")
    return "\n".join(lines)


def emit_field_index(catalog: dict) -> str:
    """C++ initializer list for unordered_set field_index()."""
    ids = []
    for sec, data in catalog["sections"].items():
        for f in data["fields"]:
            if f.get("nested"):
                continue
            ids.append(f"{sec}.{f['snake']}")
            for alias in f.get("ini_order") or []:
                if alias != f["snake"]:
                    ids.append(f"{sec}.{alias}")
    # also flat aliases under Optimizer for ssot_has_field convenience
    for a in catalog.get("flat_aliases", []):
        ids.append(f"Optimizer.{a['flat_key']}")
    ids = sorted(set(ids))
    lines = [
        "// AUTO-GENERATED from schema/eon_params.capnp — do not edit.",
        "// Include inside a std::unordered_set<std::string> initializer.",
    ]
    for i in ids:
        lines.append(f'    "{i}",')
    lines.append("")
    return "\n".join(lines)


def _scalar_fields(section: dict) -> list:
    return [f for f in section["fields"] if not f.get("nested")]


def _projected(catalog: dict) -> list[tuple[str, dict]]:
    return [
        (name, data)
        for name, data in catalog["sections"].items()
        if data.get("project")
    ]


def _const_name(section: str, field: dict) -> str:
    safe = re.sub(r"[^A-Za-z0-9]+", "_", section)
    return f"{safe}_{field['snake']}".upper()


def _ini_expr(section_expr: str, field: dict, default_expr: str, keys: list[str]) -> str:
    getter = _GETTER.get(field["type"])
    if getter is None:
        raise SystemExit(f"no ini getter for {field['capnp']} ({field['type']})")
    expr = default_expr
    for key in reversed(keys):
        expr = f'ini.{getter}({section_expr}, "{key}", {expr})'
    if field.get("cxx_cast"):
        expr = f"{field['cxx_cast']}({expr})"
    return expr


def _member_default(field: dict, by_capnp: dict) -> str:
    member = field.get("ini_fallback_member")
    if not member:
        return f"o.{field['capnp']}"
    other = by_capnp.get(member)
    if other is None:
        raise SystemExit(f"{field['capnp']} falls back to unknown field {member}")
    if other["type"] == field["type"]:
        return f"o.{member}"
    cpp = _CPP_TYPE.get(field["type"])
    if cpp is None:
        raise SystemExit(f"no C++ type for {field['type']}")
    return f"static_cast<{cpp}>(o.{member})"


def _assign(field: dict, expr: str) -> str:
    return f"o.{field['capnp']} = {expr};"


def _ini_statements(section: str, data: dict) -> list[str]:
    fields = _scalar_fields(data)
    by_capnp = {f["capnp"]: f for f in fields}
    lines = [
        f"auto &o = ParametersLoadAccess::{data['accessor']}(params);",
        f'const char *sec = "{section}";',
    ]
    ordered = [f for f in fields if not f.get("ini_fallback_member")]
    ordered.extend(f for f in fields if f.get("ini_fallback_member"))
    for field in ordered:
        keys = field.get("ini_order") or [field["snake"]]
        expr = _ini_expr("sec", field, _member_default(field, by_capnp), keys)
        lines.append(_assign(field, expr))
    overlays = data.get("overlays") or []
    if any(str(item.get("when", "")).startswith("backend:") for item in overlays):
        lines.append("const std::string be = toLowerCase(o.backend);")
    for overlay in overlays:
        chunk = []
        section_expr = f'"{overlay["section"]}"'
        if "map" in overlay:
            pairs = [(item["member"], [item["key"]]) for item in overlay["map"]]
        else:
            pairs = []
            for member in overlay["fields"]:
                field = by_capnp[member]
                keys = (
                    field.get("overlay_order")
                    or field.get("ini_order")
                    or [field["snake"]]
                )
                pairs.append((member, keys))
        for member, keys in pairs:
            field = by_capnp[member]
            expr = _ini_expr(section_expr, field, f"o.{field['capnp']}", keys)
            chunk.append(_assign(field, expr))
        when = overlay.get("when")
        if when:
            names = when.split(":", 1)[1].split(",")
            cond = " || ".join(f'be == "{name}"' for name in names)
            lines.append(f"if ({cond}) {{")
            lines.extend(f"  {row}" for row in chunk)
            lines.append("}")
        else:
            lines.extend(chunk)
    return lines


def emit_ini(catalog: dict) -> str:
    lines = [
        "// AUTO-GENERATED from schema/eon_params.capnp — do not edit.",
        "// Regenerate: python tools/params_ssot/codegen.py",
        "// Included inside namespace eonc::config, after toLowerCase.",
        "inline void project_ssot_ini(INIReader &ini, Parameters &params) {",
    ]
    for section, data in _projected(catalog):
        if "ini" not in data["project"]:
            continue
        body = _ini_statements(section, data)
        when = data.get("ini_when")
        if when:
            enum = when.split("=", 1)[1]
            lines.append(
                "  if (params.potential_options().potential == "
                f"PotType::{enum}) {{"
            )
            lines.extend(f"    {row}" for row in body)
            lines.append("  }")
        else:
            lines.extend(f"  {row}" for row in body)
    lines.append("}")
    lines.append("")
    return "\n".join(lines)


def _json_arm(field: dict, key: str, keyword: str) -> str:
    member = field["capnp"]
    return (
        f'{keyword} (s.contains("{key}")) '
        f'o.{member} = s.at("{key}").get<decltype(o.{member})>();'
    )


def emit_json(catalog: dict) -> str:
    lines = [
        "// AUTO-GENERATED from schema/eon_params.capnp — do not edit.",
        "// Regenerate: python tools/params_ssot/codegen.py",
        "// Included inside namespace eonc::config, after JSON_OPT.",
        "inline void project_ssot_json_write(json &j, const Parameters &p) {",
    ]
    for section, data in _projected(catalog):
        if "json" not in data["project"]:
            continue
        json_name = data.get("json_name") or section
        fields = _scalar_fields(data)
        lines.append("  {")
        lines.append(
            f"    const auto &o = ParametersLoadAccess::{data['accessor']}(p);"
        )
        lines.append(f'    j["{json_name}"] = {{')
        for index, field in enumerate(fields):
            comma = "," if index + 1 < len(fields) else ""
            lines.append(
                f'        {{"{field["snake"]}", o.{field["capnp"]}}}{comma}'
            )
        lines.append("    };")
        lines.append("  }")
    lines.append("}")
    lines.append("")
    lines.append("inline void project_ssot_json_read(const json &j, Parameters &p) {")
    for section, data in _projected(catalog):
        if "json" not in data["project"]:
            continue
        json_name = data.get("json_name") or section
        fields = _scalar_fields(data)
        by_capnp = {field["capnp"]: field for field in fields}
        ordered = [f for f in fields if not f.get("ini_fallback_member")]
        ordered.extend(f for f in fields if f.get("ini_fallback_member"))
        lines.append(f'  if (j.contains("{json_name}")) {{')
        lines.append(f'    auto &s = j.at("{json_name}");')
        lines.append(f"    auto &o = ParametersLoadAccess::{data['accessor']}(p);")
        for field in ordered:
            keys = field.get("ini_order") or [field["snake"]]
            for index, key in enumerate(keys):
                lines.append(
                    "    " + _json_arm(field, key, "if" if index == 0 else "else if")
                )
            if field.get("ini_fallback_member"):
                lines.append(
                    f"    else o.{field['capnp']} = "
                    f"{_member_default(field, by_capnp)};"
                )
        lines.append("  }")
    lines.append("}")
    lines.append("")
    return "\n".join(lines)


def emit_apply(catalog: dict) -> str:
    lines = [
        "// AUTO-GENERATED from schema/eon_params.capnp — do not edit.",
        "// Regenerate: python tools/params_ssot/codegen.py",
        "// Included at the end of apply_ssot_defaults. GD is in scope.",
    ]
    for section, data in _projected(catalog):
        lines.append("  {")
        lines.append(f"    auto &o = ParametersLoadAccess::{data['accessor']}(p);")
        for field in _scalar_fields(data):
            cname = _const_name(section, field)
            if field["type"] == "Text":
                lines.append(
                    f"    o.{field['capnp']} = std::string(GD::{cname});"
                )
            else:
                lines.append(f"    o.{field['capnp']} = GD::{cname};")
        lines.append("  }")
    if not _projected(catalog):
        lines.append("  // No struct is tagged project.")
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    text = CAPNP.read_text(encoding="utf-8")
    structs = parse_capnp(text)
    if "MainOptions" not in structs:
        print("error: failed to parse MainOptions from", CAPNP, file=sys.stderr)
        return 1
    catalog = build_catalog(structs)

    out_json_tree = REPO / "schema" / "eon_params_catalog.json"
    out_json_pkg = (
        REPO
        / "packages"
        / "eon-schema"
        / "src"
        / "eon_schema"
        / "ssot"
        / "eon_params_catalog.json"
    )
    out_py = REPO / "eon" / "_params_ssot_catalog.py"
    out_h = REPO / "include" / "eon" / "generated" / "ParametersSSOTDefaults.h"
    out_idx = REPO / "include" / "eon" / "generated" / "ParametersSSOTFieldIndex.inc"
    out_ini = REPO / "include" / "eon" / "generated" / "ParametersSSOTIni.inc"
    out_json_inc = REPO / "include" / "eon" / "generated" / "ParametersSSOTJson.inc"
    out_apply = REPO / "include" / "eon" / "generated" / "ParametersSSOTApply.inc"
    out_h.parent.mkdir(parents=True, exist_ok=True)
    out_json_pkg.parent.mkdir(parents=True, exist_ok=True)

    catalog_json = json.dumps(catalog, indent=2, sort_keys=True) + "\n"
    out_json_tree.write_text(catalog_json)
    out_json_pkg.write_text(catalog_json)
    out_py.write_text(emit_python(catalog))
    out_h.write_text(emit_cpp_header(catalog))
    out_idx.write_text(emit_field_index(catalog))
    out_ini.write_text(emit_ini(catalog))
    out_json_inc.write_text(emit_json(catalog))
    out_apply.write_text(emit_apply(catalog))
    print("read", CAPNP.relative_to(REPO))
    print("wrote", out_json_tree.relative_to(REPO), "(SSoT catalog)")
    print("wrote", out_json_pkg.relative_to(REPO), "(eon-schema vendored)")
    print("wrote", out_py.relative_to(REPO))
    print("wrote", out_h.relative_to(REPO))
    print("wrote", out_idx.relative_to(REPO))
    print("wrote", out_ini.relative_to(REPO))
    print("wrote", out_json_inc.relative_to(REPO))
    print("wrote", out_apply.relative_to(REPO))
    nfields = sum(
        len([f for f in s["fields"] if not f.get("nested")])
        for s in catalog["sections"].values()
    )
    print(f"sections={len(catalog['sections'])} scalar_fields={nfields} aliases={len(catalog['flat_aliases'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
