#!/usr/bin/env python3
"""Generate consumer projections from Cap'n Proto L0 SSoT.

Authoring SSoT: schema/eon_params.capnp (monorepo root).

The PyPI split ``eon-schema`` vendors a copy under
``packages/eon-schema/src/eon_schema/ssot/`` for standalone installs;
codegen refreshes that catalog JSON in lockstep. Fat release tarballs
(``git archive`` of the full monorepo) still ship everything for
conda-forge ``eon-feedstock`` — layout cleanup is fine; the release
contract is one fat archive + optional split PyPI packages.

Outputs (relative to repo root):
  schema/eon_params_catalog.json
  packages/eon-schema/.../ssot/eon_params_catalog.json  (vendored for split)
  eon/_params_ssot_catalog.py
  include/eon/generated/ParametersSSOTDefaults.h
  include/eon/generated/ParametersSSOTFieldIndex.inc
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
    "Optimizer.CG.no_overshooting": "cg_no_overshooting",
    "Optimizer.CG.knock_out_max_move": "cg_knock_out_max_move",
    "Optimizer.CG.line_search": "cg_line_search",
    "Optimizer.CG.line_converged": "cg_line_converged",
    "Optimizer.CG.line_search_max_iter": "cg_max_iter_line_search",
    "Optimizer.CG.max_iter_before_reset": "cg_max_iter_before_reset",
    "Optimizer.Quickmin.steepest_descent": "qm_steepest_descent",
    "Optimizer.SD.alpha": "sd_alpha",
    "Optimizer.SD.two_point": "sd_two_point",
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


def parse_capnp(text: str) -> dict:
    structs: dict[str, list] = {}
    lines = []
    for line in text.splitlines():
        if "#" in line and '"' not in line.split("#")[0]:
            line = line.split("#")[0]
        lines.append(line)
    body = "\n".join(lines)
    struct_re = re.compile(r"struct\s+(\w+)\s*\{([^}]*)\}", re.S)
    field_re = re.compile(
        r"(\w+)\s+@(\d+)\s*:\s*([A-Za-z0-9_.]+)(?:\s*=\s*([^;]+))?\s*;",
        re.S,
    )
    for m in struct_re.finditer(body):
        sname = m.group(1)
        fields = []
        for fm in field_re.finditer(m.group(2)):
            fname, ord_, typ, default = fm.groups()
            nested = typ[0].isupper() and typ.endswith("Options") and default is None
            fields.append(
                {
                    "name": fname,
                    "ordinal": int(ord_),
                    "type": typ,
                    "default": None if nested else parse_default(default, typ),
                    "snake": to_snake(fname),
                    **({"nested": True} if nested else {}),
                }
            )
        structs[sname] = fields
    return structs


def build_catalog(structs: dict) -> dict:
    sections = {}
    for sname, fields in structs.items():
        if sname == "EonParameters":
            continue
        sec = SECTION_MAP.get(sname)
        if not sec:
            continue
        sections[sec] = {
            "struct": sname,
            "fields": [
                {
                    "snake": f["snake"],
                    "capnp": f["name"],
                    "ordinal": f["ordinal"],
                    "type": f["type"],
                    "default": f["default"],
                    **({"nested": True} if f.get("nested") else {}),
                }
                for f in fields
            ],
        }

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
    out_h.parent.mkdir(parents=True, exist_ok=True)
    out_json_pkg.parent.mkdir(parents=True, exist_ok=True)

    catalog_json = json.dumps(catalog, indent=2, sort_keys=True) + "\n"
    out_json_tree.write_text(catalog_json)
    out_json_pkg.write_text(catalog_json)
    out_py.write_text(emit_python(catalog))
    out_h.write_text(emit_cpp_header(catalog))
    out_idx.write_text(emit_field_index(catalog))
    print("read", CAPNP.relative_to(REPO))
    print("wrote", out_json_tree.relative_to(REPO), "(SSoT catalog)")
    print("wrote", out_json_pkg.relative_to(REPO), "(eon-schema vendored)")
    print("wrote", out_py.relative_to(REPO))
    print("wrote", out_h.relative_to(REPO))
    print("wrote", out_idx.relative_to(REPO))
    nfields = sum(
        len([f for f in s["fields"] if not f.get("nested")])
        for s in catalog["sections"].values()
    )
    print(f"sections={len(catalog['sections'])} scalar_fields={nfields} aliases={len(catalog['flat_aliases'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
