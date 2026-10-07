"""The server catalog lists the path-integral keys."""

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _fields(section: dict) -> dict:
    return {field["snake"]: field["default"] for field in section["fields"]}


def test_catalog_lists_dynamics_and_instanton_path_keys():
    catalog = json.loads(
        (ROOT / "schema" / "eon_params_catalog.json").read_text(encoding="utf-8")
    )
    vendored = json.loads(
        (
            ROOT
            / "packages"
            / "eon-schema"
            / "src"
            / "eon_schema"
            / "ssot"
            / "eon_params_catalog.json"
        ).read_text(encoding="utf-8")
    )
    dynamics = _fields(catalog["sections"]["Dynamics"])
    instanton = _fields(catalog["sections"]["Instanton"])
    assert dynamics["path_beads"] == 8
    assert dynamics["path_springs"] == "trotter"
    assert dynamics["path_pile_tau"] == 100.0
    assert dynamics["path_seed"] == 1
    assert instanton["springs"] == "trotter"
    assert _fields(vendored["sections"]["Dynamics"]) == dynamics
    assert _fields(vendored["sections"]["Instanton"]) == instanton
    py = (ROOT / "eon" / "_params_ssot_catalog.py").read_text(encoding="utf-8")
    assert '"path_beads"' in py
    assert '"springs"' in py
