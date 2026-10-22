"""Basin hopping reads displacement_algorithm, default standard."""

from __future__ import annotations

import json
from pathlib import Path

import yaml

from eon_schema.config.models import BasinHoppingConfig

ROOT = Path(__file__).resolve().parents[1]


def test_displacement_algorithm_matches_the_client():
    header = (ROOT / "include" / "eon" / "ParametersOptions.h").read_text()
    assert 'displacement_algorithm{"standard"}' in header
    ini = (ROOT / "client" / "ParametersINI.cpp").read_text()
    assert '"displacement_algorithm"' in ini
    loaded = yaml.load(
        (ROOT / "eon" / "config.yaml").read_text(), Loader=yaml.BaseLoader
    )
    options = loaded["Basin Hopping"]["options"]
    assert "max_displacement_algorithm" not in options
    assert options["displacement_algorithm"]["default"] == "standard"
    assert BasinHoppingConfig().displacement_algorithm == "standard"
    schema = json.loads(
        (ROOT / "docs" / "source" / "_extra" / "schema" / "basin-hopping.json").read_text()
    )
    assert "max_displacement_algorithm" not in schema["properties"]
    assert schema["properties"]["displacement_algorithm"]["default"] == "standard"
    index = json.loads(
        (ROOT / "docs" / "source" / "_extra" / "schema" / "index.json").read_text()
    )
    keys = index["sections"]["Basin Hopping"]["keys"]
    assert "displacement_algorithm" in keys
    assert "max_displacement_algorithm" not in keys
