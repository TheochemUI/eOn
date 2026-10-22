"""Positive confinement uses the client key confine_positive_min_move."""

from __future__ import annotations

import json
from pathlib import Path

import yaml

from eon_schema.config.models import SaddleSearchConfig

ROOT = Path(__file__).resolve().parents[1]


def test_confine_positive_min_move_matches_the_client():
    ini = (ROOT / "client" / "ParametersINI.cpp").read_text()
    assert '"confine_positive_min_move"' in ini
    loaded = yaml.load(
        (ROOT / "eon" / "config.yaml").read_text(), Loader=yaml.BaseLoader
    )
    options = loaded["Saddle Search"]["options"]
    assert "confine_positive_min_force" not in options
    assert options["confine_positive_min_move"]["default"] == "0.5"
    assert SaddleSearchConfig().confine_positive_min_move == 0.5
    schema = json.loads(
        (ROOT / "docs" / "source" / "_extra" / "schema" / "saddle-search.json").read_text()
    )
    assert "confine_positive_min_force" not in schema["properties"]
    prop = schema["properties"]["confine_positive_min_move"]
    assert prop["default"] == 0.5
    assert "force component" in prop["description"]
    index = json.loads(
        (ROOT / "docs" / "source" / "_extra" / "schema" / "index.json").read_text()
    )
    keys = index["sections"]["Saddle Search"]["keys"]
    assert "confine_positive_min_move" in keys
    assert "confine_positive_min_force" not in keys
