"""The omitted jump_max matches the client default of 10."""

from __future__ import annotations

import json
from pathlib import Path

import yaml

from eon.config import ConfigClass
from eon_schema.config.models import BasinHoppingConfig

ROOT = Path(__file__).resolve().parents[1]


def _server_default(section, key):
    for item in ConfigClass().format:
        if item.name != section:
            continue
        for option in item.keys:
            if option.name == key:
                return option.default
    raise AssertionError("%s %s" % (section, key))


def test_jump_max_default_matches_the_client():
    header = (ROOT / "include" / "eon" / "ParametersOptions.h").read_text()
    assert "jump_max{10}" in header
    assert BasinHoppingConfig().jump_max == 10
    loaded = yaml.load(
        (ROOT / "eon" / "config.yaml").read_text(), Loader=yaml.BaseLoader
    )
    assert loaded["Basin Hopping"]["options"]["jump_max"]["default"] == "10"
    assert _server_default("Basin Hopping", "jump_max") == "10"
    schema = json.loads(
        (ROOT / "docs" / "source" / "_extra" / "schema" / "basin-hopping.json").read_text()
    )
    assert schema["properties"]["jump_max"]["default"] == 10
