"""The omitted displace_radius matches the client default of 4 angstroms."""

from __future__ import annotations

import json
from pathlib import Path

import yaml

from eon.config import ConfigClass
from eon_schema.config.models import SaddleSearchConfig

ROOT = Path(__file__).resolve().parents[1]


def test_displace_radius_default_matches_the_client():
    header = (ROOT / "include" / "eon" / "ParametersOptions.h").read_text()
    assert "displace_radius{4.0}" in header
    assert SaddleSearchConfig().displace_radius == 4.0
    loaded = yaml.load(
        (ROOT / "eon" / "config.yaml").read_text(), Loader=yaml.BaseLoader
    )
    assert loaded["Saddle Search"]["options"]["displace_radius"]["default"] == "4.0"
    schema = json.loads(
        (ROOT / "docs" / "source" / "_extra" / "schema" / "saddle-search.json").read_text()
    )
    assert schema["properties"]["displace_radius"]["default"] == 4.0


def test_omitted_displace_radius_is_four(tmp_path):
    path = tmp_path / "config.ini"
    path.write_text("[Main]\njob = akmc\n")
    cfg = ConfigClass()
    cfg.init(str(path))
    assert cfg.disp_radius == 4.0
