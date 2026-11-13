"""GPR dimer angles match the client, and the classical dimer stays 5."""

from __future__ import annotations

import json
from pathlib import Path

import yaml

from eon.config import ConfigClass
from eon_schema.config.models import DimerConfig, GPRDimerConfig

ROOT = Path(__file__).resolve().parents[1]


def _server_default(section, key):
    for item in ConfigClass().format:
        if item.name != section:
            continue
        for option in item.keys:
            if option.name == key:
                return option.default
    raise AssertionError("%s %s" % (section, key))


def test_gpr_angles_match_the_client():
    header = (ROOT / "include" / "eon" / "ParametersOptions.h").read_text()
    assert "converged_angle{0.08}" in header
    assert "relax_conv_angle{0.001}" in header
    assert "converged_angle{5.0}" in header
    gpr = GPRDimerConfig()
    assert gpr.converged_angle == 0.08
    assert gpr.relaxation_converged_angle == 0.001
    assert DimerConfig().converged_angle == 5.0
    loaded = yaml.load(
        (ROOT / "eon" / "config.yaml").read_text(), Loader=yaml.BaseLoader
    )
    options = loaded["GPR Dimer"]["options"]
    assert options["converged_angle"]["default"] == "0.08"
    assert options["relaxation_converged_angle"]["default"] == "0.001"
    assert _server_default("GPR Dimer", "converged_angle") == "0.08"
    assert _server_default("GPR Dimer", "relaxation_converged_angle") == "0.001"
    schema = json.loads(
        (ROOT / "docs" / "source" / "_extra" / "schema" / "gpr-dimer.json").read_text()
    )
    assert schema["properties"]["converged_angle"]["default"] == 0.08
    assert schema["properties"]["relaxation_converged_angle"]["default"] == 0.001
