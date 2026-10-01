"""RgpotPot keys in config.yaml and the pydantic model are one set."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml
from eon_schema.config import RgpotPot

REPO = Path(__file__).resolve().parents[3]
YAML_PATH = REPO / "eon" / "config.yaml"


def _normalize(value: object) -> str:
    if isinstance(value, str):
        lowered = value.lower()
        if lowered in {"true", "false"}:
            return lowered
        try:
            return format(float(value), ".15g")
        except ValueError:
            return value
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int | float):
        return format(value, ".15g")
    return str(value)


def test_rgpot_pot_keys_match_config_yaml() -> None:
    """A key in only one of the two lists fails this test."""
    if not YAML_PATH.is_file():
        pytest.skip("monorepo config.yaml is not next to this package")
    raw = yaml.load(YAML_PATH.read_text(encoding="utf-8"), Loader=yaml.BaseLoader)
    options = raw["RgpotPot"]["options"]
    yaml_keys = set(options)
    model_keys = set(RgpotPot.model_fields)
    only_yaml = sorted(yaml_keys - model_keys)
    only_model = sorted(model_keys - yaml_keys)
    assert not only_yaml and not only_model, (
        f"only in config.yaml: {only_yaml}; only in RgpotPot: {only_model}"
    )
    for key in sorted(yaml_keys):
        info = RgpotPot.model_fields[key]
        assert _normalize(options[key]["default"]) == _normalize(info.default), key


def test_cpmd_aliases_load_and_are_not_written() -> None:
    """A dump carries one cutoff and one functional key per section."""
    from eon_schema.config import Cpmd

    rg = RgpotPot(cutOffRy=60.0, cpmd_functional="PBE")
    assert rg.cutOffRy == 60.0
    dumped = rg.model_dump()
    assert "cutoff_ry" in dumped and "functional" in dumped
    assert not {"cutOffRy", "cpmd_cut_off_ry", "cpmd_functional"} & set(dumped)

    cpmd = Cpmd(cutoff_ry=60.0)
    assert cpmd.cutoff_ry == 60.0
    dumped = cpmd.model_dump()
    assert "cutOffRy" in dumped and "functional" in dumped
    assert not {"cutoff_ry", "cpmd_cut_off_ry", "cpmd_functional"} & set(dumped)
