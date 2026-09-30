"""AmselConfig fields are the [amsel] keys in eon/config.yaml."""
from __future__ import annotations

from pathlib import Path

import yaml
from eon_schema.config import AmselConfig

REPO = Path(__file__).resolve().parents[1]
YAML_PATH = REPO / "eon" / "config.yaml"

_AMSEL_KEYS = (
    "discover_decide",
    "e_min_init",
    "e_min_step",
    "e_min_floor",
    "cv_threshold",
    "on_error",
)


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


def test_amsel_config_fields_match_config_yaml() -> None:
    raw = yaml.load(YAML_PATH.read_text(encoding="utf-8"), Loader=yaml.BaseLoader)
    options = raw["amsel"]["options"]
    assert tuple(options) == _AMSEL_KEYS
    model_keys = set(AmselConfig.model_fields)
    assert model_keys == set(options)
    for key in _AMSEL_KEYS:
        info = AmselConfig.model_fields[key]
        assert _normalize(options[key]["default"]) == _normalize(info.default), key
    assert set(options["on_error"]["values"]) == {
        "fallback_single",
        "unavailable_mcamc",
        "raise",
    }


def test_eon_schema_reexports_amsel_config() -> None:
    from eon.schema import AmselConfig as FromEon

    assert FromEon is AmselConfig
