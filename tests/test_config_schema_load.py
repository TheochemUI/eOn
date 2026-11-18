"""``ConfigClass`` reads ``config.yaml`` with the C parser where there is one.

The schema file is parsed on every server start. libyaml's loader gives the
same documents as the pure Python one at a fraction of the time, so the
schema read must come out the same either way, and must not spend its time
in the Python scanner when the C loader is present.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest
import yaml

_TESTS = Path(__file__).resolve().parent
if str(_TESTS) not in sys.path:
    sys.path.insert(0, str(_TESTS))

from driver_oracles import count_events  # noqa: E402

from eon import config as config_module  # noqa: E402

SCHEMA = Path(config_module.__file__).parent / "config.yaml"


def schema_of(config):
    return [
        (
            section.name,
            [(k.name, k.kind, k.default, list(k.values)) for k in section.keys],
        )
        for section in config.format
    ]


def reference_schema():
    y = yaml.load(SCHEMA.read_text(), Loader=yaml.BaseLoader)
    out = []
    for name in y:
        keys = []
        for key, attr in y[name]["options"].items():
            keys.append((key, attr["kind"], attr["default"], list(attr.get("values", []))))
        out.append((name, keys))
    return out


def test_the_schema_is_what_the_pure_python_loader_reads():
    assert schema_of(config_module.ConfigClass()) == reference_schema()


@pytest.mark.skipif(
    not hasattr(yaml, "CBaseLoader"), reason="PyYAML was built without libyaml"
)
def test_the_python_scanner_is_not_used_when_libyaml_is_present():
    with count_events() as counts:
        config_module.ConfigClass()
    assert counts["py:Scanner.check_token"] == 0
