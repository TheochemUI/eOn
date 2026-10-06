"""The RgpotPot guide and the schema literal use one spelling."""

from __future__ import annotations

import configparser
import io
import re
from pathlib import Path

from eon_schema.config.models import PotentialConfig, RgpotPot

ROOT = Path(__file__).resolve().parents[3]
PAGE = ROOT / "docs" / "source" / "user_guide" / "rgpot_pot.md"


def ini_blocks(text: str) -> list[str]:
    blocks: list[str] = []
    marker = "```{code-block} ini\n"
    rest = text
    while True:
        start = rest.find(marker)
        if start < 0:
            break
        rest = rest[start + len(marker) :]
        end = rest.find("```")
        if end < 0:
            raise AssertionError("rgpot_pot.md ini block is not closed")
        blocks.append(rest[:end])
        rest = rest[end + 3 :]
    return blocks


def test_guide_cpmd_ini_uses_the_schema_literal():
    text = PAGE.read_text(encoding="utf-8")
    assignments = re.findall(r"(?m)^potential\s*=\s*(\S+)\s*$", text)
    assert assignments, "rgpot_pot.md shows no potential assignment"
    assert set(assignments) == {"rgpot"}
    assert "`RGPOT` and `rgpot`" not in text

    cpmd = [block for block in ini_blocks(text) if "backend = cpmdc" in block]
    assert cpmd, "rgpot_pot.md has no cpmdc ini block"
    for block in cpmd:
        parser = configparser.ConfigParser()
        parser.optionxform = str
        parser.read_file(io.StringIO(block))
        potential = PotentialConfig(potential=parser.get("Potential", "potential"))
        rgpot = RgpotPot(backend=parser.get("RgpotPot", "backend"))
        assert potential.potential == "rgpot"
        assert rgpot.backend == "cpmdc"
