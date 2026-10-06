"""The potential reference page renders RgpotPot."""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PAGE = ROOT / "docs" / "source" / "user_guide" / "potential.md"
MODELS = (
    ROOT
    / "packages"
    / "eon-schema"
    / "src"
    / "eon_schema"
    / "config"
    / "models.py"
)
DIRECTIVE = "autopydantic_model:: eon.schema.RgpotPot"


def test_potential_page_renders_rgpot_pot():
    page = PAGE.read_text(encoding="utf-8")
    supported, _, rest = page.partition("## Configuration")
    assert re.search(r"(?m)^RGPOT\s*$", supported)
    configurations = rest.split("## Potential configurations", 1)[1]
    assert DIRECTIVE in configurations
    start = MODELS.read_text(encoding="utf-8").index("class RgpotPot")
    block = MODELS.read_text(encoding="utf-8")[start:]
    block = block[: block.index("class Metatomic")]
    assert "params_path:" in block
    assert "ranks_per_image:" in block
