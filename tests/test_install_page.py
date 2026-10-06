"""The install page does not treat the forge package as a CPMD install."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PAGE = ROOT / "docs" / "source" / "install" / "index.md"


def test_forge_package_does_not_contain_the_cpmd_engine():
    text = PAGE.read_text(encoding="utf-8")
    assert "does not contain the CPMD engine" in text
    start = text.split("# Obtaining sources", 1)[0]
    assert "rgpot_cpmd_blyp" not in start
