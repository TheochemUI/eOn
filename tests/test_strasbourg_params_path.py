"""The guide and the CPMD example name params_path as the Strasbourg message."""

import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GUIDE = ROOT / "docs" / "source" / "user_guide" / "rgpot_pot.md"
EXAMPLE = ROOT / "examples" / "rgpot_cpmd_blyp" / "config.ini"
PAGE = "https://eondocs.org/user_guide/rgpot_pot.html"


def test_guide_and_example_name_the_strasbourg_message():
    guide = GUIDE.read_text(encoding="utf-8")
    example = EXAMPLE.read_text(encoding="utf-8")
    assert "params_path" in guide
    assert "Strasbourg message" in guide
    assert "params_path =" in example
    assert "Strasbourg" in example


def test_published_rgpot_page_contains_params_path():
    request = urllib.request.Request(PAGE, headers={"User-Agent": "eon-url-check"})
    with urllib.request.urlopen(request, timeout=30) as response:
        assert response.status == 200
        body = response.read().decode("utf-8", errors="replace")
    assert "params_path" in body
