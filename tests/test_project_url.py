"""The project documentation URL answers, and the coarse page matches the gate."""

import re
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PROJECTS = (
    ROOT / "pyproject.toml",
    ROOT / "pyproject-pyeonclient.toml",
)
PAGE = ROOT / "docs" / "source" / "user_guide" / "coarse_graining.md"
AKMC = ROOT / "eon" / "akmc.py"
GATE = "use_mcamc is not required"


def documentation_urls() -> list[str]:
    urls = []
    for path in PROJECTS:
        match = re.search(
            r'(?m)^Documentation = "([^"]+)"',
            path.read_text(encoding="utf-8"),
        )
        assert match, path.name
        urls.append(match.group(1))
    return urls


def test_project_documentation_url_returns_200():
    for url in documentation_urls():
        request = urllib.request.Request(
            url,
            method="GET",
            headers={"User-Agent": "eon-url-check"},
        )
        with urllib.request.urlopen(request, timeout=30) as response:
            assert response.status == 200


def test_coarse_page_states_the_shipped_discover_decide_gate():
    source = AKMC.read_text(encoding="utf-8")
    assert GATE in source
    page = PAGE.read_text(encoding="utf-8")
    assert "`use_mcamc` is not required" in page
    assert "MRM" in page
    assert "FPTA" in page
    assert "runs only when `use_mcamc` is true" not in page
