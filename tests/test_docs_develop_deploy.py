"""A develop push deploys the instanton page and a live serve-mode href."""

import re
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WORKFLOW = ROOT / ".github" / "workflows" / "ci_docs.yml"
MAIN = ROOT / "docs" / "source" / "user_guide" / "main.md"
SERVE = ROOT / "docs" / "source" / "user_guide" / "serve_mode.md"
HREF = "https://rgpot.rgoswami.me/howto/integration.html"


def _get(url: str) -> tuple[int, str]:
    request = urllib.request.Request(url, headers={"User-Agent": "eon-url-check"})
    with urllib.request.urlopen(request, timeout=30) as response:
        return response.status, response.read().decode("utf-8", errors="replace")


def test_develop_deploys_and_main_names_instanton():
    workflow = WORKFLOW.read_text(encoding="utf-8")
    branches = workflow.split("branches:", 1)[1].split("tags:", 1)[0]
    assert "develop" in branches
    assert re.search(r"(?m)^Instanton\s*$", MAIN.read_text(encoding="utf-8"))
    serve = SERVE.read_text(encoding="utf-8")
    assert HREF in serve
    assert "integration_guide.html" not in serve


def test_live_instanton_rgpot_and_serve_href():
    status, _ = _get("https://eondocs.org/user_guide/instanton.html")
    assert status == 200
    status, body = _get("https://eondocs.org/user_guide/rgpot_pot.html")
    assert status == 200
    assert "ranks_per_image" in body
    status, _ = _get(HREF)
    assert status == 200
