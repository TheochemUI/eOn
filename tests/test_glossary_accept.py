"""Glossary short forms stay on the repo accept list."""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GLOSSARY = ROOT / "docs" / "source" / "glossary.md"
ACCEPT = ROOT / ".proseguard" / "accept.txt"
INDEX = ROOT / "docs" / "source" / "index.md"


def _accepted() -> set[str]:
    terms = set()
    for line in ACCEPT.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            terms.add(line)
    return terms


def _fence() -> str:
    text = GLOSSARY.read_text(encoding="utf-8")
    return text.split("```{glossary}", 1)[1].split("```", 1)[0]


def test_index_lists_the_glossary():
    index = INDEX.read_text(encoding="utf-8")
    assert re.search(r"(?m)^glossary\s*$", index)


def test_glossary_short_forms_are_on_the_accept_list():
    fence = _fence()
    assert "PI-QTST" in fence
    short = sorted(set(re.findall(r"\b[A-Z][A-Z0-9]{1,}\b", fence)))
    missing = [word for word in short if word not in _accepted()]
    assert missing == []
