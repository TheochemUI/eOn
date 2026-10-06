"""The KDB page is the definition of the three match numbers."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PAGE = ROOT / "docs" / "source" / "user_guide" / "kdb.md"
MODELS = (
    ROOT
    / "packages"
    / "eon-schema"
    / "src"
    / "eon_schema"
    / "config"
    / "models.py"
)


def test_kdb_page_names_aselite_package_and_defines_the_match_numbers():
    page = PAGE.read_text(encoding="utf-8")
    assert "kdb_nf" in page and "neighbor fudge" in page
    assert "kdb_dc" in page and "angstrom" in page
    assert "kdb_mac" in page and "cosine" in page
    assert "aselite" in page
    assert "tsase" in page
    assert "tsase" in page[page.index("aselite") - 80 : page.index("aselite") + 80]

    models = MODELS.read_text(encoding="utf-8")
    block = models[models.index("class KDBConfig") : models.index("class RecyclingConfig")]
    assert "not sure" not in block.lower()
    assert "Neighbor fudge" in block

    extra = []
    for path in (ROOT / "docs" / "source").rglob("*.md"):
        if path.resolve() == PAGE.resolve():
            continue
        if "neighbor fudge" in path.read_text(encoding="utf-8").lower():
            extra.append(path.relative_to(ROOT).as_posix())
    for path in (ROOT / "eon").glob("*.py"):
        if "neighbor fudge" in path.read_text(encoding="utf-8").lower():
            extra.append(path.relative_to(ROOT).as_posix())
    assert extra == []
