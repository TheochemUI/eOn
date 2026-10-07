"""The guide says the env var replaces the ini key, and the group size matches the header."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GUIDE = ROOT / "docs" / "source" / "user_guide" / "rgpot_pot.md"
HEADER = ROOT / "include" / "eon" / "ParametersOptions.h"
MODELS = (
    ROOT
    / "packages"
    / "eon-schema"
    / "src"
    / "eon_schema"
    / "config"
    / "models.py"
)


def test_rgpot_params_path_replaces_the_ini_key():
    text = GUIDE.read_text(encoding="utf-8")
    assert "RGPOT_PARAMS_PATH" in text
    assert "replaces the ini key" in text


def _header_ranks_comment() -> str:
    lines = HEADER.read_text(encoding="utf-8").splitlines()
    index = next(i for i, line in enumerate(lines) if "ranks_per_image{0}" in line)
    words = []
    cursor = index - 1
    while cursor >= 0 and lines[cursor].strip().startswith("//"):
        words.append(lines[cursor].split("//", 1)[1].strip())
        cursor -= 1
    return " ".join(reversed(words))


def test_ranks_description_matches_the_header():
    comment = _header_ranks_comment()
    models = MODELS.read_text(encoding="utf-8")
    start = models.index("ranks_per_image: int")
    chunk = " ".join(models[start : start + 500].split())
    for phrase in (
        "spread NEB images over the groups",
        "0 keeps one session on every rank",
    ):
        assert phrase in comment
        assert phrase in chunk
