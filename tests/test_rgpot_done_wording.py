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


def test_ranks_description_matches_the_header():
    header = HEADER.read_text(encoding="utf-8")
    models = MODELS.read_text(encoding="utf-8")
    start = models.index("ranks_per_image: int")
    chunk = models[start : start + 500]
    for phrase in (
        "spread NEB images over the groups",
        "0 keeps one session on every rank",
    ):
        assert phrase in header
        assert phrase in chunk
