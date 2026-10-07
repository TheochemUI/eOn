"""The optimizer is not what keeps the instanton path unused."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_the_guide_names_what_keeps_the_path_unused():
    guide = (ROOT / "docs" / "source" / "user_guide" / "instanton.md").read_text(
        encoding="utf-8"
    )
    assert "The optimizer is not what keeps the path unused." in guide
    assert "fluctuation prefactor" in guide
    assert "collapses onto the saddle" in guide
