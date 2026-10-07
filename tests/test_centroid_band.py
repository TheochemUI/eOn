"""A ring on the solid-state band is one centroid and one cell."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_the_solid_state_band_keeps_one_cell_per_image():
    guide = (ROOT / "docs" / "source" / "user_guide" / "neb.md").read_text(
        encoding="utf-8"
    )
    assert "A ring on this band is one centroid and one cell." in guide
    assert "the ring springs stay inside the image" in guide
    assert "not a step of this band" in guide
    band = (ROOT / "client" / "NudgedElasticBand.cpp").read_text(encoding="utf-8")
    assert "The band does not store a cell on each bead." in band
