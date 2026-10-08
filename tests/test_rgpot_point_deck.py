"""The point config names the checked-in CPMD message."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_config_names_the_checked_in_message():
    text = (ROOT / "examples" / "rgpot-point" / "config.ini").read_text(
        encoding="utf-8"
    )
    assert "job = point" in text
    assert "potential = RGPOT" in text
    assert "params_path = messages/si3n4-isomer1.params.bin" in text
    guide = (ROOT / "docs" / "source" / "user_guide" / "rgpot_pot.md").read_text(
        encoding="utf-8"
    )
    assert "examples/rgpot-point/config.ini" in guide
    assert "SYMMETRY" in guide
    assert "NEWCODE" in guide
    assert "GC-CUTOFF" in guide
    assert "KLEINMAN-BYLANDER" in guide
