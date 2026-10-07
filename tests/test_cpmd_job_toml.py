"""The rgpot guide shows one TOML for the in-process and file routes."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_guide_shows_one_toml_for_both_routes():
    text = (ROOT / "docs" / "source" / "user_guide" / "rgpot_pot.md").read_text(
        encoding="utf-8"
    )
    assert "```toml" in text
    fence = text.split("```toml", 1)[1].split("```", 1)[0]
    assert 'functional = "BLYP"' in fence
    assert "cutoff_ry = 70.0" in fence
    assert "gc_cutoff = 1.0e-7" in fence
    assert "job_toml.py" in text
    assert "params_path" in text
    assert "cpmdc_params_render_input_deck" in text
    assert "si3n4-isomer1.toml" in text
    assert "asin-tls.toml" in text
