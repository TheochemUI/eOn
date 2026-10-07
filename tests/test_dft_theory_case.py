"""NWChem DFT input blocks accept any case of theory=dft."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_dft_theory_comparison_is_case_insensitive():
    source = (
        ROOT / "client" / "potentials" / "Rgpot" / "RGPotEngine.cpp"
    ).read_text(encoding="utf-8")
    assert "nwchemDftInputBlock(opt.theory" in source
    assert 'theoryKey == "dft"' in source
