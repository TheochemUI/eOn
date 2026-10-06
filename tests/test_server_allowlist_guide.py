"""Guide samples stay inside the server allowlist, and name the climb gate."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
AKMC = ROOT / "docs" / "source" / "user_guide" / "akmc.md"
MINIMIZATION = ROOT / "docs" / "source" / "user_guide" / "minimization.md"
NEB_MD = ROOT / "docs" / "source" / "user_guide" / "neb.md"
NEB_CPP = ROOT / "client" / "NudgedElasticBand.cpp"
NEB_TEST = ROOT / "client" / "unit_tests" / "NEBRegressionTest.cpp"


def allowlist():
    from eon.config import ConfigClass

    return {
        section.name: {key.name: key for key in section.keys}
        for section in ConfigClass().format
    }


def test_akmc_sample_uses_the_section_the_server_reads():
    sections = allowlist()
    assert sections["AKMC"]["confidence_scheme"].default == "old"
    page = AKMC.read_text(encoding="utf-8")
    assert "[AKMC]" in page
    assert "confidence_scheme = old" in page
    assert "keys under `[akmc]` are ignored" in page


def test_rejected_sample_keys_are_absent_from_the_allowlist():
    sections = allowlist()
    assert "Refine" not in sections
    assert "parallel" not in sections["Main"]
    neb = sections["Nudged Elastic Band"]
    for key in (
        "climbing_image_band_slack",
        "match_endpoints",
        "match_method",
        "ci_after_rel",
    ):
        assert key not in neb
    opt = sections["Optimizer"]["opt_method"]
    assert opt.default == "cg"
    assert "qm" in opt.values
    assert "quickmin" not in opt.values

    minimization = MINIMIZATION.read_text(encoding="utf-8")
    assert "The default `opt_method` is `cg`." in minimization
    assert "`quickmin` is not a token." in minimization
    assert "server loads exits" in minimization
    neb_page = NEB_MD.read_text(encoding="utf-8")
    assert "Neither key is in `eon/config.yaml`" in neb_page
    assert "server exits on the unknown option" in neb_page
    assert "`[Main] parallel` is a client key" in neb_page


def test_guide_states_the_climbing_image_endpoint_gate():
    source = NEB_CPP.read_text(encoding="utf-8")
    assert "const bool climb = ci_active && maxEnergy > endpointEnergy;" in source
    assert (
        "const bool climb =\n"
        "      ci_active && enthalpy[static_cast<size_t>(highest)] > endpointEnthalpy;"
        in source
    )
    regression = NEB_TEST.read_text(encoding="utf-8")
    assert 'NEB does not climb below the endpoint energies' in regression
    page = NEB_MD.read_text(encoding="utf-8")
    assert "strictly above the higher endpoint" in page
    assert "A fixed-cell band compares potential energy." in page
    assert "A solid-state band compares enthalpy." in page
    assert "A monotonic band keeps the spring" in page
