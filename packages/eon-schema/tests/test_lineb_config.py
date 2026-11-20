"""LinebConfig round-trip."""

from eon_schema.config.models import LinebConfig


def test_lineb_defaults():
    cfg = LinebConfig()
    assert cfg.beads == 0
    assert cfg.molecular is False
    assert cfg.spring == 0.1
    assert cfg.endpoint_scale == 0.2


def test_lineb_json_round_trip():
    configured = LinebConfig(
        initial_path="neb.con",
        energy=0.25,
        beads=32,
        molecular=True,
        opt_tol=1e-4,
        perp_tol=2e-4,
    )
    restored = LinebConfig.model_validate_json(configured.model_dump_json())
    assert restored.initial_path == "neb.con"
    assert restored.energy == 0.25
    assert restored.beads == 32
    assert restored.molecular is True
    assert restored.opt_tol == 1e-4
    assert restored.perp_tol == 2e-4
