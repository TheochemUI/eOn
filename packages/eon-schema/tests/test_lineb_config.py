"""LinebConfig round-trip and beads validation."""

import pytest
from pydantic import ValidationError

from eon_schema.config.models import LinebConfig


def test_lineb_defaults():
    cfg = LinebConfig()
    assert cfg.beads == 0
    assert cfg.molecular is False
    assert cfg.interpolated_force_clock is False
    assert cfg.spring == 0.1
    assert cfg.endpoint_scale == 0.2


def test_lineb_json_round_trip():
    configured = LinebConfig(
        initial_path="neb.con",
        energy=0.25,
        beads=32,
        endpoint_kappa=2.5,
        molecular=True,
        interpolated_force_clock=True,
        interpolated_hessian=True,
        image_hessian_step=0.005,
        opt_tol=1e-4,
        perp_tol=2e-4,
    )
    restored = LinebConfig.model_validate_json(configured.model_dump_json())
    assert restored.initial_path == "neb.con"
    assert restored.energy == 0.25
    assert restored.endpoint_kappa == 2.5
    assert restored.beads == 32
    assert restored.molecular is True
    assert restored.interpolated_force_clock is True
    assert restored.interpolated_hessian is True
    assert restored.image_hessian_step == 0.005
    assert restored.opt_tol == 1e-4
    assert restored.perp_tol == 2e-4


def test_lineb_beads_one_refused():
    with pytest.raises(ValidationError, match="beads must be 0"):
        LinebConfig(beads=1)


def test_lineb_beads_two_ok():
    assert LinebConfig(beads=2).beads == 2


def test_image_displacement_requires_image_hessians():
    with pytest.raises(ValidationError, match="requires interpolated_hessian"):
        LinebConfig(image_hessian_step=0.005)
