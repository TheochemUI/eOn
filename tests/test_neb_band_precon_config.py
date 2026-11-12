"""Server allowlist and schema for the rgsaddle Exp / ODE12r band keys."""

from eon.config import ConfigClass
from eon_schema.config.models import NudgedElasticBandConfig


def _neb_keys():
    for section in ConfigClass().format:
        if section.name == "Nudged Elastic Band":
            return {key.name: key for key in section.keys}
    raise AssertionError("Nudged Elastic Band is missing from config.yaml")


def test_band_precon_defaults_are_off_and_listed():
    keys = _neb_keys()
    model = NudgedElasticBandConfig()
    expected = {
        "band_precon": "none",
        "band_precon_A": 0.0,
        "band_precon_mu": 0.0,
        "band_precon_r_cut_factor": 0.0,
        "band_precon_c_stab": 0.0,
        "band_stepper": "auto",
    }
    for name, default in expected.items():
        assert name in keys
        assert keys[name].default == str(default) or keys[name].default == default
        assert getattr(model, name) == default
    assert keys["band_precon"].values == ["none", "exp"]
    assert keys["band_stepper"].values == ["auto", "ode12r"]
    assert model.band_precon == "none"
    assert model.band_stepper == "auto"


def test_band_precon_model_accepts_the_exp_ode12r_selection():
    tuned = NudgedElasticBandConfig(
        band_precon="exp",
        band_precon_A=3.0,
        band_precon_mu=1.0,
        band_precon_r_cut_factor=2.0,
        band_precon_c_stab=0.1,
        band_stepper="ode12r",
    )
    assert tuned.band_precon == "exp"
    assert tuned.band_precon_A == 3.0
    assert tuned.band_precon_mu == 1.0
    assert tuned.band_precon_r_cut_factor == 2.0
    assert tuned.band_precon_c_stab == 0.1
    assert tuned.band_stepper == "ode12r"
