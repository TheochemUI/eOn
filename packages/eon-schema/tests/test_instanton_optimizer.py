"""Rate optimizer selection and saved-reference paths retain their values."""

import struct

import pytest
from pydantic import ValidationError

from eon_schema.config.models import InstantonConfig


def test_instanton_optimizer_defaults_and_roundtrip():
    defaults = InstantonConfig()
    assert defaults.newton_limit == (1 << (8 * struct.calcsize("l") - 1)) - 1
    assert defaults.reactant_hessian_filename == ""
    assert defaults.saddle_hessian_filename == ""
    configured = InstantonConfig(
        newton_limit=0,
        reactant_hessian_filename="reactant-hessian.dat",
        saddle_hessian_filename="saddle-hessian.dat",
    )
    restored = InstantonConfig.model_validate_json(configured.model_dump_json())
    assert restored.newton_limit == 0
    assert restored.reactant_hessian_filename == "reactant-hessian.dat"
    assert restored.saddle_hessian_filename == "saddle-hessian.dat"
    assert InstantonConfig(newton_limit=2048).newton_limit == 2048


@pytest.mark.parametrize("limit", [-1, -0.5])
def test_instanton_optimizer_rejects_negative_limit(limit):
    with pytest.raises(ValidationError):
        InstantonConfig(newton_limit=limit)
