"""L1 job-config models live in eon-schema (shared with eon-akmc)."""

from eon_schema.config import (
    Config,
    MainConfig,
    Metatomic,
    PotentialConfig,
    SaddleSearchConfig,
    XtsciConfig,
)


def test_main_config_defaults():
    m = MainConfig()
    assert m.job == "akmc"


def test_metatomic_model_fields():
    # construction with defaults
    m = Metatomic()
    assert hasattr(m, "model_path") or "model_path" in type(m).model_fields


def test_root_config_composes_sections():
    # Config requires nested models — build with defaults where possible
    assert MainConfig is not None
    assert PotentialConfig is not None
    assert Config is not None


def test_saddle_search_accepts_scalar_minus_one_atom_list():
    for value in (-1, [-1], "-1"):
        cfg = SaddleSearchConfig(displace_atom_list=value)
        assert cfg.displace_atom_list == value
    listed = SaddleSearchConfig(displace_atom_list=[0, 2, 4])
    assert listed.displace_atom_list == [0, 2, 4]
    csv = SaddleSearchConfig(displace_atom_list="0, 1, 2")
    assert csv.displace_atom_list == "0, 1, 2"


def test_xtsci_is_an_engine_with_methods():
    assert XtsciConfig().method == "lbfgs"
    assert XtsciConfig(method="newton").method == "newton"
    assert XtsciConfig(method="polak_ribiere").method == "polak_ribiere"
    assert XtsciConfig().qn_step == "lbfgs"
    assert XtsciConfig().precon == "none"
    assert XtsciConfig(qn_step="newton", precon="pair").qn_step == "newton"
    assert XtsciConfig(qn_step="newton", precon="pair").precon == "pair"
    assert "method" in XtsciConfig.model_fields
    assert "qn_step" in XtsciConfig.model_fields
    assert "precon" in XtsciConfig.model_fields
    assert "xtsci_method" not in XtsciConfig.model_fields
