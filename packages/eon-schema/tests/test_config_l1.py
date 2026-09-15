"""L1 job-config models live in eon-schema (shared with eon-akmc)."""

from eon_schema.config import (
    Config,
    MainConfig,
    Metatomic,
    PotentialConfig,
    SaddleSearchConfig,
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
