"""Python Parameters and the schema INI expose the RgpotPot CPMD keys."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BIND = ROOT / "client" / "python" / "bind" / "bind_parameters.cpp"
SCHEMA = ROOT / "packages" / "eon-schema" / "src"

PROPERTIES = (
    "rgpot_scf_type",
    "rgpot_functional",
    "rgpot_cutoff_ry",
    "rgpot_charge",
    "rgpot_multiplicity",
    "rgpot_engine_root",
    "rgpot_scratch_dir",
    "rgpot_input_block",
    "rgpot_permanent_dir",
    "rgpot_params_path",
    "rgpot_ranks_per_image",
    "rgpot_xtb_paramset",
    "rgpot_xtb_accuracy",
    "rgpot_xtb_electronic_temperature",
    "rgpot_xtb_max_iterations",
    "rgpot_xtb_charge",
    "rgpot_xtb_uhf",
)


def test_parameters_binding_names_the_cpmd_keys():
    text = BIND.read_text(encoding="utf-8")
    for name in PROPERTIES:
        assert f'"{name}"' in text


def test_write_models_ini_emits_params_path_and_ranks(tmp_path):
    import sys

    sys.path.insert(0, str(SCHEMA))
    from eon_schema.config.ini import write_models_ini
    from eon_schema.config.models import RgpotPot

    path = tmp_path / "rgpot.ini"
    write_models_ini(
        path,
        RgpotPot(params_path="/data/si3n4.bin", ranks_per_image=6),
    )
    text = path.read_text(encoding="utf-8")
    assert "params_path = /data/si3n4.bin" in text
    assert "ranks_per_image = 6" in text
