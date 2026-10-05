"""config.ini JSON Schema and the llms.txt map."""

from __future__ import annotations

import json
from pathlib import Path

from eon_schema.config.ini import INI_FIELD_ALIASES, MODEL_INI_SECTION
from eon_schema.config.jsonschema import (
    DEFAULT_SITE,
    render_llms_txt,
    section_json_schema,
    section_models,
    user_guide_pages,
    write_config_schemas,
    write_llms_txt,
)
from eon_schema.config.models import (
    BasinHoppingConfig,
    DebugConfig,
    DimerConfig,
    DistributedReplicaConfig,
    MainConfig,
    NudgedElasticBandConfig,
    OptimizerConfig,
    PathsConfig,
    PotentialConfig,
    SaddleSearchConfig,
    XtsciConfig,
)

REPO = Path(__file__).resolve().parents[3]
GUIDE = REPO / "docs" / "source" / "user_guide"
PUBLISHED = REPO / "docs" / "source" / "_extra"


def test_every_mapped_section_has_ini_keys_and_defaults():
    seen = set()
    for cls_name, section, model in section_models():
        assert cls_name in MODEL_INI_SECTION
        schema = section_json_schema(model, section)
        assert schema["title"] == section
        assert schema["additionalProperties"] is False
        assert schema["$schema"].endswith("/schema")
        props = schema["properties"]
        assert set(props) == {
            INI_FIELD_ALIASES.get(
                (section, name), name if not field.alias else field.alias
            )
            for name, field in model.model_fields.items()
        }
        seen.add(section)
    assert seen == set(MODEL_INI_SECTION.values())


def test_dimer_uses_ini_names_and_model_defaults():
    schema = section_json_schema(DimerConfig, "Dimer")
    props = schema["properties"]
    assert "dimer_improved" not in props
    assert props["improved"]["default"] is True
    assert props["max_iterations"]["default"] == 1000
    assert props["rotations_min"]["default"] == 1
    assert props["torque_min"]["default"] == 0.1
    assert props["torque_max"]["default"] == 1.0
    assert props["remove_rotation"]["default"] is False
    assert props["rotations_max"]["default"] == 10
    assert props["finite_angle"]["default"] == 0.005
    assert props["converged_angle"]["default"] == 5.0
    assert props["opt_method"]["enum"] == ["sd", "cg", "lbfgs"]
    assert props["opt_method"]["default"] == "cg"
    assert "Steepest descent" in props["opt_method"]["description"]
    assert "examples" not in schema


def test_main_job_enum_and_unfixed_random_seed():
    props = section_json_schema(MainConfig, "Main")["properties"]
    assert props["job"]["default"] == "akmc"
    assert props["job"]["enum"] == list(
        MainConfig.model_fields["job"].annotation.__args__
    )
    assert props["temperature"]["default"] == 300.0
    assert props["finite_difference"]["default"] == 0.01
    assert "default" not in props["random_seed"]
    assert "Options:" in props["job"]["description"]


def test_paths_and_potential_use_effective_defaults():
    paths = section_json_schema(PathsConfig, "Paths")["properties"]
    assert paths["jobs_out"]["default"] == PathsConfig().jobs_out
    assert paths["main_directory"]["default"] == "./"
    assert paths["potential_files"]["default"] == "potfiles"
    log = section_json_schema(PotentialConfig, "Potential")["properties"][
        "log_potential"
    ]
    assert log["default"] is False


def test_nonfinite_and_unchanged_odd_default():
    neb = section_json_schema(NudgedElasticBandConfig, "Nudged Elastic Band")
    assert neb["properties"]["ci_after"]["default"] == "inf"
    assert neb["properties"]["images"]["default"] == 5
    assert neb["properties"]["spring"]["default"] == 5.0
    hop = section_json_schema(BasinHoppingConfig, "Basin Hopping")
    assert hop["properties"]["max_displacement_algorithm"]["default"] == 0
    assert hop["properties"]["displacement_distribution"]["enum"] == [
        "gaussian",
        "uniform",
    ]


def test_alias_with_space_and_optimizer_literals():
    dist = section_json_schema(DistributedReplicaConfig, "Distributed Replica")
    assert dist["properties"]["sampling steps"]["default"] == 500
    opt = section_json_schema(OptimizerConfig, "Optimizer")["properties"]
    assert opt["opt_method"]["default"] == "cg"
    assert opt["opt_method"]["enum"][0] == "box"
    assert "xtsci" in opt["opt_method"]["enum"]
    assert opt["converged_force"]["default"] == 0.01
    xt = section_json_schema(XtsciConfig, "Xtsci")["properties"]
    assert xt["method"]["default"] == "lbfgs"
    assert "polak_ribiere" in xt["method"]["enum"]
    assert xt["qn_step"]["default"] == "lbfgs"
    debug = section_json_schema(DebugConfig, "Debug")["properties"]
    assert debug["neb_mmf_estimator"]["enum"] == ["dimer", "lanczos"]
    search = section_json_schema(SaddleSearchConfig, "Saddle Search")["properties"]
    assert search["method"]["default"] == "min_mode"
    assert "artn" in search["min_mode_method"]["enum"]


def test_write_roundtrip_matches_index(tmp_path: Path):
    index = write_config_schemas(tmp_path, site=DEFAULT_SITE)
    assert set(index["sections"]) == set(MODEL_INI_SECTION.values())
    for section, entry in index["sections"].items():
        document = json.loads((tmp_path / entry["file"]).read_text(encoding="utf-8"))
        assert document["title"] == section
        assert list(document["properties"]) == entry["keys"]
        assert document["$id"] == entry["id"]


def test_llms_lists_only_existing_user_guide_pages(tmp_path: Path):
    assert GUIDE.is_dir()
    text = render_llms_txt(GUIDE, site=DEFAULT_SITE)
    stems = {path.stem for path in GUIDE.glob("*.md")}
    listed = []
    for line in text.splitlines():
        if "/user_guide/" not in line or not line.startswith("- ["):
            continue
        url = line.split("](", 1)[1].split(")", 1)[0]
        stem = url.rsplit("/", 1)[-1].removesuffix(".html")
        listed.append(stem)
        assert stem in stems
    assert set(listed) == stems
    assert listed[0] == "index"
    assert "dimer" in listed
    assert f"{DEFAULT_SITE}/schema/dimer.json" in text
    assert f"{DEFAULT_SITE}/schema/index.json" in text
    assert "rotation_backend" not in text
    dimer_line = next(
        line for line in text.splitlines() if "/user_guide/dimer.html" in line
    )
    assert "Schema: [Dimer]" in dimer_line
    assert "Dimer method" in dimer_line
    path = write_llms_txt(tmp_path / "llms.txt", GUIDE)
    assert path.read_text(encoding="utf-8") == text


def test_user_guide_order_follows_index_toc():
    pages = [path.stem for path in user_guide_pages(GUIDE)]
    assert pages[0] == "index"
    assert pages[1:3] == ["pyeonclient", "neighbor_lists"]
    assert "stdpar" in pages
    assert len(pages) == len(list(GUIDE.glob("*.md")))


def test_published_tree_matches_generator(tmp_path: Path):
    write_config_schemas(tmp_path / "schema", site=DEFAULT_SITE)
    write_llms_txt(tmp_path / "llms.txt", GUIDE, site=DEFAULT_SITE)
    published_schema = PUBLISHED / "schema"
    assert published_schema.is_dir()
    fresh = {
        path.name: path.read_text(encoding="utf-8")
        for path in (tmp_path / "schema").iterdir()
    }
    on_disk = {
        path.name: path.read_text(encoding="utf-8")
        for path in published_schema.iterdir()
    }
    assert fresh == on_disk
    assert (PUBLISHED / "llms.txt").read_text(encoding="utf-8") == (
        tmp_path / "llms.txt"
    ).read_text(encoding="utf-8")
