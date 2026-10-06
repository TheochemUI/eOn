"""The eOn-devel bundle pulls foss and builds Cap'n Proto."""

import ast
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EASYCONFIG = ROOT / "eessi" / "eOn-devel-2026.06-GCCcore-15.2.0.eb"
README = ROOT / "eessi" / "README.md"
PAGE = ROOT / "docs" / "source" / "install" / "eessi.md"
CAPNP_SHA256 = "fa02378ad522b318916b9ad928d1372fc9abd43dd1f4f0392e50450f5c87828f"
SEPARATE_LOAD = (
    "module load foss/2026.1 eOn-devel/2026.06-GCCcore-15.2.0 "
    "CapnProto/1.4.0-GCCcore-15.2.0"
)


def _assigned(text, name):
    tree = ast.parse(text)
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        for target in node.targets:
            if isinstance(target, ast.Name) and target.id == name:
                return ast.literal_eval(node.value)
    raise AssertionError(f"{name} is not assigned")


def test_bundle_depends_on_foss_and_builds_capnproto():
    text = EASYCONFIG.read_text(encoding="utf-8")
    dependencies = _assigned(text, "dependencies")
    components = _assigned(text, "components")
    compilers = _assigned(text, "modextravars")
    assert ("foss/2026.1", "EXTERNAL_MODULE") in dependencies
    assert components[0][0] == "CapnProto"
    assert components[0][1] == "1.4.0"
    spec = components[0][2]
    assert spec["easyblock"] == "ConfigureMake"
    assert spec["checksums"] == [CAPNP_SHA256]
    assert "capnproto-c++-1.4.0.tar.gz" in spec["sources"]
    assert compilers["CC"] == "gcc"
    assert compilers["CXX"] == "g++"
    assert compilers["FC"] == "gfortran"


def test_guide_loads_the_bundle_without_a_second_module():
    for path in (README, PAGE):
        text = path.read_text(encoding="utf-8")
        assert "module load eOn-devel/2026.06-GCCcore-15.2.0" in text
        assert SEPARATE_LOAD not in text
        assert "capnp-rpc" in text
        assert "foss/2026.1" in text
