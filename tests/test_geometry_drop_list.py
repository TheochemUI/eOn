"""Geometry keeps forces, atom ids, and fixed axes, and names what it drops."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DROPS = "drops ConFrame charges, spins, magmoms, angles, headers, and specVersion"


def test_geometry_names_the_conframe_fields_it_drops():
    schema = (ROOT / "schema" / "eon_job_result.capnp").read_text(encoding="utf-8")
    wire = (ROOT / "eon" / "geometry" / "wire.py").read_text(encoding="utf-8")
    assert DROPS in schema
    assert DROPS in wire
    struct = schema[schema.find("struct Geometry") : schema.find("struct JobRequest")]
    for name in (
        "charges",
        "spins",
        "magmoms",
        "angles",
        "preboxHeader",
        "postboxHeader",
        "specVersion",
    ):
        assert name not in struct
    assert "forces @8" in struct
    assert "atomId @9" in struct
    assert "fixedAxes @10" in struct
    client = (ROOT / "eon" / "communicator_inprocess.py").read_text(encoding="utf-8")
    start = client.find("job_result = _job_result(")
    window = client[start : start + 500]
    assert "geometry_wire_lists(" in window
