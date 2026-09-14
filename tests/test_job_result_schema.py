"""Monorepo: JobResult Cap'n Proto file + results.dat adapters."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_schema_file_in_monorepo():
    p = ROOT / "schema" / "eon_job_result.capnp"
    assert p.is_file()
    text = p.read_text(encoding="utf-8")
    assert "struct JobResult" in text
    assert "struct JobRequest" in text
    assert "struct Geometry" in text
    # Flat geometry, not con text blobs
    assert "positions @0" in text
    assert ".con" not in text or "not .con text" in text


def test_adapters_importable_from_package():
    import sys

    sys.path.insert(0, str(ROOT / "packages" / "eon-schema" / "src"))
    from eon_schema.jobs import job_result_capnp_path, results_dat_to_dict

    assert job_result_capnp_path().is_file()
    d = results_dat_to_dict("0 termination_reason\ngood termination_reason_text\n")
    assert d["termination_reason"] == 0
