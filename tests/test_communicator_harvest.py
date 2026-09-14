"""Non-bundled harvest honours return_files.dat and keeps declared suffixes."""
from __future__ import annotations

from eon.communicator import harvest_job_files, read_return_files_manifest


def test_manifest_harvests_declared_json_and_log(tmp_path):
    (tmp_path / "results.dat").write_text("0 termination_reason\n")
    (tmp_path / "_potcalls.json").write_text("{}\n")
    (tmp_path / "client_quill.log").write_text("log\n")
    (tmp_path / "freqs.dat").write_text("undeclared\n")
    (tmp_path / "return_files.dat").write_text(
        "results.dat\n_potcalls.json\nclient_quill.log\nmissing.dat\n"
    )
    files = harvest_job_files(str(tmp_path))
    assert set(files) == {"results.dat", "_potcalls.json", "client_quill.log"}
    files["results.dat"].seek(0)
    assert "termination_reason" in files["results.dat"].read()


def test_missing_manifest_keeps_json_and_log(tmp_path):
    (tmp_path / "results.dat").write_text("0 termination_reason\n")
    (tmp_path / "_potcalls.json").write_text("{}\n")
    (tmp_path / "client_quill.log").write_text("log\n")
    (tmp_path / "notes.txt").write_text("nope\n")
    files = harvest_job_files(str(tmp_path))
    assert "results.dat" in files
    assert "_potcalls.json" in files
    assert "client_quill.log" in files
    assert "notes.txt" not in files


def test_read_return_files_manifest_none_when_absent(tmp_path):
    assert read_return_files_manifest(str(tmp_path)) is None
