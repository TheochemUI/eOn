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
    assert "struct EngineCompatibility" in text
    assert "struct LandfoldArtifact" in text
    assert "landfoldArtifacts @35 :List(LandfoldArtifact);" in text
    assert "forces @8 :List(Float64);" in text
    assert "atomId @9 :List(UInt64);" in text
    assert "fixedAxes @10 :List(UInt8);" in text
    vend = (
        ROOT
        / "packages"
        / "eon-schema"
        / "src"
        / "eon_schema"
        / "jobs"
        / "eon_job_result.capnp"
    )
    assert vend.read_text(encoding="utf-8") == text
    assert "struct OptimizerProvenance" in text
    assert "struct EngineCompatibility" in text
    assert "struct EindirAbi" in text
    assert "struct RgpotIdentity" in text
    assert "optimizer @35 :OptimizerProvenance;" in text
    assert "compatibility @37 :EngineCompatibility;" in text
    assert "rgpot @38 :RgpotIdentity;" in text
    assert "termination @30" in text
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


def test_packaged_schema_matches_monorepo_schema():
    canonical = ROOT / "schema" / "eon_job_result.capnp"
    packaged = (
        ROOT
        / "packages"
        / "eon-schema"
        / "src"
        / "eon_schema"
        / "jobs"
        / "eon_job_result.capnp"
    )
    assert packaged.read_bytes().replace(b"\r\n", b"\n") == canonical.read_bytes().replace(
        b"\r\n", b"\n"
    )


def test_results_dat_adapter_legacy_optimizer_has_no_xtsci_stamp():
    import sys

    sys.path.insert(0, str(ROOT / "packages" / "eon-schema" / "src"))
    from eon_schema.jobs import job_result_scalars_from_results_dat, results_dat_to_dict

    sample = """0 termination_reason
minimization job_type
cg optimizer_backend
eon.optimizer.v1 optimizer_provenance_schema
eon.compatibility.v1 compatibility_schema
3 compatibility_readcon_spec_version
0.14.9 compatibility_readcon_min_version
lj rgpot_name
eon engine_id
3.3.1 engine_version
abc123 engine_build_identity
"""
    raw = results_dat_to_dict(sample)
    assert raw["termination_reason"] == 0
    assert raw["optimizer_backend"] == "cg"
    parsed = job_result_scalars_from_results_dat(sample)
    assert parsed["optimizer"]["backend"] == "cg"
    assert parsed["optimizer"]["schema"] == "eon.optimizer.v1"
    assert parsed["optimizer"]["xts_abi"] == {"major": 0, "minor": 0, "layout": 0}
    assert parsed["optimizer"]["has_eindir"] is False
    assert "engine_compatibility" not in parsed["compatibility"]
    assert parsed["rgpot"]["name"] == "lj"
    assert parsed["compatibility"]["engine"]["id"] == "eon"


def test_results_dat_adapter_xtsci_optimizer_and_eindir():
    import sys

    sys.path.insert(0, str(ROOT / "packages" / "eon-schema" / "src"))
    from eon_schema.jobs import job_result_scalars_from_results_dat, job_result_from_wire, job_result_to_wire

    sample = """0 termination_reason
minimization job_type
xtsci optimizer_backend
eon.optimizer.v1 optimizer_provenance_schema
eon.compatibility.v1 compatibility_schema
eon engine_id
eon.objective compatibility_engine_protocol_family
1 compatibility_engine_protocol_major
0 compatibility_engine_protocol_minor
1 compatibility_engine_abi_major
10 compatibility_engine_abi_minor
2 compatibility_engine_layout_revision
abc123 engine_build_identity
1 optimizer_xts_abi_major
10 optimizer_xts_abi_minor
2 optimizer_xts_abi_layout
1 optimizer_eindir_abi_major
0 optimizer_eindir_abi_minor
1 optimizer_eindir_objective_layout
64 optimizer_eindir_objective_size
8 optimizer_eindir_objective_align
1 optimizer_eindir_dlpack_major
0 optimizer_eindir_dlpack_minor
3 optimizer_eindir_features
eon.rgpot.v1 rgpot_schema
morse_pt rgpot_name
3.2.0 rgpot_version
3 compatibility_readcon_spec_version
0.14.9 compatibility_readcon_min_version
0.2.3 compatibility_eon_schema_min_version
1.10.4 compatibility_rgpycrumbs_min_version
1.9.17 compatibility_chemparseplot_min_version
"""
    parsed = job_result_scalars_from_results_dat(sample)
    assert parsed["optimizer"]["backend"] == "xtsci"
    assert parsed["optimizer"]["xts_abi"] == {"major": 1, "minor": 10, "layout": 2}
    assert parsed["optimizer"]["has_eindir"] is True
    assert parsed["optimizer"]["eindir"]["abi_major"] == 1
    assert parsed["optimizer"]["eindir"]["objective_size"] == 64
    assert parsed["optimizer"]["eindir"]["features"] == 3
    stamp = parsed["compatibility"]["engine_compatibility"]
    assert stamp["schema"] == "eon.compatibility.v1"
    assert stamp["engineId"] == "eon"
    assert stamp["protocolFamily"] == "eon.objective"
    assert stamp["abiMinor"] == 10
    assert stamp["readconMinVersion"] == "0.14.9"
    assert stamp["chemparseplotMinVersion"] == "1.9.17"
    assert stamp["rgpotName"] == "morse_pt"
    assert parsed["rgpot"] == {"schema": "eon.rgpot.v1", "name": "morse_pt", "version": "3.2.0"}
    wire = job_result_to_wire(parsed)
    assert wire["optimizer"]["backend"] == "xtsci"
    assert wire["optimizer"]["xtsAbiMinor"] == 10
    assert wire["optimizer"]["eindir"]["objectiveSize"] == 64
    assert wire["compatibility"]["engineId"] == "eon"
    assert wire["rgpot"]["version"] == "3.2.0"
    back = job_result_from_wire(wire)
    assert back["optimizer"]["backend"] == "xtsci"
    assert back["optimizer"]["eindir"]["features"] == 3
    assert back["rgpot"]["name"] == "morse_pt"


def test_results_dat_adapter_does_not_upgrade_unknown_compatibility_schema():
    import sys

    sys.path.insert(0, str(ROOT / "packages" / "eon-schema" / "src"))
    from eon_schema.jobs import job_result_scalars_from_results_dat

    parsed = job_result_scalars_from_results_dat(
        """eon.compatibility.v0 compatibility_schema
eon engine_id
eon.objective compatibility_engine_protocol_family
1 compatibility_engine_protocol_major
0 compatibility_engine_protocol_minor
1 compatibility_engine_abi_major
0 compatibility_engine_abi_minor
2 compatibility_engine_layout_revision
abc123 engine_build_identity
"""
    )
    assert "engine_compatibility" not in parsed["compatibility"]
