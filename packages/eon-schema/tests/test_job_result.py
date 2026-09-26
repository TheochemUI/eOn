from eon_schema.jobs import (
    dict_to_results_dat,
    job_result_capnp_path,
    job_result_dumps,
    job_result_from_wire,
    job_result_loads,
    job_result_scalars_from_results_dat,
    job_result_to_wire,
    results_dat_to_dict,
)


def test_job_result_capnp_exists():
    path = job_result_capnp_path()
    assert path.is_file()
    text = path.read_text(encoding="utf-8")
    assert "struct JobResult" in text
    assert "struct Geometry" in text
    assert "struct JobRequest" in text
    assert "positions" in text
    assert "statusCode" in text
    assert "enum TerminationCode" in text
    assert "dimerRestoredBest @21" in text
    assert "termination @30" in text
    assert "body :union" in text
    assert "struct MinimizationBody" in text
    assert "struct NEBBody" in text
    assert "struct ProcessSearchBody" in text
    assert "struct EngineCompatibility" in text
    assert "struct LandfoldArtifact" in text
    assert "landfoldArtifacts @35 :List(LandfoldArtifact);" in text
    # Body union already owns @31-@34; the artifact list is append-only.
    assert "unset @31" in text
    assert "forces @8 :List(Float64);" in text
    assert "atomId @9 :List(UInt64);" in text
    assert "fixedAxes @10 :List(UInt8);" in text
    assert "struct OptimizerProvenance" in text
    assert "optimizer @35 :OptimizerProvenance;" in text
    assert "struct EindirAbi" in text
    assert "rgpot @38 :RgpotIdentity;" in text


def test_results_dat_roundtrip_scalars():
    sample = (
        "0 termination_reason\n"
        "good termination_reason_text\n"
        "process_search job_type\n"
        "lj potential_type\n"
        "42 total_force_calls\n"
        "1.234567890123e+00 potential_energy_saddle\n"
        "1.000000000000e+00 potential_energy_reactant\n"
        "1.100000000000e+00 potential_energy_product\n"
        "2.345678901234e-01 barrier_reactant_to_product\n"
    )
    d = results_dat_to_dict(sample)
    assert d["termination_reason"] == 0
    assert d["job_type"] == "process_search"
    assert abs(d["potential_energy_saddle"] - 1.234567890123) < 1e-12
    mapped = job_result_scalars_from_results_dat(sample)
    assert mapped["status_code"] == 0
    assert mapped["force_calls"]["total"] == 42
    assert abs(mapped["barrier_reactant_to_product"] - 0.2345678901234) < 1e-10
    # round-trip text
    again = dict_to_results_dat(
        {
            "termination_reason": 0,
            "termination_reason_text": "good",
            "total_force_calls": 42,
        }
    )
    assert "0 termination_reason" in again
    assert "42 total_force_calls" in again


def test_job_result_wire_roundtrip():
    src = {
        "job_type": "minimization",
        "status_code": 0,
        "status_text": "good",
        "potential_energy": 1.25,
        "force_calls": {"total": 7, "minimization": 7},
    }
    wire = job_result_to_wire(src)
    assert wire["jobType"] == "minimization"
    assert wire["statusCode"] == 0
    assert wire["forceCalls"]["total"] == 7
    back = job_result_from_wire(wire)
    assert back["job_type"] == "minimization"
    assert back["status_code"] == 0
    assert back["force_calls"]["total"] == 7
    blob = job_result_dumps(src)
    again = job_result_loads(blob)
    assert again["job_type"] == "minimization"
    assert again["force_calls"]["total"] == 7


def test_landfold_artifact_roundtrip():
    src = {
        "job_type": "minimization",
        "status_code": 0,
        "landfold_artifacts": [
            {
                "schema": "landfold.analysis.v1",
                "source_run_id": "run-7",
                "input_digest": "sha256:abc",
                "engine_compatibility": {
                    "schema": "eon.compatibility.v1",
                    "engine_id": "eon",
                    "protocol_family": "eindir",
                    "protocol_major": 1,
                    "protocol_minor": 2,
                    "abi_major": 3,
                    "abi_minor": 4,
                    "layout_revision": 9,
                    "build_identity": "eon-schema",
                },
            }
        ],
    }
    wire = job_result_to_wire(src)
    assert wire["landfoldArtifacts"][0]["sourceRunId"] == "run-7"
    assert wire["landfoldArtifacts"][0]["engineCompatibility"]["engineId"] == "eon"
    back = job_result_from_wire(wire)
    assert back["landfold_artifacts"] == src["landfold_artifacts"]
    again = job_result_loads(job_result_dumps(src))
    arts = again["landfold_artifacts"]
    assert arts[0]["schema"] == "landfold.analysis.v1"
    assert arts[0]["source_run_id"] == "run-7"
    assert arts[0]["input_digest"] == "sha256:abc"
    compat = arts[0]["engine_compatibility"]
    assert compat["schema"] == "eon.compatibility.v1"
    assert compat["engine_id"] == "eon"
    assert compat["protocol_family"] == "eindir"
    assert int(compat["protocol_major"]) == 1
    assert int(compat["protocol_minor"]) == 2
    assert int(compat["abi_major"]) == 3
    assert int(compat["abi_minor"]) == 4
    assert int(compat["layout_revision"]) == 9
    assert compat["build_identity"] == "eon-schema"
def test_legacy_backend_omits_engine_compatibility():
    parsed = job_result_scalars_from_results_dat(
        "cg optimizer_backend\n"
        "eon.optimizer.v1 optimizer_provenance_schema\n"
        "eon.compatibility.v1 compatibility_schema\n"
        "eon engine_id\n"
    )
    assert parsed["optimizer"]["backend"] == "cg"
    assert parsed["optimizer"]["xts_abi"]["major"] == 0
    assert "engine_compatibility" not in parsed["compatibility"]


def test_xtsci_backend_round_trips_abi_fields():
    parsed = job_result_scalars_from_results_dat(
        "xtsci optimizer_backend\n"
        "eon.optimizer.v1 optimizer_provenance_schema\n"
        "eon.compatibility.v1 compatibility_schema\n"
        "eon engine_id\n"
        "eon.objective compatibility_engine_protocol_family\n"
        "1 compatibility_engine_protocol_major\n"
        "0 compatibility_engine_protocol_minor\n"
        "1 compatibility_engine_abi_major\n"
        "10 compatibility_engine_abi_minor\n"
        "2 compatibility_engine_layout_revision\n"
        "build engine_build_identity\n"
        "1 optimizer_xts_abi_major\n"
        "10 optimizer_xts_abi_minor\n"
        "2 optimizer_xts_abi_layout\n"
        "1 optimizer_eindir_abi_major\n"
        "0 optimizer_eindir_abi_minor\n"
        "1 optimizer_eindir_objective_layout\n"
        "64 optimizer_eindir_objective_size\n"
        "8 optimizer_eindir_objective_align\n"
        "1 optimizer_eindir_dlpack_major\n"
        "0 optimizer_eindir_dlpack_minor\n"
        "1 optimizer_eindir_features\n"
        "lj rgpot_name\n"
        "3.2.0 rgpot_version\n"
        "eon.rgpot.v1 rgpot_schema\n"
    )
    assert parsed["optimizer"]["xts_abi"] == {"major": 1, "minor": 10, "layout": 2}
    assert parsed["optimizer"]["eindir"]["objective_layout"] == 1
    assert parsed["compatibility"]["engine_compatibility"]["abiMinor"] == 10
    assert parsed["rgpot"]["version"] == "3.2.0"
    wire = job_result_to_wire({"optimizer": parsed["optimizer"], "rgpot": parsed["rgpot"]})
    again = job_result_from_wire(wire)
    assert again["optimizer"]["backend"] == "xtsci"
    assert again["optimizer"]["has_eindir"] is True
