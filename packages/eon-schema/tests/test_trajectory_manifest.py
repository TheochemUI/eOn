import math

import pytest

from eon_schema.jobs import (
    TrajectoryManifestError,
    geometry_digest,
    landfold_consume,
    landfold_embedding_inputs,
    landfold_fes_inputs,
    results_dat_to_dict,
    trajectory_manifest,
    trajectory_manifest_dumps,
    trajectory_manifest_loads,
)


def _frame(x, energy=None, frame_id=None):
    frame = {
        "positions": [0.0, 0.0, 0.0, x, 0.0, 0.0],
        "atomic_numbers": [1, 1],
        "box": [10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 10.0],
    }
    if energy is not None:
        frame["energy"] = energy
        frame["has_energy"] = True
    if frame_id is not None:
        frame["frame_id"] = frame_id
    return frame


def test_digest_tracks_exact_coordinate_bits():
    base = {"positions": [0.0, 0.0, 0.0]}
    nudged = {"positions": [math.nextafter(0.0, 1.0), 0.0, 0.0]}
    assert geometry_digest(base) != geometry_digest(nudged)
    assert geometry_digest(base) == geometry_digest({"positions": [0.0, 0.0, 0.0]})


def test_manifest_round_trip_and_landfold_inputs():
    manifest = trajectory_manifest(
        [_frame(1.0, energy=-0.5, frame_id="a"), _frame(1.25, energy=-0.25, frame_id="b")],
        {"potential_type": "MORSE_PT", "revision": "a" * 40},
        source_run_id="run-7",
    )
    assert manifest["schema"] == "eon.trajectory.v1"
    assert manifest["potential"]["kernel"] == "MorsePot"
    assert manifest["potential"]["provider"] == "rgpot"
    assert manifest["frames"][0]["digest"].startswith("sha256:")
    again = trajectory_manifest_loads(trajectory_manifest_dumps(manifest))
    assert again["geometry_digest"] == manifest["geometry_digest"]
    assert again["potential_digest"] == manifest["potential_digest"]
    consumed = landfold_consume(again)
    embedding = consumed["embedding"]
    assert embedding["schema"] == "landfold.embedding.input.v1"
    assert embedding["n_atoms"] == 2
    assert embedding["coordinates"][0][3] == 1.0
    assert embedding["frame_ids"] == ["a", "b"]
    assert consumed["fes"]["schema"] == "landfold.fes.input.v1"
    assert consumed["fes"]["energies"] == pytest.approx([-0.5, -0.25])
    assert consumed["potential"]["revision"] == "a" * 40


def test_order_changes_geometry_digest_not_frame_digest():
    left = trajectory_manifest(
        [_frame(1.0, energy=1.0), _frame(2.0, energy=2.0)],
        "lj",
        source_run_id="run",
    )
    right = trajectory_manifest(
        [_frame(2.0, energy=2.0), _frame(1.0, energy=1.0)],
        "lj",
        source_run_id="run",
    )
    assert left["frames"][0]["digest"] == right["frames"][1]["digest"]
    assert left["geometry_digest"] != right["geometry_digest"]


def test_results_dat_is_not_a_manifest_source():
    text = "0 termination_reason\nlj potential_type\n"
    parsed = results_dat_to_dict(text)
    assert parsed["potential_type"] == "lj"
    with pytest.raises(TrajectoryManifestError, match="results.dat"):
        trajectory_manifest_loads(text)
    with pytest.raises(TrajectoryManifestError, match="results.dat"):
        trajectory_manifest(text, "lj", source_run_id="run")


def test_non_rgpot_and_bad_revision_are_rejected():
    with pytest.raises(TrajectoryManifestError, match="rgpot"):
        trajectory_manifest([_frame(1.0)], "emt", source_run_id="run")
    with pytest.raises(TrajectoryManifestError, match="40 hex"):
        trajectory_manifest(
            [_frame(1.0, energy=1.0)],
            {"potential_type": "lj", "revision": "abc"},
            source_run_id="run",
        )


def test_rgpot_backend_alias_and_config_are_part_of_identity():
    plain = trajectory_manifest(
        [_frame(1.0, energy=1.0)],
        {"potential_type": "rgpot", "backend": "nwchem"},
        source_run_id="run",
    )
    assert plain["potential"]["backend"] == "nwchemc"
    assert plain["potential"]["kernel"] == "RgpotPot"
    other = trajectory_manifest(
        [_frame(1.0, energy=1.0)],
        {"potential_type": "zbl", "config": {"cut_inner": 2.0, "cut_global": 2.5}},
        source_run_id="run",
    )
    shifted = trajectory_manifest(
        [_frame(1.0, energy=1.0)],
        {"potential_type": "zbl", "config": {"cut_inner": 2.1, "cut_global": 2.5}},
        source_run_id="run",
    )
    assert other["potential"]["kernel"] == "ZBLPot"
    assert other["potential_digest"] != shifted["potential_digest"]
    assert other["geometry_digest"] == shifted["geometry_digest"]
    restored = trajectory_manifest_loads(trajectory_manifest_dumps(other))
    assert restored["potential"]["config"]["cut_inner"]["__f64__"] == other["potential"]["config"]["cut_inner"]["__f64__"]
    assert restored["potential_digest"] == other["potential_digest"]


def test_tampered_frame_and_incomplete_fes_are_rejected():
    manifest = trajectory_manifest(
        [_frame(1.0, energy=1.0), _frame(2.0)],
        "lj",
        source_run_id="run",
    )
    with pytest.raises(TrajectoryManifestError, match="energy"):
        landfold_fes_inputs(manifest)
    consumed = landfold_consume(manifest)
    assert consumed["fes"] is None
    assert landfold_embedding_inputs(manifest)["n_atoms"] == 2
    broken = trajectory_manifest_dumps(manifest)
    broken = broken.replace(manifest["frames"][0]["positions_hex"][3], "0x1.8p+1", 1)
    with pytest.raises(TrajectoryManifestError, match="digest"):
        trajectory_manifest_loads(broken)


def test_mixed_atom_counts_cannot_embed():
    wide = {"positions": [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0], "energy": 1.0}
    manifest = trajectory_manifest([_frame(1.0, energy=1.0), wide], "ljcluster", source_run_id="run")
    with pytest.raises(TrajectoryManifestError, match="atom count"):
        landfold_embedding_inputs(manifest)
