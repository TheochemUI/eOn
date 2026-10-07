"""A dynamics rewind past a transition keeps the bond boost."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_rewind_keeps_the_bond_boost():
    safe = (ROOT / "client" / "SafeHyperJob.cpp").read_text(encoding="utf-8")
    replica = (ROOT / "client" / "ReplicaDynamicsJob.cpp").read_text(
        encoding="utf-8"
    )
    parrep = (ROOT / "client" / "ParallelReplicaJob.cpp").read_text(
        encoding="utf-8"
    )
    case = (ROOT / "client" / "unit_tests" / "BondBoostTest.cpp").read_text(
        encoding="utf-8"
    )
    assert "assignKeepingBias" in safe
    assert "assignKeepingBias" in replica
    assert "assignKeepingBias(transitionStructure)" in parrep
    assert "a dropped bond boost leaves the next Verlet step on stale bias" in case
    assert "REQUIRE((kept.second - kept.first).norm() > 1e-4)" in case
    assert "REQUIRE(dropped.second.isApprox(dropped.first, 1e-12))" in case
