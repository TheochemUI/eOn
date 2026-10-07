"""Dynamics and replica exchange lock a seeded result."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_dynamics_and_replica_exchange_are_pinned():
    source = (ROOT / "client" / "unit_tests" / "JobIntegrationTest.cpp").read_text(
        encoding="utf-8"
    )
    assert "-23.58285447425" in source
    assert "samplingCalls == 161" in source
