"""Outcome is job type and status. Trajectories stay readcon."""

import json
from pathlib import Path

from eon.contracts.load_outcome import load_outcome, outcome_to_json

ROOT = Path(__file__).resolve().parents[1]


def test_outcome_schema_is_job_type_and_status():
    text = (ROOT / "schema" / "eon_outcome.capnp").read_text(encoding="utf-8")
    assert "struct Outcome" in text
    assert "jobType @0" in text
    assert "statusCode @1" in text
    assert "statusText @2" in text
    assert "positions" not in text
    assert "readcon" in text
    raw = outcome_to_json(
        {"jobType": "minimization", "statusCode": 0, "statusText": "good"}
    )
    back = load_outcome(raw)
    assert back == {
        "jobType": "minimization",
        "statusCode": 0,
        "statusText": "good",
    }
    assert json.loads(raw)["statusCode"] == 0
