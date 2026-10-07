"""Load an Outcome. JSON is a debug decode of schema/eon_outcome.capnp."""

from __future__ import annotations

import json
from typing import Any, Mapping

_FIELDS = ("jobType", "statusCode", "statusText")


def load_outcome(payload: str | Mapping[str, Any]) -> dict[str, Any]:
    """Return job type and status. Trajectories are not part of this record."""
    if isinstance(payload, str):
        data = json.loads(payload)
    else:
        data = dict(payload)
    if not isinstance(data, dict):
        raise TypeError("outcome payload must be an object")
    if "positions" in data or "trajectory" in data:
        raise ValueError("trajectories stay readcon")
    missing = [key for key in _FIELDS if key not in data]
    if missing:
        raise KeyError(",".join(missing))
    return {
        "jobType": str(data["jobType"]),
        "statusCode": int(data["statusCode"]),
        "statusText": str(data["statusText"]),
    }


def outcome_to_json(outcome: Mapping[str, Any]) -> str:
    """Debug JSON of the same three fields. Not a second schema."""
    return json.dumps(load_outcome(outcome), sort_keys=True)
