"""An in-process AKMC search does not create a client workdir."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _slice(text: str, start_at: str, end_at: str) -> str:
    start = text.find(start_at)
    assert start != -1, start_at
    end = text.find(end_at, start + len(start_at))
    assert end != -1, end_at
    return text[start:end]


def test_inprocess_register_skips_the_client_workdir():
    text = (ROOT / "eon" / "explorer.py").read_text(encoding="utf-8")
    register = _slice(
        text,
        "def register_results(self):",
        "def __init__(self, states, previous_state, state, superbasin=None, config=None):",
    )
    guard = register.find("if not inprocess:")
    mkdir_at = register.find("jobs_in.mkdir")
    assert guard != -1 and mkdir_at != -1 and guard < mkdir_at
    assert "tot_searches = num_registered" in register
    make = _slice(text, "def make_jobs(self):", "def register_results(self):")
    stat_at = make.find('os.stat("pos.con")')
    make_guard = make.find("if not inprocess:")
    assert make_guard != -1 and stat_at != -1 and make_guard < stat_at
