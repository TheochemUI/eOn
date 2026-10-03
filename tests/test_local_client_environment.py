"""The Local communicator starts eonclient with UCX memory events off.

An MPI-linked eonclient loads UCX in every run; its memory hooks cost about
0.4 s per process on conda's OpenMPI. Local clients have no MPI peers.
"""
from __future__ import annotations

from eon.communicator import client_environment


def test_memory_events_off_by_default():
    env = client_environment({"PATH": "/usr/bin"})
    assert env["UCX_MEM_EVENTS"] == "no"
    assert env["PATH"] == "/usr/bin"


def test_user_setting_wins():
    env = client_environment({"UCX_MEM_EVENTS": "yes"})
    assert env["UCX_MEM_EVENTS"] == "yes"


def test_caller_environment_is_not_modified():
    base = {"HOME": "/tmp"}
    client_environment(base)
    assert "UCX_MEM_EVENTS" not in base
