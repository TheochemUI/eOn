"""Continuous parallel replica on a cluster waits between passes."""

import sys
import time

import pytest

from eon.config import ConfigClass
from eon.parallelreplica import main as replica_main
from tests.test_catalog_and_mains import _ini


def test_continuous_cluster_waits_ten_seconds(tmp_path, monkeypatch):
    root, ini = _ini(tmp_path, "parallel_replica")
    monkeypatch.chdir(root)
    monkeypatch.setattr(sys, "argv", ["replica", str(ini), "-C", "-q"])
    monkeypatch.setattr("eon.locking.LockFile.aquirelock", lambda self: True)

    real_init = ConfigClass.init_from_cli

    def init_as_cluster(self, args):
        real_init(self, args)
        self.comm_type = "cluster"

    monkeypatch.setattr(ConfigClass, "init_from_cli", init_as_cluster)
    slept = []

    def fast_sleep(seconds):
        slept.append(seconds)
        raise SystemExit(7)

    monkeypatch.setattr(time, "sleep", fast_sleep)
    with pytest.raises(SystemExit) as caught:
        replica_main()
    assert caught.value.code == 7
    assert slept == [10.0]
