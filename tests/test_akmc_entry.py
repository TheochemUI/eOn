"""The aKMC entry point prints status, resets, and refuses incompatible settings."""

import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from eon import fileio as io
from eon.akmc import main
from eon.structure import Structure


class _Comm:
    def __init__(self):
        self.submitted = []

    def get_results(self, path, keep):
        return iter(())

    def queued_search_count(self):
        return 0

    def get_queue_size(self):
        return 0

    def submit_jobs(self, searches, invariants):
        self.submitted.extend(searches)

    def cancel_state(self, number):
        return 0


def _layout(tmp_path, extra=""):
    root = tmp_path / "run"
    root.mkdir()
    atoms = Structure(2)
    atoms.names = ["Pt", "Pt"]
    atoms.mass[:] = 195.084
    atoms.box = np.diag([20.0, 20.0, 20.0])
    atoms.r[1] = [2.5, 0.0, 0.0]
    io.savecon(str(root / "pos.con"), atoms)
    ini = root / "config.ini"
    ini.write_text(
        "\n".join(
            [
                "[Main]",
                "job = akmc",
                "temperature = 300",
                "[Recycling]",
                "use_recycling = false",
                "use_sb_recycling = false",
                "[Coarse Graining]",
                "use_mcamc = false",
                "use_askmc = false",
                "[Paths]",
                f"main_directory = {root}",
                f"results = {root}",
                f"states = {root / 'states'}",
                extra,
                "",
            ]
        )
    )
    return root, ini


def test_status_exits_after_reporting_the_state(tmp_path, monkeypatch, capsys):
    _root, ini = _layout(tmp_path)
    monkeypatch.setattr(sys, "argv", ["akmc", str(ini), "-s", "-q"])
    monkeypatch.setattr(
        "eon.communicator.get_communicator",
        lambda config: SimpleNamespace(get_queue_size=lambda: 3),
    )
    with pytest.raises(SystemExit) as caught:
        main()
    assert caught.value.code == 0
    out = capsys.readouterr().out
    assert "Current state: 0" in out
    assert "Searches in queue: 3" in out


def test_forced_reset_removes_the_dynamics_file(tmp_path, monkeypatch, capsys):
    root, ini = _layout(tmp_path)
    (root / "dynamics.txt").write_text("steps\n")
    monkeypatch.setattr(sys, "argv", ["akmc", str(ini), "-R", "-f", "-q"])
    with pytest.raises(SystemExit) as caught:
        main()
    assert caught.value.code == 0
    assert not (root / "dynamics.txt").is_file()
    assert capsys.readouterr().out == ""


def test_declined_restart_exits(tmp_path, monkeypatch):
    _root, ini = _layout(tmp_path)
    monkeypatch.setattr(sys, "argv", ["akmc", str(ini), "-r", "-q"])
    monkeypatch.setattr("builtins.input", lambda prompt: "n")
    with pytest.raises(SystemExit) as caught:
        main()
    assert caught.value.code == 1


def test_min_mode_with_dynamics_confidence_exits(tmp_path, monkeypatch):
    _root, ini = _layout(
        tmp_path,
        "\n".join(
            [
                "[AKMC]",
                "confidence_scheme = dynamics",
                "[Saddle Search]",
                "method = min_mode",
            ]
        ),
    )
    monkeypatch.setattr(sys, "argv", ["akmc", str(ini), "-q"])
    with pytest.raises(SystemExit) as caught:
        main()
    assert caught.value.code == 1


def test_one_unlocked_pass_returns(tmp_path, monkeypatch):
    root, ini = _layout(tmp_path)
    monkeypatch.chdir(root)
    monkeypatch.setattr(sys, "argv", ["akmc", str(ini), "-q", "-n"])
    comm = _Comm()
    monkeypatch.setattr("eon.communicator.get_communicator", lambda config: comm)
    assert main() is None
