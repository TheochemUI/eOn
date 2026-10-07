"""Dynamics movies read reactant frames from ingest_directory."""

from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _body(text: str, name: str) -> str:
    start = text.find(f"def {name}(")
    assert start != -1, name
    end = text.find("\ndef ", start + 1)
    assert end != -1
    return text[start:end]


def test_dynamics_reads_frames_from_ingest_directory(tmp_path: Path):
    movie = (ROOT / "eon" / "movie.py").read_text(encoding="utf-8")
    dynamics = _body(movie, "dynamics")
    assert "frame_texts_for_paths(" in dynamics
    assert "get_reactant(" not in dynamics
    corpus = (ROOT / "eon" / "concorpus.py").read_text(encoding="utf-8")
    index = _body(corpus, "_index_tree")
    assert "ingest_directory(" in index
    read = _body(corpus, "frame_texts_for_paths")
    assert "readonly=True" in read
    assert "get_frame_texts(" in read

    pytest.importorskip("readcon_db")
    import numpy as np

    from eon import fileio as io
    from eon.movie import dynamics as load_dynamics
    from eon.structure import Structure

    def atoms(x: float) -> Structure:
        structure = Structure(1)
        structure.names = ["Cu"]
        structure.mass[:] = 63.546
        structure.box = np.diag([10.0, 10.0, 10.0])
        structure.r[0] = [x, 0.0, 0.0]
        return structure

    for number, x in ((0, 1.0), (1, 2.5)):
        state = tmp_path / "states" / str(number)
        state.mkdir(parents=True)
        io.savecon(str(state / "reactant.con"), atoms(x))

    class _State:
        def __init__(self, number: int):
            self.reactant_path = tmp_path / "states" / str(number) / "reactant.con"

        def get_reactant(self):
            raise AssertionError("dynamics must not read reactant.con itself")

    class _States:
        def get_num_states(self):
            return 2

        def get_state(self, number: int):
            return _State(number)

    frames = load_dynamics(tmp_path, _States(), unique=True)
    assert frames[0].r[0, 0] == pytest.approx(1.0)
    assert frames[1].r[0, 0] == pytest.approx(2.5)
