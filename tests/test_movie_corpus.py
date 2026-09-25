"""Dynamics movies read reactant frames from a readcon-db corpus."""

from __future__ import annotations

from pathlib import Path

from eon.concorpus import directories_with_con, frame_texts_for_paths
from eon.movie import dynamics


class _Corpus:
    """Stand-in for readcon_db.ConCorpus keyed by the corpus directory."""

    blobs: dict[str, dict[tuple[int, int], str]] = {}

    def __init__(self, path, readonly=False):
        self.path = path
        self.readonly = readonly

    def ingest_directory(self, directory, start_traj_id=1):
        rows = []
        traj_id = start_traj_id
        store = self.blobs.setdefault(self.path, {})
        for path in sorted(Path(directory).glob("*.con")):
            store[(traj_id, 0)] = path.read_text()
            rows.append((traj_id, 1, str(path)))
            traj_id += 1
        return rows

    def get_frame_texts(self, keys):
        assert self.readonly
        store = self.blobs[self.path]
        return [store[tuple(key)] for key in keys]


def test_directories_with_con_is_sorted_and_not_recursive_files(tmp_path: Path):
    (tmp_path / "states" / "1").mkdir(parents=True)
    (tmp_path / "states" / "0").mkdir()
    (tmp_path / "states" / "0" / "reactant.con").write_text("a")
    (tmp_path / "states" / "1" / "reactant.con").write_text("b")
    (tmp_path / "states" / "0" / "procdata").mkdir()
    (tmp_path / "states" / "0" / "procdata" / "saddle_0.con").write_text("c")
    (tmp_path / "notes.txt").write_text("no")
    found = directories_with_con(tmp_path)
    assert found == [
        tmp_path / "states" / "0",
        tmp_path / "states" / "0" / "procdata",
        tmp_path / "states" / "1",
    ]


def test_frame_texts_follow_ingest_order(tmp_path: Path):
    state0 = tmp_path / "states" / "0"
    state1 = tmp_path / "states" / "1"
    state0.mkdir(parents=True)
    state1.mkdir()
    (state0 / "reactant.con").write_text("reactant-zero")
    (state1 / "reactant.con").write_text("reactant-one")
    texts = frame_texts_for_paths(
        tmp_path,
        [state1 / "reactant.con", state0 / "reactant.con"],
        corpus_cls=_Corpus,
    )
    assert texts == ["reactant-one", "reactant-zero"]


def test_dynamics_uses_corpus_frame_text(monkeypatch, tmp_path: Path):
    seen = {}

    def fake_texts(tree, paths):
        seen["tree"] = tree
        seen["paths"] = list(paths)
        return ["FRAME-0", "FRAME-1"]

    def fake_loadcon(filein, reset=True):
        return filein.read()

    monkeypatch.setattr("eon.movie.frame_texts_for_paths", fake_texts)
    monkeypatch.setattr("eon.movie.io.loadcon", fake_loadcon)

    class _State:
        def __init__(self, number):
            self.reactant_path = tmp_path / "states" / str(number) / "reactant.con"

        def get_reactant(self):
            raise AssertionError("dynamics must not read reactant.con itself")

    class _States:
        def get_num_states(self):
            return 2

        def get_state(self, number):
            return _State(number)

    assert dynamics(tmp_path, _States(), unique=True) == ["FRAME-0", "FRAME-1"]
    assert seen["tree"] == tmp_path
    assert [p.name for p in seen["paths"]] == ["reactant.con", "reactant.con"]
