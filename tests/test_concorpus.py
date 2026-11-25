"""Con path I/O stores the frame in a readcon-db corpus."""

from __future__ import annotations

from pathlib import Path

import pytest
import readcon

readcon_db = pytest.importorskip("readcon_db")

from eon import fileio as io
from eon.concorpus import corpus_dir, stored_frame_text


def _one_frame_fixture() -> Path:
    for path in sorted(Path("tests/data").rglob("*.con")):
        try:
            frames = readcon.read_con(str(path))
        except OSError:
            continue
        if len(frames) == 1:
            return path
    pytest.skip("no single-frame con fixture")


def _stored(path: Path) -> str:
    stored = stored_frame_text(path)
    assert stored
    return stored


def test_savecon_mirrors_into_readcon_db(tmp_path: Path):
    fixture = _one_frame_fixture()
    (tmp_path / "config.ini").write_text("[Main]\n")
    src = tmp_path / "in.con"
    src.write_text(fixture.read_text())
    atoms = io.loadcon(str(src))
    dest = tmp_path / "reactant.con"
    io.savecon(str(dest), atoms)

    stored = _stored(dest)
    original = readcon.read_con(str(dest))[0]
    mirrored = readcon.read_con_string(stored)[0]
    assert len(mirrored.atoms) == len(original.atoms)
    assert corpus_dir(dest) == tmp_path / "readcon.db"


def test_loadcon_mirrors_an_existing_file(tmp_path: Path):
    fixture = _one_frame_fixture()
    dest = tmp_path / "copied.con"
    dest.write_text(fixture.read_text())
    io.loadcon(str(dest))
    stored = _stored(dest)
    assert len(readcon.read_con_string(stored)[0].atoms) == len(
        readcon.read_con(str(dest))[0].atoms
    )


def _cu_atom(x: float):
    import numpy as np

    from eon.structure import Structure

    structure = Structure(1)
    structure.names = ["Cu"]
    structure.mass[:] = [63.546]
    structure.box = np.eye(3) * 10.0
    structure.r[0] = [x, 2.0, 3.0]
    return structure


def test_a_second_savecon_replaces_the_stored_text(tmp_path: Path):
    (tmp_path / "config.ini").write_text("[Main]\n")
    dest = tmp_path / "pos.con"
    io.savecon(str(dest), _cu_atom(1.0))
    first = dest.read_text()
    assert _stored(dest) == first
    io.savecon(str(dest), _cu_atom(1.5))
    second = dest.read_text()
    assert second != first
    assert _stored(dest) == second


def test_a_second_store_of_one_path_returns_the_new_text(tmp_path: Path):
    from eon.concorpus import load_frame_text, store_frame_text

    (tmp_path / "config.ini").write_text("[Main]\n")
    dest = tmp_path / "saddle.con"
    io.savecon(str(dest), _cu_atom(1.0))
    first_key = store_frame_text(dest, dest.read_text())
    io.savecon(str(dest), _cu_atom(1.5))
    second = dest.read_text()
    second_key = store_frame_text(dest, second)
    assert first_key is not None and second_key is not None
    assert second_key != first_key
    assert load_frame_text(corpus_dir(dest), *second_key) == second


def test_a_compressed_append_stays_out_of_the_corpus(tmp_path: Path):
    import numpy as np

    from eon.concorpus import flush_corpus_appends
    from eon.structure import Structure

    def _scattered(seed: int) -> Structure:
        structure = Structure(32)
        structure.names = ["Cu"] * 32
        structure.mass[:] = 63.546
        structure.box = np.eye(3) * 10.0
        structure.r[:] = np.random.default_rng(seed).uniform(0.0, 10.0, (32, 3))
        return structure

    (tmp_path / "config.ini").write_text("[Main]\n")
    dest = tmp_path / "movie.con.gz"
    io.savecon(str(dest), _scattered(0))
    io.savecon(str(dest), _scattered(1), "a")
    flush_corpus_appends()
    frame_text = readcon.write_con_string([_scattered(1).to_conframe()])
    assert dest.stat().st_size > len(frame_text.encode())
    assert len(readcon.read_con(str(dest))) == 2
    assert not (tmp_path / "readcon.db").exists()


def _xs(texts: list[str]) -> list[float]:
    return [readcon.read_con_string(text)[0].atoms[0].x for text in texts]


def _disk_xs(path: Path) -> list[float]:
    return [frame.atoms[0].x for frame in readcon.read_con(str(path))]


def test_the_generation_holds_the_initial_frame_and_the_appends(tmp_path: Path):
    from eon.concorpus import generation_frame_texts

    (tmp_path / "config.ini").write_text("[Main]\n")
    dest = tmp_path / "movie.con"
    io.savecon(str(dest), _cu_atom(1.0))
    io.savecon(str(dest), _cu_atom(2.0), "a")
    assert _disk_xs(dest) == [1.0, 2.0]
    assert _xs(generation_frame_texts(dest)) == [1.0, 2.0]


def test_an_overwrite_with_the_same_first_frame_starts_a_new_generation(
    tmp_path: Path,
):
    from eon.concorpus import generation_frame_texts, generation_id, load_frame_text

    (tmp_path / "config.ini").write_text("[Main]\n")
    dest = tmp_path / "movie.con"
    io.savecon(str(dest), _cu_atom(1.0))
    first = dest.read_text()
    io.savecon(str(dest), _cu_atom(2.0), "a")
    io.savecon(str(dest), _cu_atom(1.0))
    io.savecon(str(dest), _cu_atom(3.0), "a")
    assert _disk_xs(dest) == [1.0, 3.0]
    assert _xs(generation_frame_texts(dest)) == [1.0, 3.0]
    old = generation_id(str(dest.resolve()), first)
    assert _xs([load_frame_text(corpus_dir(dest), old, i) for i in (0, 1)]) == [
        1.0,
        2.0,
    ]
    assert load_frame_text(corpus_dir(dest), old, 2) is None


def test_a_repeated_overwrite_without_appends_reuses_the_generation(tmp_path: Path):
    from eon.concorpus import generation_id, start_generation

    (tmp_path / "config.ini").write_text("[Main]\n")
    dest = tmp_path / "pos.con"
    io.savecon(str(dest), _cu_atom(1.0))
    text = dest.read_text()
    assert start_generation(dest, text) == generation_id(str(dest.resolve()), text)


def _held_generation(dest: Path) -> int | None:
    from eon import concorpus

    return concorpus._generations.get(str(dest.resolve()))


@pytest.mark.parametrize(
    "xs",
    [(1.0, 2.0, 3.0), (1.0, 2.123456789, 3.987654321), (1.0, 1.25, 1.5, 1.75)],
)
def test_an_append_in_a_new_process_extends_the_generation_on_disk(
    tmp_path: Path, monkeypatch, xs
):
    from eon import concorpus
    from eon.concorpus import generation_frame_texts, generation_id

    (tmp_path / "config.ini").write_text("[Main]\n")
    dest = tmp_path / "movie.con"
    io.savecon(str(dest), _cu_atom(xs[0]))
    first = dest.read_text()
    for x in xs[1:-1]:
        io.savecon(str(dest), _cu_atom(x), "a")
    concorpus.flush_corpus_appends()
    monkeypatch.setattr(concorpus, "_generations", {})
    io.savecon(str(dest), _cu_atom(xs[-1]), "a")
    assert _held_generation(dest) == generation_id(str(dest.resolve()), first)
    assert _xs(generation_frame_texts(dest)) == pytest.approx(list(xs), abs=1e-6)


def test_a_resume_matches_frames_readcon_db_wrote_with_six_decimals(
    tmp_path: Path, monkeypatch
):
    from eon import concorpus
    from eon.concorpus import generation_id, load_frame_text

    (tmp_path / "config.ini").write_text("[Main]\n")
    dest = tmp_path / "movie.con"
    io.savecon(str(dest), _cu_atom(1.0))
    first = dest.read_text()
    io.savecon(str(dest), _cu_atom(2.123456789), "a")
    concorpus.flush_corpus_appends()
    tid = generation_id(str(dest.resolve()), first)
    stored = load_frame_text(corpus_dir(dest), tid, 1)
    disk = readcon.write_con_string(readcon.read_con(str(dest))[1:2])
    assert "2.123457" in stored and "2.123456789" in disk
    monkeypatch.setattr(concorpus, "_generations", {})
    io.savecon(str(dest), _cu_atom(3.0), "a")
    assert _held_generation(dest) == tid
    assert load_frame_text(corpus_dir(dest), generation_id(
        str(dest.resolve()), first, 1), 0) is None


def test_a_file_shorter_than_its_generation_starts_a_new_one(
    tmp_path: Path, monkeypatch
):
    from eon import concorpus
    from eon.concorpus import generation_frame_texts, generation_id, load_frame_text

    (tmp_path / "config.ini").write_text("[Main]\n")
    dest = tmp_path / "movie.con"
    io.savecon(str(dest), _cu_atom(1.0))
    first = dest.read_text()
    for x in (2.0, 3.0, 4.0):
        io.savecon(str(dest), _cu_atom(x), "a")
    concorpus.flush_corpus_appends()
    monkeypatch.setattr(concorpus, "_generations", {})
    dest.write_text(
        readcon.write_con_string(
            [_cu_atom(1.0).to_conframe(), _cu_atom(2.0).to_conframe()]
        )
    )
    io.savecon(str(dest), _cu_atom(5.0), "a")
    resolved = str(dest.resolve())
    assert _held_generation(dest) == generation_id(resolved, first, 1)
    assert _xs(generation_frame_texts(dest)) == [1.0, 2.0, 5.0]
    old = generation_id(resolved, first)
    assert _xs([load_frame_text(corpus_dir(dest), old, i) for i in range(4)]) == [
        1.0,
        2.0,
        3.0,
        4.0,
    ]


def test_a_new_process_does_not_extend_a_generation_the_file_left(
    tmp_path: Path, monkeypatch
):
    from eon import concorpus
    from eon.concorpus import generation_frame_texts

    (tmp_path / "config.ini").write_text("[Main]\n")
    dest = tmp_path / "movie.con"
    io.savecon(str(dest), _cu_atom(1.0))
    io.savecon(str(dest), _cu_atom(2.0), "a")
    concorpus.flush_corpus_appends()
    monkeypatch.setattr(concorpus, "_generations", {})
    dest.write_text(
        readcon.write_con_string(
            [_cu_atom(1.0).to_conframe(), _cu_atom(5.0).to_conframe()]
        )
    )
    io.savecon(str(dest), _cu_atom(6.0), "a")
    assert _xs(generation_frame_texts(dest)) == [1.0, 5.0, 6.0]


def test_a_changed_middle_frame_starts_a_new_generation(
    tmp_path: Path, monkeypatch
):
    from eon import concorpus
    from eon.concorpus import generation_frame_texts, generation_id, load_frame_text

    (tmp_path / "config.ini").write_text("[Main]\n")
    dest = tmp_path / "movie.con"
    # A whole-file seed stores each frame's own text, so every stored
    # frame here is round-trip text rather than six-decimal text.
    texts = [
        readcon.write_con_string([_cu_atom(x).to_conframe()])
        for x in (1.0, 2.0, 3.0, 4.0, 6.0)
    ]
    dest.write_text("".join(texts))
    resolved = str(dest.resolve())
    db = concorpus._corpus(corpus_dir(dest))
    assert db is not None
    old = generation_id(resolved, texts[0])
    db.append_trajectory_str(old, "".join(texts), source=resolved)
    monkeypatch.setattr(concorpus, "_generations", {})
    replaced = list(texts)
    replaced[1] = readcon.write_con_string([_cu_atom(5.0).to_conframe()])
    dest.write_text("".join(replaced))
    io.savecon(str(dest), _cu_atom(7.0), "a")
    assert _disk_xs(dest) == [1.0, 5.0, 3.0, 4.0, 6.0, 7.0]
    assert _xs(generation_frame_texts(dest)) == [1.0, 5.0, 3.0, 4.0, 6.0, 7.0]
    assert _xs([load_frame_text(corpus_dir(dest), old, i) for i in range(5)]) == [
        1.0,
        2.0,
        3.0,
        4.0,
        6.0,
    ]
    assert load_frame_text(corpus_dir(dest), old, 5) is None
