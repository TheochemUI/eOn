"""The catalog saddle frame keeps the process mode in displacements."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_insert_writes_the_mode_into_the_saddle_frame(tmp_path):
    src = (ROOT / "eon" / "process_catalog.py").read_text(encoding="utf-8")
    start = src.find("def insert(")
    end = src.find("\ndef ", start + 1)
    body = src[start:end]
    call = body.find("saddle_text = _saddle_text_with_mode(")
    store = body.find("store_frame_text(")
    assert call != -1 and store != -1 and call < store

    import numpy as np

    from eon import fileio as io
    from eon import process_catalog as catalog
    from eon.structure import Structure, structure_order

    atoms = Structure(2)
    atoms.names = ["Pt", "Pt"]
    atoms.mass[:] = 195.084
    atoms.box = np.diag([20.0, 20.0, 20.0])
    atoms.r[0] = [1.0, 2.0, 3.0]
    atoms.r[1] = [2.4, 2.0, 3.0]
    mode = np.array([[0.2, -0.1, 0.0], [0.0, 0.3, -0.2]], dtype=float)
    path = tmp_path / "saddle.con"
    io.savecon(str(path), atoms)
    text = catalog._saddle_text_with_mode(path.read_text(encoding="utf-8"), mode)
    assert "displacements" in text
    import readcon

    frame = readcon.read_con_string(text)[0]
    assert frame.disp is not None
    ids = np.array([atom.atom_id for atom in frame.atoms], dtype=np.uint64)
    assert np.allclose(np.asarray(frame.disp)[structure_order(ids)], mode)
