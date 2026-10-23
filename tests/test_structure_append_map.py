"""Structure.append grows the file-order map with the geometry."""

from __future__ import annotations

from eon.structure import Structure


def test_append_grows_file_to_struct():
    structure = Structure(1)
    structure.append([1.0, 0.0, 0.0], 1, "H", 1.0)
    assert list(structure.file_to_struct) == [0, 1]
    assert len(structure.file_to_struct) == len(structure)
