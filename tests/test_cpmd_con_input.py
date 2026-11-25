"""Preserve geometry while removing stored evaluations from CPMD point inputs."""
import importlib.util
from pathlib import Path

import pytest
import readcon


SOURCE = Path(__file__).resolve().parents[1] / "examples/ext_pot_cpmd/con_to_input.py"
SPEC = importlib.util.spec_from_file_location("cpmd_con_input", SOURCE)
CON_INPUT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CON_INPUT)


def cached_frame(offset=0.0):
    atoms = [readcon.Atom("Si", 1.0+offset, 2.0, 3.0, atom_id=11, mass=28.0855,
                          fixed=[False, True, False], fx=1.0, fy=2.0, fz=3.0, energy=4.0),
             readcon.Atom("N", 4.0+offset, 5.0, 6.0, atom_id=21, mass=14.0067,
                          fixed=[True, True, True], fx=4.0, fy=5.0, fz=6.0, energy=7.0)]
    frame = readcon.ConFrame(cell=[12.0, 13.0, 14.0], angles=[90.0]*3, atoms=atoms)
    frame.set_energy(-123.45)
    return frame


def test_cached_evaluations_are_removed_with_exact_geometry(tmp_path):
    source, output = tmp_path/"cached.con", tmp_path/"pos.con"
    readcon.write_con(str(source), [cached_frame()], precision=17)
    source_bytes = source.read_bytes()
    original = readcon.read_first_frame(str(source))
    assert original.energy == -123.45 and original.has_forces
    assert all(atom.has_energy for atom in original.atoms)
    record = CON_INPUT.prepare(source, output)
    clean = readcon.read_first_frame(str(output))
    assert clean.energy is None and not clean.has_forces
    assert all(not atom.has_energy for atom in clean.atoms)
    assert clean.cell == original.cell and clean.angles == original.angles
    for old, new in zip(original.atoms, clean.atoms):
        assert (new.symbol, new.atom_id, new.mass, new.fixed, new.x, new.y, new.z) == (
            old.symbol, old.atom_id, old.mass, old.fixed, old.x, old.y, old.z)
    assert source.read_bytes() == source_bytes
    assert record["input_had_energy"] and record["input_had_forces"]


def test_explicit_frame_selection(tmp_path):
    source = tmp_path/"trajectory.con"
    readcon.write_con(str(source), [cached_frame(), cached_frame(0.5)], precision=17)
    record = CON_INPUT.prepare(source, tmp_path/"last.con", -1)
    assert record["frame"] == 1 and record["frames"] == 2
    assert readcon.read_first_frame(str(tmp_path/"last.con")).atoms[0].x == 1.5


def test_out_of_range_frame_writes_no_input(tmp_path):
    source, output = tmp_path/"source.con", tmp_path/"pos.con"
    readcon.write_con(str(source), [cached_frame()])
    with pytest.raises(ValueError, match="outside"):
        CON_INPUT.prepare(source, output, 2)
    assert not output.exists()


def test_existing_input_is_preserved(tmp_path):
    source, output = tmp_path/"source.con", tmp_path/"pos.con"
    readcon.write_con(str(source), [cached_frame()])
    output.write_text("preserve existing input\n")
    with pytest.raises(FileExistsError, match="overwrite"):
        CON_INPUT.prepare(source, output)
    assert output.read_text() == "preserve existing input\n"
