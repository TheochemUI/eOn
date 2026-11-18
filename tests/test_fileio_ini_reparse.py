"""``fileio.ini`` parses a file when its contents change, not on every get.

A state's info file is read through ``ini.get`` several hundred times in
one server cycle. The parsed values must stay what a fresh read of the
file gives, including after another writer changes the file in place.
"""

from __future__ import annotations

import configparser
import sys
from pathlib import Path

import pytest

_TESTS = Path(__file__).resolve().parent
if str(_TESTS) not in sys.path:
    sys.path.insert(0, str(_TESTS))

from driver_oracles import count_events  # noqa: E402

from eon import fileio as io  # noqa: E402


def parses(counts):
    return counts["py:RawConfigParser._read"]


def test_repeated_gets_parse_the_file_once(tmp_path):
    path = tmp_path / "info"
    path.write_text("[MetaData]\nreactant energy = -10318.5\nprevious state = 3\n")
    info = io.ini(str(path))
    with count_events() as counts:
        for _ in range(50):
            assert info.get("MetaData", "reactant energy") == -10318.5
            assert info.get("MetaData", "previous state") == 3
    assert parses(counts) == 1


def test_a_file_rewritten_in_place_with_the_same_size_is_seen(tmp_path):
    path = tmp_path / "info"
    path.write_text("[MetaData]\nvalue = 1.5\n")
    info = io.ini(str(path))
    assert info.get("MetaData", "value") == 1.5
    path.write_text("[MetaData]\nvalue = 2.5\n")
    assert info.get("MetaData", "value") == 2.5
    path.write_text("[MetaData]\nvalue = 9.5\n")
    assert info.get("MetaData", "value") == 9.5


def test_a_file_replaced_by_another_writer_is_seen(tmp_path):
    path = tmp_path / "info"
    path.write_text("[A]\nx = 1\n")
    info = io.ini(str(path))
    assert info.get("A", "x") == 1
    other = io.ini(str(path))
    other.set("A", "x", 2)
    other.set("B", "y", "hello")
    assert info.get("A", "x") == 2
    assert info.get("B", "y") == "hello"


def test_set_is_visible_to_the_next_get_and_costs_one_parse_at_most(tmp_path):
    path = tmp_path / "info"
    info = io.ini(str(path))
    assert info.get("MetaData", "reactant energy", None) is None
    info.set("MetaData", "reactant energy", -3.25)
    with count_events() as counts:
        for _ in range(20):
            assert info.get("MetaData", "reactant energy") == -3.25
    assert parses(counts) <= 1
    assert io.ini(str(path)).get("MetaData", "reactant energy") == -3.25


def test_a_file_that_appears_later_is_seen_and_a_missing_one_gives_the_default(
    tmp_path,
):
    path = tmp_path / "info"
    info = io.ini(str(path))
    assert info.get("A", "x", "fallback") == "fallback"
    assert info.get("A", "x", "fallback") == "fallback"
    path.write_text("[A]\nx = 7\n")
    assert info.get("A", "x", "fallback") == 7
    with pytest.raises(NameError):
        info.get("A", "missing")


def test_a_file_that_does_not_parse_does_not_poison_the_next_read(tmp_path):
    path = tmp_path / "info"
    path.write_text("no section header\n")
    info = io.ini(str(path))
    with pytest.raises(configparser.MissingSectionHeaderError):
        info.get("A", "x", "fallback")
    assert info.get("A", "x", "fallback") == "fallback"
    path.write_text("[A]\nx = 4\n")
    assert info.get("A", "x") == 4


def test_a_later_file_overrides_an_earlier_one(tmp_path):
    first, second = tmp_path / "a.ini", tmp_path / "b.ini"
    first.write_text("[S]\nk = 1\nonly_first = yes\n")
    second.write_text("[S]\nk = 2\n")
    info = io.ini([str(first), str(second)])
    assert info.get("S", "k") == 2
    assert info.get("S", "only_first") == "yes"
    second.write_text("[S]\nk = 3\n")
    assert info.get("S", "k") == 3


def test_values_keep_their_types(tmp_path):
    path = tmp_path / "info"
    path.write_text("[T]\ni = 12\nf = 1.5e-3\nyes = True\nno = false\ns = text\n")
    info = io.ini(str(path))
    got = [info.get("T", k) for k in ("i", "f", "yes", "no", "s")]
    assert got == [12, 1.5e-3, True, False, "text"]
    assert [type(v) for v in got] == [int, float, bool, bool, str]
