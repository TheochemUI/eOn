"""Tests for the displacement atom list script feature."""

import logging
import os
from pathlib import Path
from unittest import mock

import pytest

from eon import _utils as utl


class TestScriptConfig:
    def test_from_eon_config(self, tmp_path):
        """Verify ScriptConfig.from_eon_config creates correct paths."""
        config = mock.Mock()
        config.displace_atom_kmc_state_script = "my_script.py"
        config.path_scratch = str(tmp_path / "scratch")
        config.path_root = str(tmp_path)

        sc = utl.ScriptConfig.from_eon_config(config)

        # script_path should be resolved against root since it's relative
        assert sc.script_path == tmp_path / "my_script.py"
        assert sc.root_path == tmp_path

    def test_resolves_relative_paths(self, tmp_path):
        """Verify __post_init__ resolves relative paths against root_path."""
        sc = utl.ScriptConfig(
            script_path=Path("scripts/run.py"),
            scratch_path=Path("scratch"),
            root_path=tmp_path,
        )
        assert sc.script_path == tmp_path / "scripts" / "run.py"
        assert sc.scratch_path == tmp_path / "scratch"

    def test_absolute_paths_unchanged(self, tmp_path):
        """Verify absolute paths are not modified by __post_init__."""
        abs_script = tmp_path / "abs_script.py"
        abs_scratch = tmp_path / "abs_scratch"
        sc = utl.ScriptConfig(
            script_path=abs_script,
            scratch_path=abs_scratch,
            root_path=Path("/some/other/root"),
        )
        assert sc.script_path == abs_script
        assert sc.scratch_path == abs_scratch


class TestGenIdsFromCon:
    def test_success(self, tmp_path):
        """Mock subprocess to return atom indices, verify string output."""
        script = tmp_path / "script.py"
        script.write_text("pass")
        scratch = tmp_path / "scratch"
        scratch.mkdir()

        sconf = utl.ScriptConfig(
            script_path=script,
            scratch_path=scratch,
            root_path=tmp_path,
        )

        fake_proc = mock.Mock()
        fake_proc.stdout = "  1 3 5  \n"
        reactant = mock.Mock()
        lgr = logging.getLogger("test")

        with mock.patch("eon._utils.subprocess.run", return_value=fake_proc):
            with mock.patch("eon._utils.eio.savecon"):
                result = utl.gen_ids_from_con(sconf, reactant, lgr)

        assert result == "1 3 5"
        assert isinstance(result, str)

    def test_missing_script(self, tmp_path):
        """Verify returns '' when script doesn't exist."""
        sconf = utl.ScriptConfig(
            script_path=tmp_path / "nonexistent.py",
            scratch_path=tmp_path / "scratch",
            root_path=tmp_path,
        )
        lgr = logging.getLogger("test")

        result = utl.gen_ids_from_con(sconf, mock.Mock(), lgr)
        assert result == ""


class TestStateGetDisplacementAtomList:
    def _make_state(self, tmp_path):
        """Create a minimal State with a pre-existing directory."""
        from eon import fileio as io
        from eon.state import State

        state_dir = tmp_path / "state_0"
        state_dir.mkdir()
        (state_dir / "procdata").mkdir()

        # Create a minimal State without triggering the disk-init branch
        # by ensuring the directory already exists.
        statelist = mock.Mock()
        config = mock.Mock()
        state = object.__new__(State)
        state.config = config
        state.statelist = statelist
        state.path = str(state_dir)
        state.number = 0
        state.procs = None
        state.proc_repeat_count = None
        state.procdata_path = os.path.join(state.path, "procdata")
        state.reactant_path = os.path.join(state.path, "reactant.con")
        state.proctable_path = os.path.join(state.path, "processtable")
        state.search_result_path = os.path.join(state.path, "search_results.txt")
        state.tar_path = os.path.join(state.path, "procdata.tar")
        state.info = io.ini(os.path.join(state.path, "info"))
        return state

    def test_caches_result(self, tmp_path):
        """Verify the result is cached in state.info and second call reads from cache."""
        state = self._make_state(tmp_path)
        state.get_reactant = mock.Mock(return_value=mock.Mock())

        config = mock.Mock()
        config.displace_atom_kmc_state_script = "script.py"
        config.path_scratch = str(tmp_path / "scratch")
        config.path_root = str(tmp_path)

        with mock.patch("eon._utils.gen_ids_from_con", return_value="0 2 4") as mock_gen:
            result1 = state.get_displacement_atom_list(config)
            assert result1 == "0 2 4"
            assert mock_gen.call_count == 1

            # Second call should read from cache
            result2 = state.get_displacement_atom_list(config)
            assert result2 == "0 2 4"
            # gen_ids_from_con should NOT be called again
            assert mock_gen.call_count == 1

    def test_reads_cache(self, tmp_path):
        """Verify that when info already has the value, no script is run."""
        state = self._make_state(tmp_path)
        # Pre-populate the cache
        state.info.set("Saddle Search", "displace_atom_list", "10 20 30")

        config = mock.Mock()

        with mock.patch("eon._utils.gen_ids_from_con") as mock_gen:
            result = state.get_displacement_atom_list(config)
            assert result == "10 20 30"
            mock_gen.assert_not_called()

    def test_empty_script_output(self, tmp_path):
        """Verify returns '' when script produces no output."""
        state = self._make_state(tmp_path)
        state.get_reactant = mock.Mock(return_value=mock.Mock())

        config = mock.Mock()
        config.displace_atom_kmc_state_script = "script.py"
        config.path_scratch = str(tmp_path / "scratch")
        config.path_root = str(tmp_path)

        with mock.patch("eon._utils.gen_ids_from_con", return_value=""):
            result = state.get_displacement_atom_list(config)
            assert result == ""
            assert isinstance(result, str)
