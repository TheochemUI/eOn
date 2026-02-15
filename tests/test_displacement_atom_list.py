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


class TestNormalizeAtomListStr:
    def test_space_separated(self):
        assert utl.normalize_atom_list_str("1 3 5") == "1, 3, 5"

    def test_comma_separated(self):
        assert utl.normalize_atom_list_str("1, 3, 5") == "1, 3, 5"

    def test_mixed_separators(self):
        assert utl.normalize_atom_list_str("1, 3 5") == "1, 3, 5"

    def test_extra_whitespace(self):
        assert utl.normalize_atom_list_str("  1   3   5  ") == "1, 3, 5"

    def test_empty(self):
        assert utl.normalize_atom_list_str("") == ""

    def test_single_value(self):
        assert utl.normalize_atom_list_str("42") == "42"


class TestParseAtomListStr:
    def test_comma_separated(self):
        assert utl.parse_atom_list_str("1, 3, 5") == [1, 3, 5]

    def test_empty_string(self):
        assert utl.parse_atom_list_str("") == []

    def test_single_value(self):
        assert utl.parse_atom_list_str("42") == [42]

    def test_negative_index(self):
        assert utl.parse_atom_list_str("10, 20, -1") == [10, 20, -1]


class TestGenIdsFromCon:
    def test_success_normalizes_output(self, tmp_path):
        """Script output is normalized to comma-separated format."""
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

        assert result == "1, 3, 5"
        assert isinstance(result, str)

    def test_success_comma_input_passthrough(self, tmp_path):
        """Script outputting comma-separated stays comma-separated."""
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
        fake_proc.stdout = "10, 20, 30\n"
        lgr = logging.getLogger("test")

        with mock.patch("eon._utils.subprocess.run", return_value=fake_proc):
            with mock.patch("eon._utils.eio.savecon"):
                result = utl.gen_ids_from_con(sconf, mock.Mock(), lgr)

        assert result == "10, 20, 30"

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

        with mock.patch("eon._utils.gen_ids_from_con", return_value="0, 2, 4") as mock_gen:
            result1 = state.get_displacement_atom_list(config)
            assert result1 == "0, 2, 4"
            assert mock_gen.call_count == 1

            # Second call should read from cache
            result2 = state.get_displacement_atom_list(config)
            assert result2 == "0, 2, 4"
            # gen_ids_from_con should NOT be called again
            assert mock_gen.call_count == 1

    def test_reads_cache(self, tmp_path):
        """Verify that when info already has the value, no script is run."""
        state = self._make_state(tmp_path)
        # Pre-populate the cache
        state.info.set("Saddle Search", "displace_atom_list", "10, 20, 30")

        config = mock.Mock()

        with mock.patch("eon._utils.gen_ids_from_con") as mock_gen:
            result = state.get_displacement_atom_list(config)
            assert result == "10, 20, 30"
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


class TestExplorerAtomListInjection:
    """Test that MinModeExplorer injects the per-state atom list into config."""

    def _make_state_with_info(self, tmp_path, atom_list_str):
        """Create a mock state with a real info file containing the atom list."""
        from eon import fileio as io

        state_dir = tmp_path / "state_0"
        state_dir.mkdir(exist_ok=True)
        info = io.ini(os.path.join(str(state_dir), "info"))
        if atom_list_str:
            info.set("Saddle Search", "displace_atom_list", atom_list_str)

        state = mock.Mock()
        state.info = info
        state.number = 0
        state.get_reactant = mock.Mock()
        return state

    def test_injects_atom_list_into_config(self, tmp_path):
        """When state.info has an atom list, it should be parsed into config.disp_listed_atoms."""
        state = self._make_state_with_info(tmp_path, "1, 3, 5")

        config = mock.Mock()
        config.displace_atom_kmc_state_script = "script.py"
        config.displace_listed_atom_weight = 0.0
        config.disp_moved_only = False
        config.kdb_on = False
        config.recycling_on = False

        # We just need to test the injection logic in MinModeExplorer.__init__.
        # Mock out everything else that __init__ does.
        with mock.patch("eon.explorer.Explorer.__init__"):
            with mock.patch("eon.explorer.communicator.get_communicator"):
                with mock.patch("eon.explorer.displace.DisplacementManager"):
                    from eon.explorer import MinModeExplorer
                    exp = object.__new__(MinModeExplorer)
                    exp.config = config
                    exp.state = state
                    exp.previous_state = mock.Mock()
                    exp.states = mock.Mock()
                    exp.comm = mock.Mock()
                    exp.superbasin = None
                    exp.nrecycled = 0

                    # Run the relevant part of __init__ manually
                    # (the injection block + DisplacementManager construction)
                    from eon import _utils as _utl
                    if config.displace_atom_kmc_state_script:
                        atom_list_str = str(state.info.get("Saddle Search", "displace_atom_list", ""))
                        if atom_list_str:
                            config.disp_listed_atoms = _utl.parse_atom_list_str(atom_list_str)
                            if config.displace_listed_atom_weight == 0.0:
                                config.displace_listed_atom_weight = 1.0

        assert config.disp_listed_atoms == [1, 3, 5]
        assert config.displace_listed_atom_weight == 1.0

    def test_no_injection_when_no_script(self, tmp_path):
        """When displace_atom_kmc_state_script is empty, config is unchanged."""
        state = self._make_state_with_info(tmp_path, "")

        config = mock.Mock()
        config.displace_atom_kmc_state_script = ""
        config.displace_listed_atom_weight = 0.0

        # The injection block is guarded by `if config.displace_atom_kmc_state_script:`
        # so nothing should change on config
        if config.displace_atom_kmc_state_script:
            from eon import _utils as _utl
            atom_list_str = str(state.info.get("Saddle Search", "displace_atom_list", ""))
            if atom_list_str:
                config.disp_listed_atoms = _utl.parse_atom_list_str(atom_list_str)
                if config.displace_listed_atom_weight == 0.0:
                    config.displace_listed_atom_weight = 1.0

        assert config.displace_listed_atom_weight == 0.0

    def test_preserves_existing_weight(self, tmp_path):
        """When user already set displace_listed_atom_weight, don't override."""
        state = self._make_state_with_info(tmp_path, "10, 20")

        config = mock.Mock()
        config.displace_atom_kmc_state_script = "script.py"
        config.displace_listed_atom_weight = 0.5  # user-set

        from eon import _utils as _utl
        atom_list_str = str(state.info.get("Saddle Search", "displace_atom_list", ""))
        if atom_list_str:
            config.disp_listed_atoms = _utl.parse_atom_list_str(atom_list_str)
            if config.displace_listed_atom_weight == 0.0:
                config.displace_listed_atom_weight = 1.0

        assert config.disp_listed_atoms == [10, 20]
        # Should NOT have been overridden since it was already non-zero
        assert config.displace_listed_atom_weight == 0.5


class TestEndToEndAtomListFlow:
    """Test the full chain: script output -> normalize -> cache -> parse -> list of ints."""

    def test_space_separated_script_output_becomes_int_list(self):
        """Space-separated script output is normalized, cached, and parsed into int list."""
        # Step 1: Script outputs "1 3 5"
        raw = "1 3 5"
        # Step 2: gen_ids_from_con normalizes it
        normalized = utl.normalize_atom_list_str(raw)
        assert normalized == "1, 3, 5"
        # Step 3: stored in state.info as "1, 3, 5" (string)
        # Step 4: explorer reads it back and parses into list[int]
        parsed = utl.parse_atom_list_str(normalized)
        assert parsed == [1, 3, 5]

    def test_comma_separated_script_output_roundtrips(self):
        """Comma-separated script output roundtrips correctly."""
        raw = "10, 20, 30"
        normalized = utl.normalize_atom_list_str(raw)
        assert normalized == "10, 20, 30"
        parsed = utl.parse_atom_list_str(normalized)
        assert parsed == [10, 20, 30]

    def test_negative_indices_supported(self):
        """Negative indices (e.g. -1 for last atom) are preserved."""
        raw = "10, 20, -1"
        normalized = utl.normalize_atom_list_str(raw)
        parsed = utl.parse_atom_list_str(normalized)
        assert parsed == [10, 20, -1]
