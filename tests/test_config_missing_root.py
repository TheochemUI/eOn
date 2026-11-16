"""A missing main_directory still reaches the confirmation prompt."""

from __future__ import annotations

from eon.config import ConfigClass


def test_missing_main_directory_still_prompts(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    missing = tmp_path / "no-such-root"
    (tmp_path / "config.ini").write_text(
        "[Main]\njob = akmc\n[Paths]\nmain_directory = %s\n" % missing
    )
    seen = []

    def answer(prompt):
        seen.append(prompt)
        return "y"

    monkeypatch.setattr("builtins.input", answer)
    cfg = ConfigClass()
    cfg.init()
    assert seen
    assert cfg.path_root == str(missing)
