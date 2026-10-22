"""Letter spellings of booleans must load."""

from __future__ import annotations

import pytest

from eon.config import ConfigClass


def _load(tmp_path, raw):
    path = tmp_path / "config.ini"
    path.write_text("[Main]\ncheckpoint = %s\n" % raw)
    cfg = ConfigClass()
    cfg.init(str(path))
    return cfg


@pytest.mark.parametrize("raw", ["T", "t"])
def test_letter_true_sets_checkpoint(tmp_path, raw):
    cfg = _load(tmp_path, raw)
    assert cfg.main_checkpoint is True


@pytest.mark.parametrize("raw", ["F", "f"])
def test_letter_false_clears_checkpoint(tmp_path, raw):
    cfg = _load(tmp_path, raw)
    assert cfg.main_checkpoint is False
