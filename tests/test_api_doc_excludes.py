"""API docs do not publish the version module or the cassowary test script."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONF = ROOT / "docs" / "source" / "conf.py"
OHTST = ROOT / "include" / "eon" / "OHTSTJob.h"


def test_autodoc_excludes_version_and_test_script():
    text = CONF.read_text(encoding="utf-8")
    block = text.split('"exclude_files"', 1)[1].split("]", 1)[0]
    assert '"version.py"' in block
    assert '"test.py"' in block
    assert not (ROOT / "eon" / "test.py").is_file()


def test_ohtst_file_tag_is_outside_the_class_comment():
    text = OHTST.read_text(encoding="utf-8")
    before_class = text.split("class OHTSTJob", 1)[0]
    class_comment = before_class.rsplit("/**", 1)[-1]
    assert "@file" not in class_comment
    assert "@file" in before_class
