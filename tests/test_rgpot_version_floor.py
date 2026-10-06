"""An installed rgpot 3.3.0 misses the calculator-group version floor."""

import os
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MESON = ROOT / "client" / "meson.build"
GUIDE = ROOT / "docs" / "source" / "user_guide" / "rgpot_pot.md"


def rgpot_pkgconfig_floor(text: str) -> str:
    match = re.search(
        r"dependency\(\s*'rgpot',\s*version:\s*'>=([^']+)'",
        text,
    )
    assert match, "client/meson.build has no rgpot pkg-config floor"
    return match.group(1)


def test_rgpot_3_3_0_fails_the_pkg_config_floor(tmp_path):
    floor = rgpot_pkgconfig_floor(MESON.read_text(encoding="utf-8"))
    pc = tmp_path / "rgpot.pc"
    pc.write_text("Name: rgpot\nDescription: rgpot\nVersion: 3.3.0\n", encoding="utf-8")
    env = os.environ.copy()
    env["PKG_CONFIG_PATH"] = str(tmp_path)
    result = subprocess.run(
        ["pkg-config", "--atleast-version", floor, "rgpot"],
        env=env,
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0


def test_guide_names_the_first_calculator_group_export():
    page = GUIDE.read_text(encoding="utf-8")
    assert "cpmdc_bind_calculator" in page
    assert "rgpot::bindCalculators" in page
    assert "3.3.0" in page
    assert "2.5.0" not in page
