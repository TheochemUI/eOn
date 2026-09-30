"""with_gprd=auto survives a failed gpr_optim fetch, and a boolean build dir can reconfigure."""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import textwrap
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MIGRATE = ROOT / "scripts" / "migrate_with_gprd_option.py"


def _meson() -> str:
    meson = shutil.which("meson")
    if meson is None:
        pytest.fail("meson is not on PATH")
    return meson


def _run(cmd: list[str], cwd: Path) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    env["GIT_TERMINAL_PROMPT"] = "0"
    env.setdefault("GIT_CONFIG_GLOBAL", "/dev/null")
    env.setdefault("GIT_CONFIG_SYSTEM", "/dev/null")
    return subprocess.run(
        cmd,
        cwd=cwd,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
        timeout=90,
    )


def test_client_passes_the_feature_as_required():
    meson = (ROOT / "client" / "meson.build").read_text()
    match = re.search(r"subproject\(\s*'gpr_optim'.*?\)", meson, re.S)
    assert match is not None
    call = match.group(0)
    assert "required:" in call
    assert "_gprd_opt" in call
    assert "libgprd_proj.found()" in meson
    assert "with_gprd=auto: subproject gpr_optim did not configure" in meson


def _write_parent(root: Path, option_text: str, required: str) -> None:
    (root / "meson_options.txt").write_text(option_text)
    (root / "meson.build").write_text(
        textwrap.dedent(
            f"""\
            project('gprd-auto-fallback', version: '0')
            opt = get_option('with_gprd')
            sub = subproject('gpr_optim', required: {required})
            if sub.found()
                error('gpr_optim configured')
            endif
            """
        )
    )
    nested = root / "subprojects" / "gpr_optim"
    nested.mkdir(parents=True)
    (nested / "meson.build").write_text(
        textwrap.dedent(
            """\
            project('gpr_optim', version: '0')
            dependency('minimage', fallback: ['minimage', 'minimage_dep'])
            """
        )
    )
    wrap = nested / "subprojects"
    wrap.mkdir()
    (wrap / "minimage.wrap").write_text(
        textwrap.dedent(
            """\
            [wrap-git]
            directory = minimage
            url = https://127.0.0.1:9/minimage.git
            revision = 7ce993f86a1f00caab57d4aa23a579551d2b5c28
            depth = 1

            [provide]
            minimage = minimage_dep
            """
        )
    )


def test_auto_configure_survives_a_failed_minimage_fetch(tmp_path: Path):
    _write_parent(
        tmp_path,
        "option('with_gprd', type: 'feature', value: 'auto')\n",
        "opt",
    )
    result = _run([_meson(), "setup", str(tmp_path / "bb"), str(tmp_path)], tmp_path)
    assert result.returncode == 0, result.stdout
    assert "gpr_optim configured" not in result.stdout
    lowered = result.stdout.lower()
    assert "buildable" in lowered or "disabling" in lowered


def test_enabled_configure_stops_when_minimage_cannot_be_fetched(tmp_path: Path):
    _write_parent(
        tmp_path,
        "option('with_gprd', type: 'feature', value: 'enabled')\n",
        "opt",
    )
    result = _run(
        [_meson(), "setup", str(tmp_path / "bb"), str(tmp_path), "-Dwith_gprd=enabled"],
        tmp_path,
    )
    assert result.returncode != 0, result.stdout
    assert "7ce993f" in result.stdout or "buildable" in result.stdout or "Git" in result.stdout


@pytest.mark.parametrize(
    ("stored", "mapped"),
    [("true", "enabled"), ("false", "disabled")],
)
def test_boolean_build_dir_reconfigures_after_migration(
    tmp_path: Path, stored: str, mapped: str
):
    src = tmp_path / "src"
    src.mkdir()
    (src / "meson.build").write_text("project('bool-gprd', version: '0')\n")
    (src / "meson_options.txt").write_text(
        f"option('with_gprd', type: 'boolean', value: {stored})\n"
    )
    build = tmp_path / "bb"
    setup = _run([_meson(), "setup", str(build), str(src)], src)
    assert setup.returncode == 0, setup.stdout

    (src / "meson_options.txt").write_text(
        "option('with_gprd', type: 'feature', value: 'auto',\n"
        "       description: 'GP dimer')\n"
    )
    bare = _run([_meson(), "setup", "--reconfigure", str(build)], src)
    assert bare.returncode != 0, bare.stdout
    assert "not boolean" in bare.stdout

    migrated = _run([_meson_python(), str(MIGRATE), str(build)], src)
    assert migrated.returncode == 0, migrated.stdout
    assert mapped in migrated.stdout

    again = _run([_meson(), "setup", "--reconfigure", str(build)], src)
    assert again.returncode == 0, again.stdout
    intro = _run([_meson(), "introspect", "--buildoptions", str(build)], src)
    assert intro.returncode == 0, intro.stdout
    options = json.loads(intro.stdout)
    found = [opt for opt in options if opt.get("name") == "with_gprd"]
    assert found and found[0]["value"] == mapped


def _meson_python() -> str:
    import sys

    try:
        import mesonbuild  # noqa: F401

        return sys.executable
    except ImportError:
        meson = Path(_meson())
        first = meson.read_text(encoding="utf-8", errors="replace").splitlines()[0]
        if first.startswith("#!"):
            return first[2:].strip().split()[0]
        return sys.executable
