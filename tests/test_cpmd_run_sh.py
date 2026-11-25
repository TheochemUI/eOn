"""The CPMD point launcher picks its converter without Bash 4 syntax."""
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

RUN_SH = Path(__file__).resolve().parents[1] / "examples/ext_pot_cpmd/run.sh"

pytestmark = pytest.mark.skipif(
    sys.platform == "win32", reason="POSIX shell launcher"
)


def _suffix_block() -> str:
    match = re.search(r"^case .*?^esac$", RUN_SH.read_text(), re.M | re.S)
    assert match, "run.sh selects its converter outside a case block"
    return match.group(0)


@pytest.mark.parametrize(
    ("name", "converter"),
    [
        ("IN.XYZ", "xyz_to_con.py"),
        ("in.xyz", "xyz_to_con.py"),
        ("In.XyZ", "xyz_to_con.py"),
        ("in.con", "con_to_input.py"),
        ("xyz.con", "con_to_input.py"),
        ("in.xyz.con", "con_to_input.py"),
    ],
)
def test_the_suffix_block_runs_under_posix_sh(name, converter):
    # /bin/sh is dash on Debian and Bash 3.2 in POSIX mode on macOS; both
    # reject ${input,,} with "Bad substitution".
    script = f'set -eu\ninput="$1"\n{_suffix_block()}\nprintf %s "$converter"\n'
    out = subprocess.run(
        ["/bin/sh", "-c", script, "sh", f"/some dir/{name}"],
        capture_output=True,
        text=True,
    )
    assert out.returncode == 0, out.stderr
    assert out.stdout == converter


def _stub(path: Path, body: str) -> None:
    path.write_text(f"#!/bin/sh\n{body}\n")
    path.chmod(0o755)


@pytest.mark.skipif(shutil.which("bash") is None, reason="needs bash")
@pytest.mark.parametrize(
    ("name", "converter"),
    [("IN.XYZ", "xyz_to_con.py"), ("in.con", "con_to_input.py")],
)
def test_run_sh_hands_the_input_to_its_converter(tmp_path, name, converter):
    example = tmp_path / "example"
    example.mkdir()
    shutil.copy(RUN_SH, example / "run.sh")
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    uv_log = tmp_path / "uv.log"
    _stub(bin_dir / "uv", 'printf "%s\\n" "${3##*/}" >> "$UV_LOG"\n: > "$5"')
    _stub(bin_dir / "eonclient", "printf '%s potential_energy\\n' -1.5 > results.dat")
    _stub(bin_dir / "cpmd.x", "exit 0")
    env_file = tmp_path / ".env"
    env_file.write_text(
        f"CPMD_BIN={bin_dir / 'cpmd.x'}\n"
        "PP_LIBRARY_PATH=/pp\nCPMD_PP='14:Si.psp'\nCPMD_MPI_LAUNCH=env\n"
        f"EON_CLIENT={bin_dir / 'eonclient'}\n"
    )
    source = tmp_path / name
    source.write_text("geometry\n")
    env = dict(
        os.environ,
        PATH=f"{bin_dir}{os.pathsep}{os.environ['PATH']}",
        EON_ENV_FILE=str(env_file),
        UV_LOG=str(uv_log),
    )
    out = subprocess.run(
        ["bash", str(example / "run.sh"), str(source)],
        capture_output=True,
        text=True,
        env=env,
    )
    assert out.returncode == 0, out.stderr
    assert uv_log.read_text().split() == [converter]
    assert "Energy: -1.5000000000 eV" in out.stdout
