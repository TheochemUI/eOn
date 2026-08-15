"""Packaging: PyPI metadata + real wheel build path (torch-style, not conda-forge)."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PYPROJECT = ROOT / "pyproject-pyeonclient.toml"
BUILD_SCRIPT = ROOT / "scripts" / "pyeonclient_build_wheel.sh"


def _load_toml():
    try:
        import tomllib
    except ImportError:  # pragma: no cover
        import tomli as tomllib  # type: ignore
    return tomllib.loads(PYPROJECT.read_text())


def _have_compiler() -> bool:
    return shutil.which("c++") is not None or shutil.which("g++") is not None


def test_pyproject_identity_and_extras():
    assert PYPROJECT.is_file()
    data = _load_toml()
    proj = data["project"]
    assert proj["name"] == "pyeonclient"
    assert proj["version"] >= "0.3.0"
    deps = proj["dependencies"]
    assert any("numpy" in d for d in deps)
    assert any("readcon" in d for d in deps)
    opt = proj["optional-dependencies"]
    assert "metatomic" in opt
    assert any("torch" in d for d in opt["metatomic"])
    assert any("metatomic" in d for d in opt["metatomic"])
    assert "ase" in opt
    setup = data["tool"]["meson-python"]["args"]["setup"]
    assert "-Dwith_pyeonclient=true" in setup
    assert "-Dwith_rgpot=true" in setup
    assert "-Dwith_metatomic=false" in setup
    assert "-Dinstall_eon_server=false" in setup
    assert "-Db_lto=false" in setup


def test_build_wheel_script_variants_disambiguated():
    assert BUILD_SCRIPT.is_file()
    text = BUILD_SCRIPT.read_text()
    assert "PYEONCLIENT_VARIANT" in text
    assert "with_metatomic=true" in text
    # metatomic wheels must not share the bare PyPI version basename
    assert "+metatomic" in text
    assert "base wheel must not use +metatomic" in text or "local version" in text
    # portable PyPI: vendor capnp / strip host RPATH (torch-style)
    assert "pyeonclient_repair_wheel" in text
    repair = ROOT / "scripts" / "pyeonclient_repair_wheel.sh"
    assert repair.is_file()
    rtext = repair.read_text()
    assert "patchelf" in rtext
    assert "libcapnp" in rtext or "vendor" in rtext
    # A host ~/.cargo/config.toml naming an absolute linker reaches the
    # readcon-core wrap through cargo and the conda cc rejects it, so the
    # script exports RUSTFLAGS rather than inheriting the host's.
    assert "export RUSTFLAGS=" in text


def test_workflow_publish_excludes_metatomic_local_version():
    wf = ROOT / ".github" / "workflows" / "pyeonclient-wheels.yml"
    assert wf.is_file()
    text = wf.read_text()
    assert "+metatomic" in text
    assert "refusing to publish metatomic" in text or "exclude +metatomic" in text
    assert "Tag metatomic wheels" in text


def test_runtime_feature_flags_when_built():
    pyec = pytest.importorskip("pyeonclient")
    assert hasattr(pyec, "built_with_metatomic")
    assert hasattr(pyec, "built_with_rgpot")
    assert isinstance(pyec.built_with_metatomic(), bool)
    assert isinstance(pyec.built_with_rgpot(), bool)
    data = _load_toml()
    # local builds may be 0.3.0; installed wheel may match
    assert pyec.__version__.split("+")[0] == data["project"]["version"].split("+")[0]


@pytest.mark.skipif(
    not _have_compiler() or os.environ.get("PYEONCLIENT_SKIP_WHEEL_BUILD") == "1",
    reason="no C++ compiler or PYEONCLIENT_SKIP_WHEEL_BUILD=1",
)
def test_base_wheel_build_and_import(tmp_path):
    """Drive the real build backend for the base (RGPOT, no metatomic) wheel.

    Installs into an isolated venv and checks:
    - built_with_rgpot() is True
    - built_with_metatomic() is False
    - extension .so has no libtorch NEEDED (thin base product)
    """
    if not BUILD_SCRIPT.is_file():
        pytest.skip("build script missing")

    # Require build tools
    for tool in ("meson", "ninja"):
        if shutil.which(tool) is None:
            # may be available via python -m
            pass

    env = os.environ.copy()
    env["PYEONCLIENT_VARIANT"] = "base"
    # build in-repo dist/
    dist_before = set((ROOT / "dist").glob("*.whl")) if (ROOT / "dist").is_dir() else set()
    proc = subprocess.run(
        ["bash", str(BUILD_SCRIPT)],
        cwd=str(ROOT),
        env=env,
        capture_output=True,
        text=True,
        timeout=int(os.environ.get("PYEONCLIENT_WHEEL_BUILD_TIMEOUT", "3600")),
    )
    # always keep log for auditors
    log_path = Path(
        os.environ.get(
            "PYEONCLIENT_WHEEL_BUILD_LOG",
            str(tmp_path / "base-wheel-build.log"),
        )
    )
    log_path.write_text(
        f"exit={proc.returncode}\n--- stdout ---\n{proc.stdout}\n--- stderr ---\n{proc.stderr}\n"
    )
    if proc.returncode != 0:
        # surface tail of log for pytest
        pytest.fail(
            f"base wheel build failed rc={proc.returncode}\n"
            f"{(proc.stdout + proc.stderr)[-4000:]}"
        )

    whls = sorted((ROOT / "dist").glob("pyeonclient-*.whl"))
    assert whls, "no wheel in dist/"
    whl = whls[-1]
    assert "+metatomic" not in whl.name, whl.name

    # The build script produces the wheel; vendoring the non-system shared
    # libraries into it is a separate step, which is why the script itself
    # reports "repair may have been skipped". The probe below strips
    # LD_LIBRARY_PATH to check the wheel stands on its own, so it has to run
    # against the repaired artifact, the one that ships. This mirrors
    # .github/workflows/pyeonclient-wheels.yml, which runs auditwheel and then
    # this script.
    repair = ROOT / "scripts" / "pyeonclient_repair_wheel.sh"
    if repair.is_file():
        # The script walks NEEDED entries to a fixpoint, but resolves each one
        # through ldd, which searches the default loader path. capnp's own
        # dependencies (libkj, libkj-async) live in the environment's lib
        # directory, so without a search path it reports "cannot resolve
        # NEEDED libkj-async" and leaves them out of the wheel.
        repair_env = os.environ.copy()
        prefix = os.environ.get("CONDA_PREFIX")
        if prefix:
            search = os.pathsep.join(
                p for p in (
                    repair_env.get("PYEONCLIENT_VENDOR_SEARCH_PATHS"),
                    os.path.join(prefix, "lib"),
                ) if p
            )
            repair_env["PYEONCLIENT_VENDOR_SEARCH_PATHS"] = search
        rep = subprocess.run(
            ["bash", str(repair), str(whl)],
            cwd=str(ROOT),
            env=repair_env,
            capture_output=True,
            text=True,
        )
        (tmp_path / "wheel-repair.log").write_text(rep.stdout + "\n" + rep.stderr)
        if rep.returncode != 0:
            pytest.fail(
                f"wheel repair failed rc={rep.returncode}\n"
                f"{(rep.stdout + rep.stderr)[-4000:]}"
            )
        whls = sorted((ROOT / "dist").glob("pyeonclient-*.whl"))
        whl = whls[-1]

    venv = tmp_path / "venv"
    subprocess.run([sys.executable, "-m", "venv", str(venv)], check=True)
    pip = venv / "bin" / "pip"
    py = venv / "bin" / "python"
    subprocess.run(
        [str(pip), "install", "-q", "--upgrade", "pip"],
        check=True,
    )
    # wheel may need numpy/readcon
    subprocess.run(
        [str(pip), "install", "-q", str(whl), "numpy", "readcon"],
        check=True,
        capture_output=True,
        text=True,
    )
    probe = subprocess.run(
        [
            str(py),
            "-c",
            "import pyeonclient as p, pathlib, subprocess, sys\n"
            "print('version', p.__version__)\n"
            "print('rgpot', p.built_with_rgpot())\n"
            "print('metatomic', p.built_with_metatomic())\n"
            "assert p.built_with_rgpot() is True\n"
            "assert p.built_with_metatomic() is False\n"
            "so = pathlib.Path(p.__file__).parent\n"
            "cores = list(so.glob('_core*.so')) + list(so.glob('**/*_core*.so'))\n"
            "print('cores', cores)\n"
            "assert cores, 'missing _core extension'\n"
            "out = subprocess.check_output(['readelf', '-d', str(cores[0])], text=True)\n"
            "print(out)\n"
            "assert 'libtorch' not in out and 'libc10' not in out\n"
            # portable: no build-host absolute RUNPATH/RPATH
            "assert '/home/' not in out and '/pixi/envs' not in out\n"
            # if capnp is NEEDED it must be loadable without host LD_LIBRARY_PATH
            "import os\n"
            "env = {k: v for k, v in os.environ.items() if k not in ('LD_LIBRARY_PATH', 'LIBRARY_PATH')}\n"
            "ldd = subprocess.check_output(['ldd', str(cores[0])], text=True, env=env)\n"
            "print(ldd)\n"
            "assert 'not found' not in ldd, ldd\n"
            "print('BASE_WHEEL_INSTALL_OK')\n",
        ],
        # The probe exists to test the wheel installed in this venv. A
        # PYTHONPATH inherited from the caller shadows it with the source
        # tree, whose extension resolves its libraries from a build
        # directory, so the probe would report on the wrong package.
        env={
            k: v for k, v in os.environ.items()
            if k not in ("PYTHONPATH", "LD_LIBRARY_PATH", "LIBRARY_PATH")
        },
        capture_output=True,
        text=True,
        check=False,
    )
    (tmp_path / "install-probe.log").write_text(
        probe.stdout + "\n" + probe.stderr
    )
    assert probe.returncode == 0, probe.stdout + probe.stderr
    assert "BASE_WHEEL_INSTALL_OK" in probe.stdout
