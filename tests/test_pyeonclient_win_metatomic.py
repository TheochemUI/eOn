"""pyeonclient wheels include a win-64 metatomic build."""

from pathlib import Path

WORKFLOW = (
    Path(__file__).resolve().parents[1]
    / ".github"
    / "workflows"
    / "pyeonclient-wheels.yml"
)


def test_win64_metatomic_wheel_is_in_the_matrix():
    text = WORKFLOW.read_text(encoding="utf-8")
    assert "name: abi3-win\n" in text
    assert "name: abi3-win-metatomic\n" in text
    assert "os: windows-2022" in text
    assert 'cibw_build: "cp312-win_amd64"' in text
    assert "name: abi3-manylinux-metatomic" in text
    assert "CIBW_BEFORE_BUILD_LINUX" in text
    assert "CIBW_BEFORE_BUILD_WINDOWS" in text
    assert "pyeonclient_cibw_before_windows.py ${{ matrix.variant }}" in text
    before = (
        Path(__file__).resolve().parents[1]
        / "scripts"
        / "pyeonclient_cibw_before_windows.py"
    ).read_text(encoding="utf-8")
    assert "numpy" in before
    assert "meson>=1.3" in before
    assert "<1.11" not in before
    pyproject = (
        Path(__file__).resolve().parents[1] / "pyproject-pyeonclient.toml"
    ).read_text(encoding="utf-8")
    assert 'meson>=1.3"' in pyproject
    assert "<1.11" not in pyproject
    assert 'variant == "metatomic"' in before
    win = text.split("name: abi3-win-metatomic", 1)[1].split("steps:", 1)[0]
    assert "variant: metatomic" in win
    assert "windows-2022" in win


def test_windows_pots_are_static_and_flang_runtime_is_a_file():
    root = Path(__file__).resolve().parents[1]
    extpot = (root / "client" / "potentials" / "ExtPot" / "meson.build").read_text(
        encoding="utf-8"
    )
    lammps = (root / "client" / "potentials" / "LAMMPS" / "meson.build").read_text(
        encoding="utf-8"
    )
    bindings = (root / "client" / "python" / "meson.build").read_text(encoding="utf-8")
    client = (root / "client" / "meson.build").read_text(encoding="utf-8")
    for meson in (extpot, lammps):
        assert "if is_windows" in meson
        assert "static_library(" in meson
        assert "get_option('default_library') != 'shared'" in meson
    assert "static_library(\n        'eonclib'," in client
    assert "'-Wno-unused-parameter'" not in bindings
    wrap = (root / "subprojects" / "rgpot.wrap").read_text(encoding="utf-8")
    assert "rgpot-windows-flang-rt.patch" in wrap
    patch = (
        root / "subprojects" / "packagefiles" / "rgpot-windows-flang-rt.patch"
    ).read_text(encoding="utf-8")
    assert "FortranRuntime.dynamic.lib" in patch
    assert "FortranDecimal.dynamic.lib" in patch
    assert "find_library never sees the runtime" in patch
    whole = (
        root / "subprojects" / "packagefiles" / "rgpot-windows-link-whole.patch"
    ).read_text(encoding="utf-8")
    assert "rgpot-windows-link-whole.patch" in wrap
    assert "link_whole: [rgpot_core] + _linkto + _rpc_whole" in whole
    assert "_readcon_link_dep" in client
    assert "dependencies: _deps + [_readcon_link_dep]" in client
    pyproject = (root / "pyproject-pyeonclient.toml").read_text(encoding="utf-8")
    assert "allow-windows-internal-shared-libs = true" in pyproject
    init = (
        root / "client" / "python" / "pyeonclient" / "__init__.py"
    ).read_text(encoding="utf-8")
    assert "os.add_dll_directory" in init
    assert ".pyeonclient.mesonpy.libs" in init
    assert "metatensor_torch" in init
    assert "_dirs_with_dlls" in init
    assert "output: 'readcon_db.dll'" in bindings
