"""GCC and Clang flags that cl.exe rejects stay behind a compiler check."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PY_MESON = ROOT / "client" / "python" / "meson.build"
CLIENT_MESON = ROOT / "client" / "meson.build"


def test_pyeonclient_does_not_pass_gcc_warning_flags_raw():
    text = PY_MESON.read_text(encoding="utf-8")
    assert "cppc.get_supported_arguments(['-Wno-unused-parameter'])" in text
    assert text.count("'-Wno-unused-parameter'") == 1
    assert "cppc.get_argument_syntax() == 'gcc'" in text
    assert "'-isystem' + _nb_include" in text
    assert "'-isystem' + _nb_robin" in text
    assert "'/I' + _nb_include" in text
    assert "'/I' + _nb_robin" in text


def test_client_keeps_gnu_as_needed_off_msvc_and_macos():
    text = CLIENT_MESON.read_text(encoding="utf-8")
    gate = (
        "_readcon_db_as_needed = cppc.get_argument_syntax() == 'gcc' "
        "and host_system != 'darwin'"
    )
    assert gate in text
    assert text.count("if _readcon_db_as_needed") == 2
    assert "'-Wl,--no-as-needed'" in text
    assert "'-Wl,--as-needed'" in text
    assert "if host_system != 'darwin'" not in text
    # The client warning set is probed, so MSVC does not see -Wno-*.
    assert "cppc.get_supported_arguments(" in text
    assert "'-Wno-unused-parameter'" in text
    # Torch includes use -isystem only for gcc argument syntax.
    assert "_torch_isystem = cppc.get_argument_syntax() == 'gcc'" in text
    assert "cppc.get_supported_arguments(['-w'])" in text
