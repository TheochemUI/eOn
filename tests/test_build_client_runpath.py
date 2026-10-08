"""The checked-in client build rejects a /home/ entry in the run path."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_build_client_script_rejects_a_home_runpath():
    script = ROOT / "scripts" / "build-client.sh"
    text = script.read_text(encoding="utf-8")
    assert "--check" in text
    assert "--sanitize" in text
    assert "readelf" in text
    assert "patchelf" in text
    assert "'/home/'" in text or '"/home/"' in text
    assert "eonclient" in text
    assert "libcpmdc.so" in text
    assert "$HOME/var/scratch" not in text
    assert "-march=native" not in text
    assert "/var/lib/builds" not in text
    readme = (ROOT / "eessi" / "README.md").read_text(encoding="utf-8")
    assert "scripts/build-client.sh" in readme
