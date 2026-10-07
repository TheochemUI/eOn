"""Runtime owns its loaders with make_unique."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_runtime_owns_loaders_with_make_unique():
    source = (ROOT / "client" / "Runtime.cpp").read_text(encoding="utf-8")
    for name in (
        "IRAResource",
        "ARTnResource",
        "PluginLoader",
        "MetatomicLoader",
    ):
        assert f"new {name}" not in source
        assert f"std::make_unique<{name}>" in source
