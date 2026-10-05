"""The relax-engine schema unit is published as .cpp for MSVC cl."""

import subprocess
import sys
from pathlib import Path


def test_cpp_replaces_the_capnp_c_plus_plus_unit(tmp_path: Path):
    schema = tmp_path / "schema" / "eon_relax_engine.capnp"
    schema.parent.mkdir()
    schema.write_text("# schema stand-in\n")
    out = tmp_path / "out"
    header = out / "eon_relax_engine.capnp.h"
    source = out / "eon_relax_engine.capnp.cpp"
    fake = tmp_path / "capnp"
    fake.write_text(
        "\n".join(
            [
                "#!/usr/bin/env python3",
                "import sys",
                "from pathlib import Path",
                "outdir = None",
                "schema = None",
                "for arg in sys.argv[1:]:",
                "    if arg.startswith('-oc++:'):",
                "        outdir = Path(arg.split(':', 1)[1])",
                "    elif not arg.startswith('-'):",
                "        schema = Path(arg)",
                "assert outdir is not None and schema is not None",
                "outdir.mkdir(parents=True, exist_ok=True)",
                "(outdir / (schema.name + '.h')).write_text('header')",
                "(outdir / (schema.name + '.c++')).write_text('unit')",
                "",
            ]
        )
    )
    fake.chmod(0o755)
    subprocess.run(
        [
            sys.executable,
            str(Path(__file__).resolve().parents[1] / "scripts" / "capnp_cpp_unit.py"),
            "--capnp",
            str(fake),
            "--import-path",
            str(tmp_path / "include"),
            "--src-prefix",
            str(schema.parent),
            "--schema",
            str(schema),
            "--header",
            str(header),
            "--source",
            str(source),
        ],
        check=True,
    )
    assert header.read_text() == "header"
    assert source.read_text() == "unit"
    assert not (out / "eon_relax_engine.capnp.c++").exists()
