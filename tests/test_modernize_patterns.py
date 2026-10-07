"""Registry lifetime, C destroy ownership, and the remaining casts."""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CLIENT = ROOT / "client"
VOID_CAST = re.compile(r"(?<![\w])\(void\)(?!\s*\*)")


def _text(rel: str) -> str:
    return (ROOT / rel).read_text(encoding="utf-8")


def test_process_registry_stays_on_the_heap():
    source = _text("client/PotRegistry.cpp")
    assert "new PotRegistry" not in source
    assert "std::make_unique<PotRegistry>().release()" in source
    assert "function-local static would already be destroyed" in source


def test_c_destroy_takes_ownership():
    relax = _text("client/relax/RelaxEngineAbi.cpp")
    assert "std::unique_ptr<EonRelaxEngine>" in relax
    assert "delete eng" not in relax
    assert "new (std::nothrow)" not in relax
    mta = _text("client/potentials/Metatomic/MetatomicCAbi.cpp")
    assert "std::unique_ptr<EonMtaPot>" in mta
    assert "delete pot" not in mta
    engine = _text("client/potentials/Metatomic/MetatomicEngineAbi.cpp")
    assert "std::unique_ptr<RgpotMtaPot>" in engine
    assert "delete pot" not in engine
    meson = _text("client/meson.build")
    pot = meson.split("metatomic_pot = library(", 1)[1].split("potentials +=", 1)[0]
    assert "potentials/Metatomic/MetatomicCAbi.cpp" in pot


def test_startup_and_rgpot_lines_do_not_flush_every_newline():
    assert "std::endl" not in _text("client/CommandLine.cpp")
    assert "std::endl" not in _text("client/potentials/Rgpot/RgpotPot.cpp")


def test_python_bindings_only_placement_new():
    for path in (CLIENT / "python" / "bind").glob("*.cpp"):
        for line in path.read_text(encoding="utf-8").splitlines():
            code = line.split("//", 1)[0]
            if "new " not in code:
                continue
            assert "new (self)" in code, f"{path.name}: {line.strip()}"


def test_no_discard_c_cast_outside_thirdparty():
    hits = []
    for path in CLIENT.rglob("*"):
        if path.suffix not in {".cpp", ".h", ".hpp"}:
            continue
        if "thirdparty" in path.parts:
            continue
        for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            if VOID_CAST.search(line):
                hits.append(f"{path.relative_to(ROOT)}:{number}")
    assert hits == []
