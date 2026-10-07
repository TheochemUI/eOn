"""Matter and Potential stay off the parameter header."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INC = ROOT / "include" / "eon"


def quoted_includes(text: str) -> list[str]:
    names = []
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped.startswith("#include") or '"' not in stripped:
            continue
        names.append(stripped.split('"', 2)[1])
    return names


def closure(start: str) -> set[str]:
    seen: set[str] = set()
    stack = [start]
    while stack:
        name = stack.pop()
        if name in seen:
            continue
        path = INC / name
        if not path.is_file():
            continue
        seen.add(name)
        for child in quoted_includes(path.read_text(encoding="utf-8")):
            stack.append(child)
    return seen


def test_matter_closure_omits_parameters():
    seen = closure("Matter.h")
    assert "Parameters.h" not in seen
    assert "ParametersOptions.h" not in seen


def test_potential_closure_omits_parameters():
    seen = closure("Potential.h")
    assert "Parameters.h" not in seen
    assert "ParametersOptions.h" not in seen


def test_job_closure_includes_parameters():
    seen = closure("Job.h")
    assert "Parameters.h" in seen
    assert "ParametersOptions.h" in seen
