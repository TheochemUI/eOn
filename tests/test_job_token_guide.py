"""The guide's job spellings, basin-hopping factor, and communicator path."""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HEADER = ROOT / "include" / "eon" / "BaseStructures.h"
INDEX = ROOT / "docs" / "source" / "user_guide" / "index.md"
BASIN_CPP = ROOT / "client" / "BasinHoppingJob.cpp"
BASIN_MD = ROOT / "docs" / "source" / "user_guide" / "basin_hopping.md"
COMM = ROOT / "docs" / "source" / "user_guide" / "communicator.md"
PYPROJECT = ROOT / "pyproject.toml"


def job_type_spellings(header: str) -> set[str]:
    """Ini spellings magic_enum accepts for JobType, case folded."""
    start = header.index("enum class JobType")
    body = header[start : header.index("};", start)]
    names = re.findall(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*(?:=|,)", body, re.M)
    return {name.lower() for name in names if name != "Unknown"}


def test_guide_job_tokens_match_the_client_enum():
    names = job_type_spellings(HEADER.read_text(encoding="utf-8"))
    guide = INDEX.read_text(encoding="utf-8")
    assert "safe_hyperdynamics" in names
    assert "hyperdynamics" not in names
    assert "finite_difference" in names
    assert "finite_differences" not in names
    assert "dynamics" in names
    assert "molecular_dynamics" not in names
    assert "safe_hyperdynamics" in guide
    assert "finite_differences" in guide
    assert "finite_difference.md" in guide
    assert "molecular_dynamics" in guide
    assert "(use `dynamics`)" in guide
    assert "hyperdynamics" in guide


def test_guide_acceptance_matches_metropolis_probability():
    source = BASIN_CPP.read_text(encoding="utf-8")
    start = source.index("BasinHoppingJob::metropolisProbability")
    function = source[start : source.index("\n}", start)]
    assert "return std::exp(-de / (kB * temperature));" in function
    page = BASIN_MD.read_text(encoding="utf-8")
    assert "exp(-de/(kB*temperature))" in page
    assert "quenching_steps" in page


def test_communicator_script_path_is_the_sge_directory():
    page = COMM.read_text(encoding="utf-8")
    assert 'script_path = "/home/user/eon/tools/clusters/sge"' in page
    assert "sge6.2" not in page
    assert (ROOT / "tools" / "clusters" / "sge").is_dir()
    assert 'eon-server = "eon.server:main"' in PYPROJECT.read_text(encoding="utf-8")
    assert "eon-server" in page
