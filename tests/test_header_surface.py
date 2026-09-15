"""Installed public headers do not gate types on WITH_* or include mpi.h."""

from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
EON = ROOT / "include" / "eon"


@pytest.mark.parametrize(
    "rel,banned",
    [
        ("Parameters.h", ("mpi.h", "EONMPI")),
        ("EigenmodeStrategy.h", ("#ifdef WITH_GPRD", "AtomicGPDimer.h", "variant")),
        ("GPSurrogateJob.h", ("WITH_CATLEARN", "CatLearnPot.h")),
        ("ARTnSaddleSearch.h", ("WITH_ARTN", "ARTnResource.h")),
        ("Potential.h", ("Parameters.h",)),
        ("Matter.h", ("Parameters.h",)),
    ],
)
def test_core_headers_have_stable_surface(rel, banned):
    text = (EON / rel).read_text()
    for token in banned:
        assert token not in text, f"{rel} still contains {token}"


def test_parameters_mpi_is_optional_header():
    assert (EON / "ParametersMpi.h").is_file()
    text = (EON / "ParametersMpi.h").read_text()
    assert "mpi.h" in text
