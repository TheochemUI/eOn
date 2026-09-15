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


def test_parameters_h_is_the_aggregate():
    p = (EON / "Parameters.h").read_text()
    assert '#include "ParametersOptions.h"' in p
    assert "struct neb_options_t" not in p
    assert "using neb_options_t" in p
    assert "private:" in p
    assert "main_options_t main_options;" not in p
    assert "main_options_t main_options_{};" in p


def test_parameters_options_defines_job_and_pot_structs():
    o = (EON / "ParametersOptions.h").read_text()
    for name in (
        "potential_options_t",
        "neb_options_t",
        "dimer_options_t",
        "metatomic_options_t",
        "oh_tst_options_t",
    ):
        assert f"struct {name}" in o or f"using {name}" in o
