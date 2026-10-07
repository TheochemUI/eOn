"""A library-route success whose log segfaults is not a passing dimer."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ENGINE = ROOT / "client" / "potentials" / "Rgpot" / "RGPotEngine.cpp"
POT = ROOT / "client" / "potentials" / "Rgpot" / "RgpotPot.cpp"


def library_route_accepted(results: str, log: str) -> bool:
    if "Segmentation fault" in log:
        return False
    return (
        "barrier_reactant_to_product" in results
        and "barrier_product_to_reactant" in results
    )


def test_success_paired_with_segfault_fails():
    results = (
        "0 termination_reason\n"
        "Success termination_reason_text\n"
        "0.3007 barrier_reactant_to_product\n"
        "0.3007 barrier_product_to_reactant\n"
    )
    assert library_route_accepted(results, "Real time: 1.0 seconds\n")
    assert not library_route_accepted(results, "Segmentation fault\n")


def test_hard_exit_is_rearmed_after_the_run():
    engine = ENGINE.read_text(encoding="utf-8")
    start = engine.index("void RGPotEngine::armGroupedExit()")
    body = engine[start : engine.index("\nvoid ", start + 1)]
    assert "call_once" not in body
    assert "::on_exit(api->hard_exit, nullptr);" in body
    worker = POT.read_text(encoding="utf-8")
    stop = worker.index("if (hdr[0] == kStop)")
    exit_at = worker.index("std::exit(0);", stop)
    window = worker[stop:exit_at]
    assert "armGroupedExit()" in window
