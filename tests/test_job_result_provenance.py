"""Job envelopes carry optimizer provenance taken from the parameters."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

JOBS = (
    "client/HessianJob.cpp",
    "client/NudgedElasticBandJob.cpp",
    "client/PrefactorJob.cpp",
    "client/MonteCarloJob.cpp",
    "client/FiniteDifferenceJob.cpp",
    "client/MinimizationJob.cpp",
    "client/DynamicsJob.cpp",
    "client/BasinHoppingJob.cpp",
    "client/PointJob.cpp",
    "client/ProcessSearchJob.cpp",
)


def _body(text: str, signature: str, nxt: str) -> str:
    start = text.find(signature)
    assert start != -1, signature
    end = text.find(nxt, start + len(signature))
    assert end != -1, nxt
    return text[start:end]


def test_jobs_fill_optimizer_provenance():
    header = (ROOT / "include/eon/JobResult.h").read_text(encoding="utf-8")
    assert "JobResultProvenance provenanceForJob" in header
    search = _body(
        header,
        "std::string processSearchString() const {",
        "std::string toString() const {",
    )
    assert "provenance.text()" in search

    filled = "env.provenance = provenanceForJob(params);"
    for rel in JOBS:
        text = (ROOT / rel).read_text(encoding="utf-8")
        assert text.count(filled) >= 1, rel
    instanton = (ROOT / "client/InstantonJob.cpp").read_text(encoding="utf-8")
    assert instanton.count(filled) >= 2

    src = (ROOT / "client/JobResultProvenance.cpp").read_text(encoding="utf-8")
    assert "xts_abi_stamp()" in src
    assert "OptType::XTSCI" in src
    assert "WITH_XTSCI" in src
    assert "WITH_RGPOT" in src
    assert "rgpot_pin" in src
    assert "'JobResultProvenance.cpp'" in (
        ROOT / "client/meson.build"
    ).read_text(encoding="utf-8")

    client = (ROOT / "client/ClientEON.cpp").read_text(encoding="utf-8")
    guard = client.find("have_backend")
    write = client.find("provenance.text()")
    assert "provenanceForJob(parameters)" in client
    assert guard != -1 and write != -1 and guard < write

    case = (ROOT / "client/unit_tests/JobResultProvenanceTest.cpp").read_text(
        encoding="utf-8"
    )
    assert "a job envelope takes optimizer provenance from Parameters" in case
    assert "provenanceForJob(params)" in case
    assert 'REQUIRE(text.find("cg optimizer_backend\\n")' in case
    assert 'REQUIRE(xts_text.find("xtsci optimizer_backend\\n")' in case
    assert "optimizer_xts_abi_major" in case
