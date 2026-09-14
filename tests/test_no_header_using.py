"""Job headers in this slice must not inject using eonc:: into the global ns."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
JOB_HEADERS = [
    ROOT / "include" / "eon" / "TestJob.h",
    ROOT / "include" / "eon" / "MinimizationJob.h",
    ROOT / "include" / "eon" / "PointJob.h",
    ROOT / "include" / "eon" / "DynamicsJob.h",
    ROOT / "include" / "eon" / "HessianJob.h",
    ROOT / "include" / "eon" / "FiniteDifferenceJob.h",
]


def test_job_headers_have_no_file_scope_using():
    offenders: list[str] = []
    for path in JOB_HEADERS:
        for i, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            stripped = line.strip()
            if stripped.startswith("using eonc::") or stripped.startswith(
                "using namespace "
            ):
                offenders.append(f"{path.name}:{i}:{stripped}")
    assert not offenders, "file-scope using remains:\n" + "\n".join(offenders)
