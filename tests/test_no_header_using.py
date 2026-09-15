"""Job headers in this slice must not inject using eonc:: into the global ns."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
# Job.h still injects using eonc::Job (wide call sites). All other *Job.h
# must stay clean. Optimizer impl headers are also clean; Optimizer.h
# still has the factory alias.
JOB_HEADERS = sorted(
    p
    for p in (ROOT / "include" / "eon").glob("*Job.h")
    if p.name != "Job.h"
) + [
    ROOT / "include" / "eon" / name
    for name in (
        "ConjugateGradients.h",
        "FIRE.h",
        "LBFGS.h",
        "Quickmin.h",
        "SteepestDescent.h",
    )
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
